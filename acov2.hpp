
#ifndef ACO_HPP
#define ACO_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper.hpp"
#include "local_search.hpp"
#include "nn.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Data structures 
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef unordered_set<int> uset;
typedef map<int, bool> mapb;
typedef pair<int, int> pi;
#define F first
#define S second
#define PB push_back
#define MP make_pair
#define all(a) (a).begin(), (a).end() // Aplicar a toda la estructura, e.g. sort(all(a))


vi determine_candidate_list(int n, vvi dist_matrix, uset unvisited, int candidate_list_size, int current_city){
    uset unvisited_copy = unvisited;
    vi candidate_list(candidate_list_size);
    for (int i = 0; i < candidate_list_size; i++){
        int min_dist = INF;
        int min_dist_city = -1;
        uset::iterator it;
        for (it = unvisited_copy.begin(); it != unvisited_copy.end(); ++it){
            if (dist_matrix[current_city][*it] < min_dist){
                min_dist = dist_matrix[current_city][*it];
                min_dist_city = *it;
            }
        }
        candidate_list[i] = min_dist_city;
        unvisited_copy.erase(min_dist_city);
    }
    if (unvisited.size() != n - 1) cout << "ERROR: unvisited size is not n - 1" << endl;
    return candidate_list;
}

vvi determine_initial_candidate_lists(int n, vvi dist_matrix, int candidate_list_size){
    vvi candidate_lists(n);
    uset unvisited;
    for (int i = 0; i < n; i++){
        unvisited.insert(i);
    }
    for (int i = 0; i < n; i++){
        unvisited.erase(i);
        candidate_lists[i] = determine_candidate_list(n, dist_matrix, unvisited, candidate_list_size, i);
        unvisited.insert(i);
    }
    return candidate_lists;
}


int choose_with_list(int n, vvi dist_matrix, vector<vector<float>> pheromone_matrix, int current_city, vi candidate_list, vector<bool> visited, int cl_size, int Q0, float beta){
    // Choose next city
    // 95% of the time, choose the best edge according to pheromone (exploitation)
    // 5% of the time, choose a random edge based on probability (exploration)
    int next_city = -1;
    if (rand() % 100 < Q0){         // Exploitation
        // Choose best edge
        float best_transition = -1.0f;
        for (int i = 0; i < cl_size; i++){
            int city = candidate_list[i];
            if (visited[city]) continue;
            float edge_transition = pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
            if (edge_transition > best_transition){
                best_transition = edge_transition;
                next_city = city;
            }
        }
    } else {                     // Exploration
        vector<float> probabilities(n, 0);    // Probabiliy of visiting city l from the current city
        float denominator = 0.0f;
        for (int i = 0; i < cl_size; i++){
            int city = candidate_list[i];
            if (visited[city]) continue;
            denominator += pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
        }

        for (int i = 0; i < cl_size; i++){
            int city = candidate_list[i];
            if (visited[city]) continue;
            probabilities[city] = pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta) / denominator;
        }
        
        // Create a distribution based on the probabilities
        discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
        // Choose a city based on the distribution
        next_city = distribution(rng);
    }
    return next_city;
}



int choose_with_uset(int n, vvi dist_matrix, vector<vector<float>> pheromone_matrix, int current_city, vector<bool> visited, int Q0, float beta){
    // Choose next city
    // 95% of the time, choose the best edge according to pheromone (exploitation)
    // 5% of the time, choose a random edge based on probability (exploration)
    int next_city = -1;
    if (rand() % 100 < Q0){         // Exploitation
        // Choose best edge
        float best_transition = -1.0f;
        for (int i = 0; i < n; i++){
            if (visited[i]) continue;
            float edge_transition = pheromone_matrix[current_city][i] * pow(1.0 / dist_matrix[current_city][i], beta);
            if (edge_transition > best_transition){
                best_transition = edge_transition;
                next_city = i;
            }
        }
    } else {                     // Exploration
        vector<float> probabilities(n, 0);    // Probabiliy of visiting city l from the current city
        float denominator = 0.0f;
        for (int i = 0; i < n; i++){
            if (visited[i]) continue;
            denominator += pheromone_matrix[current_city][i] * pow(1.0 / dist_matrix[current_city][i], beta);
        }
        
        for (int i = 0; i < n; i++){
            if (visited[i]) continue;
            probabilities[i] = pheromone_matrix[current_city][i] * pow(1.0 / dist_matrix[current_city][i], beta) / denominator;
        }
        
        // Create a distribution based on the probabilities
        discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
        // Choose a city based on the distribution
        next_city = distribution(rng);
    }
    return next_city;
}

vi ACO_candidate(Problem * problem)
{
    float start_time = (float) clock();

    int n = problem->n;
    int numIterations = 0;         // End condition for the main loop
    int candidate_list_size = 15;  // Size of the candidate list
    float alpha = 2.0f;
    float phi = 0.1f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.1f;               // Evaporation rate
    int Q0 = 95;                    // Probability of choosing the best solution
    int m = 10;                     // Number of ants

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    //float tau0 = 1.0f / (n * length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1000 * n / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    float tau0 = 1 / (n * (float) length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 100.0f;
    float tau;

    // Pheromone matrix
    vector<vector<float>> pheromone_matrix(n, vector<float>(n, tau0));
    vvi candidate_lists = determine_initial_candidate_lists(n, problem->dist_matrix, candidate_list_size);


    // m = 10 from original paper
    // Ants and their paths
    vector<vi> ants(m, vi(n, 0));

    int best_ant_length = INF;      // Global histórico
    vi best_ant_sol;                // Global histórico    

    int next_city;               

    while (((float) clock() - start_time) / CLOCKS_PER_SEC < 180){         // Run for 3 minutes
        vector<vector<bool>> visited = vector<vector<bool>>(m, vector<bool>(n, false)); // Unordered set of visited cities, one set for each ant

        // Put each ant in a random initial city
        for (int k = 0; k < m; k++){        // For each ant
            int start_city = rand() % n;
            ants[k][0] = start_city;
            visited[k][start_city] = true;     // Remove initial city from unvisited
        }

        // Unordered set of visited cities, one set for each ant
        

        for (int j = 0; j < n - 1; j++){     // For each city
            vi candidate_list = candidate_lists[j];
            for (int k = 0; k < m; k++){        // For each ant
                vector<bool> ant_visited = visited[k];
                int current_city = ants[k][j];

                // Choose next city
                // 95% of the time, choose the best edge according to pheromone (exploitation)
                // 5% of the time, choose a random edge based on probability (exploration)
                next_city = -1;
                float thereAreCandidates = false;
                for (int l = 0; l < candidate_list_size; l++){
                    if (!ant_visited[candidate_list[l]]){
                        thereAreCandidates = true;
                        break;
                    }
                }

                if (thereAreCandidates){
                    next_city = choose_with_list(n, problem->dist_matrix, pheromone_matrix, current_city, candidate_list, ant_visited, candidate_list_size, Q0, beta);          
                } else {
                    next_city = choose_with_uset(n, problem->dist_matrix, pheromone_matrix, current_city, visited[k], Q0, beta);    
                }

                
                // Update ant's path
                ants[k][j + 1] = next_city;
                // Remove next city from unvisited
                visited[k][next_city] = true;

                // Update pheromone matrix
                tau = (1 - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
            }
        }
        
        int random_ant = rand() % m;
        for (int k = 0; k < m; k++){
            // Update pheromone from last city to first city
            tau = (1 - rho) * pheromone_matrix[ants[k][n - 1]][ants[k][0]] + rho * tau0;
            pheromone_matrix[ants[k][n - 1]][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][ants[k][n - 1]] = tau;

            if (k == random_ant){
                //ants[k] = two_five_opt(problem, ants[k]);
                //ants[k] = two_opt_greedy(problem, ants[k]);
            }

            // Calculate cost of each ant's path
            // TODO: improve updating it in the algorithm?
            int ant_length = compute_length(problem, ants[k]);
            if (ant_length < best_ant_length){
                best_ant_length = ant_length;
                best_ant_sol = ants[k];
                /*
                for (int i = 0; i < n; i++){
                    cout << best_ant_sol[i] << " ";
                }
                */
                cout << "Best ant length: " << best_ant_length << endl;
            }
        }

        best_ant_sol = two_five_opt(problem, best_ant_sol);
        best_ant_sol = two_opt_greedy(problem, best_ant_sol);
        //cout << "Best ant length: " << best_ant_length << endl;

        // Update pheromone matrix globally
        // Evaporate for every edge
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                pheromone_matrix[i][j] *= (1 - phi);
            }
        }

        // Update pheromone matrix with best ant's path
        float constant = phi / best_ant_length;
        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            pheromone_matrix[node_1][node_2] += constant;
            pheromone_matrix[node_2][node_1] += constant;
        }
        pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] += constant;
        pheromone_matrix[best_ant_sol[0]][best_ant_sol[n - 1]] += constant;
        numIterations++;
    }
    return best_ant_sol;
}

#endif /* ACO_HPP */
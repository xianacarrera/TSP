
#ifndef ACO_HPP
#define ACO_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper.hpp"
#include "nn.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Data structures 
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef unordered_set<int> uset;
typedef pair<int, int> pi;
#define F first
#define S second
#define PB push_back
#define MP make_pair
#define all(a) (a).begin(), (a).end() // Aplicar a toda la estructura, e.g. sort(all(a))


vi dynamic_ACO(Problem * problem){
    int n = problem->n;
    int numIterations = 1250;       // End condition for the main loop
    float alpha = 0.1f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.1f;               // Evaporation rate
    int Q0 = 95;                    // Probability of choosing the best solution

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    float tau0 = 1.0f / (n * compute_length(problem, NN(problem, rand() % n))); // Initial pheromone value
    //float tau0 = 0.1f;
    float tau;

    // Pheromone matrix
    vector<vector<float>> pheromone_matrix(n, vector<float>(n, tau0));


    // m = 10 from original paper
    // Ants and their paths
    vector<vi> ants(10, vi(n, 0));

    int best_ant_length = INF;      // Global histórico
    vi best_ant_sol;                // Global histórico    
    int worst_ant_length = 0;       // Global histórico
    vi worst_ant_sol;               // Global histórico

    int next_city;               

    for (int i = 0; i < numIterations; i++){         // While (!end_condition)
        // Put each ant in a random initial city
        for (int k = 0; k < 10; k++){        // For each ant
            ants[k][0] = rand() % n;
        }

        // Unordered set of unvisited cities, one set for each ant
        vector<uset> unvisited(10);
        for (int k = 0; k < 10; k++){
            for (int j = 0; j < n; j++){
                unvisited[k].insert(j);
            }
            unvisited[k].erase(ants[k][0]);     // Remove initial city from unvisited
        }

        for (int k = 0; k < 10; k++){        // For each ant
            uset ant_unvisited = unvisited[k];
            for (int j = 0; j < n - 1; j++){     // For each city
                int current_city = ants[k][j];


                // Choose next city
                // 95% of the time, choose the best edge according to pheromone (exploitation)
                // 5% of the time, choose a random edge based on probability (exploration)
                next_city = 0;
                if (rand() % 100 < Q0){         // Exploitation
                    // Choose best edge
                    float best_transition = -1.0f;
                    uset::iterator it;
                    for (it = ant_unvisited.begin(); it != ant_unvisited.end(); ++it) {
                        // edge_transition = tau * (1 / cost)^beta
                        float edge_transition = pheromone_matrix[current_city][*it] * pow(1.0 / problem->dist_matrix[current_city][*it], beta);
                        if (edge_transition > best_transition){
                            best_transition = edge_transition;
                            next_city = *it;
                        }
                    }
                } else {                     // Exploration
                    vector<float> probabilities(n, 0);    // Probabiliy of visiting city l from the current city
                    float denominator = 0.0f;
                    uset::iterator l;
                    for (l = ant_unvisited.begin(); l != ant_unvisited.end(); ++l) {
                        denominator += pheromone_matrix[current_city][*l] * pow(1.0 / problem->dist_matrix[current_city][*l], beta);
                    }
                    
                    uset::iterator l2;
                    for (l2 = ant_unvisited.begin(); l2 != ant_unvisited.end(); ++l2) {
                        probabilities[*l2] = pheromone_matrix[current_city][*l2] * pow(1.0 / problem->dist_matrix[current_city][*l2], beta) / denominator;
                    }
                    
                    // Create a distribution based on the probabilities
                    discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
                    // Choose a city based on the distribution
                    next_city = distribution(rng);
                }
                
                // Update ant's path
                ants[k][j + 1] = next_city;
                // Remove next city from unvisited
                ant_unvisited.erase(next_city);

                // Update pheromone matrix
                tau = pheromone_matrix[current_city][next_city];
                tau = (1 - rho * tau) * tau + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
            }
            // Update pheromone from last city to first city
            tau = pheromone_matrix[next_city][ants[k][0]];
            tau = (1 - rho * tau) * tau + rho * tau0;
            pheromone_matrix[next_city][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][next_city] = tau;
            
            // Calculate cost of each ant's path
            // TODO: improve updating it in the algorithm?
            int ant_length = compute_length(problem, ants[k]);
            if (ant_length < best_ant_length){
                best_ant_length = ant_length;
                best_ant_sol = ants[k];
            } 
            if (ant_length > worst_ant_length){
                worst_ant_length = ant_length;
                worst_ant_sol = ants[k];
            }
        }
        // Update pheromone matrix globally for best tour
        float constant = alpha / best_ant_length;
        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            tau = pheromone_matrix[node_1][node_2];
            tau = (1 - alpha * tau) * tau + constant;
            pheromone_matrix[node_1][node_2] = tau;
            pheromone_matrix[node_2][node_1] = tau;
        }
        //?????????
        tau = pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]];
        tau = (1 - alpha * tau) * tau + constant;
        pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] = tau;
        pheromone_matrix[best_ant_sol[0]][best_ant_sol[n - 1]] = tau;


        // Update pheromone matrix globally worst tour
        constant = alpha / worst_ant_length;
        for (int j = 0; j < n - 1; j++){
            int node_1 = worst_ant_sol[j];
            int node_2 = worst_ant_sol[j + 1];
            tau = pheromone_matrix[node_1][node_2];
            tau = (1 - alpha * tau) * tau + constant;
            pheromone_matrix[node_1][node_2] = tau;
            pheromone_matrix[node_2][node_1] = tau;
        }
        //?????????
        tau = (1 - alpha) * pheromone_matrix[worst_ant_sol[n - 1]][worst_ant_sol[0]] + constant;
        pheromone_matrix[worst_ant_sol[n - 1]][worst_ant_sol[0]] = tau;
        pheromone_matrix[worst_ant_sol[0]][worst_ant_sol[n - 1]] = tau;
    }
    return best_ant_sol;
}



vi determine_candidate_list(int n, vvi dist_matrix, uset unvisited, int candidate_list_size, int current_city){
    uset unvisited_copy(all(unvisited));
    vi candidate_list(candidate_list_size);
    for (int i = 0; i < candidate_list_size; i++){
        int min_dist = INT_MAX;
        int min_dist_city = -1;
        uset::iterator it;
        for (it = unvisited_copy.begin(); it != unvisited_copy.end(); ++it){
            if (dist_matrix[current_city][*it] < min_dist){
                min_dist = dist_matrix[current_city][*it];
                min_dist_city = *it;
            }
        }
        candidate_list.push_back(min_dist_city);
        unvisited_copy.erase(min_dist_city);
    }
    return candidate_list;
}

vvi determine_initial_candidate_lists(int n, vvi dist_matrix, int candidate_list_size){
    vvi candidate_lists(n, vi(candidate_list_size));
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


pi choose_with_list(int n, vvi dist_matrix, vector<vector<float>> pheromone_matrix, int current_city, vi candidate_list, vector<bool> cl_visited, int cl_size, int Q0, float beta){
    // Choose next city
    // 95% of the time, choose the best edge according to pheromone (exploitation)
    // 5% of the time, choose a random edge based on probability (exploration)
    pi info;
    if (rand() % 100 < Q0){         // Exploitation
        // Choose best edge
        float best_transition = -1.0f;
        for (int i = 0; i < cl_size; i++){
            if (cl_visited[i]) continue;
            int city = candidate_list[i];
            float edge_transition = pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
            if (edge_transition > best_transition){
                best_transition = edge_transition;
                info.F = city;
                info.S = i;
            }
        }
    } else {                     // Exploration
        vector<float> probabilities(n, 0);    // Probabiliy of visiting city l from the current city
        float denominator = 0.0f;
        for (int i = 0; i < cl_size; i++){
            if (cl_visited[i]) continue;
            int city = candidate_list[i];
            denominator += pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
        }

        for (int i = 0; i < cl_size; i++){
            if (cl_visited[i]) continue;
            int city = candidate_list[i];
            probabilities[city] = pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta) / denominator;
        }
        
        // Create a distribution based on the probabilities
        discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
        // Choose a city based on the distribution
        info.F = distribution(rng);
        info.S = find(candidate_list.begin(), candidate_list.end(), info.F) - candidate_list.begin();
    }
    return info;
}



int choose_with_uset(int n, vvi dist_matrix, vector<vector<float>> pheromone_matrix, int current_city, uset unvisited, int Q0, float beta){
    // Choose next city
    // 95% of the time, choose the best edge according to pheromone (exploitation)
    // 5% of the time, choose a random edge based on probability (exploration)
    int next_city = 0;
    if (rand() % 100 < Q0){         // Exploitation
        // Choose best edge
        float best_transition = -1.0f;
        uset::iterator it;
        for (it = unvisited.begin(); it != unvisited.end(); ++it) {
            // edge_transition = tau * (1 / cost)^beta
            float edge_transition = pheromone_matrix[current_city][*it] * pow(1.0 / dist_matrix[current_city][*it], beta);
            if (edge_transition > best_transition){
                best_transition = edge_transition;
                next_city = *it;
            }
        }
    } else {                     // Exploration
        vector<float> probabilities(n, 0);    // Probabiliy of visiting city l from the current city
        float denominator = 0.0f;
        uset::iterator l;
        for (l = unvisited.begin(); l != unvisited.end(); ++l) {
            denominator += pheromone_matrix[current_city][*l] * pow(1.0 / dist_matrix[current_city][*l], beta);
        }
        
        uset::iterator l2;
        for (l2 = unvisited.begin(); l2 != unvisited.end(); ++l2) {
            probabilities[*l2] = pheromone_matrix[current_city][*l2] * pow(1.0 / dist_matrix[current_city][*l2], beta) / denominator;
        }
        
        // Create a distribution based on the probabilities
        discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
        // Choose a city based on the distribution
        next_city = distribution(rng);
    }
    return next_city;
}


/*
  Candidate list [2] is a strategy that tries to 
improve the performance of an algorithm. It is a 
static data structure that lists a limited number of 
preferred closest cities to be visited ordered by 
increasing distance. In ACS, the decision process to 
move from city r to city s is first need to consider 
those preferred cities listed in the candidate list. Only 
when an ant cannot find suitable city to choose then 
the decision to choose a city will consider those 
which are outside the candidate list.
*/

vi ACO_candidate_list(Problem * problem)
{
    int n = problem->n;
    int candidate_list_size = 15;
    int numIterations = 1250;       // End condition for the main loop
    float alpha = 0.1f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.1f;               // Evaporation rate
    int Q0 = 95;                    // Probability of choosing the best solution

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    float tau0 = 1.0f / (n * compute_length(problem, NN(problem, rand() % n))); // Initial pheromone value
    //float tau0 = 0.1f;
    float tau;

    // Pheromone matrix
    vector<vector<float>> pheromone_matrix(n, vector<float>(n, tau0));
    vvi candidate_lists = determine_initial_candidate_lists(n, problem->dist_matrix, candidate_list_size);

    // m = 10 from original paper
    // Ants and their paths
    vector<vi> ants(10, vi(n, 0));

    int best_ant_length = INF;      // Global histórico
    vi best_ant_sol;                // Global histórico    

    int next_city;               

    for (int i = 0; i < numIterations; i++){         // While (!end_condition)

        vector<umap> unvisited(10);

        for (int k = 0; k < 10; k++){        // For each ant
            vi candidate_lists_indexes = vi(n, 0);

            ants[k][0] = rand() % n;         // Set initial city
            
            for (int j = 0; j < n; j++){
                unvisited[k][j] = false;
            }
            unvisited[k][j] = true;     // Remove initial city from unvisited

            if (unvisited[k].size() != n - 1){
                cout << "Error: unvisited size is not n - 1" << endl;
                exit(1);
            }

            for (int j = 0; j < n - 1; j++){     // For each city
                int current_city = ants[k][j];

                int next_city;
                if (candidate_lists_indexes[current_city] == candidate_list_size){
                    next_city = choose_with_uset(n, problem->dist_matrix, pheromone_matrix, current_city, unvisited[k], Q0, beta);
                } else {
                    vi candidate_list = candidate_lists[current_city];
                    auto info = choose_with_list(n, problem->dist_matrix, pheromone_matrix, current_city, candidate_list, candidate_lists_visited[current_city], candidate_list_size, Q0, beta);
                    next_city = info.F;
                    candidate_lists_indexes[current_city]++;
                    candidate_lists_visited[current_city][info.S] = true;
                }


                // Update ant's path
                ants[k][j + 1] = next_city;
                // Remove next city from unvisited
                unvisited[k].erase(next_city);

                // Update pheromone matrix
                tau = (1 - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
            }
            // Update pheromone from last city to first city
            tau = (1 - rho) * pheromone_matrix[next_city][ants[k][0]] + rho * tau0;
            pheromone_matrix[next_city][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][next_city] = tau;
            
            // Calculate cost of each ant's path
            // TODO: improve updating it in the algorithm?
            int ant_length = compute_length(problem, ants[k]);
            if (ant_length < best_ant_length){
                best_ant_length = ant_length;
                best_ant_sol = ants[k];
            }
        }
        // Update pheromone matrix globally

        float constant = alpha / best_ant_length;
        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            tau = (1 - alpha) * pheromone_matrix[node_1][node_2] + constant;
            pheromone_matrix[node_1][node_2] = tau;
            pheromone_matrix[node_2][node_1] = tau;
        }
        //?????????
        tau = (1 - alpha) * pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] + constant;
        pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] = tau;
        pheromone_matrix[best_ant_sol[0]][best_ant_sol[n - 1]] = tau;
    }
    return best_ant_sol;
}

#endif /* ACO_HPP */
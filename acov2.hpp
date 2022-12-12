
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



queue<int> determine_ESACO_candidate_list(int n, vvi dist_matrix, uset unvisited, int candidate_list_size, int current_city){
    uset unvisited_copy = unvisited;
    queue<int> candidate_list;
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
        candidate_list.push(min_dist_city);
        unvisited_copy.erase(min_dist_city);
    }
    if (unvisited.size() != n - 1) cout << "ERROR: unvisited size is not n - 1" << endl;
    return candidate_list;
}

vector<queue<int>> determine_ESACO_initial_candidate_lists(int n, vvi dist_matrix, int candidate_list_size){
    vector<queue<int>> candidate_lists(n);
    uset unvisited;
    for (int i = 0; i < n; i++){
        unvisited.insert(i);
    }
    for (int i = 0; i < n; i++){
        unvisited.erase(i);
        candidate_lists[i] = determine_ESACO_candidate_list(n, dist_matrix, unvisited, candidate_list_size, i);
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
        for (int i = 0; i < cl_size; i++){    // Do not choose the closest city
            int city = candidate_list[i];
            if (current_city == city) continue;
            if (visited[city]) continue;
            float edge_transition = pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
            if (edge_transition > best_transition){
                best_transition = edge_transition;
                next_city = city;
            }
        }
    } else {                     // Exploration
        float denominator = 0.0f;
        for (int i = 0; i < cl_size; i++){
            int city = candidate_list[i];
            if (current_city == city) continue;
            if (visited[city]) continue;
            denominator += pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
        }

        // Get a random number between 0 and 1
        float random = (float)rand() / (float) (RAND_MAX);
        float cum_prob = 0.0f;
        for (int i = 0; i < cl_size; i++){
            int city = candidate_list[i];
            if (current_city == city) continue;
            if (visited[city]) continue;
            cum_prob += pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta) / denominator;
            if (random <= cum_prob){
                return city;
            }
        }
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
            if (current_city == i) continue;
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
            if (current_city == i) continue;
            if (visited[i]) continue;
            denominator += pheromone_matrix[current_city][i] * pow(1.0 / dist_matrix[current_city][i], beta);
        }
        
        // Get a random number between 0 and 1
        float random = (float)rand() / (float) (RAND_MAX);
        float cum_prob = 0.0f;
        for (int i = 0; i < n; i++){
            if (current_city == i) continue;
            if (visited[i]) continue;
            cum_prob += pheromone_matrix[current_city][i] * pow(1.0 / dist_matrix[current_city][i], beta) / denominator;
            if (random <= cum_prob){
                return i;
            }
        }
    }
    return next_city;
}

vi ACO_candidate(Problem * problem, float alpha, float beta, float rho, int Q0)
{
    float start_time = (float) clock();

    int n = problem->n;
    int numIterations = 0;         // End condition for the main loop
    int candidate_list_size = n / 4;  // Size of the candidate list
    /*
    float alpha = 2.0f;
    float phi = 0.1f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.1f;               // Evaporation rate
    int Q0 = 95;                    // Probability of choosing the best solution
    */
    int m = 10;                     // Number of ants

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    //float tau0 = 1.0f / (n * length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1000 * n / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    float tau0 = 1.0f / (n * (float) length_NN(problem, rand() % n)); // Initial pheromone value
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
            // Get a random number between 0 and n - 1
            int start_city = (rand() / RAND_MAX) * n;
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
                for (int l = 1; l < candidate_list_size; l++){    // Do not consider the closest city
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
                pheromone_matrix[i][j] *= (1 - alpha);
            }
        }

        // Update pheromone matrix with best ant's path
        float constant = alpha / best_ant_length;
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

vvi ESACO_update_candidate_lists(vvi old_candidate_lists, vector<vector<float>> pheromone_matrix, int n, vi best_ant_sol){
    // Update candidate lists
    vvi new_candidate_lists(n);
    for (int i = 0; i < n; i++){
        vi candidate_list = old_candidate_lists[i];

        int best_pheromone = -1;
        int best_neighbour = -1;
        for (int j = 0; j < n; j++){
            if (i == j) continue;
            if (pheromone_matrix[i][j] > best_pheromone){
                best_pheromone = pheromone_matrix[i][j];
                best_neighbour = j;
            }
        }

        // Remove last node of candidates_list
        candidate_list.pop_back();
        // Insert best neighbour at the beginning of the list
        candidate_list.insert(candidate_list.begin(), best_neighbour);
        // Remove last node of candidates_list
        candidate_list.pop_back();
        // Insert successor of i in best tour as first in the candidates_list
        auto idx = find(best_ant_sol.begin(), best_ant_sol.end(), i) + 1 - best_ant_sol.begin();
        if (idx == n) idx = 0;
        int successor = best_ant_sol[idx];
        candidate_list.insert(candidate_list.begin(), successor);

        new_candidate_lists[i] = candidate_list;
    }
    return new_candidate_lists;
}


/*
 * 1. Candidate lists as queues
 */
vi ESACO(Problem * problem, float phi, float beta, float rho, int Q0)
{
    float start_time = (float) clock();

    int n = problem->n;
    int numIterations = 0;         // End condition for the main loop
    int candidate_list_size = n / 4;  // Size of the candidate list
    /*
    float phi = 0.15f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.15f;               // Evaporation rate
    int Q0 = 95;                    // Probability of choosing the best solution
    */
    int m = 10;                     // Number of ants

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    float tau0 = n / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1000 * n / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1 / (n * (float) length_NN(problem, rand() % n)); // Initial pheromone value
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

    //while ((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){         // Run for 3 minutes
    while((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){
        vector<vector<bool>> visited = vector<vector<bool>>(m, vector<bool>(n, false)); // Unordered set of visited cities, one set for each ant

        // Put each ant in a random initial city
        for (int k = 0; k < m; k++){        // For each ant
            int start_city = rand() % n;
            ants[k][0] = start_city;
            visited[k][start_city] = true;     // Remove initial city from unvisited
        }

        // Unordered set of visited cities, one set for each ant
        

        for (int j = 0; j < n - 1; j++){     // For each city
            for (int k = 0; k < m; k++){        // For each ant
                vector<bool> ant_visited = visited[k];
                int current_city = ants[k][j];
                vi candidate_list = candidate_lists[current_city];

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
                
                for(int i = 0; i < j + 1; i++) {
                    if (ants[k][i] == ants[k][j + 1]) {
                        cout << "ERROR: Ant " << k << " has visited city " << ants[k][i] << " twice" << endl;
                        cout << "j" << j << endl;
                        cout << "Ant's path: ";
                        for (int i = 0; i < n; i++) {
                            cout << ants[k][i] << " ";
                        }
                        cout << thereAreCandidates << endl;
                        exit(1);
                    }
                }
                
                // Update pheromone matrix
                tau = (1 - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
            }
        }
        
        for (int k = 0; k < m; k++){
            // Update pheromone from last city to first city
            tau = (1 - rho) * pheromone_matrix[next_city][ants[k][0]] + rho * tau0;
            pheromone_matrix[next_city][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][next_city] = tau;
            
            ants[k] = two_opt_greedy(problem, ants[k]);

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

        if (best_ant_length < 1.1 * problem->best_known){
            best_ant_sol = two_five_opt(problem, best_ant_sol);
        } else {
            best_ant_sol = two_opt_simplified(problem, best_ant_sol);
        }
        int old_ant_length = best_ant_length;
        best_ant_length = compute_length(problem, best_ant_sol);
        if (best_ant_length < old_ant_length){
            cout << "Best ant length: " << best_ant_length << endl;
        }

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

        candidate_lists = ESACO_update_candidate_lists(candidate_lists, pheromone_matrix, n, best_ant_sol);
    }
    return best_ant_sol;
}



vi ACO_v2_with_candidate(Problem * problem, float alpha, float beta, float rho, int Q0)
{
    float start_time = (float) clock();

    int n = problem->n;
    int numIterations = 0;       // End condition for the main loop
    int candidate_list_size = 25;
    /*
    float alpha = 2.0f;
    float phi = 0.1f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.1f;               // Evaporation rate
    int Q0 = 95;                    // Probability of choosing the best solution
    */
    alpha = 0.10f;
    rho = 0.1f;
    beta = 2.5f;                     // Number of ants
    Q0 = 98;
    int m = 12;                     // Number of ants

    vi best_ant_sol;                // Global histórico    
    int best_ant_length = INF;      // Global histórico

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    //float tau0 = 1.0f / (n * length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1000 * n / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    float tau0 = 1.0f / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 0.01f;
    //float tau0 = 1.0f / best_ant_length;
    //float tau0 = 100.0f;
    float tau;

    // Pheromone matrix
    vector<vector<float>> pheromone_matrix(n, vector<float>(n, tau0));


    // m = 10 from original paper
    // Ants and their paths
    vector<vi> ants(m, vi(n, 0));



    int next_city;              
    
    vvi candidate_lists = determine_initial_candidate_lists(n, problem->dist_matrix, candidate_list_size);
 
    while((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){
        vector<vector<bool>> visited = vector<vector<bool>>(m, vector<bool>(n, false)); // Unordered set of visited cities, one set for each ant
        vi ants_lengths(m, 0);      // Length of each ant's path

        // Put each ant in a random initial city
        for (int k = 0; k < m; k++){        // For each ant
            int start_city = rand() % n;
            ants[k][0] = start_city;
            visited[k][start_city] = true;     // Remove initial city from unvisited
        }

        for (int j = 0; j < n - 1; j++){     // For each city
            for (int k = 0; k < m; k++){        // For each ant
                vector<bool> ant_visited = visited[k];
                int current_city = ants[k][j];
                vi candidate_list = candidate_lists[current_city];

                // Choose next city
                // 95% of the time, choose the best edge according to pheromone (exploitation)
                // 5% of the time, choose a random edge based on probability (exploration)
                next_city = -1;
                float thereAreCandidates = false;
                for (int l = 1; l < candidate_list_size; l++){      // Do not use the closest city
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
                ants_lengths[k] += problem->dist_matrix[current_city][next_city];
                /*
                for(int i = 0; i < j + 1; i++) {
                    if (ants[k][i] == ants[k][j + 1]) {
                        cout << "ERROR: Ant " << k << " has visited city " << ants[k][i] << " twice" << endl;
                        cout << "j" << j << endl;
                        cout << "Ant's path: ";
                        for (int i = 0; i < n; i++) {
                            cout << ants[k][i] << " ";
                        }
                        cout << thereAreCandidates << endl;
                        exit(1);
                    }
                }*/
                
                // Update pheromone matrix
                tau = (1 - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
            }
        }


 //       if ((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180) return best_ant_sol;
        
        
        for (int k = 0; k < m; k++){
            // Update pheromone from last city to first city
            tau = (1 - rho) * pheromone_matrix[ants[k][n-1]][ants[k][0]] + rho * tau0;
            pheromone_matrix[ants[k][n-1]][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][ants[k][n-1]] = tau;
            ants_lengths[k] += problem->dist_matrix[ants[k][n-1]][ants[k][0]];

            if (ants_lengths[k] < best_ant_length){
                best_ant_length = ants_lengths[k];
                best_ant_sol = ants[k];
                /*
                for (int i = 0; i < n; i++){
                    cout << best_ant_sol[i] << " ";
                }
                */
                cout << "Best ant length 1: " << best_ant_length << endl;
            }
        }

        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC >= 180) return best_ant_sol;
        int random_ant = rand() % m;
        cout << "old length: " << ants_lengths[random_ant] << endl;
        //cout << "computed length " << compute_length(problem, ants[random_ant]) << endl;
        pair<vi, int> new_ant_sol_and_length = two_five_opt_2(problem, ants[random_ant]);
        ants[random_ant] = new_ant_sol_and_length.first;
        ants_lengths[random_ant] += new_ant_sol_and_length.second;
        cout << "new length: " << ants_lengths[random_ant] << endl;
        //cout << "computed length " << compute_length(problem, ants[random_ant]) << endl;

        if (ants_lengths[random_ant] < best_ant_length){
            best_ant_length = ants_lengths[random_ant];
            best_ant_sol = ants[random_ant];
            cout << "Best ant length: " << best_ant_length << endl;
        }

        pair<vi, int> best_pair;
        if (numIterations % 3 == 0){
            best_pair = two_five_opt_2(problem, best_ant_sol);
        } else {
            best_pair = two_opt_simplified_2(problem, best_ant_sol);
        }
        best_ant_sol = best_pair.first;
        best_ant_length += best_pair.second;
        cout << "Best ant length 2: " << best_ant_length << endl;

        // Update pheromone matrix globally
        // Evaporate for every edge
        
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                pheromone_matrix[i][j] *= (1 - alpha);
            }
        }

        // Update pheromone matrix with best ant's path
        float constant = alpha / best_ant_length;
        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            pheromone_matrix[node_1][node_2] += constant;
            pheromone_matrix[node_2][node_1] += constant;
        }
        /*
               float constant = alpha / best_ant_length;

        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            tau = (1 - alpha) * pheromone_matrix[node_1][node_2] + constant;
            pheromone_matrix[node_1][node_2] = tau;
            pheromone_matrix[node_2][node_1] = tau;
        }*/
        tau = (1 - alpha) * pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] + constant;
        pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] = tau;
        pheromone_matrix[best_ant_sol[0]][best_ant_sol[n - 1]] = tau;
        numIterations++;
        
    }
    return best_ant_sol;
}


vi ESACO_v2(Problem * problem, float alpha, float beta, float rho, int Q0)
{
    float start_time = (float) clock();

    int n = problem->n;
    int numIterations = 0;       // End condition for the main loop
    int candidate_list_size = 40;
    /*
    float alpha = 2.0f;
    float phi = 0.1f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.1f;               // Evaporation rate
    int Q0 = 95;                    // Probability of choosing the best solution
    *//*
    alpha = 0.10f;
    rho = 0.1f;
    beta = 2.5f;                     // Number of ants
    Q0 = 98;
    */
    int m = 10;                     // Number of ants

    vi best_ant_sol = best_NN(problem);                // Global histórico    
    int best_ant_length = compute_length(problem, best_ant_sol);      // Global histórico

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    //float tau0 = 1.0f / (n * length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1000 * n / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1.0f / (n * (float) length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 0.01f;
    float tau0 = 1.0f / (n * best_ant_length);
    //float tau0 = 100.0f;
    float tau;

    // Pheromone matrix
    vector<vector<float>> pheromone_matrix(n, vector<float>(n, tau0));


    // m = 10 from original paper
    // Ants and their paths
    vector<vi> ants(m, vi(n, 0));



    int next_city;              
    
    vvi candidate_lists = determine_initial_candidate_lists(n, problem->dist_matrix, candidate_list_size);
 
    while((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){
        vector<vector<bool>> visited = vector<vector<bool>>(m, vector<bool>(n, false)); // Unordered set of visited cities, one set for each ant
        vi ants_lengths(m, 0);      // Length of each ant's path

        // Put each ant in a random initial city
        for (int k = 0; k < m; k++){        // For each ant
            int start_city = rand() % n;
            ants[k][0] = start_city;
            visited[k][start_city] = true;     // Remove initial city from unvisited
        }

        for (int j = 0; j < n - 1; j++){     // For each city
            for (int k = 0; k < m; k++){        // For each ant
                vector<bool> ant_visited = visited[k];
                int current_city = ants[k][j];
                vi candidate_list = candidate_lists[current_city];

                // Choose next city
                // 95% of the time, choose the best edge according to pheromone (exploitation)
                // 5% of the time, choose a random edge based on probability (exploration)
                next_city = -1;
                
                float thereAreCandidates = false;
                for (int l = 0; l < candidate_list_size; l++){      // Do not use the closest city
                    if (!ant_visited[candidate_list[l]]){
                        thereAreCandidates = true;
                        break;
                    }
                }


                if (thereAreCandidates){
                    next_city = choose_with_list(n, problem->dist_matrix, pheromone_matrix, current_city, candidate_list, ant_visited, candidate_list_size, Q0, beta);          
                } else {
                    next_city = choose_with_uset(n, problem->dist_matrix, pheromone_matrix, current_city, ant_visited, Q0, beta);    
                }
                
                // Update ant's path
                ants[k][j + 1] = next_city;
                // Remove next city from unvisited
                visited[k][next_city] = true;
                ants_lengths[k] += problem->dist_matrix[current_city][next_city];
                /*
                for(int i = 0; i < j + 1; i++) {
                    if (ants[k][i] == ants[k][j + 1]) {
                        cout << "ERROR: Ant " << k << " has visited city " << ants[k][i] << " twice" << endl;
                        cout << "j" << j << endl;
                        cout << "Ant's path: ";
                        for (int i = 0; i < n; i++) {
                            cout << ants[k][i] << " ";
                        }
                        cout << thereAreCandidates << endl;
                        exit(1);
                    }
                }*/
                
                // Update pheromone matrix
                tau = (1.0f - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
                
            }
        }


 //       if ((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180) return best_ant_sol;
        
        int best_iter_len = INF;
        vi best_iter_sol = vi(n, 0);
        for (int k = 0; k < m; k++){
            // Update pheromone from last city to first city
            tau = (1.0f - rho) * pheromone_matrix[ants[k][n-1]][ants[k][0]] + rho * tau0;
            pheromone_matrix[ants[k][n-1]][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][ants[k][n-1]] = tau;
            ants_lengths[k] += problem->dist_matrix[ants[k][n-1]][ants[k][0]];

            if (ants_lengths[k] < best_iter_len){
                best_iter_len = ants_lengths[k];
                best_iter_sol = ants[k];
                /*
                for (int i = 0; i < n; i++){
                    cout << best_ant_sol[i] << " ";
                }
                */
                //cout << "Best ant length 1: " << best_ant_length << endl;
            }
            
        }

//        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180) return best_ant_sol;
        vi improvement(n, 0);
        int improvement_length = 0;
        if (numIterations % 4 == 1 || numIterations % 4 == 3){
            //cout << "old length: " << ants_lengths[random_ant] << endl;
            //cout << "computed length " << compute_length(problem, ants[random_ant]) << endl;
            pair<vi, int> new_ant_sol_and_length = two_five_opt_2(problem, best_iter_sol);
            improvement = new_ant_sol_and_length.first;
            improvement_length = best_iter_len + new_ant_sol_and_length.second;
            //cout << "new length: " << ants_lengths[random_ant] << endl;
            //cout << "computed length " << compute_length(problem, ants[random_ant]) << endl;

            if (improvement_length > best_iter_len){
                improvement_length = best_iter_len;
                improvement = best_iter_sol;
                //cout << "Best ant length: " << best_ant_length << endl;
            }
        }

        pair<vi, int> best_pair;
        if (numIterations % 4 == 0){
            best_pair = two_five_opt_2(problem, best_iter_sol);
            improvement = best_pair.first;
            improvement_length = best_iter_len + best_pair.second;

        } else if (numIterations % 4 == 2){
            best_pair = two_opt_greedy_2(problem, best_ant_sol);
            improvement = best_pair.first;
            improvement_length = best_ant_length + best_pair.second;
        }

        if (improvement_length < best_ant_length){
            best_ant_length = improvement_length;
            best_ant_sol = improvement;
        }

        // Update pheromone matrix globally
        // Evaporate for every edge
        /*
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                pheromone_matrix[i][j] *= (1 - alpha);
            }
        }

        // Update pheromone matrix with best ant's path
        float constant = alpha / best_ant_length;
        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            pheromone_matrix[node_1][node_2] += constant;
            pheromone_matrix[node_2][node_1] += constant;
        }*/
        
               float constant = alpha / best_ant_length;

        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            tau = (1 - alpha) * pheromone_matrix[node_1][node_2] + constant;
            pheromone_matrix[node_1][node_2] = tau;
            pheromone_matrix[node_2][node_1] = tau;
        }
        tau = (1 - alpha) * pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] + constant;
        pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] = tau;
        pheromone_matrix[best_ant_sol[0]][best_ant_sol[n - 1]] = tau;
        numIterations++;
        cout << "Best ant length: " << best_ant_length << endl;
        
        //candidate_lists = ESACO_update_candidate_lists(candidate_lists, pheromone_matrix, n, best_ant_sol);
    }
    return best_ant_sol;
}




#endif /* ACO_HPP */
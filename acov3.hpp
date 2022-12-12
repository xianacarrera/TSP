
#ifndef ACO_HPP
#define ACO_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper_v3.hpp"
#include "local_search_v3.hpp"
#include "nn_v3.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Data structures 
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vector<float>> vvf;
typedef unordered_set<int> uset;
typedef map<int, bool> mapb;
typedef pair<int, int> pi;
#define F first
#define S second
#define PB push_back
#define MP make_pair
#define all(a) (a).begin(), (a).end() // Aplicar a toda la estructura, e.g. sort(all(a))


void determine_candidate_list(vi &candidate_list, int n, vvf &dist_matrix, uset &unvisited, int candidate_list_size, int current_city){
    uset unvisited_copy = unvisited;
    for (int i = 0; i < candidate_list_size; i++){
        float min_dist = INF;
        float min_dist_city = -1;
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
}

void determine_initial_candidate_lists(vvi &candidate_lists, int n, vvf &dist_matrix, int candidate_list_size){
    uset unvisited;
    for (int i = 0; i < n; i++){
        unvisited.insert(i);
    }
    for (int i = 0; i < n; i++){
        unvisited.erase(i);
        determine_candidate_list(candidate_lists[i], n, dist_matrix, unvisited, candidate_list_size, i);
        unvisited.insert(i);
    }
}


int choose_with_list(int n, vvf &dist_matrix, vvf &pheromone_matrix, int current_city, vi &candidate_list, vector<bool> &visited, int cl_size, int Q0, float beta, int ignore_best){
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
        for (int i = ignore_best; i < cl_size; i++){
            int city = candidate_list[i];
            if (current_city == city) continue;
            if (visited[city]) continue;
            denominator += pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
        }

        // Get a random number between 0 and 1
        float random = (float)rand() / (float) (RAND_MAX);
        float cum_prob = 0.0f;
        for (int i = ignore_best; i < cl_size; i++){
            int city = candidate_list[i];
            if (current_city == city) continue;
            if (visited[city]) continue;
            cum_prob += pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta) / denominator;
            if (random <= cum_prob){
                next_city = city;
                break;
            }
        }
        if (next_city == -1){
            for (int i = 0; i < n; i++){
                if (current_city == i) continue;
                if (visited[i]) continue;
                next_city = i;
                break;
            }
        }
    }
    return next_city;
}



int choose_with_uset(int n, vvf &dist_matrix, vvf &pheromone_matrix, int current_city, vector<bool> &visited, int Q0, float beta, int ignore_best){
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
        float denominator = 0.0f;
        int min_dist = INF;
        int closest_city = -1;
        for (int i = 0; i < n; i++){
            if (current_city == i) continue;
            if (visited[i]) continue;
            int dist = dist_matrix[current_city][i];
            if (dist < min_dist){
                min_dist = dist;
                closest_city = i;
            }
            denominator += pheromone_matrix[current_city][i] * pow(1.0 / dist, beta);
        }
        
        // Get a random number between 0 and 1
        float random = (float)rand() / (float) (RAND_MAX);
        float cum_prob = 0.0f;
        for (int i = 0; i < n; i++){
            if (current_city == i) continue;
            if (visited[i]) continue;
            cum_prob += pheromone_matrix[current_city][i] * pow(1.0 / dist_matrix[current_city][i], beta) / denominator;
            if (ignore_best && i == closest_city) continue;
            if (random <= cum_prob){
                next_city = i;
                break;
            }
        }
        if (next_city == -1){
            for (int i = 0; i < n; i++){
                if (current_city == i) continue;
                if (visited[i]) continue;
                next_city = i;
                break;
            }
        }

        if (next_city == -1){
            cout << "ERROR: next_city is -1" << endl;
        }
    }
    return next_city;
}

void ESACO_update_candidate_lists(vvi &candidate_lists, vvf &pheromone_matrix, int n, vi &best_ant_sol){
    // Update candidate lists
    for (int i = 0; i < n; i++){

        float best_pheromone = -1;
        int best_neighbour = -1;
        for (int j = 0; j < n; j++){
            if (i == j) continue;
            if (pheromone_matrix[i][j] > best_pheromone){
                best_pheromone = pheromone_matrix[i][j];
                best_neighbour = j;
            }
        }

        // Remove last node of candidates_list
        candidate_lists[i].pop_back();
        // Insert best neighbour at the beginning of the list
        candidate_lists[i].insert(candidate_lists[i].begin(), best_neighbour);
        // Remove last node of candidates_list
        candidate_lists[i].pop_back();
        // Insert successor of i in best tour as first in the candidates_list
        auto idx = find(best_ant_sol.begin(), best_ant_sol.end(), i) + 1 - best_ant_sol.begin();
        if (idx == n) idx = 0;
        int successor = best_ant_sol[idx];
        candidate_lists[i].insert(candidate_lists[i].begin(), successor);
    }
}


vi ESACO_v3(int n, vvf &dist_matrix, float alpha, float beta, float rho, int Q0, bool esaco, int m, int opt_tau0, int stopOnIter, int n0, int n1)
{
    float start_time = (float) clock();

    int numIterations = 0;       // End condition for the main loop
    int candidate_list_size = min(40, n / 4);
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
//int m = 10;                     // Number of ants

    vi best_ant_sol;                // Global histórico    
    float best_ant_length = INF;      // Global histórico

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    //float tau0 = 1.0f / (n * length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1000 * n / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = n /  ( (float) length_NN(n, dist_matrix, rand() % n)); // Initial pheromone value
    //float tau0 = 0.5f;
    //float tau0 = 1.0f / (n * best_ant_length);
    //float tau0 = 100.0f;
    float tau0, tau;

    if (opt_tau0 == 0){
        tau0 = n /  ( (float) length_NN(n, dist_matrix, rand() % n));
    } else if (opt_tau0 == 1){
        tau0 = 1.0f / (n * (float) length_NN(n, dist_matrix, rand() % n));
    } else {
        cout << "ERROR: Invalid opt_tau0 value\n";
    }

    // Pheromone matrix
    vvf pheromone_matrix(n, vector<float>(n, tau0));

    // m = 10 from original paper
    // Ants and their paths
    vvi ants(m, vi(n, 0));



    int next_city;              
    vvi candidate_lists = vvi(n, vi(candidate_list_size, 0));
    determine_initial_candidate_lists(candidate_lists, n, dist_matrix, candidate_list_size);
 
    while(numIterations < stopOnIter && (((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){
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
                    next_city = choose_with_list(n, dist_matrix, pheromone_matrix, current_city, candidate_list, ant_visited, candidate_list_size, Q0, beta, opt_tau0);          
                    /*
                    if (next_city == -1){
                        cout << "ERROR: next_city is -1 IN LIST\n";
                        exit(1);
                    }*/
                } else {
                    next_city = choose_with_uset(n, dist_matrix, pheromone_matrix, current_city, ant_visited, Q0, beta, opt_tau0);    
                    /*
                    if (next_city == -1){
                        cout << "ERROR: next_city is -1 IN SET\n";
                        exit(1);
                    }*/
                }
                
                // Update ant's path
                ants[k][j + 1] = next_city;
                // Remove next city from unvisited
                visited[k][next_city] = true;
                ants_lengths[k] += dist_matrix[current_city][next_city];
                
                // Update pheromone matrix
                tau = (1.0f - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
                
            }
        }


        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;
        
        int best_iter_len = INF;
        vi best_iter_sol;
        for (int k = 0; k < m; k++){
            // Update pheromone from last city to first city
            tau = (1.0f - rho) * pheromone_matrix[ants[k][n-1]][ants[k][0]] + rho * tau0;
            pheromone_matrix[ants[k][n-1]][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][ants[k][n-1]] = tau;
            ants_lengths[k] += dist_matrix[ants[k][n-1]][ants[k][0]];
/*
            if (!verify_sol(n, ants[k])) {
                cout << "NOT A VALID SOLUTION IN ANT " << "\n";
                exit(1);
            }
*/
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

        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;
 
        if (numIterations % 4 == 0 || numIterations % 4 == 2){
            best_iter_len += two_five_opt_2(n, dist_matrix, best_iter_sol);
            /*
            if (!verify_sol(n, best_iter_sol)) {
                cout << "NOT A VALID SOLUTION IN BEST ITER SOL 2.5 " << "\n";
                exit(1);
            }
            */
            if (best_iter_len < best_ant_length){
                best_ant_length = best_iter_len;
                best_ant_sol = best_iter_sol;
            }
        }

if(numIterations > n0){
        if (numIterations % 4 == 1 && numIterations >= n1){
            best_ant_length += two_five_opt_2(n, dist_matrix, best_ant_sol);
            /*
            if (!verify_sol(n, best_ant_sol)) {
                cout << "NOT A VALID SOLUTION IN BEST ANT SOL 2.5 " << "\n";
                exit(1);
            }
            */

        } else if (numIterations % 4 == 3 || (numIterations < n1 && numIterations % 4 == 1)){
            best_ant_length += two_opt_greedy_2(n, dist_matrix, best_ant_sol);
            /*
            if (!verify_sol(n, best_ant_sol)) {
                cout << "NOT A VALID SOLUTION IN BEST ANT LENGTH GREEDY " << "\n";
                exit(1);
            }
            */
        }
}

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
       // cout << "Best ant length: " << best_ant_length << endl;

       if (esaco)  ESACO_update_candidate_lists(candidate_lists, pheromone_matrix, n, best_ant_sol);
       //if (numIterations % 150 == 0) cout << "Iteration: " << numIterations << " Best ant length: " << best_ant_length << endl;
        
                numIterations++;

    }
    return best_ant_sol;
}


vi ESACO_for_RAT(int n, vvf &dist_matrix, float alpha, float beta, float rho, int Q0, bool esaco, int m, int opt_tau0, int stopOnIter, int n0, int n1)
{
    float start_time = (float) clock();

    int numIterations = 0;       // End condition for the main loop
    int candidate_list_size = min(40, n / 4);
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
//int m = 10;                     // Number of ants


    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    //float tau0 = 1.0f / (n * length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 1000 * n / ((float) length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = n /  ( (float) length_NN(n, dist_matrix, rand() % n)); // Initial pheromone value
    //float tau0 = 0.5f;
    //float tau0 = 1.0f / (n * best_ant_length);
    //float tau0 = 100.0f;
    float tau0, tau;

    vi best_ant_sol = best_NN(n, dist_matrix);                // Global histórico    

    int best_ant_length = compute_length(n, dist_matrix, best_ant_sol);      // Global histórico


    float improvement, improvement2;
    while(improvement || improvement2){
        improvement = two_five_opt_2(n, dist_matrix, best_ant_sol);
    best_ant_length += improvement;
        improvement2 = two_opt_greedy_2(n, dist_matrix, best_ant_sol);
    best_ant_length += improvement2;

    }


    tau0 = n /  ( (float) length_NN(n, dist_matrix, rand() % n));

    // Pheromone matrix
    vvf pheromone_matrix(n, vector<float>(n, tau0));

    // m = 10 from original paper
    // Ants and their paths
    vvi ants(m, vi(n, 0));



    int next_city;              
    vvi candidate_lists = vvi(n, vi(candidate_list_size, 0));
    determine_initial_candidate_lists(candidate_lists, n, dist_matrix, candidate_list_size);
 
    while(numIterations < stopOnIter && (((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){
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
                    next_city = choose_with_list(n, dist_matrix, pheromone_matrix, current_city, candidate_list, ant_visited, candidate_list_size, Q0, beta, !esaco);          
                    /*
                    if (next_city == -1){
                        cout << "ERROR: next_city is -1 IN LIST\n";
                        exit(1);
                    }*/
                } else {
                    next_city = choose_with_uset(n, dist_matrix, pheromone_matrix, current_city, ant_visited, Q0, beta, !esaco);    
                    /*
                    if (next_city == -1){
                        cout << "ERROR: next_city is -1 IN SET\n";
                        exit(1);
                    }*/
                }
                
                // Update ant's path
                ants[k][j + 1] = next_city;
                // Remove next city from unvisited
                visited[k][next_city] = true;
                ants_lengths[k] += dist_matrix[current_city][next_city];
                
                // Update pheromone matrix
                tau = (1.0f - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
                
            }
        }


        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;
        
        int best_iter_len = INF;
        vi best_iter_sol;
        for (int k = 0; k < m; k++){
            // Update pheromone from last city to first city
            tau = (1.0f - rho) * pheromone_matrix[ants[k][n-1]][ants[k][0]] + rho * tau0;
            pheromone_matrix[ants[k][n-1]][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][ants[k][n-1]] = tau;
            ants_lengths[k] += dist_matrix[ants[k][n-1]][ants[k][0]];
/*
            if (!verify_sol(n, ants[k])) {
                cout << "NOT A VALID SOLUTION IN ANT " << "\n";
                exit(1);
            }
*/
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

        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;
 
        if (numIterations % 4 == 0 || numIterations % 4 == 2){
            best_iter_len += two_five_opt_2(n, dist_matrix, best_iter_sol);
            /*
            if (!verify_sol(n, best_iter_sol)) {
                cout << "NOT A VALID SOLUTION IN BEST ITER SOL 2.5 " << "\n";
                exit(1);
            }
            */
            if (best_iter_len < best_ant_length){
                best_ant_length = best_iter_len;
                best_ant_sol = best_iter_sol;
            }
        }

if(numIterations > n0){
        if (numIterations % 4 == 1 && numIterations >= n1){
            best_ant_length += two_five_opt_2(n, dist_matrix, best_ant_sol);
            /*
            if (!verify_sol(n, best_ant_sol)) {
                cout << "NOT A VALID SOLUTION IN BEST ANT SOL 2.5 " << "\n";
                exit(1);
            }
            */

        } else if (numIterations % 4 == 3 || (numIterations < n1 && numIterations % 4 == 1)){
            best_ant_length += two_opt_greedy_2(n, dist_matrix, best_ant_sol);
            /*
            if (!verify_sol(n, best_ant_sol)) {
                cout << "NOT A VALID SOLUTION IN BEST ANT LENGTH GREEDY " << "\n";
                exit(1);
            }
            */
        }
}

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
       // cout << "Best ant length: " << best_ant_length << endl;

       if (esaco)  ESACO_update_candidate_lists(candidate_lists, pheromone_matrix, n, best_ant_sol);
       //if (numIterations % 150 == 0) cout << "Iteration: " << numIterations << " Best ant length: " << best_ant_length << endl;
        
                numIterations++;

    }
    return best_ant_sol;
}

#endif /* ACO_H */
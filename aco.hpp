/*
 * Xiana Carrera Alonso
 * 17th AI Cup
 * 2022
 * 
 * Ant Colony Optimization algorithms
 */

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

// Abbreviations for data structures 
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



// Function to create candidate lists for each city (the closest cities in order)
void determine_candidate_list(vi &candidate_list, int n, vvf &dist_matrix, uset &unvisited, int candidate_list_size, int current_city){
    uset unvisited_copy = unvisited;            // Set of all cities except current city
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
        candidate_list[i] = min_dist_city;    // Add closest city to candidate list
        unvisited_copy.erase(min_dist_city);  // Remove closest city from unvisited
    }
    if (unvisited.size() != n - 1) cout << "ERROR: unvisited size is not n - 1" << endl;            // Sanity check
}

// Function to create candidate lists for all cities
void determine_initial_candidate_lists(vvi &candidate_lists, int n, vvf &dist_matrix, int candidate_list_size){
    uset unvisited;
    for (int i = 0; i < n; i++){
        unvisited.insert(i);        // Set of all cities
    }
    for (int i = 0; i < n; i++){
        unvisited.erase(i);         // Set of all cities except current city
        // Determine candidate list for current city (i)   
        determine_candidate_list(candidate_lists[i], n, dist_matrix, unvisited, candidate_list_size, i);
        unvisited.insert(i);        // Add current city back to unvisited
    }
}


// Choose the next city using its candidate list
int choose_with_list(int n, vvf &dist_matrix, vvf &pheromone_matrix, int current_city, vi &candidate_list, vector<bool> &visited, int cl_size, int Q0, float beta, int ignore_best){
    // Choose next city
    // Q0% of the time, choose the best edge according to pheromone (exploitation) -> choose the first unvisited city in the candidate list
    // (100 - Q0)% of the time, choose a random edge based on probability (exploration) -> choose a random unvisited city in the candidate list
    int next_city = -1;
    if (rand() % 100 < Q0){         // Exploitation: choose the best edge
        float best_transition = -1.0f;

        for (int i = 0; i < cl_size; i++){    
            int city = candidate_list[i];
            if (current_city == city) continue;         // Don't choose the current city
            if (visited[city]) continue;                // Don't choose a city that has already been visited

            // Calculate transition probability using pheromone and distance
            float edge_transition = pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
            if (edge_transition > best_transition){
                best_transition = edge_transition;
                next_city = city;
            }
        }
    } else {                     // Exploration: choose a random edge
        float denominator = 0.0f;

        // Only the unvisited cities will end up having probability > 0
        for (int i = ignore_best; i < cl_size; i++){              // If ignore_best is true, the first city in the candidate list will be ignored (i.e. the best city)
            int city = candidate_list[i];
            if (current_city == city) continue;
            if (visited[city]) continue;
            denominator += pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta);
        }

        // Get a random number between 0 and 1
        float random = (float) rand() / (float) (RAND_MAX);
        float cum_prob = 0.0f;              // Cumulative probability
        for (int i = ignore_best; i < cl_size; i++){            
            int city = candidate_list[i];
            if (current_city == city) continue;
            if (visited[city]) continue;

            // For the unvisited cities, calculate their weighted probability and add it to the cumulative probability
            cum_prob += pheromone_matrix[current_city][city] * pow(1.0 / dist_matrix[current_city][city], beta) / denominator;
            if (random <= cum_prob){        // If the random number is less than the cumulative probability, choose that city
                next_city = city;
                break;
            }
        }
        if (next_city == -1){               // If no city was chosen, choose the first unvisited city in the candidate list (can happen with float rounding errors)
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


// If the candidate list has no unvisited cities, choose the next city using the entire set of unvisited cities
int choose_with_uset(int n, vvf &dist_matrix, vvf &pheromone_matrix, int current_city, vector<bool> &visited, int Q0, float beta, int ignore_best){
    // The algorithm performs the same as choose_with_list, but uses a set instead of a list
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
    }
    return next_city;
}


// The ESACO algorithm updates candidate lists by removing the two last nodes and inserting in the front the best neighbour (according to pheromone) and the city that 
// follows the current one in the best solution found so far 
void ESACO_update_candidate_lists(vvi &candidate_lists, vvf &pheromone_matrix, int n, vi &best_ant_sol){
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
        // Insert best neighbour according to pheromone at the beginning of the list
        candidate_lists[i].insert(candidate_lists[i].begin(), best_neighbour);
        // Remove last node of candidates_list
        candidate_lists[i].pop_back();
        // Insert successor of i in best tour as first in the candidates_list
        auto idx = find(best_ant_sol.begin(), best_ant_sol.end(), i) + 1 - best_ant_sol.begin();
        if (idx == n) idx = 0;          // If i is the last city in the best tour, the successor is the first city
        int successor = best_ant_sol[idx];
        candidate_lists[i].insert(candidate_lists[i].begin(), successor);
    }
}


vi ESACO(int n, vvf &dist_matrix, float alpha, float beta, float rho, int Q0, bool esaco, int m, int opt_tau0, int stopOnIter, int n0, int n1)
{
    float start_time = (float) clock();

    int numIterations = 0;                                           // Iteration counter 
    int candidate_list_size = min(40, n / 4);                        // Size of the candidate list for each node
    float tau0;                                                      // Initial pheromone value                          
    float tau;                                                       // Pheromone value for an edge
    vi best_ant_sol;                                                 // Best solution encountered  
    float best_ant_length = INF;                                     // Length of the best solution encountered
    vvi ants(m, vi(n, 0));                                           // Ants and their paths
    int next_city;                                                   // Next city to visit
    vvi candidate_lists = vvi(n, vi(candidate_list_size, 0));        // Candidate lists for each node


    // Initialize pheromone matrix according to the chosen option
    if (opt_tau0 == 0){
        tau0 = n /  ( (float) length_NN(n, dist_matrix, rand() % n));
    } else if (opt_tau0 == 1){
        tau0 = 1.0f / (n * (float) length_NN(n, dist_matrix, rand() % n));
    } else {
        cout << "ERROR: Invalid opt_tau0 value\n";
    }

    // Pheromone matrix
    vvf pheromone_matrix(n, vector<float>(n, tau0));

    // Initialize candidate lists with the closest cities to each node               
    determine_initial_candidate_lists(candidate_lists, n, dist_matrix, candidate_list_size);
 
    // The main loop can either be stopped after a certain number of iterations or after 3 minutes have passed
    while(numIterations < stopOnIter && (((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){
        vector<vector<bool>> visited = vector<vector<bool>>(m, vector<bool>(n, false));    // List of unvisited cities for each ant
        vi ants_lengths(m, 0);      // Length of each ant's path

        // Put each ant in a random initial city
        for (int k = 0; k < m; k++){        
            int start_city = rand() % n;
            ants[k][0] = start_city;
            visited[k][start_city] = true;     // Remove initial city from unvisited
        }

        for (int j = 0; j < n - 1; j++){        // For each city
            for (int k = 0; k < m; k++){        // For each ant
                vector<bool> ant_visited = visited[k];              // Unvisited cities for ant k (stored locally for faster access)
                int current_city = ants[k][j];                      // Current city of ant k (stored locally for faster access)
                vi candidate_list = candidate_lists[current_city];  // Candidate list for ant k (stored locally for faster access)

                // Choose next city
                next_city = -1;
                
                float thereAreCandidates = false;
                for (int l = 0; l < candidate_list_size; l++){       // Check if there are any unvisited cities in the candidate list
                    if (!ant_visited[candidate_list[l]]){
                        thereAreCandidates = true;                   // If there are, choose with list
                        break;
                    }
                }


                if (thereAreCandidates){
                    next_city = choose_with_list(n, dist_matrix, pheromone_matrix, current_city, candidate_list, ant_visited, candidate_list_size, Q0, beta, opt_tau0);          
                } else {                                // If there are no unvisited cities in the candidate list, choose with the total set of unvisited cities
                    next_city = choose_with_uset(n, dist_matrix, pheromone_matrix, current_city, ant_visited, Q0, beta, opt_tau0);    
                }
                
              
                ants[k][j + 1] = next_city;                 // Add next city to ant's path
                visited[k][next_city] = true;               // Remove next city from unvisited
                ants_lengths[k] += dist_matrix[current_city][next_city];          // Update ant's path length
                
                // Apply local pheromone update
                tau = (1.0f - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
                
            }
        }

        // Check again that the time limit has not been exceeded
        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;
        
        vi best_iter_sol;               // Best solution of this iteration
        int best_iter_len = INF;        // Length of best solution of this iteration

        // Update the last edge locally and also search for the best ant in the current iteration
        for (int k = 0; k < m; k++){
            // Update pheromone from last city to first city
            tau = (1.0f - rho) * pheromone_matrix[ants[k][n-1]][ants[k][0]] + rho * tau0;
            pheromone_matrix[ants[k][n-1]][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][ants[k][n-1]] = tau;
            ants_lengths[k] += dist_matrix[ants[k][n-1]][ants[k][0]];

            // Update best ant for the iteration
            if (ants_lengths[k] < best_iter_len){
                best_iter_len = ants_lengths[k];
                best_iter_sol = ants[k];
            }
            
        }

        // Check again that the time limit has not been exceeded
        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;
 
        if (numIterations % 4 == 0 || numIterations % 4 == 2){          // 50% of the time, update the iteration's best ant with 2.5-opt
            best_iter_len += two_five_opt_2(n, dist_matrix, best_iter_sol);
            if (best_iter_len < best_ant_length){
                best_ant_length = best_iter_len;
                best_ant_sol = best_iter_sol;
            }
        }

        if(numIterations > n0){             // Do not start improving the globally best ant solution until the first n0 iterations
            if (numIterations % 4 == 1 && numIterations >= n1){            // Do not start using 2.5-opt until the first n1 iterations
                best_ant_length += two_five_opt_2(n, dist_matrix, best_ant_sol);       // 2.5-opt to the global best
            } else if (numIterations % 4 == 3 || (numIterations < n1 && numIterations % 4 == 1)){
                best_ant_length += two_opt_greedy_2(n, dist_matrix, best_ant_sol);      // 2-opt greedy to the global best
            }
        }

        // Check again that the time limit has not been exceeded
        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;

        // Update pheromone matrix globally

        // Evaporate for every edge
        
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                pheromone_matrix[i][j] *= (1 - alpha);
            }
        }

        // Update pheromone matrix to the edges in the best ant's path
        float constant = alpha / best_ant_length;
        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            pheromone_matrix[node_1][node_2] += constant;
            pheromone_matrix[node_2][node_1] += constant;
        }
        
        // Update the last edge
        tau = (1 - alpha) * pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] + constant;
        pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] = tau;
        pheromone_matrix[best_ant_sol[0]][best_ant_sol[n - 1]] = tau;

        // If ESACO is used, update the candidate lists according to the pheromone matrix and the best ant solution
        if (esaco)  ESACO_update_candidate_lists(candidate_lists, pheromone_matrix, n, best_ant_sol);
        numIterations++;
    }
    return best_ant_sol;
}


vi ESACO_with_NN(int n, vvf &dist_matrix, float alpha, float beta, float rho, int Q0, bool esaco, int m, int opt_tau0, int stopOnIter, int n0, int n1)
{

    /* This version of the algorithm is intented for problems with a large number of nodes, where the ACO does not do a sufficiently large number of
       iterations to find a good solution. As it takes a long time to start the convergence, the final solution sometimes showed a big gap.

       To solve this problem, instead of starting with a random solution, we start with the best NN solution improved with 2.5-opt and 2-opt until
       no improvement is possible. This actually takes very little time, and the solution is already very good, so the ants are able to follow a 
       good path from the beginning. This way, the algorithm converges faster and the final solution is better.

       To account for local minima, all parameters that concern exploration are increased. 
    */

    float start_time = (float) clock();

    int candidate_list_size = n / 4;                                 // Size of the candidate list for each node
    float tau0;                                                      // Initial pheromone value
    float tau;                                                       // Pheromone value for an edge
    int numIterations = 0;                                           // Number of iterations
    float improvement, improvement2;                                 // Improvement in the length of the solution for 2.5-opt and 2-opt
    vvi ants(m, vi(n, 0));                                           // Ants and their paths
    int next_city;                                                   // Next city to visit
    vvi candidate_lists = vvi(n, vi(candidate_list_size, 0));        // Candidate lists for each node


    // Initialize the best ant solution with the best NN solution
    vi best_ant_sol = best_NN(n, dist_matrix);
    int best_ant_length = compute_length(n, dist_matrix, best_ant_sol);     // Length of the best ant solution


    while (improvement || improvement2){          // Improve the best ant solution with 2.5-opt and 2-opt until no improvement is possible
        improvement = two_five_opt_2(n, dist_matrix, best_ant_sol);
        best_ant_length += improvement;          // Update the length of the best ant solution
        improvement2 = two_opt_greedy_2(n, dist_matrix, best_ant_sol);
        best_ant_length += improvement2;         // Update the length of the best ant solution
    }

    // Initialize the pheromone matrix with a higher value than usual
    tau0 = n /  ( (float) length_NN(n, dist_matrix, rand() % n));
    // Search the closest nodes for each city and add them to the candidate list
    determine_initial_candidate_lists(candidate_lists, n, dist_matrix, candidate_list_size);

    vvf pheromone_matrix(n, vector<float>(n, tau0));                 // Pheromone matrix 

    while((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){             // Run for 3 minutes
        vector<vector<bool>> visited = vector<vector<bool>>(m, vector<bool>(n, false));         // Unordered set of visited cities, one set for each ant
        vi ants_lengths(m, 0);       // Length of each ant's path

        // Put each ant in a random initial city
        for (int k = 0; k < m; k++){       
            int start_city = rand() % n;
            ants[k][0] = start_city;
            visited[k][start_city] = true;     // Remove initial city from unvisited
        }

        for (int j = 0; j < n - 1; j++){        // For each city
            for (int k = 0; k < m; k++){        // For each ant
                vector<bool> ant_visited = visited[k];      // Unvisited cities for the current ant (stored locally for faster access)
                int current_city = ants[k][j];              // Current city of the ant
                vi candidate_list = candidate_lists[current_city];     // Candidate list for the current city (stored locally for faster access)


                next_city = -1;                
                float thereAreCandidates = false;
                for (int l = 0; l < candidate_list_size; l++){          // Check if there are any unvisited cities in the candidate list
                    if (!ant_visited[candidate_list[l]]){
                        thereAreCandidates = true;                      // If there are, we'll choose a city from the candidate list
                        break;
                    }
                }

                if (thereAreCandidates){
                    next_city = choose_with_list(n, dist_matrix, pheromone_matrix, current_city, candidate_list, ant_visited, candidate_list_size, Q0, beta, !esaco);          
                } else {        // If there are no unvisited cities in the candidate list, we'll choose a city from the total set of unvisited cities
                    next_city = choose_with_uset(n, dist_matrix, pheromone_matrix, current_city, ant_visited, Q0, beta, !esaco);    
                }
                
                ants[k][j + 1] = next_city;                   // Update the path of the ant
                visited[k][next_city] = true;                 // Remove the city from the unvisited set
                ants_lengths[k] += dist_matrix[current_city][next_city];        // Update the length of the ant's path
                
                // Update pheromone matrix locally (i.e., make the edge less attractive for the other ants)
                tau = (1.0f - rho) * pheromone_matrix[current_city][next_city] + rho * tau0;
                pheromone_matrix[current_city][next_city] = tau;
                pheromone_matrix[next_city][current_city] = tau;
                
            }
        }

        // Check again that the time limit has not been exceeded
        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;
        
        vi best_iter_sol;                   // Best ant solution in the current iteration
        int best_iter_len = INF;            // Length of the best ant solution in the current iteration

        // Update the last edge locally and also search for the best ant in the current iteration
        for (int k = 0; k < m; k++){        
            tau = (1.0f - rho) * pheromone_matrix[ants[k][n-1]][ants[k][0]] + rho * tau0;
            pheromone_matrix[ants[k][n-1]][ants[k][0]] = tau;
            pheromone_matrix[ants[k][0]][ants[k][n-1]] = tau;
            ants_lengths[k] += dist_matrix[ants[k][n-1]][ants[k][0]];

            // Update the best ant solution in the current iteration
            if (ants_lengths[k] < best_iter_len){
                best_iter_len = ants_lengths[k];
                best_iter_sol = ants[k];
            }
            
        }

        // Check again that the time limit has not been exceeded
        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;
 
        if (numIterations % 4 == 0 || numIterations % 4 == 2){      // Apply 2.5-opt to the best ant solution in the current iteration 50% of the time
            best_iter_len += two_five_opt_2(n, dist_matrix, best_iter_sol);
            if (best_iter_len < best_ant_length){
                best_ant_length = best_iter_len;
                best_ant_sol = best_iter_sol;
            }
        }

        if(numIterations > n0){         // Do not start improving the globally best ant solution until the first n0 iterations
            if (numIterations % 4 == 1 && numIterations >= n1){     // Apply 2.5-opt to the global best
                best_ant_length += two_five_opt_2(n, dist_matrix, best_ant_sol);
            } else if (numIterations % 4 == 3 || (numIterations < n1 && numIterations % 4 == 1)){
                best_ant_length += two_opt_greedy_2(n, dist_matrix, best_ant_sol);      // Apply 2-opt to the global best
            }
        }

        // Check again that the time limit has not been exceeded
        if ((((float) clock()) - start_time) / CLOCKS_PER_SEC > 180) return best_ant_sol;

        // Update pheromone matrix globally

        // Evaporate for every edge
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                pheromone_matrix[i][j] *= (1 - alpha);
            }
        }

        // Update pheromone matrix to the edges in the best ant's path
        float constant = alpha / best_ant_length;
        for (int j = 0; j < n - 1; j++){
            int node_1 = best_ant_sol[j];
            int node_2 = best_ant_sol[j + 1];
            pheromone_matrix[node_1][node_2] += constant;
            pheromone_matrix[node_2][node_1] += constant;
        }

        // Update the last edge
        tau = (1 - alpha) * pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] + constant;
        pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] = tau;
        pheromone_matrix[best_ant_sol[0]][best_ant_sol[n - 1]] = tau;

        // If ESACO is used, update the candidate lists according to the pheromone matrix and the best ant solution
        if (esaco)  ESACO_update_candidate_lists(candidate_lists, pheromone_matrix, n, best_ant_sol);        
        numIterations++;
    }
    return best_ant_sol;
}

#endif /* ACO_H */
/*
 * Xiana Carrera Alonso
 * 17th AI Cup
 * 2022
 * 
 * 2.5-opt and 2-opt local search
 */


#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Abbreviations for data structures
typedef vector<int> vi;
typedef vector<vector<float>> vvf;


// Applies 2-opt local search to the solution greedily and returns the gain
float two_opt_greedy_2(int n, vvf &dist_matrix, vi &solution){
    int best_i;                 // Best i for the 2-opt move
    int best_j;                 // Best j for the 2-opt move
    float best_gain;            // Best gain encountered in the iteration
    float global_gain = 0;      // Global gain of the 2-opt move
    int niter = 0;              // Number of iterations

    do{
        best_gain = 0;        // Reset the best gain
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = solution[i - 1];
                int b = solution[i];
                int c = solution[j + 1];
                int d = solution[j];
                // Compute the difference between the old and the new path
                float gain = dist_matrix[a][d] + dist_matrix[b][c] - dist_matrix[a][b] - dist_matrix[c][d];
                if (gain < best_gain){      // If the current length will be shorter, update the best gain
                    best_gain = gain;
                    best_i = i;
                    best_j = j; 
                }
            }
        }   
        if (best_gain < 0){         // If the best gain is negative (there's a gain), apply the move
            reverse(solution.begin() + best_i, solution.begin() + best_j + 1);      // Reverse the subpath
            global_gain += best_gain;
            break;
        }
        niter++;
    } while (best_gain && niter < 100);             // Don't do more than 100 iterations to avoid taking too much time
    return global_gain;
}


// Applies 2-opt local search stopping as soon as a gain is found
float two_opt_simplified_2(int n, vvf &dist_matrix, vi &solution){
    bool improved;
    int best_i;
    int best_j;
    float gain = 0;
    int n_iter = 0;

    do{
        improved = false;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = solution[i - 1];
                int b = solution[i];
                int c = solution[j + 1];
                int d = solution[j];
                if (dist_matrix[a][d] + dist_matrix[b][c] < dist_matrix[a][b] - dist_matrix[c][d]){
                    improved = true;        // If the current length will be shorter, stop
                    gain += dist_matrix[a][d] + dist_matrix[b][c] - dist_matrix[a][b] - dist_matrix[c][d];
                    best_i = i;
                    best_j = j; 
                    break;
                }
            }
            if (improved){      // Apply the exchange and stop
                reverse(solution.begin() + best_i, solution.begin() + best_j + 1);
                break;
            }
        }   
        n_iter++;
    } while (improved && n_iter < 100);         // Don't do more than 100 iterations to avoid taking too much time
    return gain;
}

// Applies 2.5-opt local search to the solution and returns the global gain
float two_five_opt_2(int n, vvf &dist_matrix, vi &solution){
    int best_i;                 // Best i for the 2.5-opt move
    int best_j;                 // Best j for the 2.5-opt move
    float best_gain = 0;        // Best gain encountered in the iteration
    float global_gain = 0;      // Total gain
    int n_iter = 0;             // Number of iterations
    short opt;                  // Option for the 2.5-opt move

    do {
        best_gain = 0;
        opt = 0;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = solution[i - 1];
                int b = solution[i];
                int c = solution[i + 1];
                int d = solution[j - 1];
                int e = solution[j];
                int f = solution[j + 1];

                // 2-opt
                float gain1 = dist_matrix[a][e] + dist_matrix[b][f] - dist_matrix[a][b] - dist_matrix[e][f];

                // node shift 1
                float gain2 = dist_matrix[a][c] + dist_matrix[b][f] + dist_matrix[b][e] - dist_matrix[a][b] - dist_matrix[b][c] - dist_matrix[e][f];
            
                // node shift 2
                float gain3 = dist_matrix[a][e] + dist_matrix[b][e] + dist_matrix[d][f] - dist_matrix[a][b] - dist_matrix[d][e] - dist_matrix[e][f];

                // Only consider the move if it's better than the best gain in the iteration
                if (gain1 >= best_gain && gain2 >= best_gain && gain3 >= best_gain) continue;

                best_i = i;
                best_j = j;

                if (gain1 < gain2 && gain1 < gain3){            // gain1 is the best
                    best_gain = gain1;
                    opt = 1;
                } else if (gain2 < gain1 && gain2 < gain3){     // gain2 is the best
                    best_gain = gain2;
                    opt = 2;
                } else if (gain3 < gain1 && gain3 < gain2) {    // gain3 is the best
                    best_gain = gain3;
                    opt = 3;
                } else {
                    opt = 0;        // Tie
                }
            }
        }

        // Apply the best gain
        vi copy_sol, copy_sol2;
        switch(opt){
            case 1:         // 2-opt
                global_gain += best_gain;       // Update the global gain
                reverse(solution.begin() + best_i, solution.begin() + best_j + 1);      
                break;
            case 2:     // Node shift 1
                global_gain += best_gain;
                copy_sol = solution;        // Copy the solution

                // Move node best_i to best_j
                for (int i = best_i; i < best_j; i++){
                    solution[i] = copy_sol[i + 1];
                }
                solution[best_j] = copy_sol[best_i];        
    
                break;
            case 3:     // Node shift 2
                global_gain += best_gain;
                copy_sol2 = solution;
                
                // Move node best_j to best_i
                solution[best_i] = copy_sol2[best_j];
                for (int i = best_i + 1; i <= best_j; i++){
                    solution[i] = copy_sol2[i - 1];
                }
           
                break;
            default:
                break;
        }
        n_iter++;
    } while (opt && n_iter < 100);          // Don't do more than 100 iterations to avoid taking too much time
    return global_gain;
}


#endif // LOCAL_SEARCH_HPP
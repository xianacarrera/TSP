
#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper_v3.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Data structures 
typedef vector<int> vi;
typedef vector<vector<float>> vvf;



float two_opt_greedy_2(int n, vvf &dist_matrix, vi &solution){
    int best_i;
    int best_j;
    float best_gain;
    float global_gain = 0;

    int niter = 0;

    do{
        best_gain = 0;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = solution[i - 1];
                int b = solution[i];
                int c = solution[j + 1];
                int d = solution[j];
                float gain = dist_matrix[a][d] + dist_matrix[b][c] - dist_matrix[a][b] - dist_matrix[c][d];
                if (gain < best_gain){
                    best_gain = gain;
                    best_i = i;
                    best_j = j; 
                }
            }
        }   
        if (best_gain < 0){
            // Exchange the edges from the best_gain
            reverse(solution.begin() + best_i, solution.begin() + best_j + 1);
            global_gain += best_gain;
            break;
        }
        niter++;
    } while (best_gain && niter < 100);
    return global_gain;
}


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
                    improved = true;
                    gain += dist_matrix[a][d] + dist_matrix[b][c] - dist_matrix[a][b] - dist_matrix[c][d];
                    best_i = i;
                    best_j = j; 
                    break;
                }
            }
            if (improved){
                // Exchange the edges from the best_gain
                reverse(solution.begin() + best_i, solution.begin() + best_j + 1);
                break;
            }
        }   
        n_iter++;
    } while (improved && n_iter < 100);
    return gain;
}


float two_five_opt_2(int n, vvf &dist_matrix, vi &solution){
    int best_i;
    int best_j;
    float best_gain = 0;
    float global_gain = 0;
    int n_iter = 0;
    int start_length = compute_length(n, dist_matrix, solution);
    short opt;
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

                // Index of the minimum gain
                if (gain1 >= best_gain && gain2 >= best_gain && gain3 >= best_gain) continue;

                best_i = i;
                best_j = j;

                if (gain1 < gain2 && gain1 < gain3){
                    best_gain = gain1;
                    opt = 1;
                } else if (gain2 < gain1 && gain2 < gain3){
                    best_gain = gain2;
                    opt = 2;
                } else if (gain3 < gain1 && gain3 < gain2) {
                    best_gain = gain3;
                    opt = 3;
                } else {
                    opt = 0;
                }
            }
        }

        // Apply the best gain
        vi copy_sol, copy_sol2;
        switch(opt){
            case 1:
                global_gain += best_gain;
                reverse(solution.begin() + best_i, solution.begin() + best_j + 1);
                break;
            case 2:     // oldIndex = best_i < newIndex = best_j + 1
                global_gain += best_gain;
                /* https://stackoverflow.com/questions/45447361/how-to-move-certain-elements-of-stdvector-to-a-new-index-within-the-vector */
                //rotate(solution.begin() + best_i, solution.begin() + best_i + 1, solution.begin() + best_j + 2);
                copy_sol = solution;

                for (int i = best_i; i < best_j; i++){
                    solution[i] = copy_sol[i + 1];
                }
                solution[best_j] = copy_sol[best_i];
    
                break;
            case 3:     // oldIndex = best_j > newIndex = best_i
                global_gain += best_gain;
                //rotate(best_solution.rend() - best_j - 1, best_solution.rend() - best_j, best_solution.rend() - best_i);

                copy_sol2 = solution;
                
                solution[best_i] = copy_sol2[best_j];
                for (int i = best_i + 1; i <= best_j; i++){
                    solution[i] = copy_sol2[i - 1];
                }
           
                /*
                                new_length36 = compute_length(problem, solution);

                if (new_length36 != global_gain + start_length){
                    cout << "Error 222222: " << new_length36 << " != " << global_gain + start_length << endl;
                    exit(1);
                }*/
                break;
            default:
                break;
        }
        n_iter++;
    } while (opt && n_iter < 100);
    return global_gain;
}


#endif // LOCAL_SEARCH_HPP
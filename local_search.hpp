
#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Data structures 
typedef vector<int> vi;
typedef vector<vi> vvi;


vi two_opt_greedy(Problem * problem, vi solution){
    int n = problem->n;
    vvi dist_matrix = problem->dist_matrix;

    vi best_solution = solution;
    int best_i;
    int best_j;
    int best_gain;

    int niter = 0;

    do{
        best_gain = 0;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = best_solution[i - 1];
                int b = best_solution[i];
                int c = best_solution[j + 1];
                int d = best_solution[j];
                int gain = dist_matrix[a][d] + dist_matrix[b][c] - dist_matrix[a][b] - dist_matrix[c][d];
                if (gain < best_gain){
                    best_gain = gain;
                    best_i = i;
                    best_j = j; 
                }
            }
        }   
        if (best_gain < 0){
            // Exchange the edges from the best_gain
            reverse(best_solution.begin() + best_i, best_solution.begin() + best_j + 1);
            break;
        }
        niter++;
    } while (best_gain);
    return best_solution;
}


pair<vi, int> two_opt_greedy_2(Problem * problem, vi solution){
    int n = problem->n;
    vvi dist_matrix = problem->dist_matrix;

    vi best_solution = solution;
    int best_i;
    int best_j;
    int best_gain;
    int global_gain = 0;

    int niter = 0;

    do{
        best_gain = 0;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = best_solution[i - 1];
                int b = best_solution[i];
                int c = best_solution[j + 1];
                int d = best_solution[j];
                int gain = dist_matrix[a][d] + dist_matrix[b][c] - dist_matrix[a][b] - dist_matrix[c][d];
                if (gain < best_gain){
                    best_gain = gain;
                    best_i = i;
                    best_j = j; 
                }
            }
        }   
        if (best_gain < 0){
            // Exchange the edges from the best_gain
            reverse(best_solution.begin() + best_i, best_solution.begin() + best_j + 1);
            global_gain += best_gain;
            break;
        }
        niter++;
    } while (best_gain && niter < 100);
    return make_pair(best_solution, global_gain);
}

vi two_opt_simplified(Problem * problem, vi solution){
    int n = problem->n;
    vvi dist_matrix = problem->dist_matrix;

    bool improved;
    vi best_solution = solution;
    int best_i;
    int best_j;

    do{
        improved = false;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = best_solution[i - 1];
                int b = best_solution[i];
                int c = best_solution[j + 1];
                int d = best_solution[j];
                if (dist_matrix[a][d] + dist_matrix[b][c] < dist_matrix[a][b] - dist_matrix[c][d]){
                    improved = true;
                    best_i = i;
                    best_j = j; 
                    break;
                }
            }
            if (improved){
                // Exchange the edges from the best_gain
                reverse(best_solution.begin() + best_i, best_solution.begin() + best_j + 1);
                break;
            }
        }   
    } while (improved);
    return best_solution;
}


pair<vi, int> two_opt_simplified_2(Problem * problem, vi solution){
    int n = problem->n;
    vvi dist_matrix = problem->dist_matrix;

    bool improved;
    vi best_solution = solution;
    int best_i;
    int best_j;
    int gain = 0;
    int n_iter = 0;

    do{
        improved = false;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = best_solution[i - 1];
                int b = best_solution[i];
                int c = best_solution[j + 1];
                int d = best_solution[j];
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
                reverse(best_solution.begin() + best_i, best_solution.begin() + best_j + 1);
                break;
            }
        }   
        n_iter++;
    } while (improved && n_iter < 100);
    return make_pair(best_solution, gain);
}



vi two_five_opt(Problem * problem, vi solution){
    int n = problem->n;
    vvi dist_matrix = problem->dist_matrix;

    vi best_solution = solution;

    int best_i;
    int best_j;
    int best_gain;
    int n_iter = 0;
    short opt = 0;
    do {
        best_gain = 0;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = best_solution[i - 1];
                int b = best_solution[i];
                int c = best_solution[i + 1];
                int d = best_solution[j - 1];
                int e = best_solution[j];
                int f = best_solution[j + 1];

                // 2-opt
                int gain1 = dist_matrix[a][e] + dist_matrix[b][f] - dist_matrix[a][b] - dist_matrix[e][f];

                // node shift 1
                int gain2 = dist_matrix[a][c] + dist_matrix[b][f] + dist_matrix[b][e] - dist_matrix[a][b] - dist_matrix[b][c] - dist_matrix[e][f];
            
                // node shift 2
                int gain3 = dist_matrix[a][e] + dist_matrix[b][e] + dist_matrix[d][f] - dist_matrix[a][b] - dist_matrix[d][e] - dist_matrix[e][f];

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
        switch(opt){
            case 1:
                reverse(best_solution.begin() + best_i, best_solution.begin() + best_j + 1);
                break;
            case 2:     // oldIndex = best_i < newIndex = best_j + 1
                /* https://stackoverflow.com/questions/45447361/how-to-move-certain-elements-of-stdvector-to-a-new-index-within-the-vector */
                rotate(best_solution.begin() + best_i, best_solution.begin() + best_i + 1, best_solution.begin() + best_j + 2);
                break;
            case 3:     // oldIndex = best_j > newIndex = best_i
                rotate(best_solution.rend() - best_j - 1, best_solution.rend() - best_j, best_solution.rend() - best_i);
                break;
            default:
                break;
        }
        n_iter++;
    } while (opt);
    return best_solution;
}


pair<vi, int> two_five_opt_2(Problem * problem, vi solution){
    int n = problem->n;
    vvi dist_matrix = problem->dist_matrix;

    vi best_solution = solution;
    int start_length = compute_length(problem, solution);

    int best_i;
    int best_j;
    int best_gain = 0;
    int global_gain = 0;
    int n_iter = 0;
    short opt = 0;
    do {
        best_gain = 0;
        for (int i = 1; i < n - 2; i++){
            for (int j = i + 1; j < n - 1; j++){
                int a = best_solution[i - 1];
                int b = best_solution[i];
                int c = best_solution[i + 1];
                int d = best_solution[j - 1];
                int e = best_solution[j];
                int f = best_solution[j + 1];

                // 2-opt
                int gain1 = dist_matrix[a][e] + dist_matrix[b][f] - dist_matrix[a][b] - dist_matrix[e][f];

                // node shift 1
                int gain2 = dist_matrix[a][c] + dist_matrix[b][f] + dist_matrix[b][e] - dist_matrix[a][b] - dist_matrix[b][c] - dist_matrix[e][f];
            
                // node shift 2
                int gain3 = dist_matrix[a][e] + dist_matrix[b][e] + dist_matrix[d][f] - dist_matrix[a][b] - dist_matrix[d][e] - dist_matrix[e][f];

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
        int new_length36;
        int new_length21;
        vi array_one, array_two, array_three, sol;
        switch(opt){
            case 1:
                global_gain += best_gain;
                reverse(best_solution.begin() + best_i, best_solution.begin() + best_j + 1);
                break;
            case 2:     // oldIndex = best_i < newIndex = best_j + 1
                global_gain += best_gain;
                /* https://stackoverflow.com/questions/45447361/how-to-move-certain-elements-of-stdvector-to-a-new-index-within-the-vector */
                //rotate(best_solution.begin() + best_i, best_solution.begin() + best_i + 1, best_solution.begin() + best_j + 2);
                array_one = {best_solution.begin(), best_solution.begin() + best_i};
                array_two = {best_solution.begin() + best_i + 1, best_solution.begin() + best_j + 1};
                array_three = {best_solution.begin() + best_j + 1, best_solution.end()};

                sol = array_one;
                sol.insert(sol.end(), array_two.begin(), array_two.end());
                sol.push_back(best_solution[best_i]);
                sol.insert(sol.end(), array_three.begin(), array_three.end());

                best_solution = sol;

/*
                new_length21 = compute_length(problem, best_solution);
                if (new_length21 != global_gain + start_length){
                    cout << "Error: " << new_length21 << " != " << global_gain + start_length << endl;
                    exit(1);
                }*/
                break;
            case 3:     // oldIndex = best_j > newIndex = best_i
                global_gain += best_gain;
                //rotate(best_solution.rend() - best_j - 1, best_solution.rend() - best_j, best_solution.rend() - best_i);

                array_one = {best_solution.begin(), best_solution.begin() + best_i};
                array_two = {best_solution.begin() + best_i, best_solution.begin() + best_j};
                array_three = {best_solution.begin() + best_j + 1, best_solution.end()};

                sol = array_one;
                sol.push_back(best_solution[best_j]);
                sol.insert(sol.end(), array_two.begin(), array_two.end());
                sol.insert(sol.end(), array_three.begin(), array_three.end());

                best_solution = sol;
                /*
                                new_length36 = compute_length(problem, best_solution);

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
    return make_pair(best_solution, global_gain);
}


#endif // LOCAL_SEARCH_HPP
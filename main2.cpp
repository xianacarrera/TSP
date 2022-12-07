#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;
int main()
{
    std::cout << "Current path is " << fs::current_path() << '\n'; // (1)
    fs::current_path(fs::temp_directory_path()); // (3)
    std::cout << "Current path is " << fs::current_path() << '\n';
}



/*
#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include <bits/stdc++.h>
#include "problem.hpp"

using namespace std;

// Es bueno desactivarlas para debuggear
#pragma GCC optimize("O3,unroll-loops")           // ofast?
#pragma GCC target("avx,avx2,fma")

// Tipos
#define ll long long   // 64 bits INT
#define ld long double // 80 bits FP 


// Constantes
#define PI 3.1415926535897932384626433832795l

// Estructuras de datos
#define ar array        
typedef vector<int> vi;
typedef pair<float, float> pf;
#define umap unordered_map
#define F first
#define S second
#define PB push_back
#define MP make_pair
#define sza(x) ((int)x.size())  
#define all(a) (a).begin(), (a).end() // Aplicar a toda la estructura, e.g. sort(all(a))

// Bucles
#define FOR(i, a, b) for(int i = a; i<= b; i++)

*/

vi ACO_LS_OLD(Problem * problem)
{
    float start_time = (float) clock();

    int n = problem->n;
    int numIterations = 1250;       // End condition for the main loop
    float alpha = 0.1f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.1f;               // Evaporation rate
    int Q0 = 95;                    // Probability of choosing the best solution

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    float tau0 = 1.0f / (n * length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 0.1f;
    float tau;

    // Pheromone matrix
    vector<vector<float>> pheromone_matrix(n, vector<float>(n, tau0));


    // m = 10 from original paper
    // Ants and their paths
    vector<vi> ants(10, vi(n, 0));

    int best_ant_length = INF;      // Global histórico
    vi best_ant_sol;                // Global histórico    

    int next_city;               

    while ((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){         // Run for 3 minutes
        bool new_ant_sol = false;

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
                new_ant_sol = true;
            }
        }
        // Update pheromone matrix globally

        if (new_ant_sol){
            two_five_opt(problem, best_ant_sol);
            two_opt_greedy(problem, best_ant_sol);
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
        pheromone_matrix[best_ant_sol[n - 1]][best_ant_sol[0]] += constant;
        pheromone_matrix[best_ant_sol[0]][best_ant_sol[n - 1]] += constant;
    }
    return best_ant_sol;
}

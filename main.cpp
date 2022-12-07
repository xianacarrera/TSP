#include <bits/stdc++.h>
#include "problem.hpp"
#include "local_search.hpp"
#include "helper.hpp"
#include "nn.hpp"
//#include "aco.hpp"
#include "acov2.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Data structures 
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef unordered_set<int> uset;
typedef pair<float, float> pf;
#define umap unordered_map
#define F first
#define S second
#define PB push_back
#define MP make_pair
#define all(a) (a).begin(), (a).end() // Aplicar a toda la estructura, e.g. sort(all(a))

const string problems[10] = {"ch130.tsp", "d198.tsp", "eil76.tsp", "fl1577.tsp", "kroA100.tsp", "lin318.tsp", "pcb442.tsp", "pr439.tsp", "rat783.tsp", "u1060.tsp"};
const string custom[1] = {"custom.tsp"};


vi ACO(Problem * problem, bool local_search)
{
    float start_time = (float) clock();

    int n = problem->n;
    int numIterations = 0;       // End condition for the main loop
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


    // m = 10 from original paper
    // Ants and their paths
    vector<vi> ants(m, vi(n, 0));

    int best_ant_length = INF;      // Global histórico
    vi best_ant_sol;                // Global histórico    

    int next_city;               

    while (((float) clock() - start_time) / CLOCKS_PER_SEC < 180){         // Run for 3 minutes
        // Put each ant in a random initial city
        for (int k = 0; k < m; k++){        // For each ant
            ants[k][0] = rand() % n;
        }

        // Unordered set of unvisited cities, one set for each ant
        vector<uset> unvisited(m);
        for (int k = 0; k < 10; k++){
            for (int j = 0; j < n; j++){
                unvisited[k].insert(j);
            }
            unvisited[k].erase(ants[k][0]);     // Remove initial city from unvisited
        }

        for (int j = 0; j < n - 1; j++){     // For each city
            for (int k = 0; k < m; k++){        // For each ant
                uset ant_unvisited = unvisited[k];
                int current_city = ants[k][j];

                // Choose next city
                // 95% of the time, choose the best edge according to pheromone (exploitation)
                // 5% of the time, choose a random edge based on probability (exploration)
                next_city = -1;
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
                    //int min_dist_to_unvisited = INF;
                    //int closest_city;
                    float denominator = 0.0f;
                    uset::iterator l;
                    for (l = ant_unvisited.begin(); l != ant_unvisited.end(); ++l) {
                        /*
                        int dist = problem->dist_matrix[current_city][*l];
                        if (dist < min_dist_to_unvisited){
                            min_dist_to_unvisited = dist;
                            closest_city = *l;
                        }
                        */
                        denominator += pheromone_matrix[current_city][*l] * pow(1.0 / problem->dist_matrix[current_city][*l], beta);
                    }
                    
                    uset::iterator l2;
                    for (l2 = ant_unvisited.begin(); l2 != ant_unvisited.end(); ++l2) {
                        probabilities[*l2] = pheromone_matrix[current_city][*l2] * pow(1.0 / problem->dist_matrix[current_city][*l2], beta) / denominator;
                    }
                    //probabilities[closest_city] = 0;    // Don't choose the closest city
                    
                    // Create a distribution based on the probabilities
                    discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
                    // Choose a city based on the distribution
                    next_city = distribution(rng);
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

        //cout << numIterations << "\n";
    }
/*
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << pheromone_matrix[i][j] << " ";
        }
    }
    cout << "\n";

    cout << "NUMBER OF ITERATIONS: " << numIterations << "\n";
    */
    return best_ant_sol;
}

vi ACO_LS(Problem * problem){
       float start_time = (float) clock();

    int n = problem->n;
    int numIterations = 1250;       // End condition for the main loop
    float phi = 0.1f;             
    float beta = 2.0f;              // Importance of distance over pheromone  
    float rho = 0.1f;               // Evaporation rate
    int Q0 = 90;                    // Probability of choosing the best solution
    int m = 10;                     // Number of ants

    // TODO: Initialize pheromone matrix
    //float tau0 = 1.0f / (n * length_best_NN(problem)); // Initial pheromone value
    float tau0 = n / (length_NN(problem, rand() % n)); // Initial pheromone value
    //float tau0 = 0.1f;
    float tau;

    // Pheromone matrix
    vector<vector<float>> pheromone_matrix(n, vector<float>(n, tau0));


    // m = 10 from original paper
    // Ants and their paths
    vector<vi> ants(m, vi(n, 0));

    int best_ant_length = INF;      // Global histórico
    vi best_ant_sol;                // Global histórico    

    int next_city;               

    while ((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){         // Run for 3 minutes
        // Put each ant in a random initial city
        for (int k = 0; k < m; k++){        // For each ant
            ants[k][0] = rand() % n;
        }

        // Unordered set of unvisited cities, one set for each ant
        vector<uset> unvisited(m);
        for (int k = 0; k < 10; k++){
            for (int j = 0; j < n; j++){
                unvisited[k].insert(j);
            }
            unvisited[k].erase(ants[k][0]);     // Remove initial city from unvisited
        }

        for (int j = 0; j < n - 1; j++){     // For each city
            for (int k = 0; k < m; k++){        // For each ant
                uset ant_unvisited = unvisited[k];
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
                unvisited[k].erase(next_city);

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
            
            // Calculate cost of each ant's path
            // TODO: improve updating it in the algorithm?
            int ant_length = compute_length(problem, ants[k]);
            if (ant_length < best_ant_length){
                best_ant_length = ant_length;
                best_ant_sol = ants[k];
            }
        }

        /* Event if the best sol hasn't changed, we still want to apply local search, as it can produce different results (because we're combining 2-5 and 2opt)*/
        two_five_opt(problem, best_ant_sol);
        two_opt_greedy(problem, best_ant_sol);

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
    }

    return best_ant_sol;
}


int main()
{
    
    // Print the seed
    int seed = 120;
    srand(120);
    rng.seed(seed);
    cout << "Seed: " << seed << endl;

    write_header("best NN + 2.5-opt + 2-opt", seed);
    for (int i = 0; i < 10; i++) {
        Problem * problem = new Problem(problems[i]);        
        //vi sol = NN(problem, rand()  % problem->n);
        //vi sol = best_NN(problem);
        //sol = two_five_opt(problem, sol);
        //sol = two_opt_greedy(problem, sol);
        vi sol = ACO_candidate(problem);
        
        if (!verify_sol(problem, sol)) {
            cout << "NOT A VALID SOLUTION FOR " << problem->name << "\n";
            exit(1);
        }
        compute_error(problem, sol);
        write_results(problem, sol);

        vvi distance_matrix = problem->dist_matrix;
        int starting_node = sol[0];
        int from_node = starting_node;

        cout << starting_node << " ";

        for (int i = 1; i < sol.size(); i++) {
            int to_node = sol[i];
            cout << to_node << ", ";
            from_node = to_node;
        }
        cout << "\n";

    }
    
   /*
    Problem * problem = new Problem(custom[0]);
    vi sol = ACO(problem);
    compute_error(problem, sol);*/
}
/*
 * Xiana Carrera Alonso
 * 17th AI Cup
 * 2022
 */

#include <bits/stdc++.h>
#include "problem.hpp"      // Problem class (includes reader)
#include "local_search.hpp" // 2-opt, 2.5-opt
#include "helper.hpp"       // Auxiliary functions
#include "nn.hpp"           // Nearest Neighbour algorithms
#include "aco.hpp"          // Ant Colony Optimization algorithms

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops") // Compilation options with loop unrolling
#pragma GCC target("avx,avx2,fma")      // Vector extensions

// Abbreviations for data structures
typedef vector<int> vi;
typedef vector<vector<float>> vvf;
typedef unordered_set<int> uset;
typedef pair<float, float> pf;
#define umap unordered_map
#define F first
#define S second
#define PB push_back
#define MP make_pair

// Names of the problems
const string files[10] = {
    "ch130.tsp",
    "d198.tsp",
    "eil76.tsp",
    "fl1577.tsp",
    "kroA100.tsp",
    "lin318.tsp",
    "pcb442.tsp",
    "pr439.tsp", 
    "rat783.tsp", 
    "u1060.tsp"
};


// Parameters for the ACO algorithms
// As problems are very different in nature, each one performs better with custom paramters (though there are subgroups with show similar characteristics)

// Alpha (global pheromone update)
float param0[10] = {0.6, 0.09, 0.35, 0.09, 0.15, 0.15, 0.09, 0.1, 0.15, 0.25};
// Beta
float param1[10] = {10, 2.5, 3, 1.8, 2.0, 2, 2.5, 2.6, 2.8, 2.7};
// Rho (local pheromone update)
float param2[10] = {0.4, 0.09, 0.08, 0.19, 0.09, 0.09, 0.09, 0.35, 0.35, 0.35};
// Q0 (exploration/exploitation)
int param3[10] = {97, 97, 97, 99, 98, 98, 98, 98, 98, 99};
// Use ESACO instead of ACO?
bool param4[10] = {false, false, true, true, true, true, false, true, true, true};
// m (number of ants)
int param5[10] = {10, 10, 10, 15, 10, 10, 15, 15, 20, 20};
// Option for initializing tau0 (starting pheromone)
int param6[10] = {1, 0, 0, 0, 0, 0, 1, 1, 0, 1};
// Maximum number of iterations (only used in testing, set to a high number in the final version to make use of the 3 minutes time limit)
int param7[10] = {1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000};
// Number of iterations without using 2.5-opt for the best tour
int param8[10] = {100, 100, 100, 100, 100, 150, 50, 70, 100, 50};



int main(int argc, char **argv)
{
    int seed, problem_idx;
    if (argc < 3)
    {
        /*
        cout << "Usage: ./main seed problem_idx" << endl;
        return 0;
        */
       seed = 52;
       problem_idx = 0;
    }
    else
    {
        seed = atoi(argv[1]);
        problem_idx = atoi(argv[2]);
    }

    srand(seed);

    Problem *problem = new Problem(files[problem_idx]);         // Read the problem
    vi sol;
    if (problem_idx == 8 || problem_idx == 9){          // Rat783 and U1060 are too big for a standard ESACO
        sol = ESACO_with_NN(problem->n, problem->dist_matrix, param0[problem_idx], param1[problem_idx], param2[problem_idx], param3[problem_idx], param4[problem_idx], param5[problem_idx], param6[problem_idx], param7[problem_idx], param8[problem_idx], 2 * param8[problem_idx]);
    } else {
        sol = ESACO(problem->n, problem->dist_matrix, param0[problem_idx], param1[problem_idx], param2[problem_idx], param3[problem_idx], param4[problem_idx], param5[problem_idx], param6[problem_idx], param7[problem_idx], param8[problem_idx], 2 * param8[problem_idx]);
    }

    if (!check_solution(problem->n, sol)){
        cout << "NOT A VALID SOLUTION FOR " << problem->name << "\n";
        exit(1);
    }

    float error = compute_error(problem, sol);                    // Compute the error and print it

    write_results("results.csv", seed, problem->name, error);     // Write the results to a file

    return 0;
}
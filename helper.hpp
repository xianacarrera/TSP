/*
 * Xiana Carrera Alonso
 * 17th AI Cup
 * 2022
 * 
 * Auxiliary functions
 */


#ifndef HELPER_HPP
#define HELPER_HPP

#include <bits/stdc++.h>
#include "problem.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Abbreviations for data structures 
typedef vector<int> vi;
typedef vector<vector<float>> vvf;
#define all(a) (a).begin(), (a).end()   // Apply to all elements of a vector


// Return a random list of n nodes
vi random_solution(Problem * problem){
    vi solution;
    for (int i = 0; i < problem->n; i++) solution.push_back(i);

    random_shuffle(all(solution));
    return solution;
}

// Get the total distance of a path
int compute_length(int n, vvf &distance_matrix, vi solution) {
    int length = 0;
    int starting_node = solution[0];
    int from_node = starting_node;      

    for (int i = 1; i < n; i++) {
        int to_node = solution[i];
        length += distance_matrix[from_node][to_node];
        from_node = to_node;        // Update the node
    }

    // Add the distance from the last node to the starting node
    length += distance_matrix[starting_node][from_node];

    return length;
}

// Compute the error (gap) of a solution
float compute_error(Problem * problem, vi solution){
    int length = compute_length(problem->n, problem->dist_matrix, solution);
    float error = 100.0f * (length - problem->best_known) / (float) problem->best_known;
    cout << problem->name << "\n";
    cout << "\tBest known: " << problem->best_known << "\n";
    cout << "\tSolution: " << length << "\n";
    cout << "\tError: " << error << "%\n";
    return error;
}

// Check if a solution is valid by
bool check_solution(int n, vi sol){
    for (int i = 0; i < n; i++) {
        bool cityFound = false;
        for (auto j : sol){
            if (i == j){
                cityFound = true;
            }
        }
        if (!cityFound) return false;
    }

    return true;
}


// Write the results to a file
void write_results(string filename, int seed, string name, float error){
    ofstream file;
    file.open(filename, ios::app);
    file << name << "; " << seed << "; " << error << "\n";
    file.close();
}

// Write the parameters used to a file
void write_params(string filename, float p0, float p1, float p2, int p3, bool p4, int p5, int p6, int p7, int p8){
    ofstream file;
    file.open(filename, ios::app);
    file << p0 << ";";
    file << p1 << ";";
    file << p2 << ";";
    file << p3 << ";";
    file << p4 << ";";
    file << p5 << ";";
    file << p6 << ";";
    file << p7 << ";";
    file << p8 << "\n\n";
    file.close();
}

#endif // HELPER_HPP
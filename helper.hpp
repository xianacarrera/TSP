#ifndef HELPER_HPP
#define HELPER_HPP

#include <bits/stdc++.h>
#include "problem.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

// Data structures 
typedef vector<int> vi;
typedef vector<vi> vvi;
#define all(a) (a).begin(), (a).end() // Aplicar a toda la estructura, e.g. sort(all(a))


vi random_solution(Problem * problem)
{
    vi solution;
    for (int i = 0; i < problem->n; i++) solution.push_back(i);

    random_shuffle(all(solution));
    return solution;
}

/*
int compute_length(Problem * problem, vi solution){
    int length = 0;
    for (int i = 0; i < problem->n; i++){
        length += problem->dist_matrix[solution[i]][solution[i + 1]];
    }
    length += problem->dist_matrix[solution[problem->n - 1]][solution[0]];
//    cout << problem->name << " length: " << length << "\n";
    return length;
}
*/

int compute_length(Problem * problem, vi solution) {
    int length = 0;
    vvi distance_matrix = problem->dist_matrix;
    int starting_node = solution[0];
    int from_node = starting_node;

    for (int i = 1; i < solution.size(); i++) {
        int to_node = solution[i];
        length += distance_matrix[from_node][to_node];
        from_node = to_node;
    }

    length += distance_matrix[starting_node][from_node];

    return length;
}

void compute_error(Problem * problem, vi solution){
    int length = compute_length(problem, solution);
    cout << problem->name << "\n";
    cout << "\tBest known: " << problem->best_known << "\n";
    cout << "\tSolution: " << length << "\n";
    cout << "\tError: " << 100.0f * (length - problem->best_known) / (float) problem->best_known << "%\n";
}


bool verify_sol(Problem * problem, vi sol){
    int sum = 0;

    for (int i = 0; i < problem->n; i++) {
        bool isIn = false;
        for (auto j : sol){
            if (j == i) isIn = true;
        }
        if (!isIn) return false;
        sum += 1;
    }

    return sum == problem->n;
}

void write_header(string filename, int seed){
    ofstream file;
    file.open(filename, ios::app);
    //file << method << "\n";
    file << seed << "\n";
    file.close();
}

void write_constants(string filename, float phi, float beta, float rho, int Q0){
    ofstream file;
    file.open(filename, ios::app);
    file << phi << ";";
    file << beta << ";";
    file << rho << ";";
    file << Q0 << "\n";
    file.close();
}

void write_results(string filename, Problem * problem, vi solution){
    ofstream file;
    file.open(filename, ios::app);
    int length = compute_length(problem, solution);
    file << problem->name << ";";
    file << problem->best_known << ";";
    file << length << ";";
    file << 100.0f * (length - problem->best_known) / (float) problem->best_known << "\n";
    file.close();
}

#endif // HELPER_HPP
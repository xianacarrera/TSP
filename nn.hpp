/*
 * Xiana Carrera Alonso
 * 17th AI Cup
 * 2022
 * 
 * Nearest Neighbours algorithm
 */


#ifndef NN_HPP
#define NN_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper.hpp"

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

using namespace std;

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


// NN starting from a random node. Returns the solution 
vi NN(int n, vvf &dist_matrix, int start_node){
    vi solution(n);
    vector<bool> visited(n, false);         // At first, all nodes are unvisited
    visited[start_node] = true;             // Except the starting node

    int last_node = start_node;     
    solution[0] = start_node;
    
    for (int i = 0; i < n - 1; i++){        // i is the number of unvisited nodes
        float min_dist = INF;
        int best_node;

        for (int j = 0; j < n; j++) {       // Find the closest node
            if (i == j || visited[j]) continue;         // If the node is the current node or already visited, skip it
            float node_dist = dist_matrix[last_node][j];
            if (node_dist < min_dist){
                min_dist = node_dist;
                best_node = j;
            }
        }

        solution[i + 1] = best_node;        // Add the closest node to the solution
        last_node = best_node;              // Update the last node
        visited[best_node] = true;          // Mark the node as visited
    }
    return solution;
}

// Function that only returns the length of the NN solution, without keeping track of the solution itself
float length_NN(int n, vvf &dist_matrix, int start_node){
    vector<bool> visited(n, false);
    visited[start_node] = true;
    float length = 0;

    int last_node = start_node;
    
    for (int i = 0; i < n - 1; i++){        
        float min_dist = INF;
        int best_node;

        for (int j = 0; j < n; j++) {
            if (i == j || visited[j]) continue;
            float node_dist = dist_matrix[last_node][j];
            if (node_dist < min_dist){
                min_dist = node_dist;
                best_node = j;
            }
        }

        last_node = best_node;
        visited[best_node] = true;
        length += min_dist;
    }
    return length;
}

// Best NN: returns the best solution of NN starting from all nodes
vi best_NN(int n, vvf &dist_matrix){
    float length_best_sol = INF;
    vi best_sol;
    for (int i = 0; i < n; i++){
        vi sol = NN(n, dist_matrix, i);
        int length_sol = compute_length(n, dist_matrix, sol);

        if (length_sol < length_best_sol){      // If the solution is better than the best solution, update the best solution
            length_best_sol = length_sol;
            best_sol = sol;
        }
    }
    return best_sol;
}

// Best NN: returns the length of the best solution of NN starting from all nodes. Does not keep track of the solution itself
float length_best_NN(int n, vvf &dist_matrix){
    float length_best_sol = INF;
    for (int i = 0; i < n; i++){
        float length_sol = length_NN(n, dist_matrix, i);
        if (length_sol < length_best_sol){
            length_best_sol = length_sol;
        }
    }
    return length_best_sol;
}

#endif /* NN_HPP */
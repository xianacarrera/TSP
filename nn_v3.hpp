#ifndef NN_HPP
#define NN_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper_v3.hpp"

#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

using namespace std;

// Data structures 
typedef vector<int> vi;
typedef vector<vector<float>> vvf;
typedef unordered_set<int> uset;
typedef pair<float, float> pf;
#define umap unordered_map
#define F first
#define S second
#define PB push_back
#define MP make_pair
#define all(a) (a).begin(), (a).end() // Aplicar a toda la estructura, e.g. sort(all(a))


vi NN(int n, vvf &dist_matrix, int start_node){
    vi solution(n);
    vector<bool> visited(n, false);
    visited[start_node] = true;

    int last_node = start_node;
    solution[0] = start_node;
    
    for (int i = 0; i < n - 1; i++){        // i is the number of unvisited nodes
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

        solution[i + 1] = best_node;
        last_node = best_node;
        visited[best_node] = true;
    }
    return solution;
}

float length_NN(int n, vvf &dist_matrix, int start_node){
    vector<bool> visited(n, false);
    visited[start_node] = true;
    float length = 0;

    int last_node = start_node;
    
    for (int i = 0; i < n - 1; i++){        // i is the number of unvisited nodes
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

vi best_NN(int n, vvf &dist_matrix)
{
    float length_best_sol = INF;
    vi best_sol;
    for (int i = 0; i < n; i++)
    {
        vi sol = NN(n, dist_matrix, i);
        int length_sol = compute_length(n, dist_matrix, sol);
        if (length_sol < length_best_sol)
        {
            length_best_sol = length_sol;
            best_sol = sol;
        }
    }
    return best_sol;
}

float length_best_NN(int n, vvf &dist_matrix)
{
    float length_best_sol = INF;
    for (int i = 0; i < n; i++)
    {
        vi sol = NN(n, dist_matrix, i);
        float length_sol = compute_length(n, dist_matrix, sol);
        if (length_sol < length_best_sol)
        {
            length_best_sol = length_sol;
        }
    }
    return length_best_sol;
}

#endif /* NN_HPP */
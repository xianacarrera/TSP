#ifndef NN_HPP
#define NN_HPP

#include <bits/stdc++.h>
#include "problem.hpp"
#include "helper.hpp"

using namespace std;

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


vi NN(Problem * problem, int start_node){
    float n = problem->n;
    vvi dist_matrix = problem->dist_matrix;

    vi solution(n);
    uset unvisited;
    for (int i = 0; i < n; i++){
        // Initialization in O(n) instead of O(nlog n)
        //unvisited.emplace_hint(unvisited.end(), i);
        unvisited.insert(i);
    }
    unvisited.erase(start_node);

    int last_node = start_node;
    solution[0] = start_node;
    
    for (int i = n - 1; i > 0; i--){        // i is the number of unvisited nodes
        int min_dist = INF;
        int best_node;

        uset::iterator it;
        for (it = unvisited.begin(); it != unvisited.end(); ++it) {
            int node_dist = dist_matrix[last_node][*it];
            if (node_dist < min_dist){
                min_dist = node_dist;
                best_node = *it;
            }
        }

        solution[n - i] = best_node;
        last_node = best_node;
        unvisited.erase(best_node);
    }
    return solution;
}

int length_NN(Problem * problem, int start_node){
    float n = problem->n;
    vvi dist_matrix = problem->dist_matrix;

    uset unvisited;
    for (int i = 0; i < n; i++){
        // Initialization in O(n) instead of O(nlog n)
        //unvisited.emplace_hint(unvisited.end(), i);
        unvisited.insert(i);
    }
    unvisited.erase(start_node);

    int last_node = start_node;
    int length = 0;
    
    for (int i = n - 1; i > 0; i--){        // i is the number of unvisited nodes
        int min_dist = INF;
        int best_node;

        uset::iterator it;
        for (it = unvisited.begin(); it != unvisited.end(); ++it) {
            int node_dist = dist_matrix[last_node][*it];
            if (node_dist < min_dist){
                min_dist = node_dist;
                best_node = *it;
            }
        }

        length += min_dist;
        last_node = best_node;
        unvisited.erase(best_node);
    }
    return length;
}

vi best_NN(Problem * problem)
{
    int length_best_sol = INF;
    vi best_sol;
    for (int i = 0; i < problem->n; i++)
    {
        vi sol = NN(problem, i);
        int length_sol = compute_length(problem, sol);
        if (length_sol < length_best_sol)
        {
            length_best_sol = length_sol;
            best_sol = sol;
        }
    }
    return best_sol;
}

int length_best_NN(Problem * problem)
{
    int length_best_sol = INF;
    for (int i = 0; i < problem->n; i++)
    {
        vi sol = NN(problem, i);
        int length_sol = compute_length(problem, sol);
        if (length_sol < length_best_sol)
        {
            length_best_sol = length_sol;
        }
    }
    return length_best_sol;
}

#endif /* NN_HPP */
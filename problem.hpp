#ifndef PROBLEM_H
#define PROBLEM_H

#include <bits/stdc++.h>
using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")


// Constants
const long long INF = 1e9;          // Infinity
const long double EPS = 1e-9;       // For floating-point comparisons, e.g. if(abs(a-b) < EPS)


// Data structures 
typedef vector<int> vi;
typedef vector<vi> vvi;
#define F first
#define S second

// Create a random number generator
mt19937 rng;


class Problem 
{
    public:
        string name;
        int n;
        int best_known;
        vector<pair<float, float>> points;
        vvi dist_matrix;

        Problem(string problem_name){
            this->name = problem_name;
            //read_input("C:/Users/Xia40/Desktop/Repositorios/TSP/problems/" + problem_name);
            read_input("./problems/" + problem_name);
            create_dist_matrix();
        }

    private:

        void read_input(string file_name){
            ifstream file(file_name);
            string line;
            float x, y;

            getline(file, line);      // NAME
            getline(file, line);      // TYPE
            getline(file, line);      // COMMENT

            getline(file, line);      // DIMENSION
            this->n = stoi(line.substr(11, line.length() - 11));

            getline(file, line);      // EDGE_WEIGHT_TYPE
            getline(file, line);      // BEST_KNOWN
            this->best_known = stoi(line.substr(13, line.length() - 13));
            getline(file, line);      // NODE_COORD_SECTION

            for (int i = 0; i < this->n; i++){
                getline(file, line);

                int space = line.find(" ");
                line = line.substr(space + 1, line.length() - space - 1);
                space = line.find(" ");

                x = stof(line.substr(0, space));
                y = stof(line.substr(space + 1, line.length() - space - 1));
                this->points.push_back(make_pair(x, y));
            }

            file.close();
        }

        void create_dist_matrix()
        {
            vvi dist_matrix(this->n, vi(this->n, 0));
            for (int i = 0; i < this->n; i++)
            {
                for (int j = 0; j < this->n; j++)
                {
                    dist_matrix[i][j] = round(sqrt(pow(this->points[i].F - this->points[j].F, 2) + pow(this->points[i].S - this->points[j].S, 2)));
                }
            }
            this->dist_matrix = dist_matrix;
        }
};


#endif // PROBLEM_H
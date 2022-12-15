/*
 * Xiana Carrera Alonso
 * 17th AI Cup
 * 2022
 * 
 * Problem class with parser
 */


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


// Abbreviations for data structures
typedef vector<int> vi;
typedef vector<vector<float>> vvf;
#define F first
#define S second


class Problem 
{
    public:
        string name;            // Name of the problem
        int n;                  // Number of nodes
        int best_known;         // Length of the best known solution
        vector<pair<float, float>> points;      // Coordinates of the nodes
        vvf dist_matrix;                        // Distance matrix

        Problem(string problem_name){
            this->name = problem_name;
            read_input("./problems/" + problem_name);
            create_dist_matrix();   
        }

    private:

        vector<string> split(string s, string delimiter) {
            size_t end;        // Position of the end of the token
            size_t start = 0;  // Position of the start of the token
            size_t delim_len = delimiter.length();

            string t;           // The token
            vector<string> res;     // The result

            // Find the first token
            while ((end = s.find(delimiter, start)) != string::npos) {
                t = s.substr (start, end - start);
                start = end + delim_len;        // Move the start position to the end of the token
                res.push_back(t);                  // Add the token to the result
            }

            res.push_back(s.substr(start));    // Add the last token to the result
            return res;
        }


        void read_input(string file_name){          // Read the input file line by line
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

                if (line.size() > 0 && line[0] == ' ') line = line.substr(1, line.length() - 1);     // Remove the first space if it exists
                auto splited = split(line, " ");        // Split the line by spaces

                // splited[0] is the index of the line
                x = stof(splited[1]);
                y = stof(splited[2]);

                this->points.push_back(make_pair(x, y));
            }

            file.close();
        }

        void create_dist_matrix(){
            vvf dist_matrix(this->n, vector<float>(this->n, 0.0f));
            for (int i = 0; i < this->n; i++){
                for (int j = 0; j < this->n; j++){     // Calculate the integer euclidean distance between the points
                    dist_matrix[i][j] = round(sqrt(pow(this->points[i].F - this->points[j].F, 2) + pow(this->points[i].S - this->points[j].S, 2)));
                }
            }
            this->dist_matrix = dist_matrix;
        }
};


#endif // PROBLEM_H
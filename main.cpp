#include <bits/stdc++.h>
#include "problem.hpp"
#include "local_search_v3.hpp"
#include "helper_v3.hpp"
#include "nn_v3.hpp"
//#include "aco.hpp"
#include "acov3.hpp"

using namespace std;

// Optimization flags
#pragma GCC optimize("O3,unroll-loops")         
#pragma GCC target("avx,avx2,fma")

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

const string files[10] = {"ch130.tsp", "d198.tsp", "eil76.tsp", "fl1577.tsp", "kroA100.tsp", "lin318.tsp", "pcb442.tsp", "pr439.tsp", "rat783.tsp", "u1060.tsp"};

float param0[10] = {0.1, 0.09, 0.35, 0.09, 0.09, 0.09, 0.09, 0.09, 0.19, 0.09};
float param1[10] = {2.9, 2.5, 3, 2.0, 2.0, 2.5, 2.5, 2.5, 2.0, 3};
float param2[10] = {0.14, 0.09, 0.08, 0.19, 0.09, 0.09, 0.09, 0.35, 0.35, 0.35};
int param3[10] = {97, 97, 97, 99, 98, 97, 98, 97, 98, 97};
bool param4[10] = {false, false, true, false, true, true, false, true, true, true};
int param5[10] = {10, 10, 10, 10, 10, 10, 15, 10, 20, 10};
int param6[10] = {1, 0, 0, 0, 0, 1, 1, 1, 0, 1};
int param7[10] = {5000, 5000, 5000, 400000, 5000, 5000, 400000, 400000, 400000, 400000};
int param8[10] = {100, 100, 100, 100, 100, 50, 50, 50, 100, 0};

int main(int argc, char** argv)
{
    int seed, problem_idx;
    if (argc < 3){
        cout << "Usage: ./main seed problem_idx" << endl;
        return 0;
    } else {
        seed = atoi(argv[1]);
        problem_idx = atoi(argv[2]);
    }

    srand(seed);

    //int n, vvi &dist_matrix, float alpha, float beta, float rho, int Q0, bool esaco, int m, int opt_tau0, int stopOnIter)
    Problem * problem = new Problem(files[problem_idx]);
    vi sol;
    if (problem_idx == 8 || problem_idx == 5){
        sol = ESACO_for_RAT(problem->n, problem->dist_matrix, param0[problem_idx], param1[problem_idx], param2[problem_idx], param3[problem_idx], param4[problem_idx], param5[problem_idx], param6[problem_idx], param7[problem_idx], param8[problem_idx], 2 * param8[problem_idx]);
    } else {
        sol = ESACO_v3(problem->n, problem->dist_matrix, param0[problem_idx], param1[problem_idx], param2[problem_idx], param3[problem_idx], param4[problem_idx], param5[problem_idx], param6[problem_idx], param7[problem_idx], param8[problem_idx], 2 * param8[problem_idx]);
    }

    // Print the parameters
    cout << param0[problem_idx] << " " << param1[problem_idx] << " " << param2[problem_idx] << " " << param3[problem_idx] << " " << param4[problem_idx] << " " << param5[problem_idx] << " " << param6[problem_idx] << " " << param7[problem_idx] << " " << param8[problem_idx] << "\n";

    if (!verify_sol(problem->n, sol)) {
        cout << "NOT A VALID SOLUTION FOR " << problem->name << "\n";
        exit(1);
    }
    float error = compute_error(problem, sol);
    cout << "Seed " << seed << "\n";
    cout << "Error for " << problem->name << ": " << error << "\n";
    
    write_results("errores_sin_esaco.csv", seed, problem->name, error);
    write_params("errores_sin_esaco.csv", param0[problem_idx], param1[problem_idx], param2[problem_idx], param3[problem_idx], param4[problem_idx], param5[problem_idx], param6[problem_idx], param7[problem_idx], param8[problem_idx]);
    param1[problem_idx] += 0.25;
    return 0;


    


    /*
            for (rho = 0.09; rho < 0.42; rho += 0.01){
                for (phi = 0.09; phi < 0.8; phi += 0.03){
    for (Q0 = 97; Q0 < 100; Q0 += 1){
                    for (beta = 2.5f; beta < 15.0f; beta += 1){
                    //string filename = "ESACO999_" + to_string(phi) + "_" + to_string(beta) + "_" + to_string(rho) + "_" + to_string(Q0) + ".csv";
                    string filename = "VFIN.csv";
                    write_constants(filename, phi, beta, rho, Q0);
                    write_header(filename, seed);

                    cout << phi << " " << beta << " " << rho << " " << Q0 << "\n";
                    for (int i = 0; i < 10; i++) {
                        seed += 1;
                        cout << "seed " << seed << "\n";
                        srand(seed);
                        rng.seed(seed);
                        Problem * problem = new Problem(problems[0]);        
                        
                        //float start_time = (float) clock();

                        //vi sol = NN(problem, rand()  % problem->n);
                        //vi sol = best_NN(problem);
                        //sol = two_five_opt(problem, sol);
                        //sol = two_opt_greedy(problem, sol);
                        //vi sol = ACO_candidate(problem);
                        vi sol = ESACO_v3(problem->n, problem->dist_matrix, phi, beta, rho, Q0, false);
                        
                        /*
                            while((((float) clock()) - start_time) / CLOCKS_PER_SEC < 180){
                                sol = two_five_opt_2(problem, sol).first;
                                sol = two_opt_greedy_2(problem, sol).first;
                            }


                        if (!verify_sol(problem->n, sol)) {
                            cout << "NOT A VALID SOLUTION FOR " << problem->name << "\n";
                            exit(1);
                        }
                        compute_error(problem, sol);
                        write_results(filename, problem, sol);

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
                }
            }
        }
    }

   /*
    Problem * problem = new Problem(custom[0]);
    vi sol = ACO(problem);
    compute_error(problem, sol);*/
}
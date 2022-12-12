#!/usr/bin/env bash

g++ -O3 -o main main.cpp

seed=5053  

total_runs=400

for j in $( eval echo {0..9} )
do
    ./main $seed 3
done


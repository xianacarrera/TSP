#!/usr/bin/env bash

# Xiana Carrera Alonso
# AI Cup script

g++ -O3 -o main main.cpp

seed=5053      # seed for random number generator

for j in $( eval echo {0..9} )
do
    ./main $seed $j
done


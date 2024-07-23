#!/bin/bash
# Script to spawn instances of nohup'ed fit_turnover_GLMM_reshuffle.R

# As arguments, provide min reshuff ID, max reshuff ID, and number of IDs per thread, e.g.,
# fit_turnover_GLMM_reshuffle.sh 1 400 10
# for 400 reshuffles, in 40 threads each running 10 iterations
# this spawns R processes with niceness 10, so that they are low priority. 
# because the .R script runs multithreaded up to using all cpus for reasons that are not clear

echo "You asked to run reshuffling IDs $1 to $2 by $3"

for ((min=$1;min<=$2;min=min+=$3)) # c-style for loop to work with arguments
do
    let max=$min+$3-1
    nohup nice -n 10 code/fit_turnover_GLMM_reshuffle.R $min $max > logs/fit_turnover_GLMM_reshuffle$min-$max.Rout &
    sleep 1 # wait so that reading the gzipped file in each process doesn't conflict with each other
done
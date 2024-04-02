#!/bin/bash
# Script to spawn instances of nohup'ed fit_turnover_GLMM_reshuffle.R

# for shuffle IDs 1 to 1000 in 50 threads
# for min in {1..981..20}
# do
#     let max=$min+19
#     nohup code/fit_turnover_GLMM_reshuffle.R $min $max > logs/fit_turnover_GLMM_reshuffle$min-$max.Rout &
#     sleep 1 # wait so that reading the gzipped file in each process doesn't conflict with each other
# done

# for shuffle IDs 1001 to 1400 in 20 threads
for min in {1001..1381..20}
do
    let max=$min+19
    nohup code/fit_turnover_GLMM_reshuffle.R $min $max > logs/fit_turnover_GLMM_reshuffle$min-$max.Rout &
    sleep 1 # wait so that reading the gzipped file in each process doesn't conflict with each other
done

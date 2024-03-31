#!/bin/bash
# Script to spawn 50 instances of nohup'ed turnover_GLMM_reshuffle.R
# for shuffle IDs 1 to 1000

for min in {1..981..20}
do
    let max=$min+19
    nohup code/turnover_GLMM_reshuffle.R $min $max > logs/turnover_GLMM_reshuffle$min-$max.Rout &
    sleep 1 # wait so that reading the gzipped file in each process doesn't conflict with each other
done

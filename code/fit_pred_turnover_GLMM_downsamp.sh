#!/bin/bash
# Script to spawn multiple instances of nohup'ed fit_pred_turnover_GLMM_downsamp.R
# As arguments, provide the name of model to fit, min downsamp ID, max downsamp ID, and number of IDs per thread, e.g.,
# fit_pred_turnover_GLMM_downsamp.sh modAllJtu 1 100 10

echo "You asked to run $1 with downsampling IDs $2 to $3 by $4"

for ((min=$2;min<=$3;min=min+=$4)) # c-style for loop to work with arguments
do
    max=$((min+$4-1))
    max=$((max<$3 ? max : $3)) # set max to the minimum of max or $3
    nohup code/fit_pred_turnover_GLMM_downsamp.R $1 $min $max > logs/fit_pred_turnover_GLMM_downsamp_$1_$min-$max.Rout &
    sleep 1 # wait so that reading the gzipped file in each process doesn't conflict with each other
done


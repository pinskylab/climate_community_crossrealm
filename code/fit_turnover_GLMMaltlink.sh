#!/bin/bash
# Script to spawn multiple instances of nohup'ed fit_turnover_GLMMaltlink.R
# As arguments, provide the names of models to fit, e.g.,
# fit_turnover_GLMMaltlink.sh modAllJtu modAllJbeta

echo "You asked to run $# models"

# Loop until all parameters are used up
while [ "$1" != "" ]; do
    echo "Spawning $1"
    nohup code/fit_turnover_GLMMaltlink.R $1 > logs/fit_turnover_GLMM$1.Rout &
    sleep 1 # wait so that reading the gzipped file in each process doesn't conflict with each other
    # Shift all the parameters down by one
    shift
done



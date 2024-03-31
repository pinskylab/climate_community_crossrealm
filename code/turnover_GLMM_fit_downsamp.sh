#!/bin/bash
# Script to spawn multiple instances of nohup'ed turnover_GLMM_fit_downsamp.R
# As arguments, provide (first) the downsampling ID and then the names of models to fit, e.g.,
# turnover_GLMM_fit_downsamp.sh 1 modAllJtu modAllJbeta

echo "You asked to run $(($#-1)) model(s) with downsampling ID $1"
bootID=$1
shift # shift all arguments down by one to get to the model names

# Loop until all parameters are used up
while [ "$1" != "" ]; do
    echo "Spawning $1"
    nohup code/turnover_GLMM_fit_downsamp.R $1 $bootID > logs/turnover_GLMM$1_boot$bootID.Rout &
    sleep 1 # wait so that reading the gzipped file in each process doesn't conflict with the next call
    # Shift all the arguments down by one
    shift
done



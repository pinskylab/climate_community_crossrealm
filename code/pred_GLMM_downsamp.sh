#!/bin/bash
# Script to spawn multiple instances of nohup'ed pred_GLMM_boot.R
# As arguments, provide the names of models to fit, e.g.,
# pred_GLMM.sh modsdTRealmAllJtu_boot1 modrawTsdTAllJtu_boot1

echo "You asked to run $# models"

# Loop until all parameters are used up
while [ "$1" != "" ]; do
    echo "Spawning $1"
    nohup code/pred_GLMM_boot.R $1 > logs/pred_GLMM_$1.Rout &

    # Shift all the parameters down by one
    shift
done



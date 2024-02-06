#!/bin/bash
# Script to spawn multiple instances of nohup'ed pred_GLMMlin.R
# As arguments, provide the names of models to fit, e.g.,
# pred_GLMMlin.sh modsdTRealmAllJtu modrawTsdTAllJtu

echo "You asked to run $# models"

# Loop until all parameters are used up
while [ "$1" != "" ]; do
    echo "Spawning $1"
    nohup code/pred_GLMMlin.R $1 > logs/pred_GLMM_$1.Rout &

    # Shift all the parameters down by one
    shift
done



#!/bin/bash
# Script to spawn multiple instances of nohup'ed turnover_GLMMlogduration_fit.R
# As arguments, provide the names of models to fit, e.g.,
# turnover_GLMM_logdurationfit.sh modLogDurAllJtu modLogDurAllHorn

echo "You asked to run $# models"

# Loop until all parameters are used up
while [ "$1" != "" ]; do
    echo "Spawning $1"
    nohup code/turnover_GLMMlogduration_fit.R $1 > logs/turnover_GLMM$1.Rout &

    # Shift all the parameters down by one
    shift
done



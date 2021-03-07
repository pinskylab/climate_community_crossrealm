#!/usr/bin/env Rscript

# Script to preict from glmmTMB models
# Set up to be run on the command line for one model at a time
# Argument is path to .rds file with model to run
# and path to csv of newdata
# nohup code/turnover_vs_temperature_GLMM_pred.R temp/modFullMaEnMiNPHuJtu.rds temp/newdata.csv > logs/turnover_vs_temperature_GLMMmodFullMaEnMiNPHuJtu_pred.Rout &
# (this works if code is executable, e.g., chmod u+x turnover_vs_temperature_GLMM_pred.R)
# (otherwise using nohup Rscript ...)

################
# Read command line arguments

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args)<1) stop("Have to specify a model to fit", call.=FALSE)
if (length(args)>2) stop("Have to specify only 1 model to fit and the new data", call.=FALSE)
modpath <- args[1]
datpath <- args[2]

###################################
# print basic info about the job

print(paste0('This is process #', Sys.getpid()))
print(Sys.time())
print(paste('Using model', modpath))
print(paste('Using data', datpath))


##############
# Prep

library(data.table) # for handling large datasets
library(glmmTMB) # for glmm models

###############################
# load data
mod <- readRDS(modpath)

newdat <- fread(datpath)
print(paste('Read', nrow(newdat), 'rows of new data'))

###############
# Make predictions
###############
preds <- predict(object = mod, newdata = newdat, se.fit = TRUE, type = 'response', re.form = NA)
newdat$pred <- preds$fit
newdat$pred.se <- preds$se.fit


##############
# save
##############
outpath <- paste0(gsub('.rds', '', modpath), '_preds.rds')
saveRDS(newdat, file = outpath)

print(paste('Wrote', outpath))


warnings()

print(paste('Ended', Sys.time()))
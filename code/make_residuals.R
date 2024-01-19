#!/usr/bin/Rscript --vanilla

# Script to make residuals from a glmmTMB model
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/make_residuals.R modInitAllJtu > logs/make_residuals_modInitAllJtu.Rout &
# (this works if code is executable, e.g., chmod u+x code/make_residuals.R)
# (otherwise using nohup Rscript ...)

# Read command line arguments ############

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 1)
    stop("Have to specify a model to fit", call. = FALSE)
if (length(args) > 1)
    stop("Have to specify only 1 model to fit", call. = FALSE)
modname <- args[1]
infile <- paste0('temp/', modname, '.rds')
outfile <- paste0('temp/', modname, '_residuals.rds')

# print basic info about the job ############################

print(paste('This is script make_residuals.R'))
print(paste('This is process #', Sys.getpid()))
print(Sys.time())


# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models when loading on Annotate
library(DHARMa) # for simulating residuals

# load model ###############################
mod <- readRDS(file = infile) # this will throw an error if the file doesn't exist


# load data ############################
# not sure whether this is needed

trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd
trendsall[, tsign := as.factor(tsign)]
trendsall[, duration.log := log(duration)]


# calculate initial dissimilarities
trendsall[, minduration := min(duration), by = rarefyID]
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu.sc)), by = .(rarefyID)] # initial dissimilarities
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])

# Choose dataset subsets
iallJtu <-
    trendsall[, complete.cases(
        Jtu.sc,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration.sc,
        microclim.sc,
        human_bowler.sc
    )]

iallHorn <-
    trendsall[, complete.cases(
        Horn.sc,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration.sc,
        microclim.sc,
        human_bowler.sc
    )]

trendsall[iallJtu, absLat.sc := scale(abs(rarefyID_y))]


# simulate residuals ########################
res <- simulateResiduals(mod, n = 250)
saveRDS(res, file = outfile)
    
print(paste0('saved ', outfile))
print(Sys.time())
print(warnings())



#### Plot residuals for model fit to all data
# not run here
# res <- readRDS('temp/modInitAllJtu_residuals.rds')
# res <- readRDS('temp/modrawTsdTTMERealmtsigninitAllJtu_residuals.rds')
# res <- readRDS('temp/modLogDurrawTsdTTMERealmtsigninitAllJtu_residuals.rds')
# res <- readRDS('temp/modSqInitAllJtu_residuals.rds')
# plot(res)
# plotResiduals(res, form=trendsall$duration[iallJtu], xlab = 'duration', main = '')
# testDispersion(res)
# tq <- testQuantiles(res, predictor = trendsall$duration[iallJtu], plot = F)


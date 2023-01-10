#!/usr/bin/Rscript --vanilla

# Script to fit glmmTMB models with thermal_bias
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmthermal_biasAllJtu.R modrawTsdTTRealmthermal_biasAllJtu > logs/turnover_vs_temperature_GLMMmodrawTsdTTRealmthermal_biasAllJtu.Rout &
# this works if code is executable, e.g., chmod u+x code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmthermal_biasAllJtu.R
# (otherwise using nohup Rscript ...)

# Read command line arguments ############

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 1)
    stop("Have to specify a model to fit", call. = FALSE)
if (length(args) > 1)
    stop("Have to specify only 1 model to fit", call. = FALSE)
fitmod <- args[1]
MATCHMOD <-
    FALSE # indicator to check if the argument matched a model name

# print basic info about the job ############################

print(paste('This is process #', Sys.getpid()))
print(Sys.time())



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(performance) # for R2

# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd
trendsall[, tsign := as.factor(tsign)]

# Models ############################

## Choose dataset
iallJtu <-
    trendsall[, complete.cases(
        Jtu.sc,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        thermal_bias.sc
    )]



# thermal_bias:tsign:T:sdT:REALM:duration #########################


# baseline
if (fitmod == 'modrawTsdTTRealmtsignAllJtu_thermal_biasdata') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration.sc +
            REALM:duration.sc +
            REALM:tsign:tempchange_abs.sc:duration.sc +
            REALM:tsign:tempave.sc:duration.sc +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration.sc +
            (duration.sc | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

# add thermal_bias:tempchange_abs:tsign
if (fitmod == 'modrawTsdTTRealmthermal_biassdTAllJtu_thermal_biasdata') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration.sc +
            REALM:duration.sc +
            REALM:tsign:tempchange_abs.sc:duration.sc +
            REALM:tsign:tempave.sc:duration.sc +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration.sc +
            REALM:thermal_bias.sc:tsign:duration.sc +
            REALM:thermal_bias.sc:tsign:tempchange_abs.sc:duration.sc +
            (duration.sc | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}


# print and save results ############################
if (MATCHMOD == FALSE)
    stop("Model name did not match anything", call. = FALSE)
if (MATCHMOD) {
    print(summary(mod))
    saveRDS(mod, file = paste0('temp/', fitmod, '.rds'))
    print(paste0('saved ', fitmod, '.rds'))
    print(Sys.time())
    print(warnings())
    if (!grepl('dredge', fitmod)) {
        print(performance::r2(mod)) # run if not a dredge object
    }
}

print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

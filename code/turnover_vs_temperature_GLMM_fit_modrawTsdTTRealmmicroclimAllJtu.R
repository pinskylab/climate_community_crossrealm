#!/usr/bin/Rscript --vanilla

# Script to fit glmmTMB modrawTsdTTRealmmicroclimAllJtu model
# nohup code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmmicroclimAllJtu.R XX > logs/turnover_vs_temperature_GLMMXX.Rout &
# where XX is a model name (see below)
# (this works if code is executable, e.g., chmod u+x code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmmicroclimAllJtu.R)
# (otherwise using nohup Rscript ...)

# print basic info about the job ############################

print('This is script turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmmicroclimAllJtu.R')
print(paste('This is process #', Sys.getpid()))
print(Sys.time())

# read arguments #############################
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 1)
    stop("Have to specify a model to fit", call. = FALSE)
if (length(args) > 1)
    stop("Have to specify only 1 model to fit", call. = FALSE)
fitmod <- args[1]
MATCHMOD <-
    FALSE # indicator to check if the argument matched a model name

print(paste('Fitting', fitmod))



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(performance) # for R2

# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz')

trendsall[, tsign := as.factor(tsign)]

# Models ############################

## Choose dataset
iallJtu <-
    trendsall[, complete.cases(
        Jtu.sc,
        tempchange_abs.sc,
        REALM,
        tempave_metab.sc,
        durationlog.sc,
        mass.sc,
        nspp.sc,
        seas.sc,
        microclim.sc,
        npp.sc,
        human_bowler.sc
    )]

## choose model

### with covariate by tempchange_abs only
if (fitmod == 'modrawTsdTTRealmmicroclimsdTAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave.sc:duration +
            REALM:tempave.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc:duration +
            REALM:microclim.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

### has tsign and interaction of covariate only with sdT
if (fitmod == 'modrawTsdTTRealmtsignmicroclimAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc:duration +
            REALM:microclim.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}


### with covariate by tempave:tempchange_abs
if (fitmod == 'modrawTsdTTRealmmicroclimAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave.sc:duration +
            REALM:tempave.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc:duration +
            REALM:microclim.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc:tempave.sc:duration +
            REALM:microclim.sc:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
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

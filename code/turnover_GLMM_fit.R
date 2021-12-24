#!/usr/bin/Rscript --vanilla

# Script to fit glmmTMB models
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_GLMM_fit.R modRealmAllJtu > logs/turnover_vs_temperature_GLMMmodRealmAllJtu.Rout &
# (this works if code is executable, e.g., chmod u+x turnover_GLMM_fit.R)
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

iallJne <-
    trendsall[, complete.cases(
        Jne.sc,
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

iallJbeta <-
    trendsall[, complete.cases(
        Jbeta.sc,
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

iallHorn <-
    trendsall[, complete.cases(
        Horn.sc,
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

# Trend model #################################
if (fitmod == 'modAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modAllJbeta') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJbeta), 'data points'))
    mod <- glmmTMB(
        Jbeta.sc ~ duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJbeta, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modAllJne') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJne), 'data points'))
    mod <- glmmTMB(
        Jne.sc ~ duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJne, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}



if (fitmod == 'modAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}



# REALM model #################################
if (fitmod == 'modRealmAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
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
}

print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

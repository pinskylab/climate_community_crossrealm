#!/usr/bin/Rscript --vanilla

#should be /usr/bin/Rscript --vanilla
# Script to fit glmmTMB model for modrawTsdTTRealmTaxamod2FEAllJtu
# Set up to be run on the command line
# nohup code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmTaxamod2FEAllJtu.R XX > logs/turnover_vs_temperature_GLMMXX.Rout &
# where XX is a model name (see below)
# this works if code is executable, e.g., chmod u+x code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmTaxamod2FEAllJtu.R
# (otherwise use nohup Rscript ...)

# print basic info about the job ############################

print(paste('This is script turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmTaxamod2FEAllJtu.R'))
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



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(performance) # for R2

# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz')

trendsall[, tsign := as.factor(tsign)]

# Model ############################

## Choose dataset
iallJtu <-
    trendsall[, complete.cases(
        Jtu.sc,
        tempchange_abs.sc,
        REALM,
        microclim.sc,
        human_bowler.sc
    )]


# with only interaction with tempchange_abs
if (fitmod == 'modrawTsdTTRealmTaxamod2FEsdTAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave.sc:duration +
            REALM:tempave.sc:tempchange_abs.sc:duration +
            REALM:duration:taxa_mod2 + # add taxamod2 as fixed effects
            REALM:tempchange_abs.sc:duration:taxa_mod2 + # add taxamod2 as fixed effects
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}

# with full interaction with tempave:tempchange_abs
if (fitmod == 'modrawTsdTTRealmTaxamod2FEAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave.sc:duration +
            REALM:tempave.sc:tempchange_abs.sc:duration +
            REALM:tempchange_abs.sc:duration:taxa_mod2 + # add taxamod2 as fixed effects
            REALM:tempave.sc:duration:taxa_mod2 +
            REALM:tempave.sc:tempchange_abs.sc:duration:taxa_mod2 +
            (duration | STUDY_ID / rarefyID),
            data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# print and save results ############################
if (MATCHMOD == FALSE)
    stop("Model name did not match anything", call. = FALSE)
if (MATCHMOD) {
    print(summary(mod))
    saveRDS(mod, file = paste0('temp/', fitmod, '_test.rds'))
    print(paste0('saved ', fitmod, '_test.rds'))
    print(Sys.time())
    print(warnings())
    if (!grepl('dredge', fitmod)) {
        print(performance::r2(mod)) # run if not a dredge object
    }
}

print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

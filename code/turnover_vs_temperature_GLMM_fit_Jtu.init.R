#!/usr/bin/Rscript --vanilla

# Script for sensitivity analyses: does including Jtu.init:duration alter the main results?
# nohup code/turnover_vs_temperature_GLMM_fit_Jtu.init.R XX > logs/turnover_vs_temperature_GLMMXX.Rout &
# where XX is a model name (see options below)
# (this works if code is executable, e.g., chmod u+x code/turnover_vs_temperature_GLMM_fit_Jtu.init.R)
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

print(paste('This is script turnover_vs_temperature_GLMM_fit_Jtu.init.R'))
print(paste('This is process #', Sys.getpid()))
print(Sys.time())



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models when loading on Annotate
library(performance) # for R2

# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd

trendsall[, tsign := as.factor(tsign)]

# calculate initial dissimilarities
trendsall[, minduration := min(duration), by = rarefyID]
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu.sc)), by = .(rarefyID)] # initial dissimilarities
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])

# Models ############################

## Choose dataset
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

## choose model

# Trend model #################################
if (fitmod == 'modInitAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modInitAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optimizer=optim, optArgs = list(method='BFGS'))
    )
    MATCHMOD <- TRUE
}



# REALM models #################################
if (fitmod == 'modRealmInitAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            REALM:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modRealmInitAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            REALM:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Taxa_mod models ###################################
if (fitmod == 'modTaxamod2InitAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            taxa_mod2:duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


if (fitmod == 'modTaxamod2InitAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            taxa_mod2:duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# tsign:sdT:realm #########################
# no environmental temperature, no realm
if (fitmod == 'modsdTtsignInitAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modsdTtsignInitAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


# realm, no environmental temperature
if (fitmod == 'modsdTRealmtsigninitAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modsdTRealmtsigninitAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


# environmental temperature
if (fitmod == 'modrawTsdTTRealmtsigninitAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modrawTsdTTRealmtsigninitAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# latitude instead of tempave #########################

# latitude
if (fitmod == 'modabsLatsdTabsLatRealmtsignInitAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:absLat.sc:duration +
            REALM:tsign:absLat.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# covariates #########################

### microclim
if (fitmod == 'modrawTsdTTRealmtsignmicroclimInitAllJtu') {
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
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

if (fitmod == 'modrawTsdTTRealmtsignhumanInitAllJtu') {
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            REALM:human_bowler.sc:duration +
            REALM:human_bowler.sc:tempchange_abs.sc:duration +
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

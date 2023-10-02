#!/usr/bin/Rscript --vanilla

# Script for sensitivity analyses: does including Jtu.init:duration + gainlossprop:duration alter the main results?
# nohup code/turnover_vs_temperature_GLMM_fit_Jtu.init.gainloss.R XX > logs/turnover_vs_temperature_GLMMXX.Rout &
# where XX is a model name (see options below)
# (this works if code is executable, e.g., chmod u+x code/turnover_vs_temperature_GLMM_fit_Jtu.init.gainloss.R)
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

print(paste('This is script turnover_vs_temperature_GLMM_fit_Jtu.init.gainloss.R'))
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

# add gains and losses data
load('data/biotime_blowes/bt-rarefy-collate.Rdata') # loads rarefied_alpha_medians and rarefied_beta_medians
rarefied_beta_medians <- as.data.table(rarefied_beta_medians)

# add gains and losses
trendsall <- merge(trendsall, rarefied_beta_medians[, .(rarefyID, year1 = as.numeric(YEAR1), year2 = as.numeric(YEAR2), 
                                                        gainlossprop = log10((gains+1)/(losses+1)))], 
                   by = c('rarefyID', 'year1', 'year2')) # merge in log10 proportion gains and losses (with +1 to avoid 0s). not clear why YEAR1 and YEAR2 got read as characters, but coercing to numeric seems fine.
nrow(trendsall)

# missing some gains and losses
trendsall[, sum(is.na(gainlossprop))]


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
        human_bowler.sc,
        gainlossprop
    )]

iallHorn <-
    trendsall[, complete.cases(
        Horn.sc,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration.sc,
        microclim.sc,
        human_bowler.sc,
        gainlossprop
    )]

trendsall[iallJtu, absLat.sc := scale(abs(rarefyID_y))]

## choose model

# Trend model #################################
if (fitmod == 'modInitGainLossAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            (duration | STUDY_ID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modInitGainLossAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            (duration | STUDY_ID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optimizer=optim, optArgs = list(method='BFGS'))
    )
    MATCHMOD <- TRUE
}



# REALM models #################################
if (fitmod == 'modRealmInitGainLossAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:duration +
            (duration | STUDY_ID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modRealmInitGainLossAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:duration +
            (duration | STUDY_ID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Taxa_mod models ###################################
if (fitmod == 'modTaxamod2InitGainLossAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            taxa_mod2:duration +
            (duration | STUDY_ID), 
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


if (fitmod == 'modTaxamod2InitGainLossAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            taxa_mod2:duration +
            (duration | STUDY_ID), 
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# tsign:sdT:realm #########################
# no environmental temperature, no realm
if (fitmod == 'modsdTtsignInitGainLossAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modsdTtsignInitGainLossAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


# realm, no environmental temperature
if (fitmod == 'modsdTRealmtsignInitGainLossAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modsdTRealmtsignInitGainLossAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


# environmental temperature
if (fitmod == 'modrawTsdTTRealmtsignInitGainLossAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modrawTsdTTRealmtsignInitGainLossAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
        data = trendsall[iallHorn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# latitude instead of tempave #########################

# latitude
if (fitmod == 'modabsLatsdTabsLatRealmtsignInitGainLossAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:absLat.sc:duration +
            REALM:tsign:absLat.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# covariates #########################

### microclim
if (fitmod == 'modrawTsdTTRealmtsignmicroclimInitGainLossAllJtu') {
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc:duration +
            REALM:microclim.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

if (fitmod == 'modrawTsdTTRealmtsignhumanInitGainLossAllJtu') {
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            Jtu.init:gainlossprop:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            REALM:human_bowler.sc:duration +
            REALM:human_bowler.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID),
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

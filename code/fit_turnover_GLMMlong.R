#!/usr/bin/Rscript --vanilla

# Script to fit glmmTMB models only to time-series >= 7 years
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/fit_turnover_GLMMlong.R modInitLongJtu > logs/turnover_GLMMlong_modInitLongJtu.Rout &
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

print(paste('This is script fit_turnover_GLMMlong.R'))
print(paste('This is process #', Sys.getpid()))
print(Sys.time())



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB) # for ME models when loading on Annotate
library(performance) # for R2

# load data ############################

# Turnover and covariates
trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd

trendsall[, tsign := as.factor(tsign)]

# calculate initial dissimilarities
trendsall[, minduration := min(duration), by = rarefyID]
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu.sc)), by = .(rarefyID)] # initial dissimilarities
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])

#calculate time series length
trendsall[, maxduration := max(duration), by = rarefyID]

# Models ############################

## Choose dataset
iallJtu <-
    trendsall[, maxduration >= 7 & complete.cases(
        Jtu,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration.sc,
        microclim.sc,
        human_bowler.sc
    )]

iallJtumarterr <-
    trendsall[, REALM %in% c('Terrestrial', 'Marine') & maxduration >= 7 & complete.cases(
        Jtu,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration,
        microclim.sc,
        human_bowler.sc
    )]

trendsall[iallJtu, absLat.sc := scale(abs(rarefyID_y))] # scale here so that only the included datasets are used

## choose model

# Baseline trend (Year) models #################################
if (fitmod == 'modYearRealmLongJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}



# Realm x Year models #################################
if (fitmod == 'modRealmxYearLongJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


# Taxa x Year models ###################################
if (fitmod == 'modTaxaxYearLongJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            taxa_mod2*duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}



# Tchange x Year x Realm models #########################
if (fitmod == 'modTchangexYearxRealmLongJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


# Tchange x Tave x Year x Realm models #########################
if (fitmod == 'modTchangexTavexYearxRealmLongJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modTchangexTavexYearxRealmLongJtu_marterr') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtumarterr), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtumarterr, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}



# Environmental covariates #########################
### microclim
if (fitmod == 'modMicroclimLongJtu') {
    print(paste(sum(iallJtumarterr), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            REALM*microclim.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtumarterr, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

### human
if (fitmod == 'modHumanLongJtu') {
    print(paste(sum(iallJtumarterr), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            REALM*human_bowler.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtumarterr, ],
        family = ordbeta(link = 'logit'),
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

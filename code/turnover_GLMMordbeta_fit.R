#!/usr/bin/Rscript --vanilla

# Script to fit glmmTMB models with ordbeta errors
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_GLMMordbeta_fit.R modOBInitAllJtu > logs/turnover_GLMMmodOBInitAllJtu.Rout &
# (this works if code is executable, e.g., chmod u+x code/turnover_GLMMordbeta_fit.R)
# (otherwise using nohup Rscript ...)
# Note: this requires a newer version of glmmTMB, e.g, 1.1.8 (installed on Annotate2 but not Annotate)

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

print(paste('This is script turnover_GLMMordbeta_fit.R'))
print(paste('This is process #', Sys.getpid()))
print(Sys.time())



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB) # for ME ordbeta models
library(performance) # for R2

# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd

trendsall[, tsign := as.factor(tsign)]
trendsall[, duration.log := log(duration)]

# calculate initial dissimilarities
trendsall[, minduration := min(duration), by = rarefyID]
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu)), by = .(rarefyID)] # initial dissimilarities
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])

# Models ############################

## Choose dataset
iallJtu <-
    trendsall[, complete.cases(
        Jtu,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration,
        microclim.sc,
        human_bowler.sc
    )]

iallHorn <-
    trendsall[, complete.cases(
        Horn,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration,
        microclim.sc,
        human_bowler.sc
    )]

trendsall[iallJtu, absLat.sc := scale(abs(rarefyID_y))]

print(warnings())

## choose model

# Baseline trend (null) model with log time term #################################
if (fitmod == 'modOBInitAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBInitAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            Jtu.init:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optimizer=optim, optArgs = list(method='BFGS'))
    )
    MATCHMOD <- TRUE
}


# Realm models #################################
if (fitmod == 'modOBRealmInitAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBRealmInitAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            Jtu.init:duration +
            REALM:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBRRealmInitAllJtu') { # adds realm main effect
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM +
            REALM:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBRRealmInitAllHorn') { # adds realm main effect
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            Jtu.init:duration +
            REALM +
            REALM:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Taxa models ###################################
if (fitmod == 'modOBTaxamod2InitAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            taxa_mod2:duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


if (fitmod == 'modOBTaxamod2InitAllHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            Jtu.init:duration +
            taxa_mod2:duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBTTaxamod2InitAllJtu') { # adds taxa main effect
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            taxa_mod2 +
            taxa_mod2:duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


if (fitmod == 'modOBTTaxamod2InitAllHorn') { # adds taxa main effect
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            Jtu.init:duration +
            taxa_mod2 +
            taxa_mod2:duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# Tchange models with main effects #########################
if (fitmod == 'modOBsdTMERealmtsigninitAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBsdTMERealmtsigninitAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBsdTMERtsRealmtsigninitAllJtu') { # also tsign and REALM
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            REALM +
            tsign +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBsdTMERtsRealmtsigninitAllHorn') { # also tsign and REALM
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            REALM +
            tsign +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Tchange x Tave models with main effects #########################
if (fitmod == 'modOBrawTsdTTMERealmtsigninitAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERealmtsigninitAllHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsigninitAllJtu') { # also tsign and REALM
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            REALM + 
            tsign + 
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modOBrawTsdTTMERtsRealmtsigninitAllHorn') {# also tsign and REALM
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ duration +
            REALM +
            tsign +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Tchange x latitude with main effects #########################
if (fitmod == 'modOBabsLatsdTabsLatMERealmtsignInitAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:absLat.sc +
            REALM:tsign:absLat.sc:duration +
            REALM:tsign:absLat.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBabsLatsdTabsLatMERtsRealmtsignInitAllJtu') { # also REALM and tsign
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            REALM +
            tsign +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:absLat.sc +
            REALM:tsign:absLat.sc:duration +
            REALM:tsign:absLat.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Environmental covariates with main effects #########################
### microclim
if (fitmod == 'modOBrawTsdTTMERealmtsignmicroclimInitAllJtu') {
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc +
            REALM:microclim.sc:duration +
            REALM:microclim.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

### human
if (fitmod == 'modOBrawTsdTTMERealmtsignhumanInitAllJtu') {
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            REALM:human_bowler.sc +
            REALM:human_bowler.sc:duration +
            REALM:human_bowler.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

### microclim, also realm and tsign main effects
if (fitmod == 'modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu') {
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            REALM +
            tsign +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc +
            REALM:microclim.sc:duration +
            REALM:microclim.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

### human, also realm and tsign main effects
if (fitmod == 'modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu') {
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            REALM +
            tsign +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave.sc +
            REALM:tsign:tempave.sc:duration +
            REALM:tsign:tempave.sc:tempchange_abs.sc:duration +
            REALM:human_bowler.sc +
            REALM:human_bowler.sc:duration +
            REALM:human_bowler.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
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

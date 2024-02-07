#!/usr/bin/Rscript --vanilla

# Script to fit glmmTMB models with alternative link functions
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_GLMMlin_fit.R modOBInitAllJtu_lin > logs/turnover_GLMMmodOBInitAllJtu_lin.Rout &
# (this works if code is executable, e.g., chmod u+x code/turnover_GLMMlin_fit.R)
# or use the shell script to spawn multiple model fits at once:
# code/turnover_GLMMlin_fit.sh modOBRInitAllJtu_lin modOBRInitAllHorn_lin 

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

print(paste('This is script turnover_GLMMlin_fit.R'))
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

# calculate initial dissimilarities
trendsall[, minduration := min(duration), by = rarefyID]
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu)), by = .(rarefyID)] # initial dissimilarities
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])

# add gains and losses data
load('data/biotime_blowes/bt-rarefy-collate.Rdata') # loads rarefied_alpha_medians and rarefied_beta_medians
rarefied_beta_medians <- as.data.table(rarefied_beta_medians)

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



# Baseline trend (null) model  #################################
## Gaussian
if (fitmod == 'modOBRInitAllJtu_lin') { # has realm main effects
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = gaussian(link='identity'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


## Probit
if (fitmod == 'modOBRInitAllJtu_probit') { # has realm main effects
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link='probit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

## cloglog
if (fitmod == 'modOBRInitAllJtu_cloglog') { # has realm main effects
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ duration +
            Jtu.init:duration +
            REALM +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link='cloglog'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

## log
if (fitmod == 'modOBRInitAllJtu_log') { # has realm main effects
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            REALM +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = lognormal(link='log'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Realm x Year models #################################
if (fitmod == 'modOBRRealmInitAllJtu_lin') { # adds realm main effect
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
        family = gaussian(link='identity'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBRRealmInitAllJtu_probit') { # adds realm main effect
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
        family = ordbeta(link='probit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBRRealmInitAllJtu_cloglog') { # adds realm main effect
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
        family = ordbeta(link='cloglog'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBRRealmInitAllJtu_log') { # adds realm main effect
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            REALM +
            REALM:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = lognormal(link='log'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Taxa models ###################################
if (fitmod == 'modOBTTaxamod2InitAllJtu_lin') { # adds taxa main effect
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
        family = gaussian(link='identity'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBTTaxamod2InitAllJtu_probit') { # adds taxa main effect
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
        family = ordbeta(link='probit'),
        dispformula = ~ REALM#,
    ) # add dispersion formula
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBTTaxamod2InitAllJtu_cloglog') { # adds taxa main effect
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
        family = ordbeta(link='cloglog'),
        dispformula = ~ REALM#,
    ) # add dispersion formula
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBTTaxamod2InitAllJtu_log') { # adds taxa main effect
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init:duration +
            taxa_mod2 +
            taxa_mod2:duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallJtu, ],
        family = lognormal(link='log'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# Tchange models with main effects #########################
if (fitmod == 'modOBsdTMERtsRealmtsigninitAllJtu_lin') { # also tsign and REALM
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
        family = gaussian(link='identity'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBsdTMERtsRealmtsigninitAllJtu_probit') { # also tsign and REALM
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
        family = ordbeta(link='probit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBsdTMERtsRealmtsigninitAllJtu_cloglog') { # also tsign and REALM
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
        family = ordbeta(link='cloglog'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBsdTMERtsRealmtsigninitAllJtu_log') { # also tsign and REALM
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM +
            tsign +
            Jtu.init:duration +
            REALM:tsign:tempchange_abs.sc +
            REALM:tsign:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = lognormal(link='log'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Tchange x Tave models with main effects #########################
if (fitmod == 'modOBrawTsdTTMERtsRealmtsigninitAllJtu_lin') { # also tsign and REALM
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
        family = gaussian(link='identity'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsigninitAllJtu_probit') { # also tsign and REALM
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
        family = ordbeta(link='probit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsigninitAllJtu_cloglog') { # also tsign and REALM
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
        family = ordbeta(link='cloglog'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsigninitAllJtu_log') { # also tsign and REALM
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
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
        family = lognormal(link='log'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


# Environmental covariates with main effects #########################
### microclim, also realm and tsign main effects
if (fitmod == 'modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu_lin') {
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
        family = gaussian(link='identity'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu_probit') {
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
        family = ordbeta(link='probit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu_cloglog') {
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
        family = ordbeta(link='cloglog'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu_log') {
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
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
        family = lognormal(link='log'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

### human, also realm and tsign main effects
if (fitmod == 'modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu_lin') {
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
        family = gaussian(link='identity'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu_probit') {
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
        family = ordbeta(link='probit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu_cloglog') {
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
        family = ordbeta(link='cloglog'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}

if (fitmod == 'modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu_log') {
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
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
        family = lognormal(link='log'),
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

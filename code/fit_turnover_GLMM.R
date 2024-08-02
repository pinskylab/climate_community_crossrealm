#!/usr/bin/Rscript --vanilla

# Script to fit glmmTMB models with ordbeta errors
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/fit_turnover_GLMM.R modOBInitAllJtu > logs/fit_turnover_GLMMmodOBInitAllJtu.Rout &
# (this works if code is executable, e.g., chmod u+x code/fit_turnover_GLMM.R)
# or use the shell script to spawn multiple model fits at once:
# code/fit_turnover_GLMM.sh modOBRInitAllJtu modOBRInitAllHorn 
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

print(paste('This is script turnover_GLMM_fit.R'))
print(paste('This is process #', Sys.getpid()))
print(Sys.time())



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB) # for ME ordbeta models
library(performance) # for R2

# load data ############################

# Turnover and covariates
trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd

trendsall[, tsign := as.factor(tsign)]

# calculate initial dissimilarities
trendsall[, minduration := min(duration), by = rarefyID]
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu)), by = .(rarefyID)] # initial dissimilarities
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])

# add gains and losses data
load('data/biotime_blowes/bt-rarefy-collate.Rdata') # loads rarefied_alpha_medians and rarefied_beta_medians. From 03_collate_resamps_cluster.R
rarefied_beta_medians <- as.data.table(rarefied_beta_medians)

# add gains and losses
trendsall <- merge(trendsall, rarefied_beta_medians[, .(rarefyID, year1 = as.numeric(YEAR1), year2 = as.numeric(YEAR2), 
                                                        gainlossprop = log10((gains+1)/(losses+1)))], 
                   by = c('rarefyID', 'year1', 'year2')) # merge in log10 proportion gains and losses (with +1 to avoid 0s). not clear why YEAR1 and YEAR2 got read as characters, but coercing to numeric seems fine.
nrow(trendsall)

# missing some gains and losses
trendsall[, sum(is.na(gainlossprop))]

#calculate time series length
trendsall[, maxduration := max(duration), by = rarefyID]

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

iGLallJtu <-
    trendsall[, complete.cases(
        Jtu,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration,
        microclim.sc,
        human_bowler.sc,
        gainlossprop
    )]

iallJtumarterr <-
    trendsall[, REALM %in% c('Terrestrial', 'Marine') & complete.cases(
        Jtu,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration,
        microclim.sc,
        human_bowler.sc
    )]

iallJtuLong <-
    trendsall[, maxduration >= 7 & complete.cases(
        Jtu,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration,
        microclim.sc,
        human_bowler.sc
    )]

trendsall[iallJtu, absLat.sc := scale(abs(rarefyID_y))]

print(warnings())



# Baseline trend (Year) models  #################################
if (fitmod == 'modYearRealmJtu') { # has realm main effect
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
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

if (fitmod == 'modYearRealmHorn') { # has realm main effect
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ Jtu.init*duration +
            REALM +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optimizer=optim, optArgs = list(method='BFGS'))
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modYearRealmJtu_marterr') { # has realm effect. fit only to marine and terrestrial
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtumarterr), 'data points'))
    print(paste(trendsall[iallJtumarterr, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtumarterr, length(unique(rarefyID))], 'time series'))
    print(paste(trendsall[iallJtumarterr, paste(unique(REALM), collapse=',')], 'realms'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtumarterr, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Realm x Year models #################################
if (fitmod == 'modRealmxYearJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
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

if (fitmod == 'modRealmxYearHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ Jtu.init*duration +
            REALM*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Taxon x Year models ###################################
if (fitmod == 'modTaxaxYearJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
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


if (fitmod == 'modTaxaxYearHorn') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ Jtu.init*duration +
            taxa_mod2*duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# (Tchange + Year) x Realm models #########################
# no interaction of Tchange x Year
if (fitmod == 'modTchangexRealmJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*duration+
            REALM*tsign*tempchange_abs.sc +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modTchangexRealmHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ Jtu.init*duration +
            REALM*duration+
            REALM*tsign*tempchange_abs.sc +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# Tchange x Year x Realm models #########################
if (fitmod == 'modTchangexYearxRealmJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modTchangexYearxRealmHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ Jtu.init*duration +
            REALM*tsign*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}

# (Tchange + Tave) x Year x Realm models #########################
# does not have Tchange x Tave x Realm x Year
if (fitmod == 'modTchangexTavexRealmJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempchange_abs.sc*duration +
            REALM*tempave.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modTchangexTavexRealmHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ Jtu.init*duration +
            REALM*tsign*tempchange_abs.sc*duration +
            REALM*tempave.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}


# Tchange x Tave x Year x Realm models #########################
if (fitmod == 'modTchangexTavexYearxRealmJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modTchangexTavexYearxRealmHorn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallHorn), 'data points'))
    print(paste(trendsall[iallHorn, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallHorn, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Horn ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallHorn, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modTchangexTavexYearxRealmJtu_marterr') { # only marine and terrestrial
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtumarterr), 'data points'))
    print(paste(trendsall[iallJtumarterr, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtumarterr, length(unique(rarefyID))], 'time series'))
    print(paste(trendsall[iallJtumarterr, paste(unique(REALM), collapse=',')], 'realms'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtumarterr, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}



# Tchange x latitude x Year x Realm #########################
if (fitmod == 'modTchangexLatxYearxRealmJtu') { # also REALM and tsign
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    print(paste(trendsall[iallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*absLat.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}

# Environmental covariate models #########################
### microclimate, only mar and terr
if (fitmod == 'modMicroclimJtu') {
    print(paste(sum(iallJtumarterr), 'data points'))
    print(paste(trendsall[iallJtumarterr, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtumarterr, length(unique(rarefyID))], 'time series'))
    print(paste(trendsall[iallJtumarterr, paste(unique(REALM), collapse=',')], 'realms'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            REALM*microclim.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtumarterr, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
        )
    MATCHMOD <- TRUE
}

### human, only mar and terr
if (fitmod == 'modHumanJtu') {
    print(paste(sum(iallJtumarterr), 'data points'))
    print(paste(trendsall[iallJtumarterr, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iallJtumarterr, length(unique(rarefyID))], 'time series'))
    print(paste(trendsall[iallJtumarterr, paste(unique(REALM), collapse=',')], 'realms'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            REALM*human_bowler.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtumarterr, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
        )
    MATCHMOD <- TRUE
}


# Gainloss models #########################
if (fitmod == 'modYearRealmGLJtu') { # Years model
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iGLallJtu), 'data points'))
    print(paste(trendsall[iGLallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iGLallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*gainlossprop*duration +
            REALM +
            (duration | STUDY_ID),
        data = trendsall[iGLallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modRealmxYearGLJtu') { # Realm x Year model
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iGLallJtu), 'data points'))
    print(paste(trendsall[iGLallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iGLallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*gainlossprop*duration +
            REALM*duration +
            (duration | STUDY_ID),
        data = trendsall[iGLallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modTaxaxYearGLJtu') { # Taxon x Year model
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iGLallJtu), 'data points'))
    print(paste(trendsall[iGLallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iGLallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*gainlossprop*duration +
            taxa_mod2*duration +
            (duration | STUDY_ID), 
        data = trendsall[iGLallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM#,
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


if (fitmod == 'modTchangexYearxRealmGLJtu') { # Tchange x Year x Realm model
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iGLallJtu), 'data points'))
    print(paste(trendsall[iGLallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iGLallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*gainlossprop*duration +
            REALM*tsign*tempchange_abs.sc*duration +
            (duration | STUDY_ID),
        data = trendsall[iGLallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modTchangexTavexYearxRealmGLJtu') { # Tchange x Tave x Year x Realm model
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iGLallJtu), 'data points'))
    print(paste(trendsall[iGLallJtu, length(unique(STUDY_ID))], 'studies'))
    print(paste(trendsall[iGLallJtu, length(unique(rarefyID))], 'time series'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*gainlossprop*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID),
        data = trendsall[iGLallJtu, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
    )
    MATCHMOD <- TRUE
}

# Long models ###########################
if (fitmod == 'modYearRealmLongJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtuLong), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtuLong, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modRealmxYearLongJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtuLong), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtuLong, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modTaxaxYearLongJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtuLong), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            taxa_mod2*duration +
            (duration | STUDY_ID / rarefyID), 
        data = trendsall[iallJtuLong, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


if (fitmod == 'modTchangexYearxRealmLongJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtuLong), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtuLong, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

if (fitmod == 'modTchangexTavexYearxRealmLongJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtuLong), 'data points'))
    mod <- glmmTMB(
        Jtu ~ Jtu.init*duration +
            REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtuLong, ],
        family = ordbeta(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)) # to address iteration limit reached without convergence
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
    if (!grepl('dredge', fitmod)) {
        print(performance::r2(mod)) # run if not a dredge object
    }
}

print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

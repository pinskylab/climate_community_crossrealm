#!/usr/bin/env Rscript

#should be /usr/bin/Rscript --vanilla
# Script to fit glmmTMB models
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_vs_temperature_GLMM_fit.R modRFgauss > logs/turnover_vs_temperature_GLMMmodRFgauss.Rout &
# (this works if code is executable, e.g., chmod u+x turnover_vs_temperature_GLMM_fit.R)
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
trendsall <- fread('output/turnover_w_covariates.csv.gz')

trendsall[, tsign := as.factor(tsign)]
trendsall[, rarefyID_y_abs := abs(rarefyID_y)]

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

## Compare variance structures with duration.sc ############################
fixed <-
    'Jtu.sc ~ duration.sc + REALM:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + mass.sc:duration.sc + nspp.sc:duration.sc + seas.sc:duration.sc + microclim.sc:duration.sc + npp.sc:duration.sc + human_bowler.sc:REALM2:duration.sc'

if (fitmod == 'modRFgauss10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFgauss10 <-
        glmmTMB(
            formula(fixed),
            data = trends10[i10, ],
            family = gaussian(),
            control = glmmTMBControl(profile = TRUE)
        )
    summary(modRFgauss10)
    saveRDS(modRFgauss10, file = 'temp/modRFgauss10.rds')
    print('saved modRFgauss10.rds')
    print(Sys.time())
    print(performance::r2(modRFgauss10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRFbeta10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFbeta10 <-
        glmmTMB(
            formula(fixed),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            control = glmmTMBControl(profile = TRUE)
        ) # add beta errors
    summary(modRFbeta10)
    saveRDS(modRFbeta10, file = 'temp/modRFbeta10.rds')
    print('saved modRFbeta10.rds')
    print(Sys.time())
    print(performance::r2(modRFbeta10))
    MATCHMOD <- TRUE
}

if (fitmod == 'modRFrID10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFrID10 <-
        glmmTMB(
            formula(paste(fixed, '+(1|rarefyID)')),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            control = glmmTMBControl(profile = TRUE)
        ) # add random effects
    summary(modRFrID10)
    saveRDS(modRFrID10, file = 'temp/modRFrID10.rds')
    print('saved modRDrID10.rds')
    print(Sys.time())
    print(performance::r2(modRFrID10))
    MATCHMOD <- TRUE
}

if (fitmod == 'modRF2lev10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRF2lev10 <-
        glmmTMB(
            formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            control = glmmTMBControl(profile = TRUE)
        ) # add nested random effects
    summary(modRF2lev10)
    saveRDS(modRF2lev10, file = 'temp/modRF2lev10.rds')
    print('saved modRF2lev10.rds')
    print(Sys.time())
    print(performance::r2(modRF2lev10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRFnestedRE10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFnestedRE10 <-
        glmmTMB(
            formula(paste(
                fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)'
            )),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            control = glmmTMBControl(profile = TRUE)
        ) # add nested random effects
    summary(modRFnestedRE10)
    saveRDS(modRFnestedRE10, file = 'temp/modRFnestedRE10.rds')
    print('saved modRFnestedRE10.rds')
    print(Sys.time())
    print(performance::r2(modRFnestedRE10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRFslopeRE10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFslopeRE10 <-
        glmmTMB(
            formula(
                paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')
            ),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            control = glmmTMBControl(profile = TRUE)
        ) # add random slopes
    summary(modRFslopeRE10)
    saveRDS(modRFslopeRE10, file = 'temp/modRFslopeRE10.rds')
    print('saved modRFslopeRE10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRFslopeRE2lev10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFslopeRE2lev10 <-
        glmmTMB(
            formula(paste(
                fixed, '+(duration.sc|STUDY_ID/rarefyID)'
            )),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            control = glmmTMBControl(profile = TRUE)
        ) # add random slopes
    summary(modRFslopeRE2lev10)
    saveRDS(modRFslopeRE2lev10, file = 'temp/modRFslopeRE2lev10.rds')
    print('saved modRFslopeRE2lev10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2lev10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRFdisp10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFdisp10 <-
        glmmTMB(
            formula(fixed),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc,
            control = glmmTMBControl(profile = TRUE)
        )
    summary(modRFdisp10)
    saveRDS(modRFdisp10, file = 'temp/modRFdisp10.rds')
    print('saved modRFdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRF2levdisp10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFdisp2levOnlyint10 <-
        glmmTMB(
            formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc,
            control = glmmTMBControl(profile =
                                         TRUE)
        )
    summary(modRFdisp2levOnlyint10)
    saveRDS(modRFdisp2levOnlyint10, file = 'temp/modRFdisp2levOnlyint10.rds')
    print('saved modRFdisp2levOnlyint10.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp2levOnlyint10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRFslopeRE2levdisp10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFslopeRE2levdisp10 <-
        glmmTMB(
            formula(paste(
                fixed, '+(duration.sc|STUDY_ID/rarefyID)'
            )),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc,
            control = glmmTMBControl(profile =
                                         TRUE)
        ) # add dispersion formula
    summary(modRFslopeRE2levdisp10)
    saveRDS(modRFslopeRE2levdisp10, file = 'temp/modRFslopeRE2levdisp10.rds')
    print('saved modRFslopeRE2levdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisp10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRFslopeREdisp10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFslopeREdisp10 <-
        glmmTMB(
            formula(
                paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')
            ),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc,
            control = glmmTMBControl(profile =
                                         TRUE)
        ) # add dispersion formula
    summary(modRFslopeREdisp10)
    saveRDS(modRFslopeREdisp10, file = 'temp/modRFslopeREdisp10.rds')
    print('saved modRFslopeREdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeREdisp10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRF2levdisprealm10') {
    # 2 level RE, slope vs. duration, disp formula by realm
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRF2levdisprealm10 <-
        glmmTMB(
            formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc + REALM,
            control = glmmTMBControl(profile =
                                         TRUE)
        ) # add dispersion formula
    summary(modRF2levdisprealm10)
    saveRDS(modRF2levdisprealm10, file = 'temp/modRF2levdisprealm10.rds')
    print('saved modRF2levdisprealm10.rds')
    print(Sys.time())
    print(performance::r2(modRF2levdisprealm10))
    MATCHMOD <- TRUE
}


if (fitmod == 'modRFslopeRE2levdisprealm10') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    modRFslopeRE2levdisprealm10 <-
        glmmTMB(
            formula(paste(
                fixed, '+(duration.sc|STUDY_ID/rarefyID)'
            )),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc + REALM,
            control = glmmTMBControl(profile =
                                         TRUE)
        ) # add dispersion formula
    summary(modRFslopeRE2levdisprealm10)
    saveRDS(modRFslopeRE2levdisprealm10, file = 'temp/modRFslopeRE2levdisprealm10.rds')
    print('saved modRFslopeRE2levdisprealm10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisprealm10))
    MATCHMOD <- TRUE
}



# REALM models #################################
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

# latitude:sdT:REALM:duration #########################
# use lat instead of T
if (fitmod == 'modLatsdTTRealmAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:rarefyID_y:duration +
            REALM:rarefyID_y:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# use abs(lat) instead of T
if (fitmod == 'modabsLatsdTTRealmAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:rarefyID_y_abs:duration +
            REALM:rarefyID_y_abs:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# T:dT/sdT:REALM INTERACTION WITH SINGLE FACTORS ########################
# T:dT:REALM:duration like Antao #########################

# use tempchange like Antao
if (fitmod == 'modTdTTRealmAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    ) #,
    #  control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# use tempchange like Antao. scale duration
if (fitmod == 'modTdTTRealmDurscAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration.sc +
            REALM:duration.sc +
            REALM:tempchange.sc:duration.sc +
            REALM:tempave_metab.sc:duration.sc +
            REALM:tempave_metab.sc:tempchange.sc:duration.sc +
            (duration.sc | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    ) #,
    # control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}


# T:sdT:REALM:duration #########################
# use tempchange_abs (otherwise like Antao)
if (fitmod == 'modTsdTTRealmAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    ) #,
    # control = glmmTMBControl(profile = TRUE)) # add dispersion formula
MATCHMOD <- TRUE
}

# use tempchange_abs and scaled duration
if (fitmod == 'modTsdTTRealmDurscAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration.sc +
            REALM:duration.sc +
            REALM:tempchange_abs.sc:duration.sc +
            REALM:tempave_metab.sc:duration.sc +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration.sc +
            (duration.sc | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# tsign:T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmtsignAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tsign:tempchange_abs.sc:duration +
            REALM:tsign:tempave_metab.sc:duration +
            REALM:tsign:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    ) #,
    #  control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# thermal_bias:T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmthermal_biasAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            REALM:thermal_bias.sc:duration +
            REALM:thermal_bias.sc:tempchange_abs.sc:duration +
            REALM:thermal_bias.sc:tempave_metab.sc:duration +
            REALM:thermal_bias.sc:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )#,
    #control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# npp:T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmnppAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            REALM:npp.sc:duration +
            REALM:npp.sc:tempchange_abs.sc:duration +
            REALM:npp.sc:tempave_metab.sc:duration +
            REALM:npp.sc:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) 
    MATCHMOD <- TRUE
}

# seas:T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmseasAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            REALM:seas.sc:duration +
            REALM:seas.sc:tempchange_abs.sc:duration +
            REALM:seas.sc:tempave_metab.sc:duration +
            REALM:seas.sc:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    ) #,
    # control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# microclim:T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmmicroclimAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc:duration +
            REALM:microclim.sc:tempchange_abs.sc:duration +
            REALM:microclim.sc:tempave_metab.sc:duration +
            REALM:microclim.sc:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# mass:T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmmassAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            REALM:mass.sc:duration +
            REALM:mass.sc:tempchange_abs.sc:duration +
            REALM:mass.sc:tempave_metab.sc:duration +
            REALM:mass.sc:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM,
        control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}


# human:T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmhumanAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            REALM:human_bowler.sc:duration +
            REALM:human_bowler.sc:tempchange_abs.sc:duration +
            REALM:human_bowler.sc:tempave_metab.sc:duration +
            REALM:human_bowler.sc:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}


# endothermfrac:T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmendothermfracAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            REALM:endothermfrac.sc:duration +
            REALM:endothermfrac.sc:tempchange_abs.sc:duration +
            REALM:endothermfrac.sc:tempave_metab.sc:duration +
            REALM:endothermfrac.sc:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM)
    MATCHMOD <- TRUE
}



# taxamod2 (as RE):T:sdT:REALM:duration #########################
if (fitmod == 'modTsdTTRealmTaxamod2REAllJtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            REALM:duration +
            REALM:tempchange_abs.sc:duration +
            REALM:tempave_metab.sc:duration +
            REALM:tempave_metab.sc:tempchange_abs.sc:duration +
            (duration | STUDY_ID / rarefyID) +
            (REALM:tempchange_abs.sc:duration |
                 taxa_mod2) + # add taxamod2 as random effects
            (REALM:tempave_metab.sc:duration | taxa_mod2) +
            (
                REALM:tempave_metab.sc:tempchange_abs.sc:duration | taxa_mod2
            ),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM#,
        #control = glmmTMBControl(profile = TRUE)
    ) # add dispersion formula
    MATCHMOD <- TRUE
}




# T+dT:REALM:duration ANTAO-STYLE DATASET #########################

# try on Antao et al. 2020-style dataset
# use tempchange, not tempchange_abs

if (fitmod == 'modAntaoTdTTRealm10') {
    # trims out freshwater and >60 or <23.5 deg lat, like Antao et al. 2020
    i10 <-
        trends10[, REALM != 'Freshwater' &
                     abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                     complete.cases(Jtu.sc,
                                    tempchange.sc,
                                    REALM,
                                    tempave_metab.sc,
                                    durationlog.sc,
                                    nspp.sc)]
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    mod <-
        glmmTMB(
            Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc +
                tempave_metab.sc:tempchange.sc:duration.sc +
                REALM:duration.sc +
                tempchange.sc:REALM:duration.sc +
                tempave_metab.sc:REALM:duration.sc +
                (duration.sc | STUDY_ID / rarefyID),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc + REALM,
            control = glmmTMBControl(profile = TRUE)
        ) # add dispersion formula
    MATCHMOD <- TRUE
}

if (fitmod == 'modAntaoTdTTRealm10Jne') {
    i10 <-
        trends10[, REALM != 'Freshwater' &
                     abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                     complete.cases(Jne.sc,
                                    tempchange.sc,
                                    REALM,
                                    tempave_metab.sc,
                                    durationlog.sc,
                                    nspp.sc)]
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    mod <-
        glmmTMB(
            Jne.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc +
                tempave_metab.sc:tempchange.sc:duration.sc +
                REALM:duration.sc +
                tempchange.sc:REALM:duration.sc +
                tempave_metab.sc:REALM:duration.sc +
                (duration.sc | STUDY_ID / rarefyID),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc + REALM,
            control = glmmTMBControl(profile = TRUE)
        ) # add dispersion formula
    MATCHMOD <- TRUE
}

if (fitmod == 'modAntaoTdTTRealm10Horn') {
    i10 <-
        trends10[, REALM != 'Freshwater' &
                     abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                     complete.cases(Horn.sc,
                                    tempchange.sc,
                                    REALM,
                                    tempave_metab.sc,
                                    durationlog.sc,
                                    nspp.sc)]
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10), 'data points'))
    mod <-
        glmmTMB(
            Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc +
                tempave_metab.sc:tempchange.sc:duration.sc +
                REALM:duration.sc +
                tempchange.sc:REALM:duration.sc +
                tempave_metab.sc:REALM:duration.sc +
                (duration.sc | STUDY_ID / rarefyID),
            data = trends10[i10, ],
            family = beta_family(link = 'logit'),
            dispformula = ~ nspp.sc + REALM,
            control = glmmTMBControl(profile = TRUE)
        ) # add dispersion formula
    MATCHMOD <- TRUE
}



# Null models ############################

if (fitmod == 'modNull10Jtu') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ 1 +
            (duration.sc | STUDY_ID / rarefyID),
        data = trends10[i10Jtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ nspp.sc + REALM,
        # add dispersion formula
        control = glmmTMBControl(profile = TRUE)
    )
    MATCHMOD <- TRUE
}


if (fitmod == 'modNull10Horn') {
    if (MATCHMOD)
        stop('Model name matched more than one model!')
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(
        Horn.sc ~ 1 +
            (duration.sc | STUDY_ID / rarefyID),
        data = trends10[i10Horn, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ nspp.sc + REALM,
        # add dispersion formula
        control = glmmTMBControl(profile = TRUE)
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

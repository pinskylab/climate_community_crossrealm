#!/usr/bin/env Rscript

# Script to fit glmmTMB models
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_vs_temperature_GLMM.R modRFgauss > logs/turnover_vs_temperature_GLMMmodRFgauss.Rout &
# (this works if code is executable, e.g., chmod u+x turnover_vs_temperature_GLMM.R)
# (otherwise using nohup Rscript ...)

# Read command line arguments ############

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args)<1) stop("Have to specify a model to fit", call.=FALSE)
if (length(args)>1) stop("Have to specify only 1 model to fit", call.=FALSE)
fitmod <- args[1]
MATCHMOD <- FALSE # indicator to check if the argument matched a model name

# print basic info about the job ############################

print(paste('This is process #', Sys.getpid()))
print(Sys.time())



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB) # for ME models

# load data ############################

# Turnover and covariates assembled by turnover_vs_temperature_prep.Rmd
trends <- fread('output/turnover_w_covariates.csv.gz')


# Models ############################

## Models to compare variance structures ############################
#fixed <- 'Jtu.sc ~ tempchange_abs.sc*REALM + tempchange_abs.sc*tempave_metab.sc + tempchange_abs.sc*duration.sc'
fixed <- 'Jtu.sc ~ tempchange_abs.sc*REALM + tempchange_abs.sc*tempave_metab.sc + tempchange_abs.sc*duration.sc + tempchange_abs.sc*mass.sc + tempchange_abs.sc*endothermfrac.sc + tempchange_abs.sc*microclim.sc + tempchange_abs.sc*npp.sc + tempchange_abs.sc*human_bowler.sc:REALM2'
i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, endothermfrac.sc, microclim.sc, npp.sc, human_bowler.sc, nspp.sc)]

if(fitmod == 'modRFgauss'){
    print(paste(sum(i), 'data points'))
    modRFgauss <- glmmTMB(formula(fixed), data = trends[i,], family = gaussian(), control = glmmTMBControl(profile=TRUE))
    summary(modRFgauss)
    saveRDS(modRFgauss, file = 'temp/modRFgauss.rds')
    print('saved modRFgauss.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFbeta'){
    print(paste(sum(i), 'data points'))
    modRFbeta <- glmmTMB(formula(fixed), data = trends[i,], family = beta_family(link = 'logit'), control = glmmTMBControl(profile=TRUE)) # add beta errors
    summary(modRFbeta)
    saveRDS(modRFbeta, file = 'temp/modRFbeta.rds')
    print('saved modRFbeta.rds')
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFrID'){
    print(paste(sum(i), 'data points'))
    modRFrID <- glmmTMB(formula(paste(fixed, '+(1|rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random effects
    summary(modRFrID)
    saveRDS(modRFrID, file = 'temp/modRFrID.rds')
    print('saved modRDrID.rds')
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2lev'){
    print(paste(sum(i), 'data points'))
    modRF2lev <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                             control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRF2lev)
    saveRDS(modRF2lev, file = 'temp/modRF2lev.rds')
    print('saved modRF2lev.rds')
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFnestedRE'){
    print(paste(sum(i), 'data points'))
    modRFnestedRE <- glmmTMB(formula(paste(fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRFnestedRE)
    saveRDS(modRFnestedRE, file = 'temp/modRFnestedRE.rds')
    print('saved modRFnestedRE.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFslopeRE'){
    print(paste(sum(i), 'data points'))
    modRFslopeRE <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE)
    saveRDS(modRFslopeRE, file = 'temp/modRFslopeRE.rds')
    print('saved modRFslopeRE.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdisp2lev'){
    print(paste(sum(i), 'data points'))
    modRFdisp2lev <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc, 
                         control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdisp2lev)
    saveRDS(modRFdisp2lev, file = 'temp/modRFdisp2lev.rds')
    print('saved modRFdisp2lev.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdisp'){
    print(paste(sum(i), 'data points'))
    modRFdisp <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdisp)
    saveRDS(modRFdisp, file = 'temp/modRFdisp.rds')
    print('saved modRFdisp.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdisp2levOnlyint'){ # 2 level RE, no slope, + dispersion formula
    print(paste(sum(i), 'data points'))
    modRFdisp2levOnlyint <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                             dispformula = ~nspp.sc, 
                             control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp2levOnlyint)
    saveRDS(modRFdisp2levOnlyint, file = 'temp/modRFdisp2levOnlyint.rds')
    print('saved modRFdisp2levOnlyint.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdurslope2lev'){ # 2 level RE, slope vs. duration
    print(paste(sum(i), 'data points'))
    modRFdurslope2lev <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                             control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurslope2lev)
    saveRDS(modRFdurslope2lev, file = 'temp/modRFdurslope2lev.rds')
    print('saved modRFdurslope2lev.rds')
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurslope2levdisp'){ # 2 level RE, slope vs. duration, disp formula
    print(paste(sum(i), 'data points'))
    modRFdurslope2levdisp <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                             dispformula = ~nspp.sc, 
                             control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurslope2levdisp)
    saveRDS(modRFdurslope2levdisp, file = 'temp/modRFdurslope2levdisp.rds')
    print('saved modRFdurslope2levdisp.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdurslope2levdisprealm'){ # 2 level RE, slope vs. duration, disp formula by realm
    print(paste(sum(i), 'data points'))
    modRFdurslope2levdisprealm <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc+REALM, 
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurslope2levdisprealm)
    saveRDS(modRFdurslope2levdisprealm, file = 'temp/modRFdurslope2levdisprealm.rds')
    print('saved modRFdurslope2levdisprealm.rds')
    MATCHMOD <- TRUE
}


# temperature-only models (all years) ############################
if(fitmod == 'modonlyTchangeJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modonlyTchangeJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) + (duration.sc|STUDY_ID/rarefyID),
                             data = trends[i,], 
                             family = beta_family(link = 'logit'), 
                             dispformula = ~nspp.sc)
    summary(modonlyTchangeJtu)
    saveRDS(modonlyTchangeJtu, file = 'temp/modonlyTchangeJtu.rds')
    print('saved modonlyTchangeJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchangeJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modonlyTchangeJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) + (duration.sc|STUDY_ID/rarefyID),
                               data = trends[i,], 
                               family = beta_family(link = 'logit'), 
                               dispformula = ~nspp.sc)
    summary(modonlyTchangeJbeta)
    saveRDS(modonlyTchangeJbeta, file = 'temp/modonlyTchangeJbeta.rds')
    print('saved modonlyTchangeJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchangeHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modonlyTchangeHorn <- glmmTMB(Horn.sc ~ abs(tempchange) + (duration.sc|STUDY_ID/rarefyID),
                              data = trends[i,], 
                              family = beta_family(link = 'logit'), 
                              dispformula = ~nspp.sc)
    summary(modonlyTchangeHorn)
    saveRDS(modonlyTchangeHorn, file = 'temp/modonlyTchangeHorn.rds')
    print('saved modonlyTchangeHorn.rds')
    MATCHMOD <- TRUE
}


# temperature-only models (1 year) ############################
if(fitmod == 'modonlyTchange1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange, nspp.sc) & duration == 1]
    print(paste(sum(i), 'data points'))
    modonlyTchange1yrJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) + (1|STUDY_ID/rarefyID),
                                 data = trends[i,], 
                                 family = beta_family(link = 'logit'), 
                                 dispformula = ~nspp.sc)
    summary(modonlyTchange1yrJtu)
    saveRDS(modonlyTchange1yrJtu, file = 'temp/modonlyTchange1yrJtu.rds')
    print('saved modonlyTchange1yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchange1yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange, nspp.sc) & duration == 1]
    print(paste(sum(i), 'data points'))
    modonlyTchange1yrJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) + (1|STUDY_ID/rarefyID),
                                   data = trends[i,], 
                                   family = beta_family(link = 'logit'), 
                                   dispformula = ~nspp.sc)
    summary(modonlyTchange1yrJbeta)
    saveRDS(modonlyTchange1yrJbeta, file = 'temp/modonlyTchange1yrJbeta.rds')
    print('saved modonlyTchange1yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchange1yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange, nspp.sc) & duration == 1]
    print(paste(sum(i), 'data points'))
    modonlyTchange1yrHorn <- glmmTMB(Horn.sc ~ abs(tempchange) + (1|STUDY_ID/rarefyID),
                                  data = trends[i,], 
                                  family = beta_family(link = 'logit'), 
                                  dispformula = ~nspp.sc)
    summary(modonlyTchange1yrHorn)
    saveRDS(modonlyTchange1yrHorn, file = 'temp/modonlyTchange1yrHorn.rds')
    print('saved modonlyTchange1yrHorn.rds')
    MATCHMOD <- TRUE
}

# temperature-only models (10 year) ############################
if(fitmod == 'modonlyTchange10yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange, nspp.sc) & duration == 10]
    print(paste(sum(i), 'data points'))
    modonlyTchange10yrJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) + (1|STUDY_ID/rarefyID),
                                    data = trends[i,], 
                                    family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc)
    summary(modonlyTchange10yrJtu)
    saveRDS(modonlyTchange10yrJtu, file = 'temp/modonlyTchange10yrJtu.rds')
    print('saved modonlyTchange10yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchange10yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange, nspp.sc) & duration == 10]
    print(paste(sum(i), 'data points'))
    modonlyTchange10yrJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) + (1|STUDY_ID/rarefyID),
                                      data = trends[i,], 
                                      family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc)
    summary(modonlyTchange10yrJbeta)
    saveRDS(modonlyTchange10yrJbeta, file = 'temp/modonlyTchange10yrJbeta.rds')
    print('saved modonlyTchange10yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchange10yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange, nspp.sc) & duration == 10]
    print(paste(sum(i), 'data points'))
    modonlyTchange10yrHorn <- glmmTMB(Horn.sc ~ abs(tempchange) + (1|STUDY_ID/rarefyID),
                                     data = trends[i,], 
                                     family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc)
    summary(modonlyTchange10yrHorn)
    saveRDS(modonlyTchange10yrHorn, file = 'temp/modonlyTchange10yrHorn.rds')
    print('saved modonlyTchange10yrHorn.rds')
    MATCHMOD <- TRUE
}

# temperature-duration models ############################
if(fitmod == 'modTDJtu'){
    i <- trends[, complete.cases(Jtu.sc, REALM, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modTDJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * duration + (duration.sc|STUDY_ID/rarefyID),
                    data = trends[i,], 
                    family = beta_family(link = 'logit'), 
                    dispformula = ~nspp.sc)
    summary(modTDJtu)
    saveRDS(modTDJtu, file = 'temp/modTDJtu.rds')
    print('saved modTDJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTDJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, REALM, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modTDJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange)*duration + (duration.sc|STUDY_ID/rarefyID),
                      data = trends[i,], 
                      family = beta_family(link = 'logit'), 
                      dispformula = ~nspp.sc)
    summary(modTDJbeta)
    saveRDS(modTDJbeta, file = 'temp/modTDJbeta.rds')
    print('saved modTDJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTDHorn'){
    i <- trends[, complete.cases(Horn.sc, REALM, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modTDHorn <- glmmTMB(Horn.sc ~ abs(tempchange)*duration + (duration.sc|STUDY_ID/rarefyID),
                     data = trends[i,], 
                     family = beta_family(link = 'logit'), 
                     dispformula = ~nspp.sc)
    summary(modTDHorn)
    saveRDS(modTDHorn, file = 'temp/modTDHorn.rds')
    print('saved modTDHorn.rds')
    MATCHMOD <- TRUE
}

# temperature-duration models by REALM ############################
if(fitmod == 'modTDrealmJtu'){
    i <- trends[, complete.cases(Jtu.sc, REALM, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modTDrealmJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * duration + abs(tempchange) * REALM + (duration.sc|STUDY_ID/rarefyID),
                         data = trends[i,], 
                         family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc+REALM)
    summary(modTDrealmJtu)
    saveRDS(modTDrealmJtu, file = 'temp/modTDrealmJtu.rds')
    print('saved modTDrealmJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTDrealmJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, REALM, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modTDrealmJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange)*duration + abs(tempchange) * REALM + (duration.sc|STUDY_ID/rarefyID),
                           data = trends[i,], 
                           family = beta_family(link = 'logit'), 
                           dispformula = ~nspp.sc+REALM)
    summary(modTDrealmJbeta)
    saveRDS(modTDrealmJbeta, file = 'temp/modTDrealmJbeta.rds')
    print('saved modTDrealmJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTDrealmHorn'){
    i <- trends[, complete.cases(Horn.sc, REALM, tempchange, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modTDrealmHorn <- glmmTMB(Horn.sc ~ abs(tempchange)*duration + abs(tempchange) * REALM + (duration.sc|STUDY_ID/rarefyID),
                          data = trends[i,], 
                          family = beta_family(link = 'logit'), 
                          dispformula = ~nspp.sc+REALM)
    summary(modTDrealmHorn)
    saveRDS(modTDrealmHorn, file = 'temp/modTDrealmHorn.rds')
    print('saved modTDrealmHorn.rds')
    MATCHMOD <- TRUE
}

# temperature-only models by REALM (1 year) ############################
if(fitmod == 'modTrealm1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange, REALM, nspp.sc) & duration == 1]
    print(paste(sum(i), 'data points'))
    modTrealm1yrJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * REALM + (1|STUDY_ID/rarefyID),
                                    data = trends[i,], 
                                    family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc+REALM)
    summary(modTrealm1yrJtu)
    saveRDS(modTrealm1yrJtu, file = 'temp/modTrealm1yrJtu.rds')
    print('saved modTrealm1yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTrealm1yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange, REALM, nspp.sc) & duration == 1]
    print(paste(sum(i), 'data points'))
    modTrealm1yrJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) * REALM + (1|STUDY_ID/rarefyID),
                                      data = trends[i,], 
                                      family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc+REALM)
    summary(modTrealm1yrJbeta)
    saveRDS(modTrealm1yrJbeta, file = 'temp/modTrealm1yrJbeta.rds')
    print('saved modTrealm1yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTrealm1yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange, REALM, nspp.sc) & duration == 1]
    print(paste(sum(i), 'data points'))
    modTrealm1yrHorn <- glmmTMB(Horn.sc ~ abs(tempchange) * REALM + (1|STUDY_ID/rarefyID),
                                     data = trends[i,], 
                                     family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc+REALM)
    summary(modTrealm1yrHorn)
    saveRDS(modTrealm1yrHorn, file = 'temp/modTrealm1yrHorn.rds')
    print('saved modTrealm1yrHorn.rds')
    MATCHMOD <- TRUE
}


# temperature-only models by REALM (10 year) ############################
if(fitmod == 'modTrealm10yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange, REALM, nspp.sc) & duration == 10]
    print(paste(sum(i), 'data points'))
    modTrealm10yrJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * REALM + (1|STUDY_ID/rarefyID),
                               data = trends[i,], 
                               family = beta_family(link = 'logit'), 
                               dispformula = ~nspp.sc+REALM)
    summary(modTrealm10yrJtu)
    saveRDS(modTrealm10yrJtu, file = 'temp/modTrealm10yrJtu.rds')
    print('saved modTrealm10yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTrealm10yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange, REALM, nspp.sc) & duration == 10]
    print(paste(sum(i), 'data points'))
    modTrealm10yrJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) * REALM + (1|STUDY_ID/rarefyID),
                                 data = trends[i,], 
                                 family = beta_family(link = 'logit'), 
                                 dispformula = ~nspp.sc+REALM)
    summary(modTrealm10yrJbeta)
    saveRDS(modTrealm10yrJbeta, file = 'temp/modTrealm10yrJbeta.rds')
    print('saved modTrealm10yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTrealm10yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange, REALM, nspp.sc) & duration == 10]
    print(paste(sum(i), 'data points'))
    modTrealm10yrHorn <- glmmTMB(Horn.sc ~ abs(tempchange) * REALM + (1|STUDY_ID/rarefyID),
                                data = trends[i,], 
                                family = beta_family(link = 'logit'), 
                                dispformula = ~nspp.sc+REALM)
    summary(modTrealm10yrHorn)
    saveRDS(modTrealm10yrHorn, file = 'temp/modTrealm10yrHorn.rds')
    print('saved modTrealm10yrHorn.rds')
    MATCHMOD <- TRUE
}


# full models with main effects ############################

if(fitmod == 'modMainTdTSeMiNPNspMaHuJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modMainTdTSeMiNPNspMaHuJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc +
                                        REALM + 
                                        tempave_metab.sc + 
                                        duration.sc +
                                        seas.sc +
                                        microclim.sc +
                                        npp.sc +
                                        nspp.sc +
                                        mass.sc +
                                        human_bowler.sc:REALM2 +
                                        (duration.sc|STUDY_ID/rarefyID), 
                                    data = trends[i,], 
                                    family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc+REALM,
                                    control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHuJtu)
    outfile = 'temp/modMainTdTSeMiNPNspMaHuJtu.rds'
    saveRDS(modMainTdTSeMiNPNspMaHuJtu, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainTdTSeMiNPNspMaHu1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modMainTdTSeMiNPNspMaHu1yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc +
                                                 REALM + 
                                                 tempave_metab.sc + 
                                                 seas.sc +
                                                 microclim.sc +
                                                 npp.sc +
                                                 nspp.sc +
                                                 mass.sc +
                                                 human_bowler.sc:REALM2 +
                                                 (1|STUDY_ID/rarefyID), 
                                             data = trends[i,], 
                                             family = beta_family(link = 'logit'), 
                                             dispformula = ~nspp.sc+REALM,
                                             control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHu1yrJtu)
    outfile = 'temp/modMainTdTSeMiNPNspMaHu1yrJtu.rds'
    saveRDS(modMainTdTSeMiNPNspMaHu1yrJtu, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainTdTSeMiNPNspMaHu5yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 5]
    print(paste(sum(i), 'data points'))
    
    modMainTdTSeMiNPNspMaHu5yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc +
                                                 REALM + 
                                                 tempave_metab.sc + 
                                                 seas.sc +
                                                 microclim.sc +
                                                 npp.sc +
                                                 nspp.sc +
                                                 mass.sc +
                                                 human_bowler.sc:REALM2 +
                                                 (1|STUDY_ID/rarefyID), 
                                             data = trends[i,], 
                                             family = beta_family(link = 'logit'), 
                                             dispformula = ~nspp.sc+REALM, 
                                             control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHu5yrJtu)
    outfile = 'temp/modMainTdTSeMiNPNspMaHu5yrJtu.rds'
    saveRDS(modMainTdTSeMiNPNspMaHu5yrJtu, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainTdTSeMiNPNspMaHu10yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    print(paste(sum(i), 'data points'))
    
    modMainTdTSeMiNPNspMaHu10yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc +
                                                  REALM + 
                                                  tempave_metab.sc + 
                                                  seas.sc +
                                                  microclim.sc +
                                                  npp.sc +
                                                  nspp.sc +
                                                  mass.sc +
                                                  human_bowler.sc:REALM2 +
                                                  (1|STUDY_ID/rarefyID), 
                                              data = trends[i,], 
                                              family = beta_family(link = 'logit'), 
                                              dispformula = ~nspp.sc+REALM) 
    #                                        control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHu10yrJtu)
    outfile = 'temp/modMainTdTSeMiNPNspMaHu10yrJtu.rds'
    saveRDS(modMainTdTSeMiNPNspMaHu10yrJtu, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}


if(fitmod == 'modMainTdTSeMiNPNspMaHuJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modMainTdTSeMiNPNspMaHuJbeta <- glmmTMB(Jbeta.sc ~ tempchange_abs.sc +
                                                REALM + 
                                                tempave_metab.sc + 
                                                duration.sc +
                                                seas.sc +
                                                microclim.sc +
                                                npp.sc +
                                                nspp.sc +
                                                mass.sc +
                                                human_bowler.sc:REALM2 +
                                                (duration.sc|STUDY_ID/rarefyID), 
                                            data = trends[i,], 
                                            family = beta_family(link = 'logit'), 
                                            dispformula = ~nspp.sc+REALM,  # add dispersion formula
                                            control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHuJbeta)
    outfile = 'temp/modMainTdTSeMiNPNspMaHuJbeta.rds'
    saveRDS(modMainTdTSeMiNPNspMaHuJbeta, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainTdTSeMiNPNspMaHu1yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modMainTdTSeMiNPNspMaHu1yrJbeta <- glmmTMB(Jbeta.sc ~ tempchange_abs.sc +
                                                   REALM + 
                                                   tempave_metab.sc + 
                                                   seas.sc +
                                                   microclim.sc +
                                                   npp.sc +
                                                   nspp.sc +
                                                   mass.sc +
                                                   human_bowler.sc:REALM2 +
                                                   (1|STUDY_ID/rarefyID), 
                                               data = trends[i,], 
                                               family = beta_family(link = 'logit'), 
                                               dispformula = ~nspp.sc+REALM,
                                               control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHu1yrJbeta)
    outfile = 'temp/modMainTdTSeMiNPNspMaHu1yrJbeta.rds'
    saveRDS(modMainTdTSeMiNPNspMaHu1yrJbeta, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainTdTSeMiNPNspMaHu10yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    print(paste(sum(i), 'data points'))
    
    modMainTdTSeMiNPNspMaHu5yrJbeta <- glmmTMB(Jbeta.sc ~ tempchange_abs.sc +
                                                   REALM + 
                                                   tempave_metab.sc + 
                                                   seas.sc +
                                                   microclim.sc +
                                                   npp.sc +
                                                   nspp.sc +
                                                   mass.sc +
                                                   human_bowler.sc:REALM2 +
                                                   (1|STUDY_ID/rarefyID), 
                                               data = trends[i,], 
                                               family = beta_family(link = 'logit'), 
                                               dispformula = ~nspp.sc+REALM, 
                                               control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHu5yrJbeta)
    outfile = 'temp/modMainTdTSeMiNPNspMaHu5yrJbeta.rds'
    saveRDS(modMainTdTSeMiNPNspMaHu5yrJbeta, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}


if(fitmod == 'modMainTdTSeMiNPNspMaHuHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modMainTdTSeMiNPNspMaHuHorn <- glmmTMB(Horn.sc ~ tempchange_abs.sc +
                                               REALM + 
                                               tempave_metab.sc + 
                                               duration.sc +
                                               seas.sc +
                                               microclim.sc +
                                               npp.sc +
                                               nspp.sc +
                                               mass.sc +
                                               human_bowler.sc:REALM2 +
                                               (duration.sc|STUDY_ID/rarefyID), 
                                           data = trends[i,], 
                                           family = beta_family(link = 'logit'), 
                                           dispformula = ~nspp.sc+REALM,  # add dispersion formula
                                           control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHuHorn)
    outfile = 'temp/modMainTdTSeMiNPNspMaHuHorn.rds'
    saveRDS(modMainTdTSeMiNPNspMaHuHorn, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainTdTSeMiNPNspMaHu1yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modMainTdTSeMiNPNspMaHu1yrHorn <- glmmTMB(Horn.sc ~ tempchange_abs.sc +
                                                  REALM + 
                                                  tempave_metab.sc + 
                                                  seas.sc +
                                                  microclim.sc +
                                                  npp.sc +
                                                  nspp.sc +
                                                  mass.sc +
                                                  human_bowler.sc:REALM2 +
                                                  (1|STUDY_ID/rarefyID), 
                                              data = trends[i,], 
                                              family = beta_family(link = 'logit'), 
                                              dispformula = ~nspp.sc+REALM)
    #                                        control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHu1yrHorn)
    outfile = 'temp/modMainTdTSeMiNPNspMaHu1yrHorn.rds'
    saveRDS(modMainTdTSeMiNPNspMaHu1yrHorn, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainTdTSeMiNPNspMaHu10yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    print(paste(sum(i), 'data points'))
    
    modMainTdTSeMiNPNspMaHu10yrHorn <- glmmTMB(Horn.sc ~ tempchange_abs.sc +
                                                   REALM + 
                                                   tempave_metab.sc + 
                                                   seas.sc +
                                                   microclim.sc +
                                                   npp.sc +
                                                   nspp.sc +
                                                   mass.sc +
                                                   human_bowler.sc:REALM2 +
                                                   (1|STUDY_ID/rarefyID), 
                                               data = trends[i,], 
                                               family = beta_family(link = 'logit'), 
                                               dispformula = ~nspp.sc+REALM) 
    #                                        control = glmmTMBControl(profile=TRUE))
    summary(modMainTdTSeMiNPNspMaHu10yrHorn)
    outfile = 'temp/modMainTdTSeMiNPNspMaHu10yrHorn.rds'
    saveRDS(modMainTdTSeMiNPNspMaHu10yrHorn, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

# full models with main effects averaged by site ############################
# only for models of a fixed temporal duration

if(fitmod == 'modMainAveMaEnMiNPHu1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    trendsave <- trends[i, .(Jtu.sc = mean(Jtu.sc),
                             tempchange_abs.sc = mean(tempchange_abs.sc),
                             REALM = unique(REALM),
                             REALM2 = unique(REALM2),
                             tempave_metab.sc = mean(tempave_metab.sc),
                             mass.sc = mean(mass.sc),
                             endothermfrac.sc = mean(endothermfrac.sc),
                             microclim.sc = mean(microclim.sc),
                             npp.sc = mean(npp.sc),
                             human_bowler.sc = mean(human_bowler.sc),
                             nspp.sc = mean(nspp.sc),
                             STUDY_ID = unique(STUDY_ID)), by = rarefyID]
    
    print(paste(nrow(trendsave), 'data points'))
    
    modMainAveMaEnMiNPHu1yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc +
                                           REALM + 
                                           tempave_metab.sc + 
                                           mass.sc +
                                           endothermfrac.sc +
                                           microclim.sc +
                                           npp.sc +
                                           human_bowler.sc:REALM2 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trendsave, 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM,
                                       control = glmmTMBControl(profile=TRUE))
    summary(modMainAveMaEnMiNPHu1yrJtu)
    outfile = 'temp/modMainAveMaEnMiNPHu1yrJtu.rds'
    saveRDS(modMainAveMaEnMiNPHu1yrJtu, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainAveMaEnMiNPHu5yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 5]
    trendsave <- trends[i, .(Jtu.sc = mean(Jtu.sc),
                             tempchange_abs.sc = mean(tempchange_abs.sc),
                             REALM = unique(REALM),
                             REALM2 = unique(REALM2),
                             tempave_metab.sc = mean(tempave_metab.sc),
                             mass.sc = mean(mass.sc),
                             endothermfrac.sc = mean(endothermfrac.sc),
                             microclim.sc = mean(microclim.sc),
                             npp.sc = mean(npp.sc),
                             human_bowler.sc = mean(human_bowler.sc),
                             nspp.sc = mean(nspp.sc),
                             STUDY_ID = unique(STUDY_ID)), by = rarefyID]
    
    print(paste(nrow(trendsave), 'data points'))
    
    modMainAveMaEnMiNPHu5yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc +
                                           REALM + 
                                           tempave_metab.sc + 
                                           mass.sc +
                                           endothermfrac.sc +
                                           microclim.sc +
                                           npp.sc +
                                           human_bowler.sc:REALM2 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trendsave, 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM, 
                                       control = glmmTMBControl(profile=TRUE))
    summary(modMainAveMaEnMiNPHu5yrJtu)
    outfile = 'temp/modMainAveMaEnMiNPHu5yrJtu.rds'
    saveRDS(modMainAveMaEnMiNPHu5yrJtu, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainAveMaEnMiNPHu10yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    trendsave <- trends[i, .(Jtu.sc = mean(Jtu.sc),
                             tempchange_abs.sc = mean(tempchange_abs.sc),
                             REALM = unique(REALM),
                             REALM2 = unique(REALM2),
                             tempave_metab.sc = mean(tempave_metab.sc),
                             mass.sc = mean(mass.sc),
                             endothermfrac.sc = mean(endothermfrac.sc),
                             microclim.sc = mean(microclim.sc),
                             npp.sc = mean(npp.sc),
                             human_bowler.sc = mean(human_bowler.sc),
                             nspp.sc = mean(nspp.sc),
                             STUDY_ID = unique(STUDY_ID)), by = rarefyID]
    
    print(paste(nrow(trendsave), 'data points'))
    
    modMainAveMaEnMiNPHu10yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc +
                                            REALM + 
                                            tempave_metab.sc + 
                                            mass.sc +
                                            endothermfrac.sc +
                                            microclim.sc +
                                            npp.sc +
                                            human_bowler.sc:REALM2 +
                                            (1|STUDY_ID/rarefyID), 
                                        data = trendsave, 
                                        family = beta_family(link = 'logit'), 
                                        dispformula = ~nspp.sc+REALM) 
    #                                        control = glmmTMBControl(profile=TRUE))
    summary(modMainAveMaEnMiNPHu10yrJtu)
    outfile = 'temp/modMainAveMaEnMiNPHu10yrJtu.rds'
    saveRDS(modMainAveMaEnMiNPHu10yrJtu, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}


if(fitmod == 'modMainAveMaEnMiNPHu1yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    trendsave <- trends[i, .(Jbeta.sc = mean(Jbeta.sc),
                             tempchange_abs.sc = mean(tempchange_abs.sc),
                             REALM = unique(REALM),
                             REALM2 = unique(REALM2),
                             tempave_metab.sc = mean(tempave_metab.sc),
                             mass.sc = mean(mass.sc),
                             endothermfrac.sc = mean(endothermfrac.sc),
                             microclim.sc = mean(microclim.sc),
                             npp.sc = mean(npp.sc),
                             human_bowler.sc = mean(human_bowler.sc),
                             nspp.sc = mean(nspp.sc),
                             STUDY_ID = unique(STUDY_ID)), by = rarefyID]
    
    print(paste(nrow(trendsave), 'data points'))
    
    modMainAveMaEnMiNPHu1yrJbeta <- glmmTMB(Jbeta.sc ~ tempchange_abs.sc +
                                             REALM + 
                                             tempave_metab.sc + 
                                             mass.sc +
                                             endothermfrac.sc +
                                             microclim.sc +
                                             npp.sc +
                                             human_bowler.sc:REALM2 +
                                             (1|STUDY_ID/rarefyID), 
                                         data = trendsave, 
                                         family = beta_family(link = 'logit'), 
                                         dispformula = ~nspp.sc+REALM,
                                         control = glmmTMBControl(profile=TRUE))
    summary(modMainAveMaEnMiNPHu1yrJbeta)
    outfile = 'temp/modMainAveMaEnMiNPHu1yrJbeta.rds'
    saveRDS(modMainAveMaEnMiNPHu1yrJbeta, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainAveMaEnMiNPHu10yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    trendsave <- trends[i, .(Jbeta.sc = mean(Jbeta.sc),
                             tempchange_abs.sc = mean(tempchange_abs.sc),
                             REALM = unique(REALM),
                             REALM2 = unique(REALM2),
                             tempave_metab.sc = mean(tempave_metab.sc),
                             mass.sc = mean(mass.sc),
                             endothermfrac.sc = mean(endothermfrac.sc),
                             microclim.sc = mean(microclim.sc),
                             npp.sc = mean(npp.sc),
                             human_bowler.sc = mean(human_bowler.sc),
                             nspp.sc = mean(nspp.sc),
                             STUDY_ID = unique(STUDY_ID)), by = rarefyID]
    
    print(paste(nrow(trendsave), 'data points'))
    
    modMainAveMaEnMiNPHu5yrJbeta <- glmmTMB(Jbeta.sc ~ tempchange_abs.sc +
                                             REALM + 
                                             tempave_metab.sc + 
                                             mass.sc +
                                             endothermfrac.sc +
                                             microclim.sc +
                                             npp.sc +
                                             human_bowler.sc:REALM2 +
                                             (1|STUDY_ID/rarefyID), 
                                         data = trendsave, 
                                         family = beta_family(link = 'logit'), 
                                         dispformula = ~nspp.sc+REALM, 
                                         control = glmmTMBControl(profile=TRUE))
    summary(modMainAveMaEnMiNPHu5yrJbeta)
    outfile = 'temp/modMainAveMaEnMiNPHu5yrJbeta.rds'
    saveRDS(modMainAveMaEnMiNPHu5yrJbeta, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}


if(fitmod == 'modMainAveMaEnMiNPHu1yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    trendsave <- trends[i, .(Horn.sc = mean(Horn.sc),
                             tempchange_abs.sc = mean(tempchange_abs.sc),
                             REALM = unique(REALM),
                             REALM2 = unique(REALM2),
                             tempave_metab.sc = mean(tempave_metab.sc),
                             mass.sc = mean(mass.sc),
                             endothermfrac.sc = mean(endothermfrac.sc),
                             microclim.sc = mean(microclim.sc),
                             npp.sc = mean(npp.sc),
                             human_bowler.sc = mean(human_bowler.sc),
                             nspp.sc = mean(nspp.sc),
                             STUDY_ID = unique(STUDY_ID)), by = rarefyID]
    
    print(paste(nrow(trendsave), 'data points'))
    
    modMainAveMaEnMiNPHu1yrHorn <- glmmTMB(Horn.sc ~ tempchange_abs.sc +
                                            REALM + 
                                            tempave_metab.sc + 
                                            mass.sc +
                                            endothermfrac.sc +
                                            microclim.sc +
                                            npp.sc +
                                            human_bowler.sc:REALM2 +
                                            (1|STUDY_ID/rarefyID), 
                                        data = trendsave, 
                                        family = beta_family(link = 'logit'), 
                                        dispformula = ~nspp.sc+REALM)
    #                                        control = glmmTMBControl(profile=TRUE))
    summary(modMainAveMaEnMiNPHu1yrHorn)
    outfile = 'temp/modMainAveMaEnMiNPHu1yrHorn.rds'
    saveRDS(modMainAveMaEnMiNPHu1yrHorn, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}

if(fitmod == 'modMainAveMaEnMiNPHu10yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    trendsave <- trends[i, .(Horn.sc = mean(Horn.sc),
                             tempchange_abs.sc = mean(tempchange_abs.sc),
                             REALM = unique(REALM),
                             REALM2 = unique(REALM2),
                             tempave_metab.sc = mean(tempave_metab.sc),
                             mass.sc = mean(mass.sc),
                             endothermfrac.sc = mean(endothermfrac.sc),
                             microclim.sc = mean(microclim.sc),
                             npp.sc = mean(npp.sc),
                             human_bowler.sc = mean(human_bowler.sc),
                             nspp.sc = mean(nspp.sc),
                             STUDY_ID = unique(STUDY_ID)), by = rarefyID]
    
    print(paste(nrow(trendsave), 'data points'))
    
    modMainAveMaEnMiNPHu10yrHorn <- glmmTMB(Horn.sc ~ tempchange_abs.sc +
                                             REALM + 
                                             tempave_metab.sc + 
                                             mass.sc +
                                             endothermfrac.sc +
                                             microclim.sc +
                                             npp.sc +
                                             human_bowler.sc:REALM2 +
                                             (1|STUDY_ID/rarefyID), 
                                         data = trendsave, 
                                         family = beta_family(link = 'logit'), 
                                         dispformula = ~nspp.sc+REALM) 
    #                                        control = glmmTMBControl(profile=TRUE))
    summary(modMainAveMaEnMiNPHu10yrHorn)
    outfile = 'temp/modMainAveMaEnMiNPHu10yrHorn.rds'
    saveRDS(modMainAveMaEnMiNPHu10yrHorn, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}



# full models with temperature change interaction ############################


if(fitmod == 'modFullmass'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    
    modFullmass <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                            tempchange_abs.sc*tempave_metab.sc + 
                            tempchange_abs.sc*duration.sc +
                            tempchange_abs.sc*mass.sc +
                            (duration.sc|STUDY_ID/rarefyID), 
                        data = trends[i,], 
                        family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modFullmass)
    saveRDS(modFullmass, file = 'temp/modFullmass.rds')
    print('saved modFullmass.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullendo'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, endothermfrac.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    
    modFullendo <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                               tempchange_abs.sc*tempave_metab.sc + 
                               tempchange_abs.sc*duration.sc +
                               tempchange_abs.sc*endothermfrac.sc +
                               (duration.sc|STUDY_ID/rarefyID), 
                           data = trends[i,], 
                           family = beta_family(link = 'logit'), 
                           dispformula = ~nspp.sc, 
                           control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modFullendo)
    saveRDS(modFullendo, file = 'temp/modFullendo.rds')
    print('saved modFullendo.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHuJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modFullMaEnMiNPHuJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                               tempchange_abs.sc*tempave_metab.sc + 
                               tempchange_abs.sc*duration.sc +
                               tempchange_abs.sc*mass.sc +
                               tempchange_abs.sc*endothermfrac.sc +
                               tempchange_abs.sc*microclim.sc +
                               tempchange_abs.sc*npp.sc +
                               tempchange_abs.sc*human_bowler.sc:REALM2 +
                               (duration.sc|STUDY_ID/rarefyID), 
                           data = trends[i,], 
                           family = beta_family(link = 'logit'), 
                           dispformula = ~nspp.sc+REALM,  # add dispersion formula
                           control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHuJtu)
    saveRDS(modFullMaEnMiNPHuJtu, file = 'temp/modFullMaEnMiNPHuJtu.rds')
    print('saved modFullMaEnMiNPHuJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHu1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modFullMaEnMiNPHu1yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                                     tempchange_abs.sc*tempave_metab.sc + 
                                     tempchange_abs.sc*mass.sc +
                                     tempchange_abs.sc*endothermfrac.sc +
                                     tempchange_abs.sc*microclim.sc +
                                     tempchange_abs.sc*npp.sc +
                                     tempchange_abs.sc*human_bowler.sc:REALM2 +
                                     (1|STUDY_ID/rarefyID), 
                                 data = trends[i,], 
                                 family = beta_family(link = 'logit'), 
                                 dispformula = ~nspp.sc+REALM,
                                 control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHu1yrJtu)
    saveRDS(modFullMaEnMiNPHu1yrJtu, file = 'temp/modFullMaEnMiNPHu1yrJtu.rds')
    print('saved modFullMaEnMiNPHu1yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHu5yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 5]
    print(paste(sum(i), 'data points'))
    
    modFullMaEnMiNPHu5yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                                            tempchange_abs.sc*tempave_metab.sc + 
                                            tempchange_abs.sc*mass.sc +
                                            tempchange_abs.sc*endothermfrac.sc +
                                            tempchange_abs.sc*microclim.sc +
                                            tempchange_abs.sc*npp.sc +
                                            tempchange_abs.sc*human_bowler.sc:REALM2 +
                                            (1|STUDY_ID/rarefyID), 
                                        data = trends[i,], 
                                        family = beta_family(link = 'logit'), 
                                        dispformula = ~nspp.sc+REALM, 
                                        control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHu5yrJtu)
    saveRDS(modFullMaEnMiNPHu5yrJtu, file = 'temp/modFullMaEnMiNPHu5yrJtu.rds')
    print('saved modFullMaEnMiNPHu5yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHu10yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    print(paste(sum(i), 'data points'))
    
    modFullMaEnMiNPHu10yrJtu <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                                           tempchange_abs.sc*tempave_metab.sc + 
                                           tempchange_abs.sc*mass.sc +
                                           tempchange_abs.sc*endothermfrac.sc +
                                           tempchange_abs.sc*microclim.sc +
                                           tempchange_abs.sc*npp.sc +
                                           tempchange_abs.sc*human_bowler.sc:REALM2 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trends[i,], 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM, 
                                       control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHu10yrJtu)
    saveRDS(modFullMaEnMiNPHu10yrJtu, file = 'temp/modFullMaEnMiNPHu10yrJtu.rds')
    print('saved modFullMaEnMiNPHu10yrJtu.rds')
    MATCHMOD <- TRUE
}


if(fitmod == 'modFullMaEnMiNPHuJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modFullMaEnMiNPHuJbeta <- glmmTMB(Jbeta.sc ~ tempchange_abs.sc*REALM + 
                                        tempchange_abs.sc*tempave_metab.sc + 
                                        tempchange_abs.sc*duration.sc +
                                        tempchange_abs.sc*mass.sc +
                                        tempchange_abs.sc*endothermfrac.sc +
                                        tempchange_abs.sc*microclim.sc +
                                        tempchange_abs.sc*npp.sc +
                                        tempchange_abs.sc*human_bowler.sc:REALM2 +
                                        (duration.sc|STUDY_ID/rarefyID), 
                                    data = trends[i,], 
                                    family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc+REALM,  # add dispersion formula
                                    control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHuJbeta)
    saveRDS(modFullMaEnMiNPHuJbeta, file = 'temp/modFullMaEnMiNPHuJbeta.rds')
    print('saved modFullMaEnMiNPHuJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHu1yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modFullMaEnMiNPHu1yrJbeta <- glmmTMB(Jbeta.sc ~ tempchange_abs.sc*REALM + 
                                           tempchange_abs.sc*tempave_metab.sc + 
                                           tempchange_abs.sc*mass.sc +
                                           tempchange_abs.sc*endothermfrac.sc +
                                           tempchange_abs.sc*microclim.sc +
                                           tempchange_abs.sc*npp.sc +
                                           tempchange_abs.sc*human_bowler.sc:REALM2 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trends[i,], 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM,
                                       control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHu1yrJbeta)
    saveRDS(modFullMaEnMiNPHu1yrJbeta, file = 'temp/modFullMaEnMiNPHu1yrJbeta.rds')
    print('saved modFullMaEnMiNPHu1yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHu5yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 5]
    print(paste(sum(i), 'data points'))
    
    modFullMaEnMiNPHu5yrJbeta <- glmmTMB(Jbeta.sc ~ tempchange_abs.sc*REALM + 
                                           tempchange_abs.sc*tempave_metab.sc + 
                                           tempchange_abs.sc*mass.sc +
                                           tempchange_abs.sc*endothermfrac.sc +
                                           tempchange_abs.sc*microclim.sc +
                                           tempchange_abs.sc*npp.sc +
                                           tempchange_abs.sc*human_bowler.sc:REALM2 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trends[i,], 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM, 
                                       control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHu5yrJbeta)
    saveRDS(modFullMaEnMiNPHu5yrJbeta, file = 'temp/modFullMaEnMiNPHu5yrJbeta.rds')
    print('saved modFullMaEnMiNPHu5yrJbeta.rds')
    MATCHMOD <- TRUE
}


if(fitmod == 'modFullMaEnMiNPHuHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modFullMaEnMiNPHuHorn <- glmmTMB(Horn.sc ~ tempchange_abs.sc*REALM + 
                                          tempchange_abs.sc*tempave_metab.sc + 
                                          tempchange_abs.sc*duration.sc +
                                          tempchange_abs.sc*mass.sc +
                                          tempchange_abs.sc*endothermfrac.sc +
                                          tempchange_abs.sc*microclim.sc +
                                          tempchange_abs.sc*npp.sc +
                                          tempchange_abs.sc*human_bowler.sc:REALM2 +
                                          (duration.sc|STUDY_ID/rarefyID), 
                                      data = trends[i,], 
                                      family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc+REALM,  # add dispersion formula
                                      control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHuHorn)
    saveRDS(modFullMaEnMiNPHuHorn, file = 'temp/modFullMaEnMiNPHuHorn.rds')
    print('saved modFullMaEnMiNPHuHorn.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHu1yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modFullMaEnMiNPHu1yrHorn <- glmmTMB(Horn.sc ~ tempchange_abs.sc*REALM + 
                                             tempchange_abs.sc*tempave_metab.sc + 
                                             tempchange_abs.sc*mass.sc +
                                             tempchange_abs.sc*endothermfrac.sc +
                                             tempchange_abs.sc*microclim.sc +
                                             tempchange_abs.sc*npp.sc +
                                             tempchange_abs.sc*human_bowler.sc:REALM2 +
                                             (1|STUDY_ID/rarefyID), 
                                         data = trends[i,], 
                                         family = beta_family(link = 'logit'), 
                                         dispformula = ~nspp.sc+REALM,
                                         control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHu1yrHorn)
    saveRDS(modFullMaEnMiNPHu1yrHorn, file = 'temp/modFullMaEnMiNPHu1yrHorn.rds')
    print('saved modFullMaEnMiNPHu1yrHorn.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHu5yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 5]
    print(paste(sum(i), 'data points'))
    
    modFullMaEnMiNPHu5yrHorn <- glmmTMB(Horn.sc ~ tempchange_abs.sc*REALM + 
                                             tempchange_abs.sc*tempave_metab.sc + 
                                             tempchange_abs.sc*mass.sc +
                                             tempchange_abs.sc*endothermfrac.sc +
                                             tempchange_abs.sc*microclim.sc +
                                             tempchange_abs.sc*npp.sc +
                                             tempchange_abs.sc*human_bowler.sc:REALM2 +
                                             (1|STUDY_ID/rarefyID), 
                                         data = trends[i,], 
                                         family = beta_family(link = 'logit'), 
                                         dispformula = ~nspp.sc+REALM, 
                                         control = glmmTMBControl(profile=TRUE))
    summary(modFullMaEnMiNPHu5yrHorn)
    saveRDS(modFullMaEnMiNPHu5yrHorn, file = 'temp/modFullMaEnMiNPHu5yrHorn.rds')
    print('saved modFullMaEnMiNPHu5yrHorn.rds')
    MATCHMOD <- TRUE
}

# full models with duration interactions ############################

if(fitmod == 'modDurIntTdTSeMiNPNspMaHu5yrJtu'){
    trends[, maxyr := max(c(year1, year2)), by = rarefyID]
    rowkeep <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) &
                    year1 > maxyr - 5 &
                    year2 > maxyr - 5] # trim to complete cases in the last 5 years
    rkeep <- trends[rowkeep, .(nyrs = length(unique(c(year1, year2)))), by = rarefyID][nyrs == 5, rarefyID] # find timeseries with exactly 5 years
    i <- trends[, rowkeep & rarefyID %in% rkeep]
    
    trends[i, tempchangeTS_abs := abs(median(tempchange/duration)), by = rarefyID] # Thiel-Sen estimator of the abs(slope): median of all pairwise differences
 
    print(paste(sum(i), 'data points'))
    modDurIntTdTSeMiNPNspMaHu5yrJtu <- glmmTMB(Jtu.sc ~ duration.sc + REALM + nspp.sc +
                                                   tempchangeTS_abs:duration.sc +
                                                   REALM : duration.sc + 
                                                   tempave_metab.sc : duration.sc + 
                                                   seas.sc : duration.sc+
                                                   microclim.sc : duration.sc+
                                                   npp.sc : duration.sc+
                                                   nspp.sc : duration.sc+
                                                   mass.sc : duration.sc +
                                                   human_bowler.sc:REALM2 : duration.sc +
                                                   (duration.sc|STUDY_ID/rarefyID), 
                                               data = trends[i,], 
                                               family = beta_family(link = 'logit'), 
                                               dispformula = ~nspp.sc+REALM,
                                               control = glmmTMBControl(profile=TRUE))
    summary(modDurIntTdTSeMiNPNspMaHu5yrJtu)
    outfile = 'temp/modDurIntTdTSeMiNPNspMaHu5yrJtu.rds'
    saveRDS(modDurIntTdTSeMiNPNspMaHu5yrJtu, file = outfile)
    print(paste('saved', outfile))
    MATCHMOD <- TRUE
}



# Null models ############################

if(fitmod == 'modNullMaEnMiNPHuJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modNullMaEnMiNPHuJtu <- glmmTMB(Jtu.sc ~ 1 +
                                        (duration.sc|STUDY_ID/rarefyID), 
                                    data = trends[i,], 
                                    family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc+REALM,  # add dispersion formula
                                    control = glmmTMBControl(profile=TRUE))
    summary(modNullMaEnMiNPHuJtu)
    saveRDS(modNullMaEnMiNPHuJtu, file = 'temp/modNullMaEnMiNPHuJtu.rds')
    print('saved modNullMaEnMiNPHuJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modNullMaEnMiNPHu1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modNullMaEnMiNPHu1yrJtu <- glmmTMB(Jtu.sc ~ 1 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trends[i,], 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM,
                                       control = glmmTMBControl(profile=TRUE))
    summary(modNullMaEnMiNPHu1yrJtu)
    saveRDS(modNullMaEnMiNPHu1yrJtu, file = 'temp/modNullMaEnMiNPHu1yrJtu.rds')
    print('saved modNullMaEnMiNPHu1yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modNullMaEnMiNPHu5yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc, nspp.sc) & 
                    duration == 5]
    print(paste(sum(i), 'data points'))
    
    modNullMaEnMiNPHu5yrJtu <- glmmTMB(Jtu.sc ~ 1 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trends[i,], 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM, 
                                       control = glmmTMBControl(profile=TRUE))
    summary(modNullMaEnMiNPHu5yrJtu)
    saveRDS(modNullMaEnMiNPHu5yrJtu, file = 'temp/modNullMaEnMiNPHu5yrJtu.rds')
    print('saved modNullMaEnMiNPHu5yrJtu.rds')
    MATCHMOD <- TRUE
}


if(fitmod == 'modNullTdTSeMiNPNspMaHuJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, 
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modNullTdTSeMiNPNspMaHuJtu <- glmmTMB(Jtu.sc ~ 1 +
                                        (duration.sc|STUDY_ID/rarefyID), 
                                    data = trends[i,], 
                                    family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc+REALM,  # add dispersion formula
                                    control = glmmTMBControl(profile=TRUE))
    summary(modNullTdTSeMiNPNspMaHuJtu)
    saveRDS(modNullTdTSeMiNPNspMaHuJtu, file = 'temp/modNullTdTSeMiNPNspMaHuJtu.rds')
    print('saved modNullTdTSeMiNPNspMaHuJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modNullTdTSeMiNPNspMaHu1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, 
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modNullTdTSeMiNPNspMaHu1yrJtu <- glmmTMB(Jtu.sc ~ 1 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trends[i,], 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM,
                                       control = glmmTMBControl(profile=TRUE))
    summary(modNullTdTSeMiNPNspMaHu1yrJtu)
    saveRDS(modNullTdTSeMiNPNspMaHu1yrJtu, file = 'temp/modNullTdTSeMiNPNspMaHu1yrJtu.rds')
    print('saved modNullTdTSeMiNPNspMaHu1yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modNullTdTSeMiNPNspMaHu10yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, 
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    print(paste(sum(i), 'data points'))
    
    modNullTdTSeMiNPNspMaHu10yrJtu <- glmmTMB(Jtu.sc ~ 1 +
                                           (1|STUDY_ID/rarefyID), 
                                       data = trends[i,], 
                                       family = beta_family(link = 'logit'), 
                                       dispformula = ~nspp.sc+REALM, 
                                       control = glmmTMBControl(profile=TRUE))
    summary(modNullTdTSeMiNPNspMaHu10yrJtu)
    saveRDS(modNullTdTSeMiNPNspMaHu10yrJtu, file = 'temp/modNullTdTSeMiNPNspMaHu10yrJtu.rds')
    print('saved modNullTdTSeMiNPNspMaHu10yrJtu.rds')
    MATCHMOD <- TRUE
}


if(fitmod == 'modNullTdTSeMiNPNspMaHuHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, 
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc)]
    print(paste(sum(i), 'data points'))
    modNullTdTSeMiNPNspMaHuHorn <- glmmTMB(Horn.sc ~ 1 +
                                              (duration.sc|STUDY_ID/rarefyID), 
                                          data = trends[i,], 
                                          family = beta_family(link = 'logit'), 
                                          dispformula = ~nspp.sc+REALM,  # add dispersion formula
                                          control = glmmTMBControl(profile=TRUE))
    summary(modNullTdTSeMiNPNspMaHuHorn)
    saveRDS(modNullTdTSeMiNPNspMaHuHorn, file = 'temp/modNullTdTSeMiNPNspMaHuHorn.rds')
    print('saved modNullTdTSeMiNPNspMaHuHorn.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modNullTdTSeMiNPNspMaHu1yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, 
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 1]
    print(paste(sum(i), 'data points'))
    
    modNullTdTSeMiNPNspMaHu1yrHorn <- glmmTMB(Horn.sc ~ 1 +
                                                 (1|STUDY_ID/rarefyID), 
                                             data = trends[i,], 
                                             family = beta_family(link = 'logit'), 
                                             dispformula = ~nspp.sc+REALM,
                                             control = glmmTMBControl(profile=TRUE))
    summary(modNullTdTSeMiNPNspMaHu1yrHorn)
    saveRDS(modNullTdTSeMiNPNspMaHu1yrHorn, file = 'temp/modNullTdTSeMiNPNspMaHu1yrHorn.rds')
    print('saved modNullTdTSeMiNPNspMaHu1yrHorn.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modNullTdTSeMiNPNspMaHu10yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, 
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) & 
                    duration == 10]
    print(paste(sum(i), 'data points'))
    
    modNullTdTSeMiNPNspMaHu10yrHorn <- glmmTMB(Horn.sc ~ 1 +
                                                  (1|STUDY_ID/rarefyID), 
                                              data = trends[i,], 
                                              family = beta_family(link = 'logit'), 
                                              dispformula = ~nspp.sc+REALM, 
                                              control = glmmTMBControl(profile=TRUE))
    summary(modNullTdTSeMiNPNspMaHu10yrHorn)
    saveRDS(modNullTdTSeMiNPNspMaHu10yrHorn, file = 'temp/modNullTdTSeMiNPNspMaHu10yrHorn.rds')
    print('saved modNullTdTSeMiNPNspMaHu10yrHorn.rds')
    MATCHMOD <- TRUE
}


# final check that something ran ############################
if(MATCHMOD == FALSE) stop("Model name did not match anything", call.=FALSE)

warnings()

print(paste('Ended', Sys.time(), sep =''))
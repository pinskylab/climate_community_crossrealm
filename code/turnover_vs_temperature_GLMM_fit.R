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
library(performance) # for R2

# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
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
    print(Sys.time())
    print(performance::r2(modRFgauss))
    
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFbeta'){
    print(paste(sum(i), 'data points'))
    modRFbeta <- glmmTMB(formula(fixed), data = trends[i,], family = beta_family(link = 'logit'), control = glmmTMBControl(profile=TRUE)) # add beta errors
    summary(modRFbeta)
    saveRDS(modRFbeta, file = 'temp/modRFbeta.rds')
    print('saved modRFbeta.rds')
    print(Sys.time())
    print(performance::r2(modRFbeta))
    
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFrID'){
    print(paste(sum(i), 'data points'))
    modRFrID <- glmmTMB(formula(paste(fixed, '+(1|rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random effects
    summary(modRFrID)
    saveRDS(modRFrID, file = 'temp/modRFrID.rds')
    print('saved modRDrID.rds')
    print(Sys.time())
    print(performance::r2(modRFrID))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2lev'){
    print(paste(sum(i), 'data points'))
    modRF2lev <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                             control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRF2lev)
    saveRDS(modRF2lev, file = 'temp/modRF2lev.rds')
    print('saved modRF2lev.rds')
    print(Sys.time())
    print(performance::r2(modRF2lev))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFnestedRE'){
    print(paste(sum(i), 'data points'))
    modRFnestedRE <- glmmTMB(formula(paste(fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRFnestedRE)
    saveRDS(modRFnestedRE, file = 'temp/modRFnestedRE.rds')
    print('saved modRFnestedRE.rds')
    print(Sys.time())
    print(performance::r2(modRFnestedRE))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFslopeRE'){
    print(paste(sum(i), 'data points'))
    modRFslopeRE <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE)
    saveRDS(modRFslopeRE, file = 'temp/modRFslopeRE.rds')
    print('saved modRFslopeRE.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE))
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
    print(Sys.time())
    print(performance::r2(modRFdisp2lev))
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
    print(Sys.time())
    print(performance::r2(modRFdisp))
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
    print(Sys.time())
    print(performance::r2(modRFdisp2levOnlyint))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdurslope2lev'){ # 2 level RE, slope vs. duration
    print(paste(sum(i), 'data points'))
    modRFdurslope2lev <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                             control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurslope2lev)
    saveRDS(modRFdurslope2lev, file = 'temp/modRFdurslope2lev.rds')
    print('saved modRFdurslope2lev.rds')
    print(Sys.time())
    print(performance::r2(modRFdurslope2lev))
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
    print(Sys.time())
    print(performance::r2(modRFdurslope2levdisp))
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
    print(Sys.time())
    print(performance::r2(modRFdurslope2levdisprealm))
    MATCHMOD <- TRUE
}


# full models with duration interactions ############################

if(fitmod == 'modDurIntTdTSeMiNPNspMaHu3yrJtu'){
    NYR = 3
    trends[, maxyr := max(c(year1, year2)), by = rarefyID]
    rowkeep <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration,
                                       seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) &
                          year1 > maxyr - NYR &
                          year2 > maxyr - NYR] # trim to complete cases in the last NYR years
    rkeep <- trends[rowkeep, .(nyrs = length(unique(c(year1, year2)))), by = rarefyID][nyrs == NYR, rarefyID] # find timeseries with exactly 5 years
    i <- trends[, rowkeep & rarefyID %in% rkeep]
    
    trends[i, tempchangeTS_abs := abs(median(tempchange/duration)), by = rarefyID] # Thiel-Sen estimator of the abs(slope): median of all pairwise differences
    
    print(paste(sum(i), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration + REALM + nspp.sc +
                                                   tempchangeTS_abs:duration +
                                                   REALM : duration + 
                                                   tempave_metab.sc : duration + 
                                                   seas.sc : duration+
                                                   microclim.sc : duration+
                                                   npp.sc : duration+
                                                   nspp.sc : duration+
                                                   mass.sc : duration +
                                                   human_bowler.sc:REALM2 : duration +
                                                   (duration|STUDY_ID/rarefyID), 
                                               data = trends[i,], 
                                               family = beta_family(link = 'logit'), 
                                               dispformula = ~nspp.sc+REALM,
                                               control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    outfile = paste0('temp/modDurIntTdTSeMiNPNspMaHu', NYR, 'yrJtu.rds')
    saveRDS(mod, file = outfile)
    print(paste('saved', outfile))
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modDurIntTdTSeMiNPNspMaHu5yrJtu'){
    trends[, maxyr := max(c(year1, year2)), by = rarefyID]
    rowkeep <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration,
                                 seas.sc, microclim.sc, npp.sc, mass.sc, human_bowler.sc, nspp.sc) &
                    year1 > maxyr - 5 &
                    year2 > maxyr - 5] # trim to complete cases in the last 5 years
    rkeep <- trends[rowkeep, .(nyrs = length(unique(c(year1, year2)))), by = rarefyID][nyrs == 5, rarefyID] # find timeseries with exactly 5 years
    i <- trends[, rowkeep & rarefyID %in% rkeep]
    
    trends[i, tempchangeTS_abs := abs(median(tempchange/duration)), by = rarefyID] # Thiel-Sen estimator of the abs(slope): median of all pairwise differences
 
    print(paste(sum(i), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration + REALM + nspp.sc +
                                                   tempchangeTS_abs:duration +
                                                   REALM : duration + 
                                                   tempave_metab.sc : duration + 
                                                   seas.sc : duration+
                                                   microclim.sc : duration+
                                                   npp.sc : duration+
                                                   nspp.sc : duration+
                                                   mass.sc : duration +
                                                   human_bowler.sc:REALM2 : duration +
                                                   (duration|STUDY_ID/rarefyID), 
                                               data = trends[i,], 
                                               family = beta_family(link = 'logit'), 
                                               dispformula = ~nspp.sc+REALM,
                                               control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    outfile = 'temp/modDurIntTdTSeMiNPNspMaHu5yrJtu.rds'
    saveRDS(mod, file = outfile)
    print(paste('saved', outfile))
    print(Sys.time())
    print(performance::r2(mod))
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
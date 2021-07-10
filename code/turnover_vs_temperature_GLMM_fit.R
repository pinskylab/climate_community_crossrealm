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
trends3 <- fread('output/turnover_w_covariates3.csv.gz')
trends5 <- fread('output/turnover_w_covariates5.csv.gz')
trends10 <- fread('output/turnover_w_covariates10.csv.gz')
trends20 <- fread('output/turnover_w_covariates20.csv.gz')


# Models ############################

## Choose dataset
i3Jtu <- trends3[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i5Jtu <- trends5[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i10Jtu <- trends10[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i20Jtu <- trends20[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

i3Jne <- trends3[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i5Jne <- trends5[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i10Jne <- trends10[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i20Jne <- trends20[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

i3Horn <- trends3[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i5Horn <- trends5[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i10Horn <- trends10[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i20Horn <- trends20[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

## Compare variance structures with duration.sc ############################
#fixed <- 'Jtu.sc ~ tempchange_abs.sc*REALM + tempchange_abs.sc*tempave_metab.sc + tempchange_abs.sc*duration.sc'
fixed <- 'Jtu.sc ~ duration.sc + REALM:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + mass.sc:duration.sc + nspp.sc:duration.sc + seas.sc:duration.sc + microclim.sc:duration.sc + npp.sc:duration.sc + human_bowler.sc:REALM2:duration.sc'

if(fitmod == 'modRFgauss3'){
    print(paste(sum(i3), 'data points'))
    modRFgauss3 <- glmmTMB(formula(fixed), data = trends3[i3,], family = gaussian(), control = glmmTMBControl(profile=TRUE))
    summary(modRFgauss3)
    saveRDS(modRFgauss3, file = 'temp/modRFgauss3.rds')
    print('saved modRFgauss3.rds')
    print(Sys.time())
    print(performance::r2(modRFgauss3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFgauss5'){
    print(paste(sum(i5), 'data points'))
    modRFgauss5 <- glmmTMB(formula(fixed), data = trends5[i5,], family = gaussian(), control = glmmTMBControl(profile=TRUE))
    summary(modRFgauss5)
    saveRDS(modRFgauss5, file = 'temp/modRFgauss5.rds')
    print('saved modRFgauss5.rds')
    print(Sys.time())
    print(performance::r2(modRFgauss5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFgauss10'){
    print(paste(sum(i10), 'data points'))
    modRFgauss10 <- glmmTMB(formula(fixed), data = trends10[i10,], family = gaussian(), control = glmmTMBControl(profile=TRUE))
    summary(modRFgauss10)
    saveRDS(modRFgauss10, file = 'temp/modRFgauss10.rds')
    print('saved modRFgauss10.rds')
    print(Sys.time())
    print(performance::r2(modRFgauss10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFgauss20'){
    print(paste(sum(i20), 'data points'))
    modRFgauss20 <- glmmTMB(formula(fixed), data = trends20[i20,], family = gaussian(), control = glmmTMBControl(profile=TRUE))
    summary(modRFgauss20)
    saveRDS(modRFgauss20, file = 'temp/modRFgauss20.rds')
    print('saved modRFgauss20.rds')
    print(Sys.time())
    print(performance::r2(modRFgauss20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFbeta3'){
    print(paste(sum(i3), 'data points'))
    modRFbeta3 <- glmmTMB(formula(fixed), data = trends3[i3,], family = beta_family(link = 'logit'), control = glmmTMBControl(profile=TRUE)) # add beta errors
    summary(modRFbeta3)
    saveRDS(modRFbeta3, file = 'temp/modRFbeta3.rds')
    print('saved modRFbeta3.rds')
    print(Sys.time())
    print(performance::r2(modRFbeta3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFbeta5'){
    print(paste(sum(i5), 'data points'))
    modRFbeta5 <- glmmTMB(formula(fixed), data = trends5[i5,], family = beta_family(link = 'logit'), control = glmmTMBControl(profile=TRUE)) # add beta errors
    summary(modRFbeta5)
    saveRDS(modRFbeta5, file = 'temp/modRFbeta5.rds')
    print('saved modRFbeta5.rds')
    print(Sys.time())
    print(performance::r2(modRFbeta5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFbeta10'){
    print(paste(sum(i10), 'data points'))
    modRFbeta10 <- glmmTMB(formula(fixed), data = trends10[i10,], family = beta_family(link = 'logit'), control = glmmTMBControl(profile=TRUE)) # add beta errors
    summary(modRFbeta10)
    saveRDS(modRFbeta10, file = 'temp/modRFbeta10.rds')
    print('saved modRFbeta10.rds')
    print(Sys.time())
    print(performance::r2(modRFbeta10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFbeta20'){
    print(paste(sum(i20), 'data points'))
    modRFbeta20 <- glmmTMB(formula(fixed), data = trends20[i20,], family = beta_family(link = 'logit'), control = glmmTMBControl(profile=TRUE)) # add beta errors
    summary(modRFbeta20)
    saveRDS(modRFbeta20, file = 'temp/modRFbeta20.rds')
    print('saved modRFbeta20.rds')
    print(Sys.time())
    print(performance::r2(modRFbeta20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFrID3'){
    print(paste(sum(i3), 'data points'))
    modRFrID3 <- glmmTMB(formula(paste(fixed, '+(1|rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random effects
    summary(modRFrID3)
    saveRDS(modRFrID3, file = 'temp/modRFrID3.rds')
    print('saved modRDrID3.rds')
    print(Sys.time())
    print(performance::r2(modRFrID3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFrID5'){
    print(paste(sum(i5), 'data points'))
    modRFrID5 <- glmmTMB(formula(paste(fixed, '+(1|rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                         control = glmmTMBControl(profile=TRUE)) # add random effects
    summary(modRFrID5)
    saveRDS(modRFrID5, file = 'temp/modRFrID5.rds')
    print('saved modRDrID5.rds')
    print(Sys.time())
    print(performance::r2(modRFrID5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFrID10'){
    print(paste(sum(i10), 'data points'))
    modRFrID10 <- glmmTMB(formula(paste(fixed, '+(1|rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                         control = glmmTMBControl(profile=TRUE)) # add random effects
    summary(modRFrID10)
    saveRDS(modRFrID10, file = 'temp/modRFrID10.rds')
    print('saved modRDrID10.rds')
    print(Sys.time())
    print(performance::r2(modRFrID10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFrID20'){
    print(paste(sum(i20), 'data points'))
    modRFrID20 <- glmmTMB(formula(paste(fixed, '+(1|rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                         control = glmmTMBControl(profile=TRUE)) # add random effects
    summary(modRFrID20)
    saveRDS(modRFrID20, file = 'temp/modRFrID20.rds')
    print('saved modRDrID20.rds')
    print(Sys.time())
    print(performance::r2(modRFrID20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRF2lev3'){
    print(paste(sum(i3), 'data points'))
    modRF2lev3 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                             control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRF2lev3)
    saveRDS(modRF2lev3, file = 'temp/modRF2lev3.rds')
    print('saved modRF2lev3.rds')
    print(Sys.time())
    print(performance::r2(modRF2lev3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2lev5'){
    print(paste(sum(i5), 'data points'))
    modRF2lev5 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                          control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRF2lev5)
    saveRDS(modRF2lev5, file = 'temp/modRF2lev5.rds')
    print('saved modRF2lev5.rds')
    print(Sys.time())
    print(performance::r2(modRF2lev5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2lev10'){
    print(paste(sum(i10), 'data points'))
    modRF2lev10 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                          control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRF2lev10)
    saveRDS(modRF2lev10, file = 'temp/modRF2lev10.rds')
    print('saved modRF2lev10.rds')
    print(Sys.time())
    print(performance::r2(modRF2lev10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2lev20'){
    print(paste(sum(i20), 'data points'))
    modRF2lev20 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                          control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRF2lev20)
    saveRDS(modRF2lev20, file = 'temp/modRF2lev20.rds')
    print('saved modRF2lev20.rds')
    print(Sys.time())
    print(performance::r2(modRF2lev20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFnestedRE3'){
    print(paste(sum(i3), 'data points'))
    modRFnestedRE3 <- glmmTMB(formula(paste(fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRFnestedRE3)
    saveRDS(modRFnestedRE3, file = 'temp/modRFnestedRE3.rds')
    print('saved modRFnestedRE3.rds')
    print(Sys.time())
    print(performance::r2(modRFnestedRE3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFnestedRE5'){
    print(paste(sum(i5), 'data points'))
    modRFnestedRE5 <- glmmTMB(formula(paste(fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                              control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRFnestedRE5)
    saveRDS(modRFnestedRE5, file = 'temp/modRFnestedRE5.rds')
    print('saved modRFnestedRE5.rds')
    print(Sys.time())
    print(performance::r2(modRFnestedRE5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFnestedRE10'){
    print(paste(sum(i10), 'data points'))
    modRFnestedRE10 <- glmmTMB(formula(paste(fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                              control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRFnestedRE10)
    saveRDS(modRFnestedRE10, file = 'temp/modRFnestedRE10.rds')
    print('saved modRFnestedRE10.rds')
    print(Sys.time())
    print(performance::r2(modRFnestedRE10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFnestedRE20'){
    print(paste(sum(i20), 'data points'))
    modRFnestedRE20 <- glmmTMB(formula(paste(fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                              control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRFnestedRE20)
    saveRDS(modRFnestedRE20, file = 'temp/modRFnestedRE20.rds')
    print('saved modRFnestedRE20.rds')
    print(Sys.time())
    print(performance::r2(modRFnestedRE20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFslopeRE3'){
    # fails to converge
    print(paste(sum(i3), 'data points'))
    modRFslopeRE3 <- glmmTMB(formula(paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'),
                   control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE3)
    saveRDS(modRFslopeRE3, file = 'temp/modRFslopeRE3.rds')
    print('saved modRFslopeRE3.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE5'){
    # fails to converge
    print(paste(sum(i5), 'data points'))
    modRFslopeRE5 <- glmmTMB(formula(paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'),
                             control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE5)
    saveRDS(modRFslopeRE5, file = 'temp/modRFslopeRE5.rds')
    print('saved modRFslopeRE5.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE10'){
    print(paste(sum(i10), 'data points'))
    modRFslopeRE10 <- glmmTMB(formula(paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                             control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE10)
    saveRDS(modRFslopeRE10, file = 'temp/modRFslopeRE10.rds')
    print('saved modRFslopeRE10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE20'){
    print(paste(sum(i20), 'data points'))
    modRFslopeRE20 <- glmmTMB(formula(paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                             control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE20)
    saveRDS(modRFslopeRE20, file = 'temp/modRFslopeRE20.rds')
    print('saved modRFslopeRE20.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFslopeRE2lev3'){
    # fails to converge
    print(paste(sum(i3), 'data points'))
    modRFslopeRE2lev3 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'),
                             control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE2lev3)
    saveRDS(modRFslopeRE2lev3, file = 'temp/modRFslopeRE2lev3.rds')
    print('saved modRFslopeRE2lev3.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2lev3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2lev5'){
    # fails to converge
    print(paste(sum(i5), 'data points'))
    modRFslopeRE2lev5 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'),
                                 control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE2lev5)
    saveRDS(modRFslopeRE2lev5, file = 'temp/modRFslopeRE2lev5.rds')
    print('saved modRFslopeRE2lev5.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2lev5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2lev10'){
    print(paste(sum(i10), 'data points'))
    modRFslopeRE2lev10 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                 control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE2lev10)
    saveRDS(modRFslopeRE2lev10, file = 'temp/modRFslopeRE2lev10.rds')
    print('saved modRFslopeRE2lev10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2lev10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2lev20'){
    print(paste(sum(i20), 'data points'))
    modRFslopeRE2lev20 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                 control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE2lev20)
    saveRDS(modRFslopeRE2lev20, file = 'temp/modRFslopeRE2lev20.rds')
    print('saved modRFslopeRE2lev20.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2lev20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdisp3'){ # no REs, + dispersion formula
    print(paste(sum(i3), 'data points'))
    modRFdisp3 <- glmmTMB(formula(fixed), data = trends3[i3,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp3)
    saveRDS(modRFdisp3, file = 'temp/modRFdisp3.rds')
    print('saved modRFdisp3.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdisp5'){
    print(paste(sum(i5), 'data points'))
    modRFdisp5 <- glmmTMB(formula(fixed), data = trends5[i5,], family = beta_family(link = 'logit'), 
                          dispformula = ~nspp.sc, 
                          control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp5)
    saveRDS(modRFdisp5, file = 'temp/modRFdisp5.rds')
    print('saved modRFdisp5.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdisp10'){
    print(paste(sum(i10), 'data points'))
    modRFdisp10 <- glmmTMB(formula(fixed), data = trends10[i10,], family = beta_family(link = 'logit'), 
                          dispformula = ~nspp.sc, 
                          control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp10)
    saveRDS(modRFdisp10, file = 'temp/modRFdisp10.rds')
    print('saved modRFdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdisp20'){
    print(paste(sum(i20), 'data points'))
    modRFdisp20 <- glmmTMB(formula(fixed), data = trends20[i20,], family = beta_family(link = 'logit'), 
                          dispformula = ~nspp.sc, 
                          control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp20)
    saveRDS(modRFdisp20, file = 'temp/modRFdisp20.rds')
    print('saved modRFdisp20.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRF2levdisp3'){ # 2 level RE, no slope, + dispersion formula
    print(paste(sum(i3), 'data points'))
    modRFdisp2levOnlyint3 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp2levOnlyint3)
    saveRDS(modRFdisp2levOnlyint3, file = 'temp/modRFdisp2levOnlyint3.rds')
    print('saved modRFdisp2levOnlyint3.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp2levOnlyint3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2levdisp5'){
    print(paste(sum(i5), 'data points'))
    modRFdisp2levOnlyint5 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp2levOnlyint5)
    saveRDS(modRFdisp2levOnlyint5, file = 'temp/modRFdisp2levOnlyint5.rds')
    print('saved modRFdisp2levOnlyint5.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp2levOnlyint5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2levdisp10'){
    print(paste(sum(i10), 'data points'))
    modRFdisp2levOnlyint10 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp2levOnlyint10)
    saveRDS(modRFdisp2levOnlyint10, file = 'temp/modRFdisp2levOnlyint10.rds')
    print('saved modRFdisp2levOnlyint10.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp2levOnlyint10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2levdisp20'){
    print(paste(sum(i20), 'data points'))
    modRFdisp2levOnlyint20 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE))
    summary(modRFdisp2levOnlyint20)
    saveRDS(modRFdisp2levOnlyint20, file = 'temp/modRFdisp2levOnlyint20.rds')
    print('saved modRFdisp2levOnlyint20.rds')
    print(Sys.time())
    print(performance::r2(modRFdisp2levOnlyint20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFslopeRE2levdisp3'){
    # singular convergence
    print(paste(sum(i3), 'data points'))
    modRFslopeRE2levdisp3 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'),
                         dispformula = ~nspp.sc,
                         control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeRE2levdisp3)
    saveRDS(modRFslopeRE2levdisp3, file = 'temp/modRFslopeRE2levdisp3.rds')
    print('saved modRFslopeRE2levdisp3.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisp3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2levdisp5'){
    # singular convergence
    print(paste(sum(i5), 'data points'))
    modRFslopeRE2levdisp5 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'),
                                     dispformula = ~nspp.sc,
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeRE2levdisp5)
    saveRDS(modRFslopeRE2levdisp5, file = 'temp/modRFslopeRE2levdisp5.rds')
    print('saved modRFslopeRE2levdisp5.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisp))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2levdisp10'){
    print(paste(sum(i10), 'data points'))
    modRFslopeRE2levdisp10 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeRE2levdisp10)
    saveRDS(modRFslopeRE2levdisp10, file = 'temp/modRFslopeRE2levdisp10.rds')
    print('saved modRFslopeRE2levdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisp10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2levdisp20'){
    print(paste(sum(i20), 'data points'))
    modRFslopeRE2levdisp20 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeRE2levdisp20)
    saveRDS(modRFslopeRE2levdisp20, file = 'temp/modRFslopeRE2levdisp20.rds')
    print('saved modRFslopeRE2levdisp20.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisp20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFslopeREdisp3'){
    # singular convergence
    print(paste(sum(i3), 'data points'))
    modRFslopeREdisp3 <- glmmTMB(formula(paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'),
                   dispformula = ~nspp.sc,
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeREdisp3)
    saveRDS(modRFslopeREdisp3, file = 'temp/modRFslopeREdisp3.rds')
    print('saved modRFslopeREdisp3.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeREdisp3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeREdisp5'){
    # singular convergence
    print(paste(sum(i5), 'data points'))
    modRFslopeREdisp5 <- glmmTMB(formula(paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'),
                                     dispformula = ~nspp.sc,
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeREdisp5)
    saveRDS(modRFslopeREdisp5, file = 'temp/modRFslopeREdisp5.rds')
    print('saved modRFslopeREdisp5.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeREdisp5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeREdisp10'){
    print(paste(sum(i10), 'data points'))
    modRFslopeREdisp10 <- glmmTMB(formula(paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeREdisp10)
    saveRDS(modRFslopeREdisp10, file = 'temp/modRFslopeREdisp10.rds')
    print('saved modRFslopeREdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeREdisp10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeREdisp20'){
    print(paste(sum(i20), 'data points'))
    modRFslopeREdisp20 <- glmmTMB(formula(paste(fixed, '+(duration.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeREdisp20)
    saveRDS(modRFslopeREdisp20, file = 'temp/modRFslopeREdisp20.rds')
    print('saved modRFslopeREdisp20.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeREdisp20))
    MATCHMOD <- TRUE
}


if(fitmod == 'modRF2levdisprealm3'){ # 2 level RE, slope vs. duration, disp formula by realm
    print(paste(sum(i3), 'data points'))
    modRF2levdisprealm3 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                                          dispformula = ~nspp.sc+REALM, 
                                          control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRF2levdisprealm3)
    saveRDS(modRF2levdisprealm3, file = 'temp/modRF2levdisprealm3.rds')
    print('saved modRF2levdisprealm3.rds')
    print(Sys.time())
    print(performance::r2(modRF2levdisprealm3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2levdisprealm5'){ # 2 level RE, slope vs. duration, disp formula by realm
    print(paste(sum(i5), 'data points'))
    modRF2levdisprealm5 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                                          dispformula = ~nspp.sc+REALM, 
                                          control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRF2levdisprealm5)
    saveRDS(modRF2levdisprealm5, file = 'temp/modRF2levdisprealm5.rds')
    print('saved modRF2levdisprealm5.rds')
    print(Sys.time())
    print(performance::r2(modRF2levdisprealm5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2levdisprealm10'){ # 2 level RE, slope vs. duration, disp formula by realm
    print(paste(sum(i10), 'data points'))
    modRF2levdisprealm10 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                           dispformula = ~nspp.sc+REALM, 
                                           control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRF2levdisprealm10)
    saveRDS(modRF2levdisprealm10, file = 'temp/modRF2levdisprealm10.rds')
    print('saved modRF2levdisprealm10.rds')
    print(Sys.time())
    print(performance::r2(modRF2levdisprealm10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRF2levdisprealm20'){ # 2 level RE, slope vs. duration, disp formula by realm
    print(paste(sum(i20), 'data points'))
    modRF2levdisprealm20 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                           dispformula = ~nspp.sc+REALM, 
                                           control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRF2levdisprealm20)
    saveRDS(modRF2levdisprealm20, file = 'temp/modRF2levdisprealm20.rds')
    print('saved modRF2levdisprealm20.rds')
    print(Sys.time())
    print(performance::r2(modRF2levdisprealm20))
    MATCHMOD <- TRUE
}


if(fitmod == 'modRFslopeRE2levdisprealm3'){ # 2 level RE, slope vs. duration, disp formula by realm
    print(paste(sum(i3), 'data points'))
    modRFslopeRE2levdisprealm3 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc+REALM, 
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeRE2levdisprealm3)
    saveRDS(modRFslopeRE2levdisprealm3, file = 'temp/modRFslopeRE2levdisprealm3.rds')
    print('saved modRFslopeRE2levdisprealm3.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisprealm3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2levdisprealm5'){
    print(paste(sum(i5), 'data points'))
    modRFslopeRE2levdisprealm5 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                                        dispformula = ~nspp.sc+REALM, 
                                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeRE2levdisprealm5)
    saveRDS(modRFslopeRE2levdisprealm5, file = 'temp/modRFslopeRE2levdisprealm5.rds')
    print('saved modRFslopeRE2levdisprealm5.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisprealm5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2levdisprealm10'){    
    print(paste(sum(i10), 'data points'))
    modRFslopeRE2levdisprealm10 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                        dispformula = ~nspp.sc+REALM, 
                                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeRE2levdisprealm10)
    saveRDS(modRFslopeRE2levdisprealm10, file = 'temp/modRFslopeRE2levdisprealm10.rds')
    print('saved modRFslopeRE2levdisprealm10.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisprealm10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFslopeRE2levdisprealm20'){
    print(paste(sum(i20), 'data points'))
    modRFslopeRE2levdisprealm20 <- glmmTMB(formula(paste(fixed, '+(duration.sc|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                        dispformula = ~nspp.sc+REALM, 
                                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFslopeRE2levdisprealm20)
    saveRDS(modRFslopeRE2levdisprealm20, file = 'temp/modRFslopeRE2levdisprealm20.rds')
    print('saved modRFslopeRE2levdisprealm20.rds')
    print(Sys.time())
    print(performance::r2(modRFslopeRE2levdisprealm20))
    MATCHMOD <- TRUE
}


## Compare variance structures with durationlog.sc ############################
fixed <- 'Jtu.sc ~ durationlog.sc + REALM:durationlog.sc + tempchange_abs.sc:durationlog.sc + tempave_metab.sc:durationlog.sc + mass.sc:durationlog.sc + nspp.sc:durationlog.sc + seas.sc:durationlog.sc + microclim.sc:durationlog.sc + npp.sc:durationlog.sc + human_bowler.sc:REALM2:durationlog.sc'

if(fitmod == 'modRFdurlog2levdisp3'){ # 2 level RE, no slope, + dispersion formula
    print(paste(sum(i3), 'data points'))
    modRFdurlog2levdisp3 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE))
    summary(modRFdurlog2levdisp3)
    saveRDS(modRFdurlog2levdisp3, file = 'temp/modRFdurlog2levdisp3.rds')
    print('saved modRFdurlog2levdisp3.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlog2levdisp3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlog2levdisp5'){
    print(paste(sum(i5), 'data points'))
    modRFdurlog2levdisp5 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc, 
                                     control = glmmTMBControl(profile=TRUE))
    summary(modRFdurlog2levdisp5)
    saveRDS(modRFdurlog2levdisp5, file = 'temp/modRFdurlog2levdisp5.rds')
    print('saved modRFdurlog2levdisp5.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlog2levdisp5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlog2levdisp10'){
    print(paste(sum(i10), 'data points'))
    modRFdurlog2levdisp10 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc, 
                                      control = glmmTMBControl(profile=TRUE))
    summary(modRFdurlog2levdisp10)
    saveRDS(modRFdurlog2levdisp10, file = 'temp/modRFdurlog2levdisp10.rds')
    print('saved modRFdurlog2levdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlog2levdisp10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlog2levdisp20'){
    print(paste(sum(i20), 'data points'))
    modRFdurlog2levdisp20 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc, 
                                      control = glmmTMBControl(profile=TRUE))
    summary(modRFdurlog2levdisp20)
    saveRDS(modRFdurlog2levdisp20, file = 'temp/modRFdurlog2levdisp20.rds')
    print('saved modRFdurlog2levdisp20.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlog2levdisp20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdurlogslopeRE2levdisp3'){
    print(paste(sum(i3), 'data points'))
    modRFdurlogslopeRE2levdisp3 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'),
                         dispformula = ~nspp.sc,
                         control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeRE2levdisp3)
    saveRDS(modRFdurlogslopeRE2levdisp3, file = 'temp/modRFdurlogslopeRE2levdisp3.rds')
    print('saved modRFdurlogslopeRE2levdisp3.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeRE2levdisp3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlogslopeRE2levdisp5'){
    print(paste(sum(i5), 'data points'))
    modRFdurlogslopeRE2levdisp5 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'),
                                     dispformula = ~nspp.sc,
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeRE2levdisp5)
    saveRDS(modRFdurlogslopeRE2levdisp5, file = 'temp/modRFdurlogslopeRE2levdisp5.rds')
    print('saved modRFdurlogslopeRE2levdisp5.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeRE2levdisp))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlogslopeRE2levdisp10'){
    print(paste(sum(i10), 'data points'))
    modRFdurlogslopeRE2levdisp10 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc, 
                                      control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeRE2levdisp10)
    saveRDS(modRFdurlogslopeRE2levdisp10, file = 'temp/modRFdurlogslopeRE2levdisp10.rds')
    print('saved modRFdurlogslopeRE2levdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeRE2levdisp10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlogslopeRE2levdisp20'){
    print(paste(sum(i20), 'data points'))
    modRFdurlogslopeRE2levdisp20 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc, 
                                      control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeRE2levdisp20)
    saveRDS(modRFdurlogslopeRE2levdisp20, file = 'temp/modRFdurlogslopeRE2levdisp20.rds')
    print('saved modRFdurlogslopeRE2levdisp20.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeRE2levdisp20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdurlogslopeREdisp3'){
    print(paste(sum(i3), 'data points'))
    modRFdurlogslopeREdisp3 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'),
                   dispformula = ~nspp.sc,
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeREdisp3)
    saveRDS(modRFdurlogslopeREdisp3, file = 'temp/modRFdurlogslopeREdisp3.rds')
    print('saved modRFdurlogslopeREdisp3.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeREdisp3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlogslopeREdisp5'){
    print(paste(sum(i5), 'data points'))
    modRFdurlogslopeREdisp5 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'),
                                     dispformula = ~nspp.sc,
                                     control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeREdisp5)
    saveRDS(modRFdurlogslopeREdisp5, file = 'temp/modRFdurlogslopeREdisp5.rds')
    print('saved modRFdurlogslopeREdisp5.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeREdisp5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlogslopeREdisp10'){
    print(paste(sum(i10), 'data points'))
    modRFdurlogslopeREdisp10 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                  dispformula = ~nspp.sc, 
                                  control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeREdisp10)
    saveRDS(modRFdurlogslopeREdisp10, file = 'temp/modRFdurlogslopeREdisp10.rds')
    print('saved modRFdurlogslopeREdisp10.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeREdisp10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlogslopeREdisp20'){
    print(paste(sum(i20), 'data points'))
    modRFdurlogslopeREdisp20 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                  dispformula = ~nspp.sc, 
                                  control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeREdisp20)
    saveRDS(modRFdurlogslopeREdisp20, file = 'temp/modRFdurlogslopeREdisp20.rds')
    print('saved modRFdurlogslopeREdisp20.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeREdisp20))
    MATCHMOD <- TRUE
}


if(fitmod == 'modRFdurlog2levdisprealm3'){ # 2 level RE, no slope, disp formula by realm
    print(paste(sum(i3), 'data points'))
    modRFdurlog2levdisprealm3 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                                   dispformula = ~nspp.sc+REALM, 
                                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlog2levdisprealm3)
    saveRDS(modRFdurlog2levdisprealm3, file = 'temp/modRFdurlog2levdisprealm3.rds')
    print('saved modRFdurlog2levdisprealm3.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlog2levdisprealm3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlog2levdisprealm5'){
    print(paste(sum(i5), 'data points'))
    modRFdurlog2levdisprealm5 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                                   dispformula = ~nspp.sc+REALM, 
                                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlog2levdisprealm5)
    saveRDS(modRFdurlog2levdisprealm5, file = 'temp/modRFdurlog2levdisprealm5.rds')
    print('saved modRFdurlog2levdisprealm5.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlog2levdisprealm5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlog2levdisprealm10'){
    print(paste(sum(i10), 'data points'))
    modRFdurlog2levdisprealm10 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc+REALM, 
                                    control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlog2levdisprealm10)
    saveRDS(modRFdurlog2levdisprealm10, file = 'temp/modRFdurlog2levdisprealm10.rds')
    print('saved modRFdurlog2levdisprealm10.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlog2levdisprealm10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlog2levdisprealm20'){
    print(paste(sum(i20), 'data points'))
    modRFdurlog2levdisprealm20 <- glmmTMB(formula(paste(fixed, '+(1|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc+REALM, 
                                    control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlog2levdisprealm20)
    saveRDS(modRFdurlog2levdisprealm20, file = 'temp/modRFdurlog2levdisprealm20.rds')
    print('saved modRFdurlog2levdisprealm20.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlog2levdisprealm20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdurlogslopeRE2levdisprealm3'){ # 2 level RE, slope vs. durationlog, disp formula by realm
    print(paste(sum(i3), 'data points'))
    modRFdurlogslopeRE2levdisprealm3 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|STUDY_ID/rarefyID)')), data = trends3[i3,], family = beta_family(link = 'logit'), 
                                          dispformula = ~nspp.sc+REALM, 
                                          control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeRE2levdisprealm3)
    saveRDS(modRFdurlogslopeRE2levdisprealm3, file = 'temp/modRFdurlogslopeRE2levdisprealm3.rds')
    print('saved modRFdurlogslopeRE2levdisprealm3.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeRE2levdisprealm3))
    MATCHMOD <- TRUE
}   
if(fitmod == 'modRFdurlogslopeRE2levdisprealm5'){
    print(paste(sum(i5), 'data points'))
    modRFdurlogslopeRE2levdisprealm5 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|STUDY_ID/rarefyID)')), data = trends5[i5,], family = beta_family(link = 'logit'), 
                                          dispformula = ~nspp.sc+REALM, 
                                          control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeRE2levdisprealm5)
    saveRDS(modRFdurlogslopeRE2levdisprealm5, file = 'temp/modRFdurlogslopeRE2levdisprealm5.rds')
    print('saved modRFdurlogslopeRE2levdisprealm5.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeRE2levdisprealm5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlogslopeRE2levdisprealm10'){
        print(paste(sum(i10), 'data points'))
    modRFdurlogslopeRE2levdisprealm10 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|STUDY_ID/rarefyID)')), data = trends10[i10,], family = beta_family(link = 'logit'), 
                                           dispformula = ~nspp.sc+REALM, 
                                           control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeRE2levdisprealm10)
    saveRDS(modRFdurlogslopeRE2levdisprealm10, file = 'temp/modRFdurlogslopeRE2levdisprealm10.rds')
    print('saved modRFdurlogslopeRE2levdisprealm10.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeRE2levdisprealm10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFdurlogslopeRE2levdisprealm20'){
    print(paste(sum(i20), 'data points'))
    modRFdurlogslopeRE2levdisprealm20 <- glmmTMB(formula(paste(fixed, '+(durationlog.sc|STUDY_ID/rarefyID)')), data = trends20[i20,], family = beta_family(link = 'logit'), 
                                           dispformula = ~nspp.sc+REALM, 
                                           control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdurlogslopeRE2levdisprealm20)
    saveRDS(modRFdurlogslopeRE2levdisprealm20, file = 'temp/modRFdurlogslopeRE2levdisprealm20.rds')
    print('saved modRFdurlogslopeRE2levdisprealm20.rds')
    print(Sys.time())
    print(performance::r2(modRFdurlogslopeRE2levdisprealm20))
    MATCHMOD <- TRUE
}



# Realm:duration ############################

if(fitmod == 'modRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + REALM:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modRealm10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + REALM:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modRealm20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + REALM:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modRealm5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + REALM:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modRealm10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + REALM:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modRealm20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + REALM:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T:duration ############################

if(fitmod == 'modT5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modT10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modT20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modT5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modT10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modT20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# dT:duration ############################

if(fitmod == 'moddT5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddT10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddT20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'moddT5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddT10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddT20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# mass:duration #########################

if(fitmod == 'modmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmass10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmass20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc +
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modmass5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmass10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc +
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmass20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc +
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modmass5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc +
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmass10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc +
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmass20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc +
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# npp:duration #########################

if(fitmod == 'modnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc +
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modnpp10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc +
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modnpp20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modnpp5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modnpp10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modnpp20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# seas:duration #########################

if(fitmod == 'modseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modseas10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modseas20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modseas5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modseas10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modseas20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# microclim:duration #########################

if(fitmod == 'modmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmicroclim10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmicroclim20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modmicroclim5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmicroclim10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modmicroclim20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# human:duration #########################

if(fitmod == 'modhuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modhuman10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modhuman20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modhuman5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modhuman10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modhuman20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T+dT:duration ############################

if(fitmod == 'modTdT3Jtu'){
    print(paste(sum(i3Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends3[i3Jtu,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdT5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdT10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdT20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modTdT5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                       data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                       dispformula = ~nspp.sc+REALM, 
                       control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdT10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdT20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+microclim:duration ############################

if(fitmod == 'modTmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTmicroclim10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTmicroclim20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modTmicroclim5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTmicroclim10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTmicroclim20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# dT+microclim:duration ############################

if(fitmod == 'moddTmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTmicroclim10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTmicroclim20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'moddTmicroclim5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTmicroclim10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTmicroclim20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T+dT+microclim:duration ############################

if(fitmod == 'modTdTmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmicroclim10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmicroclim20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modTdTmicroclim5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmicroclim10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmicroclim20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T:dT:duration ############################

if(fitmod == 'modTdTT3Jtu'){
    print(paste(sum(i3Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                       data = trends3[i3Jtu,], family = beta_family(link = 'logit'), 
                       dispformula = ~nspp.sc+REALM, 
                       control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                       data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                       dispformula = ~nspp.sc+REALM, 
                       control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTT3Jne'){
    print(paste(sum(i3Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends3[i3Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTT3Horn'){
    print(paste(sum(i3Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends3[i3Horn,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                         data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc+REALM, 
                         control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                         data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc+REALM, 
                         control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# T:dT:REALM:duration #########################

if(fitmod == 'modTdTTRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm10Jtu'){ 
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTRealm5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTRealm5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm10Horn'){ 
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# try on Antao et al. 2020-style dataset
# use tempchange, not tempchange_abs
if(fitmod == 'modAntaoTdTTRealm5'){ # trims out freshwater and >60 or <23.5 deg lat, like Antao et al. 2020
    i5 <- trends5[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                      complete.cases(Jtu.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modAntaoTdTTRealm10'){ # trims out freshwater and >60 or <23.5 deg lat, like Antao et al. 2020
    i10 <- trends10[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                      complete.cases(Jtu.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
   MATCHMOD <- TRUE
}
if(fitmod == 'modAntaoTdTTRealm20'){ # trims out freshwater and >60 or <23.5 deg lat, like Antao et al. 2020
    i20 <- trends20[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                      complete.cases(Jtu.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm20.rds')
    print('saved modAntaoTdTTRealm20.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modAntaoTdTTRealm5Jne'){ # Jne, like gains/losses of Antao et al. 2020
    i5 <- trends5[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                      complete.cases(Jne.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modAntaoTdTTRealm10Jne'){
    i10 <- trends10[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                      complete.cases(Jne.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modAntaoTdTTRealm20Jne'){
    i20 <- trends20[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                        complete.cases(Jne.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modAntaoTdTTRealm5Horn'){ # Horn
    i5 <- trends5[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                      complete.cases(Horn.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modAntaoTdTTRealm10Horn'){
    i10 <- trends10[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                        complete.cases(Horn.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modAntaoTdTTRealm20Horn'){
    i20 <- trends20[, REALM != 'Freshwater' & abs(rarefyID_y) <= 60 & abs(rarefyID_y) >= 23.5 &
                        complete.cases(Horn.sc, tempchange.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T:dT:tsign:duration #########################

if(fitmod == 'modTdTTtsign5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc + 
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTtsign5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTtsign5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T:dT+mass:duration #########################

if(fitmod == 'modTdTTmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTmass5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTmass5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T:dT+npp:duration #########################

if(fitmod == 'modTdTTnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTnpp10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTnpp20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTnpp5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTnpp10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTnpp20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T:dT+seas:duration #########################

if(fitmod == 'modTdTTseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTseas10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTseas20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTseas5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTseas10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTseas20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T:dT+microclim:duration #########################

if(fitmod == 'modTdTTmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmicroclim10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmicroclim20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTmicroclim5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmicroclim10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmicroclim20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T:dT+human:duration #########################

if(fitmod == 'modTdTThuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTThuman10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTThuman20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTThuman5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTThuman10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTThuman20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# Null models ############################

if(fitmod == 'modNull5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ 1 +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], 
                   family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM,  # add dispersion formula
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

if(fitmod == 'modNull10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ 1 +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], 
                   family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM,  # add dispersion formula
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

if(fitmod == 'modNull20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ 1 +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], 
                   family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM,  # add dispersion formula
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

if(fitmod == 'modNull5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ 1 +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], 
                   family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM,  # add dispersion formula
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

if(fitmod == 'modNull10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ 1 +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], 
                   family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM,  # add dispersion formula
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

if(fitmod == 'modNull20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ 1 +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], 
                   family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM,  # add dispersion formula
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# print and save results ############################
if(MATCHMOD == FALSE) stop("Model name did not match anything", call.=FALSE)
if(MATCHMOD){
    print(summary(mod))
    saveRDS(mod, file = paste0('temp/', fitmod, '.rds'))
    print(paste0('saved ', fitmod, '.rds'))
    print(Sys.time())
    print(performance::r2(mod))
}

print(warnings())

print(paste('Ended', Sys.time(), sep =''))
#!/usr/bin/env Rscript

# Script to fit glmmTMB models
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_vs_temperature_GLMM_fit.R modRFgauss > logs/turnover_vs_temperature_GLMMmodRFgauss.Rout &
# (this works if code is executable, e.g., chmod u+x turnover_vs_temperature_GLMM_fit.R)
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
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(performance) # for R2

# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz')
trends5 <- fread('output/turnover_w_covariates5.csv.gz')
trends10 <- fread('output/turnover_w_covariates10.csv.gz')
trends20 <- fread('output/turnover_w_covariates20.csv.gz')

trendsall[, tsign := as.factor(tsign)]
trends5[, tsign := as.factor(tsign)]
trends10[, tsign := as.factor(tsign)]
trends20[, tsign := as.factor(tsign)]

# Models ############################

## Choose dataset
iallJtu <- trendsall[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i5Jtu <- trends5[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i10Jtu <- trends10[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i20Jtu <- trends20[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

i5Jne <- trends5[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i10Jne <- trends10[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i20Jne <- trends20[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

i5Horn <- trends5[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i10Horn <- trends10[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i20Horn <- trends20[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

trends5Jtu <- trends5[i5Jtu,]
trends5Horn <- trends5[i5Horn,]

## Compare variance structures with duration.sc ############################
#fixed <- 'Jtu.sc ~ tempchange_abs.sc*REALM + tempchange_abs.sc*tempave_metab.sc + tempchange_abs.sc*duration.sc'
fixed <- 'Jtu.sc ~ duration.sc + REALM:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + mass.sc:duration.sc + nspp.sc:duration.sc + seas.sc:duration.sc + microclim.sc:duration.sc + npp.sc:duration.sc + human_bowler.sc:REALM2:duration.sc'

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


# SINGLE FACTORS ############################
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

# sdT:duration ############################

if(fitmod == 'modsdT5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdT5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
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


# T+SINGLE FACTORS ############################
# T+REALM:duration ############################

if(fitmod == 'modTRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTRealm10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTRealm20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modTRealm5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTRealm10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTRealm20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+dT:duration ############################

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

# T+sdT:duration ############################

if(fitmod == 'modTsdT5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + 
                       tempave_metab.sc:duration.sc + 
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTsdT5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + 
                       tempave_metab.sc:duration.sc + 
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T+mass:duration ############################

if(fitmod == 'modTmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTmass10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTmass20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modTmass5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTmass10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTmass20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+npp:duration ############################

if(fitmod == 'modTnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTnpp10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTnpp20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modTnpp5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTnpp10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTnpp20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+seas:duration ############################

if(fitmod == 'modTseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTseas10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTseas20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modTseas5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTseas10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTseas20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
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




# T+human:duration ############################

if(fitmod == 'modThuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modThuman10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modThuman20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modThuman5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modThuman10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modThuman20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# dT+SINGLE FACTORS ############################
# dT+REALM:duration ############################

if(fitmod == 'moddTRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTRealm10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTRealm20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'moddTRealm5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTRealm10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTRealm20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# dT+mass:duration ############################

if(fitmod == 'moddTmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTmass10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTmass20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTmass5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTmass10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTmass20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# dT+npp:duration ############################

if(fitmod == 'moddTnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempchange_abs.sc:duration.sc + 
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTnpp10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTnpp20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'moddTnpp5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTnpp10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTnpp20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# dT+seas:duration ############################

if(fitmod == 'moddTseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTseas10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTseas20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTseas5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTseas10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTseas20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
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

# dT+human:duration ############################

if(fitmod == 'moddThuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddThuman10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddThuman20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddThuman5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddThuman10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddThuman20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# sdT+SINGLE FACTORS ############################
# sdT+REALM:duration ############################

if(fitmod == 'modsdTRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTRealm5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# sdT+mass:duration ############################

if(fitmod == 'modsdTmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempchange.sc:duration.sc 
                   + (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTmass5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# sdT+npp:duration ############################

if(fitmod == 'modsdTnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempchange.sc:duration.sc + 
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTnpp5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# sdT+seas:duration ############################

if(fitmod == 'modsdTseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTseas5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# sdT+microclim:duration ############################

if(fitmod == 'modsdTmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTmicroclim5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# sdT+human:duration ############################

if(fitmod == 'modsdThuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdThuman5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# TdT+SINGLE FACTORS ############################
# T+dT+REALM:duration #########################

if(fitmod == 'modTdTRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTRealm10Jtu'){ 
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTRealm20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTRealm5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTRealm10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTRealm20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTRealm5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTRealm10Horn'){ 
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTRealm20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc +
                       tempchange_abs.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T+dT:tsign:duration #########################

if(fitmod == 'modTdTtsign5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ tempave_metab.sc:duration.sc +
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc + 
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTtsign10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ tempave_metab.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTtsign20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ tempave_metab.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTtsign5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ tempave_metab.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTtsign10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ tempave_metab.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTtsign20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ tempave_metab.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTtsign5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ tempave_metab.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTtsign10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ tempave_metab.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTtsign20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ tempave_metab.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T+dT+mass:duration #########################

if(fitmod == 'modTdTmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmass10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmass20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTmass5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmass10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmass20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTmass5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmass10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmass20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+dT+npp:duration #########################

if(fitmod == 'modTdTnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTnpp10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTnpp20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTnpp5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTnpp10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTnpp20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+dT+seas:duration #########################

if(fitmod == 'modTdTseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTseas10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTseas20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTseas5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTseas10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTseas20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T+dT+microclim:duration #########################

if(fitmod == 'modTdTmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmicroclim10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmicroclim20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTmicroclim5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTmicroclim10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTmicroclim20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+dT+human:duration #########################

if(fitmod == 'modTdThuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdThuman10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdThuman20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdThuman5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdThuman10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdThuman20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# TsdT+SINGLE FACTORS ############################
# T+sdT+REALM:duration #########################

if(fitmod == 'modTsdTRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc +
                       tempchange.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTsdTRealm5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc +
                       tempchange.sc:duration.sc + 
                       tempave_metab.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+sdT+mass:duration #########################

if(fitmod == 'modTsdTmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTsdTmass5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+sdT+npp:duration #########################

if(fitmod == 'modTsdTnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTsdTnpp5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       npp.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+sdT+seas:duration #########################

if(fitmod == 'modTsdTseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTsdTseas5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       seas.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T+sdT+microclim:duration #########################

if(fitmod == 'modTsdTmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTsdTmicroclim5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       microclim.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T+sdT+human:duration #########################

if(fitmod == 'modTsdThuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTsdThuman5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       human_bowler.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T INTERACTION WITH SINGLE FACTORS ############################
# T:REALM:duration ############################

if(fitmod == 'modTxRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + REALM:tempave_metab.sc:duration.sc + (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTxRealm10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + REALM:tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTxRealm20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + REALM:tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'modTxRealm5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + REALM:tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTxRealm10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + REALM:tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'modTxRealm20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + REALM:tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# T:dT:duration ############################

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

# T:sdT:duration ############################

if(fitmod == 'modTsdTT5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange.sc:duration.sc 
                   + tempave_metab.sc:duration.sc 
                   + tempave_metab.sc:tempchange.sc:duration.sc 
                   + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTsdTT5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange.sc:duration.sc 
                   + tempave_metab.sc:duration.sc 
                   + tempave_metab.sc:tempchange.sc:duration.sc 
                   + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}


# T:mass:duration ############################
if(fitmod == 'modTxmass5Jtu'){
    print(paste(nrow(trends5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:mass.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTxmass5Horn'){
    print(paste(nrow(trends5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:mass.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# T:npp:duration ############################
if(fitmod == 'modTxnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:npp.sc:duration.sc + (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTxnpp5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:npp.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# T:seas:duration ############################
if(fitmod == 'modTxseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:seas.sc:duration.sc + (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTxseas5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:seas.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# T:microclim:duration ############################
if(fitmod == 'modTxmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:microclim.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTxmicroclim5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:microclim.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# T:human:duration ############################
if(fitmod == 'modTxhuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:human_bowler.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modTxhuman5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:human_bowler.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# dT INTERACTION WITH SINGLE FACTORS ############################
# dT:REALM:duration ############################

if(fitmod == 'moddTxRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + REALM:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTxRealm10Jtu'){
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + REALM:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTxRealm20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + REALM:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


if(fitmod == 'moddTxRealm5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + REALM:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTxRealm10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + REALM:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTxRealm20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + REALM:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# dT:tsign:duration #########################

if(fitmod == 'moddTtsign5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc + 
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTtsign10Jtu'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    print(paste(sum(i10Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTtsign20Jtu'){
    print(paste(sum(i20Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTtsign5Jne'){
    print(paste(sum(i5Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTtsign10Jne'){
    print(paste(sum(i10Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTtsign20Jne'){
    print(paste(sum(i20Jne), 'data points'))
    mod <- glmmTMB(Jne.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Jne,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTtsign5Horn'){ # Horn
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTtsign10Horn'){
    print(paste(sum(i10Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}
if(fitmod == 'moddTtsign20Horn'){
    print(paste(sum(i20Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# dT:mass:duration ############################
if(fitmod == 'moddTxmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:mass.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTxmass5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:mass.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# dT:npp:duration ############################
if(fitmod == 'moddTxnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempchange_abs.sc:duration.sc + 
                       tempchange_abs.sc:npp.sc:duration.sc + (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTxnpp5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:npp.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# dT:seas:duration ############################
if(fitmod == 'moddTxseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:seas.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTxseas5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:seas.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# dT:microclim:duration ############################
if(fitmod == 'moddTxmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:microclim.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTxmicroclim5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:microclim.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# dT:human:duration ############################
if(fitmod == 'moddTxhuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:human_bowler.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'moddTxhuman5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange_abs.sc:duration.sc + tempchange_abs.sc:human_bowler.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# sdT INTERACTION WITH SINGLE FACTORS ############################
# sdT:REALM:duration ############################

if(fitmod == 'modsdTxRealm5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ REALM:duration.sc + REALM:tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Jtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTxRealm5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ REALM:duration.sc + REALM:tempchange.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5Horn,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# sdT:mass:duration ############################
if(fitmod == 'modsdTxmass5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + mass.sc:duration.sc 
                   + tempchange.sc:duration.sc + tempchange.sc:mass.sc:duration.sc 
                   + (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTxmass5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + mass.sc:duration.sc + tempchange.sc:duration.sc + tempchange.sc:mass.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# sdT:npp:duration ############################
if(fitmod == 'modsdTxnpp5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + npp.sc:duration.sc + tempchange.sc:duration.sc + 
                       tempchange.sc:npp.sc:duration.sc + (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTxnpp5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + npp.sc:duration.sc + tempchange.sc:duration.sc + tempchange.sc:npp.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# sdT:seas:duration ############################
if(fitmod == 'modsdTxseas5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + seas.sc:duration.sc + tempchange.sc:duration.sc 
                   + tempchange.sc:seas.sc:duration.sc + (duration.sc||STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTxseas5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + seas.sc:duration.sc + tempchange.sc:duration.sc + tempchange.sc:seas.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# sdT:microclim:duration ############################
if(fitmod == 'modsdTxmicroclim5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + microclim.sc:duration.sc + tempchange.sc:duration.sc + tempchange.sc:microclim.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTxmicroclim5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + microclim.sc:duration.sc + tempchange.sc:duration.sc + tempchange.sc:microclim.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}

# sdT:human:duration ############################
if(fitmod == 'modsdTxhuman5Jtu'){
    print(paste(sum(i5Jtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange.sc:duration.sc + tempchange.sc:human_bowler.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Jtu, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

if(fitmod == 'modsdTxhuman5Horn'){
    print(paste(sum(i5Horn), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + human_bowler.sc:duration.sc + tempchange.sc:duration.sc + tempchange.sc:human_bowler.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5Horn, family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    MATCHMOD <- TRUE
}


# T:dT/sdT:REALM INTERACTION WITH SINGLE FACTORS ########################
# T:dT:REALM:duration like Antao #########################

# use tempchange like Antao
if(fitmod == 'modTdTTRealmAllJtu'){
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration +
                       REALM:duration +
                       REALM:tempchange.sc:duration + 
                       REALM:tempave_metab.sc:duration +
                       REALM:tempave_metab.sc:tempchange.sc:duration +
                       (duration|STUDY_ID/rarefyID), 
                   data = trendsall[iallJtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# use tempchange like Antao. scale duration
if(fitmod == 'modTdTTRealmDurscAllJtu'){
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc +
                       REALM:duration.sc +
                       REALM:tempchange.sc:duration.sc + 
                       REALM:tempave_metab.sc:duration.sc +
                       REALM:tempave_metab.sc:tempchange.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trendsall[iallJtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}


# T:sdT:REALM:duration #########################
# use tempchange_abs (otherwise like Antao)
if(fitmod == 'modTsdTTRealmAllJtu'){
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration +
                       REALM:duration +
                       REALM:tempchange_abs.sc:duration + 
                       REALM:tempave_metab.sc:duration +
                       REALM:tempave_metab.sc:tempchange_abs.sc:duration +
                       (duration|STUDY_ID/rarefyID), 
                   data = trendsall[iallJtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}

# use tempchange_abs
if(fitmod == 'modTsdTTRealmDurscAllJtu'){
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc +
                       REALM:duration.sc +
                       REALM:tempchange_abs.sc:duration.sc + 
                       REALM:tempave_metab.sc:duration.sc +
                       REALM:tempave_metab.sc:tempchange_abs.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trendsall[iallJtu,], family = beta_family(link = 'logit'), 
                   dispformula = ~REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    MATCHMOD <- TRUE
}



# T+dT:REALM:duration ANTAO-STYLE DATASET #########################

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

# DREDGE FULL MODELS #########################
if(fitmod == 'dredge5Jtu'){
    require(MuMIn)
    print(paste(nrow(trends5Jtu), 'data points'))
    modfull <- glmmTMB(formula(paste(fixed, '+(duration.sc||STUDY_ID/rarefyID)')), 
                       data = trends5Jtu, 
                       family = beta_family(link = 'logit'), 
                       dispformula = ~nspp.sc+REALM, 
                       control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(modfull))
    print(Sys.time())
    print('Starting dredge')
    mod <- lapply(dredge(modfull, fixed = c('disp(nspp.sc)', 'disp(REALM)'),
                  rank = 'AIC',
                  evaluate = FALSE), eval)
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
    if(!grepl('dredge', fitmod)){
        print(performance::r2(mod)) # run if not a dredge object
    }
}

print(warnings())

print(paste('Ended', Sys.time(), sep =''))
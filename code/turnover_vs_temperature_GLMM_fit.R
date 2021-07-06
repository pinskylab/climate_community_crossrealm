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

## Compare variance structures with duration.sc ############################
#fixed <- 'Jtu.sc ~ tempchange_abs.sc*REALM + tempchange_abs.sc*tempave_metab.sc + tempchange_abs.sc*duration.sc'
fixed <- 'Jtu.sc ~ duration.sc + REALM:duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + mass.sc:duration.sc + nspp.sc:duration.sc + seas.sc:duration.sc + microclim.sc:duration.sc + npp.sc:duration.sc + human_bowler.sc:REALM2:duration.sc'
i3 <- trends3[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i5 <- trends5[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i10 <- trends10[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i20 <- trends20[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

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
i3 <- trends3[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i5 <- trends5[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i10 <- trends10[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
i20 <- trends20[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

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



# temperature, temperature change models with duration interactions ############################

if(fitmod == 'modTdT3'){
    print(paste(sum(i3), 'data points'))
    modTdT3 <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends3[i3,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modTdT3)
    saveRDS(modTdT3, file = 'temp/modTdT3.rds')
    print('saved modTdT3.rds')
    print(Sys.time())
    print(performance::r2(modTdT3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdT5'){
    print(paste(sum(i5), 'data points'))
    modTdT5 <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends5[i5,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modTdT5)
    saveRDS(modTdT5, file = 'temp/modTdT5.rds')
    print('saved modTdT5.rds')
    print(Sys.time())
    print(performance::r2(modTdT5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdT10'){
    print(paste(sum(i10), 'data points'))
    modTdT10 <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends10[i10,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modTdT10)
    saveRDS(modTdT10, file = 'temp/modTdT10.rds')
    print('saved modTdT10.rds')
    print(Sys.time())
    print(performance::r2(modTdT10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdT20'){
    print(paste(sum(i20), 'data points'))
    modTdT20 <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends20[i20,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modTdT20)
    saveRDS(modTdT20, file = 'temp/modTdT20.rds')
    print('saved modTdT20.rds')
    print(Sys.time())
    print(performance::r2(modTdT20))
    MATCHMOD <- TRUE
}

# temperature:temperature change models with duration interactions ############################

if(fitmod == 'modTdTT3'){
    print(paste(sum(i3), 'data points'))
    modTdTT3 <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                       data = trends3[i3,], family = beta_family(link = 'logit'), 
                       dispformula = ~nspp.sc+REALM, 
                       control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modTdTT3)
    saveRDS(modTdTT3, file = 'temp/modTdTT3.rds')
    print('saved modTdTT3.rds')
    print(Sys.time())
    print(performance::r2(modTdTT3))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT5'){
    print(paste(sum(i5), 'data points'))
    modTdTT5 <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                       data = trends5[i5,], family = beta_family(link = 'logit'), 
                       dispformula = ~nspp.sc+REALM, 
                       control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modTdTT5)
    saveRDS(modTdTT5, file = 'temp/modTdTT5.rds')
    print('saved modTdTT5.rds')
    print(Sys.time())
    print(performance::r2(modTdTT5))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT10'){
    print(paste(sum(i10), 'data points'))
    modTdTT10 <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends10[i10,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modTdTT10)
    saveRDS(modTdTT10, file = 'temp/modTdTT10.rds')
    print('saved modTdTT10.rds')
    print(Sys.time())
    print(performance::r2(modTdTT10))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT20'){
    print(paste(sum(i20), 'data points'))
    modTdTT20 <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends20[i20,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modTdTT20)
    saveRDS(modTdTT20, file = 'temp/modTdTT20.rds')
    print('saved modTdTT20.rds')
    print(Sys.time())
    print(performance::r2(modTdTT20))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTT3Jne'){
    i3 <- trends3[, complete.cases(Jne.sc, tempchange_abs.sc, tempave_metab.sc, duration.sc, REALM, nspp.sc)]
    print(paste(sum(i3), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends3[i3,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTT3Jne.rds')
    print('saved modTdTT3Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT5Jne'){
    i5 <- trends5[, complete.cases(Jne.sc, tempchange_abs.sc, tempave_metab.sc, duration.sc, REALM, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTT5Jne.rds')
    print('saved modTdTT5Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT10Jne'){
    i10 <- trends10[, complete.cases(Jne.sc, tempchange_abs.sc, tempave_metab.sc, duration.sc, REALM, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTT10Jne.rds')
    print('saved modTdTT10Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT20Jne'){
    i20 <- trends20[, complete.cases(Jne.sc, tempchange_abs.sc, tempave_metab.sc, duration.sc, REALM, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTT20Jne.rds')
    print('saved modTdTT20Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTT3Horn'){
    i3 <- trends3[, complete.cases(Horn.sc, tempchange_abs.sc, tempave_metab.sc, duration.sc, REALM, nspp.sc)]
    print(paste(sum(i3), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends3[i3,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTT3Horn.rds')
    print('saved modTdTT3Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT5Horn'){
    i5 <- trends5[, complete.cases(Horn.sc, tempchange_abs.sc, tempave_metab.sc, duration.sc, REALM, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                        data = trends5[i5,], family = beta_family(link = 'logit'), 
                        dispformula = ~nspp.sc+REALM, 
                        control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTT5Horn.rds')
    print('saved modTdTT5Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT10Horn'){
    i10 <- trends10[, complete.cases(Horn.sc, tempchange_abs.sc, tempave_metab.sc, duration.sc, REALM, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                         data = trends10[i10,], family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc+REALM, 
                         control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTT10Horn.rds')
    print('saved modTdTT10Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTT20Horn'){
    i20 <- trends20[, complete.cases(Horn.sc, tempchange_abs.sc, tempave_metab.sc, duration.sc, REALM, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + tempave_metab.sc:tempchange_abs.sc:duration.sc + (duration.sc|STUDY_ID/rarefyID), 
                         data = trends20[i20,], family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc+REALM, 
                         control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTT20Horn.rds')
    print('saved modTdTT20Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

# temp:tempchange:REALM with duration interaction #########################

if(fitmod == 'modTdTTRealm5'){
    i5 <- trends5[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm5.rds')
    print('saved modTdTTRealm5.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm10'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    i10 <- trends10[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm10.rds')
    print('saved modTdTTRealm10.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm20'){
    i20 <- trends20[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm20.rds')
    print('saved modTdTTRealm20.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTRealm5Jne'){
    i5 <- trends5[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm5Jne.rds')
    print('saved modTdTTRealm5Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm10Jne'){
    i10 <- trends10[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm10Jne.rds')
    print('saved modTdTTRealm10Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm20Jne'){
    i20 <- trends20[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm20Jne.rds')
    print('saved modTdTTRealm20Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTRealm5Horn'){ # Horn
    i5 <- trends5[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm5Horn.rds')
    print('saved modTdTTRealm5Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm10Horn'){
    i10 <- trends10[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm10Horn.rds')
    print('saved modTdTTRealm10Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTRealm20Horn'){
    i20 <- trends20[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       REALM:duration.sc +
                       tempchange_abs.sc:REALM:duration.sc + 
                       tempave_metab.sc:REALM:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTRealm20Horn.rds')
    print('saved modTdTTRealm20Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm5.rds')
    print('saved modAntaoTdTTRealm5.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm10.rds')
    print('saved modAntaoTdTTRealm10.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm5Jne.rds')
    print('saved modAntaoTdTTRealm5Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm10Jne.rds')
    print('saved modAntaoTdTTRealm10Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm20Jne.rds')
    print('saved modAntaoTdTTRealm20Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm5Horn.rds')
    print('saved modAntaoTdTTRealm5Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm10Horn.rds')
    print('saved modAntaoTdTTRealm10Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    print(summary(mod))
    saveRDS(mod, file = 'temp/modAntaoTdTTRealm20Horn.rds')
    print('saved modAntaoTdTTRealm20Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

# temp:tempchange:tsign with duration interaction #########################

if(fitmod == 'modTdTTtsign5'){
    i5 <- trends5[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc + 
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign5.rds')
    print('saved modTdTTtsign5.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign10'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    i10 <- trends10[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign10.rds')
    print('saved modTdTTtsign10.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign20'){
    i20 <- trends20[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign20.rds')
    print('saved modTdTTtsign20.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTtsign5Jne'){
    i5 <- trends5[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign5Jne.rds')
    print('saved modTdTTtsign5Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign10Jne'){
    i10 <- trends10[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign10Jne.rds')
    print('saved modTdTTtsign10Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign20Jne'){
    i20 <- trends20[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign20Jne.rds')
    print('saved modTdTTtsign20Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTtsign5Horn'){ # Horn
    i5 <- trends5[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign5Horn.rds')
    print('saved modTdTTtsign5Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign10Horn'){
    i10 <- trends10[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign10Horn.rds')
    print('saved modTdTTtsign10Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTtsign20Horn'){
    i20 <- trends20[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       tsign:duration.sc +
                       tempchange_abs.sc:tsign:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTtsign20Horn.rds')
    print('saved modTdTTtsign20Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}


# temp:tempchange +mass with duration interaction #########################

if(fitmod == 'modTdTTmass5'){
    i5 <- trends5[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass5.rds')
    print('saved modTdTTmass5.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass10'){ # trims out freshwater and >60 or <23.5 deg lat, like  et al. 2020
    i10 <- trends10[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass10.rds')
    print('saved modTdTTmass10.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass20'){
    i20 <- trends20[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass20.rds')
    print('saved modTdTTmass20.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTmass5Jne'){
    i5 <- trends5[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass5Jne.rds')
    print('saved modTdTTmass5Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass10Jne'){
    i10 <- trends10[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass10Jne.rds')
    print('saved modTdTTmass10Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass20Jne'){
    i20 <- trends20[, complete.cases(Jne.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Jne.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass20Jne.rds')
    print('saved modTdTTmass20Jne.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}

if(fitmod == 'modTdTTmass5Horn'){ # Horn
    i5 <- trends5[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i5), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends5[i5,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass5Horn.rds')
    print('saved modTdTTmass5Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass10Horn'){
    i10 <- trends10[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i10), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i10,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass10Horn.rds')
    print('saved modTdTTmass10Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modTdTTmass20Horn'){
    i20 <- trends20[, complete.cases(Horn.sc, tempchange_abs.sc, REALM, tempave_metab.sc, durationlog.sc, nspp.sc)]
    print(paste(sum(i20), 'data points'))
    mod <- glmmTMB(Horn.sc ~ duration.sc + tempchange_abs.sc:duration.sc + tempave_metab.sc:duration.sc + 
                       tempave_metab.sc:tempchange_abs.sc:duration.sc + 
                       mass.sc:duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i20,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    print(summary(mod))
    saveRDS(mod, file = 'temp/modTdTTmass20Horn.rds')
    print('saved modTdTTmass20Horn.rds')
    print(Sys.time())
    print(performance::r2(mod))
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
    i <- trends5[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]

    print(paste(sum(i), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + REALM + nspp.sc +
                                                   tempchange_abs.sc:duration.sc +
                                                   REALM : duration.sc + 
                                                   tempave_metab.sc : duration.sc + 
                                                   seas.sc : duration.sc +
                                                   microclim.sc : duration.sc +
                                                   npp.sc : duration.sc +
                                                   nspp.sc : duration.sc +
                                                   mass.sc : duration.sc +
                                                   human_bowler.sc:REALM2 : duration.sc +
                                                   (duration.sc|STUDY_ID/rarefyID), 
                                               data = trends5[i,], 
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
if(fitmod == 'modDurIntTdTSeMiNPNspMaHu10yrJtu'){
    i <- trends10[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
    
    print(paste(sum(i), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + REALM + nspp.sc +
                       tempchange_abs.sc:duration.sc +
                       REALM : duration.sc + 
                       tempave_metab.sc : duration.sc + 
                       seas.sc : duration.sc +
                       microclim.sc : duration.sc +
                       npp.sc : duration.sc +
                       nspp.sc : duration.sc +
                       mass.sc : duration.sc +
                       human_bowler.sc:REALM2 : duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends10[i,], 
                   family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM,
                   control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    outfile = 'temp/modDurIntTdTSeMiNPNspMaHu10yrJtu.rds'
    saveRDS(mod, file = outfile)
    print(paste('saved', outfile))
    print(Sys.time())
    print(performance::r2(mod))
    MATCHMOD <- TRUE
}
if(fitmod == 'modDurIntTdTSeMiNPNspMaHu20yrJtu'){
    i <- trends20[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, nspp.sc, seas.sc, microclim.sc, npp.sc, human_bowler.sc)]
    
    print(paste(sum(i), 'data points'))
    mod <- glmmTMB(Jtu.sc ~ duration.sc + REALM + nspp.sc +
                       tempchange_abs.sc:duration.sc +
                       REALM : duration.sc + 
                       tempave_metab.sc : duration.sc + 
                       seas.sc : duration.sc +
                       microclim.sc : duration.sc +
                       npp.sc : duration.sc +
                       nspp.sc : duration.sc +
                       mass.sc : duration.sc +
                       human_bowler.sc:REALM2 : duration.sc +
                       (duration.sc|STUDY_ID/rarefyID), 
                   data = trends20[i,], 
                   family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc+REALM,
                   control = glmmTMBControl(profile=TRUE))
    print(summary(mod))
    outfile = 'temp/modDurIntTdTSeMiNPNspMaHu20yrJtu.rds'
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
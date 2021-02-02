#!/usr/bin/env Rscript

# Script to fit glmmTMB models
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_vs_temperature_GLMM.R modRFgauss > logs/turnover_vs_temperature_GLMMmodRFgauss.Rout &
# (this works if code is executable, e.g., chmod u+x turnover_vs_temperature_GLMM.R)
# (otherwise using nohup Rscript ...)

################
# Read command line arguments

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args)<1) stop("Have to specify a model to fit", call.=FALSE)
if (length(args)>1) stop("Have to specify only 1 model to fit", call.=FALSE)
fitmod <- args[1]
MATCHMOD <- FALSE # indicator to check if the argument matched a model name

###################################
# print basic info about the job

print('This is process #')
print(Sys.getpid())
print(Sys.time())



##############
# Prep

library(methods)
library(data.table) # for handling large datasets
library(glmmTMB) # for ME models

############
# load data

# Turnover and covariates assembled by turnover_vs_temperature_prep.Rmd
trends <- fread('output/turnover_w_covariates.csv.gz')

###############
# trim data

# trim to only data with some temperature change
# important since sign of temperature change is a variable
trends[tempchange == 0, .N] # number to remove
trends <- trends[tempchange != 0, ] # also removes any NA values

#######################
# set up useful vars

# set realm order
trends[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = FALSE)]

# set up sign of temperature change
trends[, tsign := factor(sign(tempchange))]

# realm that combines Terrestrial and Freshwater, for interacting with human impact
trends[, REALM2 := REALM]
levels(trends$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")

# group Marine invertebrates/plants in with All
trends[, taxa_mod2 := taxa_mod]
trends[taxa_mod == 'Marine invertebrates/plants', taxa_mod2 := 'All']

# calculate duration
trends[, duration := year2 - year1]

#add a comparison id
trends[, compID := paste0(rarefyID, '_', year1, '_', year2)]


#######################
## Transformations

### Adjust response away from 0-1
# transformation for 2 categories. Eq. 1 in Douma & Weedon 2019 MEE
transform01 <- function(x) (x * (length(x) - 1) + 0.5) / (length(x))

trends[, Jtu.sc := transform01(Jtu)]
trends[, Jbeta.sc := transform01(Jbeta)]
trends[, Horn.sc := transform01(Horn)]


### Log-transform some variables, then center and scale. 
trends[, tempave.sc := scale(tempave)]
trends[, tempave_metab.sc := scale(tempave_metab)]
trends[, seas.sc := scale(seas)]
trends[, microclim.sc := scale(log(microclim))]
trends[, tempchange.sc := scale(tempchange, center = FALSE)] # do not center
trends[, tempchange_abs.sc := scale(abs(tempchange), center = FALSE)] # do not center, so that 0 is still 0 temperature change
trends[, mass.sc := scale(log(mass_mean_weight))]
trends[, speed.sc := scale(log(speed_mean_weight+1))]
trends[, lifespan.sc := scale(log(lifespan_mean_weight))]
trends[, consumerfrac.sc := scale(consfrac)]
trends[, endothermfrac.sc := scale(endofrac)]
trends[, nspp.sc := scale(log(Nspp))]
trends[, thermal_bias.sc := scale(thermal_bias)]
trends[, npp.sc := scale(log(npp))]
trends[, veg.sc := scale(log(veg+1))]
trends[, duration.sc := scale(log(duration))]
trends[, human_bowler.sc := scale(log(human_bowler+1)), by = REALM2] # separate scaling by realm
trends[REALM2 == 'TerrFresh', human_footprint.sc := scale(log(human_venter+1))]
trends[REALM2 == 'Marine', human_footprint.sc := scale(log(human_halpern))]

###############
# Models
###############

#########################################
## Models to compare variance structures
fixed <- 'Jtu.sc ~ tempchange_abs.sc*REALM + tempchange_abs.sc*tempave_metab.sc + tempchange_abs.sc*duration.sc'
i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc)]

if(fitmod == 'modRFgauss'){
    modRFgauss <- glmmTMB(formula(fixed), data = trends[i,], family = gaussian(), control = glmmTMBControl(profile=TRUE))
    summary(modRFgauss)
    saveRDS(modRFgauss, file = 'temp/modRFgauss.rds')
    print('saved modRFgauss.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFbeta'){
    modRFbeta <- glmmTMB(formula(fixed), data = trends[i,], family = beta_family(link = 'logit'), control = glmmTMBControl(profile=TRUE)) # add beta errors
    summary(modRFbeta)
    saveRDS(modRFbeta, file = 'temp/modRFbeta.rds')
    print('saved modRFbeta.rds')
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFrID'){
    modRFrID <- glmmTMB(formula(paste(fixed, '+(1|rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random effects
    summary(modRFrID)
    saveRDS(modRFrID, file = 'temp/modRFrID.rds')
    print('saved modRDrID.rds')
    MATCHMOD <- TRUE
}
if(fitmod == 'modRFnestedRE'){
    modRFnestedRE <- glmmTMB(formula(paste(fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRFnestedRE)
    saveRDS(modRFnestedRE, file = 'temp/modRFnestedRE.rds')
    print('saved modRFnestedRE.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFslopeRE'){
    modRFslopeRE <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE)
    saveRDS(modRFslopeRE, file = 'temp/modRFslopeRE.rds')
    print('saved modRFslopeRE.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdisp2lev'){
    modRFdisp2lev <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc, 
                         control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdisp2lev)
    saveRDS(modRFdisp2lev, file = 'temp/modRFdisp2lev.rds')
    print('saved modRFdisp2lev.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modRFdisp'){
    modRFdisp <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdisp)
    saveRDS(modRFdisp, file = 'temp/modRFdisp.rds')
    print('saved modRFdisp.rds')
    MATCHMOD <- TRUE
}

##################################
## temperature-only models (all years)
if(fitmod == 'modonlyTchangeJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange)]
    modonlyTchangeJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                             data = trends[i,], 
                             family = beta_family(link = 'logit'), 
                             dispformula = ~nspp.sc)
    summary(modonlyTchangeJtu)
    saveRDS(modonlyTchangeJtu, file = 'temp/modonlyTchangeJtu.rds')
    print('saved modonlyTchangeJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchangeJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange)]
    modonlyTchangeJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                               data = trends[i,], 
                               family = beta_family(link = 'logit'), 
                               dispformula = ~nspp.sc)
    summary(modonlyTchangeJbeta)
    saveRDS(modonlyTchangeJbeta, file = 'temp/modonlyTchangeJbeta.rds')
    print('saved modonlyTchangeJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchangeHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange)]
    modonlyTchangeHorn <- glmmTMB(Horn.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                              data = trends[i,], 
                              family = beta_family(link = 'logit'), 
                              dispformula = ~nspp.sc)
    summary(modonlyTchangeHorn)
    saveRDS(modonlyTchangeHorn, file = 'temp/modonlyTchangeHorn.rds')
    print('saved modonlyTchangeHorn.rds')
    MATCHMOD <- TRUE
}


##################################
## temperature-only models (1 year)
if(fitmod == 'modonlyTchange1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange) & duration == 1]
    modonlyTchange1yrJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                 data = trends[i,], 
                                 family = beta_family(link = 'logit'), 
                                 dispformula = ~nspp.sc)
    summary(modonlyTchange1yrJtu)
    saveRDS(modonlyTchange1yrJtu, file = 'temp/modonlyTchange1yrJtu.rds')
    print('saved modonlyTchange1yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchange1yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange) & duration == 1]
    modonlyTchange1yrJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                   data = trends[i,], 
                                   family = beta_family(link = 'logit'), 
                                   dispformula = ~nspp.sc)
    summary(modonlyTchange1yrJbeta)
    saveRDS(modonlyTchange1yrJbeta, file = 'temp/modonlyTchange1yrJbeta.rds')
    print('saved modonlyTchange1yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchange1yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange) & duration == 1]
    modonlyTchange1yrHorn <- glmmTMB(Horn.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                  data = trends[i,], 
                                  family = beta_family(link = 'logit'), 
                                  dispformula = ~nspp.sc)
    summary(modonlyTchange1yrHorn)
    saveRDS(modonlyTchange1yrHorn, file = 'temp/modonlyTchange1yrHorn.rds')
    print('saved modonlyTchange1yrHorn.rds')
    MATCHMOD <- TRUE
}

#####################################
## temperature-only models (10 year)
if(fitmod == 'modonlyTchange10yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange) & duration == 10]
    modonlyTchange10yrJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                    data = trends[i,], 
                                    family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc)
    summary(modonlyTchange10yrJtu)
    saveRDS(modonlyTchange10yrJtu, file = 'temp/modonlyTchange10yrJtu.rds')
    print('saved modonlyTchange10yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchange10yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange) & duration == 10]
    modonlyTchange10yrJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                      data = trends[i,], 
                                      family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc)
    summary(modonlyTchange10yrJbeta)
    saveRDS(modonlyTchange10yrJbeta, file = 'temp/modonlyTchange10yrJbeta.rds')
    print('saved modonlyTchange10yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modonlyTchange10yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange) & duration == 10]
    modonlyTchange10yrHorn <- glmmTMB(Horn.sc ~ abs(tempchange) + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                     data = trends[i,], 
                                     family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc)
    summary(modonlyTchange10yrHorn)
    saveRDS(modonlyTchange10yrHorn, file = 'temp/modonlyTchange10yrHorn.rds')
    print('saved modonlyTchange10yrHorn.rds')
    MATCHMOD <- TRUE
}

#################################
## temperature-duration models
if(fitmod == 'modTDJtu'){
    i <- trends[, complete.cases(Jtu.sc, REALM, tempchange)]
    modTDJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * duration + (tempchange_abs.sc|STUDY_ID/rarefyID),
                    data = trends[i,], 
                    family = beta_family(link = 'logit'), 
                    dispformula = ~nspp.sc)
    summary(modTDJtu)
    saveRDS(modTDJtu, file = 'temp/modTDJtu.rds')
    print('saved modTDJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTDJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, REALM, tempchange)]
    modTDJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange)*duration + (tempchange_abs.sc|STUDY_ID/rarefyID),
                      data = trends[i,], 
                      family = beta_family(link = 'logit'), 
                      dispformula = ~nspp.sc)
    summary(modTDJbeta)
    saveRDS(modTDJbeta, file = 'temp/modTDJbeta.rds')
    print('saved modTDJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTDHorn'){
    i <- trends[, complete.cases(Horn.sc, REALM, tempchange)]
    modTDHorn <- glmmTMB(Horn.sc ~ abs(tempchange)*duration + (tempchange_abs.sc|STUDY_ID/rarefyID),
                     data = trends[i,], 
                     family = beta_family(link = 'logit'), 
                     dispformula = ~nspp.sc)
    summary(modTDHorn)
    saveRDS(modTDHorn, file = 'temp/modTDHorn.rds')
    print('saved modTDHorn.rds')
    MATCHMOD <- TRUE
}

###########################################
## temperature-duration models by REALM
if(fitmod == 'modTDrealmJtu'){
    i <- trends[, complete.cases(Jtu.sc, REALM, tempchange)]
    modTDrealmJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * duration + abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                         data = trends[i,], 
                         family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc)
    summary(modTDrealmJtu)
    saveRDS(modTDrealmJtu, file = 'temp/modTDrealmJtu.rds')
    print('saved modTDrealmJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTDrealmJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, REALM, tempchange)]
    modTDrealmJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange)*duration + abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                           data = trends[i,], 
                           family = beta_family(link = 'logit'), 
                           dispformula = ~nspp.sc)
    summary(modTDrealmJbeta)
    saveRDS(modTDrealmJbeta, file = 'temp/modTDrealmJbeta.rds')
    print('saved modTDrealmJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTDrealmHorn'){
    i <- trends[, complete.cases(Horn.sc, REALM, tempchange)]
    modTDrealmHorn <- glmmTMB(Horn.sc ~ abs(tempchange)*duration + abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                          data = trends[i,], 
                          family = beta_family(link = 'logit'), 
                          dispformula = ~nspp.sc)
    summary(modTDrealmHorn)
    saveRDS(modTDrealmHorn, file = 'temp/modTDrealmHorn.rds')
    print('saved modTDrealmHorn.rds')
    MATCHMOD <- TRUE
}

##############################################
## temperature-only models by REALM (1 year)
if(fitmod == 'modTrealm1yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange, REALM) & duration == 1]
    modTrealm1yrJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                    data = trends[i,], 
                                    family = beta_family(link = 'logit'), 
                                    dispformula = ~nspp.sc)
    summary(modTrealm1yrJtu)
    saveRDS(modTrealm1yrJtu, file = 'temp/modTrealm1yrJtu.rds')
    print('saved modTrealm1yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTrealm1yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange, REALM) & duration == 1]
    modTrealm1yrJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                      data = trends[i,], 
                                      family = beta_family(link = 'logit'), 
                                      dispformula = ~nspp.sc)
    summary(modTrealm1yrJbeta)
    saveRDS(modTrealm1yrJbeta, file = 'temp/modTrealm1yrJbeta.rds')
    print('saved modTrealm1yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTrealm1yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange, REALM) & duration == 1]
    modTrealm1yrHorn <- glmmTMB(Horn.sc ~ abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                     data = trends[i,], 
                                     family = beta_family(link = 'logit'), 
                                     dispformula = ~nspp.sc)
    summary(modTrealm1yrHorn)
    saveRDS(modTrealm1yrHorn, file = 'temp/modTrealm1yrHorn.rds')
    print('saved modTrealm1yrHorn.rds')
    MATCHMOD <- TRUE
}


##############################################
## temperature-only models by REALM (10 year)
if(fitmod == 'modTrealm10yrJtu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange, REALM) & duration == 10]
    modTrealm10yrJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                               data = trends[i,], 
                               family = beta_family(link = 'logit'), 
                               dispformula = ~nspp.sc)
    summary(modTrealm10yrJtu)
    saveRDS(modTrealm10yrJtu, file = 'temp/modTrealm10yrJtu.rds')
    print('saved modTrealm10yrJtu.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTrealm10yrJbeta'){
    i <- trends[, complete.cases(Jbeta.sc, tempchange, REALM) & duration == 10]
    modTrealm10yrJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                 data = trends[i,], 
                                 family = beta_family(link = 'logit'), 
                                 dispformula = ~nspp.sc)
    summary(modTrealm10yrJbeta)
    saveRDS(modTrealm10yrJbeta, file = 'temp/modTrealm10yrJbeta.rds')
    print('saved modTrealm10yrJbeta.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modTrealm10yrHorn'){
    i <- trends[, complete.cases(Horn.sc, tempchange, REALM) & duration == 10]
    modTrealm10yrHorn <- glmmTMB(Horn.sc ~ abs(tempchange) * REALM + (tempchange_abs.sc|STUDY_ID/rarefyID),
                                data = trends[i,], 
                                family = beta_family(link = 'logit'), 
                                dispformula = ~nspp.sc)
    summary(modTrealm10yrHorn)
    saveRDS(modTrealm10yrHorn, file = 'temp/modTrealm10yrHorn.rds')
    print('saved modTrealm10yrHorn.rds')
    MATCHMOD <- TRUE
}

###################################
# full models


if(fitmod == 'modFullmass'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc)]

    modFullmass <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                            tempchange_abs.sc*tempave_metab.sc + 
                            tempchange_abs.sc*duration.sc +
                            tempchange_abs.sc*mass.sc +
                            (tempchange_abs.sc|STUDY_ID/rarefyID), 
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
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, endothermfrac.sc)]
    
    modFullendo <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                               tempchange_abs.sc*tempave_metab.sc + 
                               tempchange_abs.sc*duration.sc +
                               tempchange_abs.sc*endothermfrac.sc +
                               (tempchange_abs.sc|STUDY_ID/rarefyID), 
                           data = trends[i,], 
                           family = beta_family(link = 'logit'), 
                           dispformula = ~nspp.sc, 
                           control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modFullendo)
    saveRDS(modFullendo, file = 'temp/modFullendo.rds')
    print('saved modFullendo.rds')
    MATCHMOD <- TRUE
}

if(fitmod == 'modFullMaEnMiNPHu'){
    i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc, mass.sc, endothermfrac.sc,
                                 microclim.sc, npp.sc, human_bowler.sc)]
    
    modFullendo <- glmmTMB(Jtu.sc ~ tempchange_abs.sc*REALM + 
                               tempchange_abs.sc*tempave_metab.sc + 
                               tempchange_abs.sc*duration.sc +
                               tempchange_abs.sc*mass.sc +
                               tempchange_abs.sc*endothermfrac.sc +
                               tempchange_abs.sc*microclim.sc +
                               tempchange_abs.sc*npp.sc +
                               tempchange_abs.sc*human_bowler.sc:REALM2 +
                               (tempchange_abs.sc|STUDY_ID/rarefyID), 
                           data = trends[i,], 
                           family = beta_family(link = 'logit'), 
                           dispformula = ~nspp.sc, 
                           control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modFullendo)
    saveRDS(modFullendo, file = 'temp/modFullendo.rds')
    print('saved modFullendo.rds')
    MATCHMOD <- TRUE
}

####################################
# final check that something ran
if(MATCHMOD == FALSE) stop("Model name did not match anything", call.=FALSE)

warnings()

print(paste('Ended', Sys.time(), sep =''))
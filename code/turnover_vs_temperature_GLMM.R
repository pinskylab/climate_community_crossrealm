#!/usr/bin/env Rscript

# Script to fit glmmTMB models
# Set up to be run on the command line for one model at a time
# Argument is model # to run (see below), e.g.
# nohup code/turnover_vs_temperature_GLMM.R 1 > logs/turnover_vs_temperature_GLMM1.Rout &

################
# Read command line arguments

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args)<1) stop("Have to specify a model to fit", call.=FALSE)
if (length(args)>1) stop("Have to specify only 1 model to fit", call.=FALSE)
if (length(args) == 1) fitmod = as.numeric(args[1])
if(is.na(fitmod)) stop("Have to specify model to fit with a number", call.=FALSE)



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

## Models to compare variance structures
fixed <- 'Jtu.sc ~ tempchange_abs.sc*REALM + tempchange_abs.sc*tempave_metab.sc + tempchange_abs.sc*duration.sc'
i <- trends[, complete.cases(Jtu.sc, tempchange_abs.sc, REALM, tempave_metab.sc, duration.sc)]

if(fitmod == 1){
    modRFgauss <- glmmTMB(formula(fixed), data = trends[i,], family = gaussian(), control = glmmTMBControl(profile=TRUE))
    summary(modRFgauss)
    saveRDS(modRFgauss, file = 'temp/modRFgauss.rds')
    print('saved modRFgauss.rds')
}

if(fitmod == 2){
    modRFbeta <- glmmTMB(formula(fixed), data = trends[i,], family = beta_family(link = 'logit'), control = glmmTMBControl(profile=TRUE)) # add beta errors
    summary(modRFbeta)
    saveRDS(modRFbeta, file = 'temp/modRFbeta.rds')
    print('saved modRFbeta.rds')
}
if(fitmod == 3){
    modRFrID <- glmmTMB(formula(paste(fixed, '+(1|rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random effects
    summary(modRFrID)
    saveRDS(modRFrID, file = 'temp/modRFrID.rds')
    print('saved modRDrID.rds')
}
if(fitmod == 4){
    modRFnestedRE <- glmmTMB(formula(paste(fixed, '+(1|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add nested random effects
    summary(modRFnestedRE)
    saveRDS(modRFnestedRE, file = 'temp/modRFnestedRE.rds')
    print('saved modRFnestedRE.rds')
}

if(fitmod == 5){
    modRFslopeRE <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   control = glmmTMBControl(profile=TRUE)) # add random slopes
    summary(modRFslopeRE)
    saveRDS(modRFslopeRE, file = 'temp/modRFslopeRE.rds')
    print('saved modRFslopeRE.rds')
}

# skipped 6!

if(fitmod == 7){
    modRFdisp <- glmmTMB(formula(paste(fixed, '+(tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID)')), data = trends[i,], family = beta_family(link = 'logit'), 
                   dispformula = ~nspp.sc, 
                   control = glmmTMBControl(profile=TRUE)) # add dispersion formula
    summary(modRFdisp)
    saveRDS(modRFdisp, file = 'temp/modRFdisp.rds')
    print('saved modRFdisp.rds')
}

## temperature-only models
if(fitmod == 8){
    i <- trends[, complete.cases(Jtu.sc, tempchange)]
    modonlyTchangeJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                             data = trends[i,], 
                             family = beta_family(link = 'logit'), 
                             dispformula = ~nspp.sc)
    summary(modonlyTchangeJtu)
    saveRDS(modonlyTchangeJtu, file = 'temp/modonlyTchangeJtu.rds')
    print('saved modonlyTchangeJtu.rds')
}

if(fitmod == 9){
    i <- trends[, complete.cases(Jbeta.sc, tempchange)]
    modonlyTchangeJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange) + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                               data = trends[i,], 
                               family = beta_family(link = 'logit'), 
                               dispformula = ~nspp.sc)
    summary(modonlyTchangeJbeta)
    saveRDS(modonlyTchangeJbeta, file = 'temp/modonlyTchangeJbeta.rds')
    print('saved modonlyTchangeJbeta.rds')
}

if(fitmod == 10){
    i <- trends[, complete.cases(Horn.sc, tempchange)]
    modonlyTchangeHorn <- glmmTMB(Horn.sc ~ abs(tempchange) + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                              data = trends[i,], 
                              family = beta_family(link = 'logit'), 
                              dispformula = ~nspp.sc)
    summary(modonlyTchangeHorn)
    saveRDS(modonlyTchangeHorn, file = 'temp/modonlyTchangeHorn.rds')
    print('saved modonlyTchangeHorn.rds')
}

## temperature-duration models
if(fitmod == 11){
    i <- trends[, complete.cases(Jtu.sc, REALM, tempchange)]
    modTDJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * duration + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                    data = trends[i,], 
                    family = beta_family(link = 'logit'), 
                    dispformula = ~nspp.sc)
    summary(mmodTDJtu)
    saveRDS(modTDJtu, file = 'temp/modTDJtu.rds')
    print('saved modTDJtu.rds')
}

if(fitmod == 12){
    i <- trends[, complete.cases(Jbeta.sc, REALM, tempchange)]
    modTDJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange)*duration + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                      data = trends[i,], 
                      family = beta_family(link = 'logit'), 
                      dispformula = ~nspp.sc)
    summary(modTDJbeta)
    saveRDS(modTDJbeta, file = 'temp/modTDJbeta.rds')
    print('saved modTDJbeta.rds')
}

if(fitmod == 13){
    i <- trends[, complete.cases(Horn.sc, REALM, tempchange)]
    modTDHorn <- glmmTMB(Horn.sc ~ abs(tempchange)*duration + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                     data = trends[i,], 
                     family = beta_family(link = 'logit'), 
                     dispformula = ~nspp.sc)
    summary(modTDHorn)
    saveRDS(modTDHorn, file = 'temp/modTDHorn.rds')
    print('saved modTDHorn.rds')
}

## temperature-duration models by REALM
if(fitmod == 14){
    i <- trends[, complete.cases(Jtu.sc, REALM, tempchange)]
    modTDrealmJtu <- glmmTMB(Jtu.sc ~ abs(tempchange) * duration + abs(tempchange) * REALM + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                         data = trends[i,], 
                         family = beta_family(link = 'logit'), 
                         dispformula = ~nspp.sc)
    summary(modTDrealmJtu)
    saveRDS(modTDrealmJtu, file = 'temp/modTDrealmJtu.rds')
    print('saved modTDrealmJtu.rds')
}

if(fitmod == 15){
    i <- trends[, complete.cases(Jbeta.sc, REALM, tempchange)]
    modTDrealmJbeta <- glmmTMB(Jbeta.sc ~ abs(tempchange)*duration + abs(tempchange) * REALM + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                           data = trends[i,], 
                           family = beta_family(link = 'logit'), 
                           dispformula = ~nspp.sc)
    summary(modTDrealmJbeta)
    saveRDS(modTDrealmJbeta, file = 'temp/modTDrealmJbeta.rds')
    print('saved modTDrealmJbeta.rds')
}

if(fitmod == 16){
    i <- trends[, complete.cases(Horn.sc, REALM, tempchange)]
    modTDrealmHorn <- glmmTMB(Horn.sc ~ abs(tempchange)*duration + abs(tempchange) * REALM + (tempchange_abs.sc|taxa_mod2/STUDY_ID/rarefyID),
                          data = trends[i,], 
                          family = beta_family(link = 'logit'), 
                          dispformula = ~nspp.sc)
    summary(modTDrealmHorn)
    saveRDS(modTDrealmHorn, file = 'temp/modTDrealmHorn.rds')
    print('saved modTDrealmHorn.rds')
}

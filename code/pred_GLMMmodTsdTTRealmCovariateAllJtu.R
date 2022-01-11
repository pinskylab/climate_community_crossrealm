## Make predictions data.frame from some of the covariate models
## From turnover_vs_temperature_GLMM.Rmd, but set up to run in the background
## so that we can also calculate SEs
#
# run as
# nohup Rscript --vanilla code/pred_GLMMmodTsdTTRealmCovariateAllJtu.R > logs/pred_GLMMmodTsdTTRealmCovariateAllJtu.Rout &
# or run by hand to, for example, start after the predictions have been made (line 93). would need to read the predictions in by hand.


print(paste('This is process #', Sys.getpid()))
print(Sys.time())

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(here) # for relative paths
library(lavaan) # for SEM models of slopes

scaleme <- function(x, nm){
    if(!(nm %in% scalingall[,var])) stop('nm not found in scalingall')
    if(scalingall[var==nm, log]){
        x.sc <- (log(x + scalingall[var==nm, plus]) - scalingall[var == nm, center]) / scalingall[var == nm, scale]  
    } else {
        x.sc <- (x  + scalingall[var==nm, plus] - scalingall[var == nm, center]) / scalingall[var == nm, scale]
    }
    return(x.sc)
}
unscaleme <- function(x.sc, nm){
    if(!(nm %in% scalingall[,var])) stop('nm not found in scalingall')
    if(scalingall[var==nm, log]){
        x <- exp(x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center]) - scalingall[var==nm, plus]
    } else {
        x <- x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center] - scalingall[var==nm, plus]
    }
    return(x)
}

# The scaling factors
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv'))

# The model
modTsdTTRealmmicroclimAllJtu <- readRDS('temp/modTsdTTRealmmicroclimAllJtu.rds') # has microclimates
modTsdTTRealmnppAllJtu <- readRDS('temp/modTsdTTRealmnppAllJtu.rds') # has npp
modTsdTTRealmseasAllJtu <- readRDS('temp/modTsdTTRealmseasAllJtu.rds') # has seasonality
modTsdTTRealmhumanAllJtu <- readRDS('temp/modTsdTTRealmhumanAllJtu.rds') # has human impact
modTsdTTRealmmassAllJtu <- readRDS('temp/modTsdTTRealmmassAllJtu.rds') # has mass
print('models loaded')

# set up prediction frame
newdat <- data.table(expand.grid(tempave_metab = c(5, 30), tempchange_abs = seq(0, 0.6, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial'), microclim.sc = c(-2, 2)))
newdat[, ':='(npp.sc = microclim.sc, seas.sc = microclim.sc, 
               human_bowler.sc = microclim.sc, mass.sc = microclim.sc)]
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave_metab.sc := scaleme(tempave_metab, 'tempave_metab.sc')]
newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]
newdat[, microclim := unscaleme(microclim.sc, 'microclim.sc')]
newdat[, npp := unscaleme(npp.sc, 'npp.sc')]
newdat[, seas := unscaleme(seas.sc, 'seas.sc')]
newdat[, human_bowler := unscaleme(human_bowler.sc, 'human_bowler.sc')]
newdat[, mass := unscaleme(mass.sc, 'mass.sc')]

# predict
preds.microclim <- predict(modTsdTTRealmmicroclimAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.microclim <- preds.microclim$fit
newdat$Jtu.sc.microclim.se <- preds.microclim$se.fit
print('finished microclim predictions')

preds.npp <- predict(modTsdTTRealmnppAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.npp <- preds.npp$fit
newdat$Jtu.sc.npp.se <- preds.npp$se.fit
print('finished npp predictions')

preds.seas <- predict(modTsdTTRealmseasAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.seas <- preds.seas$fit
newdat$Jtu.sc.seas.se <- preds.seas$se.fit
print('finished seas predictions')

preds.human <- predict(modTsdTTRealmhumanAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.human <- preds.human$fit
newdat$Jtu.sc.human.se <- preds.human$se.fit
print('finished human predictions')

preds.mass <- predict(modTsdTTRealmmassAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.mass <- preds.mass$fit
newdat$Jtu.sc.mass.se <- preds.mass$se.fit
print('finished mass predictions')

# write out
saveRDS(newdat, file = here('temp', 'preds_TsdTTRealmCovariate.rds'))


# calculate slopes and SE of the slope using latent variables (since predictions have SE)
mods.microclim <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.microclim\nJtu.sc.microclim ~~ ', mean(Jtu.sc.microclim.se), '*Jtu.sc.microclim'), data.frame(Jtu.sc.microclim, duration)))), 
                         by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)] # takes ~30 min
print('finished microclim SEM')
mods.npp <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.npp\nJtu.sc.npp ~~ ', mean(Jtu.sc.npp.se), '*Jtu.sc.npp'), data.frame(Jtu.sc.npp, duration)))), 
                         by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)]
print('finished npp SEM')
mods.seas <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.seas\nJtu.sc.seas ~~ ', mean(Jtu.sc.seas.se), '*Jtu.sc.seas'), data.frame(Jtu.sc.seas, duration)))), 
                   by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)]
print('finished seas SEM')
mods.human <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.human\nJtu.sc.human ~~ ', mean(Jtu.sc.human.se), '*Jtu.sc.human'), data.frame(Jtu.sc.human, duration)))), 
                   by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)]
print('finished human SEM')
mods.mass <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.mass\nJtu.sc.mass ~~ ', mean(Jtu.sc.mass.se), '*Jtu.sc.mass'), data.frame(Jtu.sc.mass, duration)))), 
                   by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)]
print('finished mass SEM')


# extract slopes and SEs
slopes.microclim <- mods.microclim[, parameterEstimates(mod[[1]]), 
                         by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                       .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_microclim = est, slope_microclim.se = se)]
slopes.npp <- mods.npp[, parameterEstimates(mod[[1]]), 
                         by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                       .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_npp = est, slope_npp.se = se)]
slopes.seas <- mods.seas[, parameterEstimates(mod[[1]]), 
                       by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                     .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_seas = est, slope_seas.se = se)]
slopes.human <- mods.human[, parameterEstimates(mod[[1]]), 
                       by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                     .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_human = est, slope_human.se = se)]
slopes.mass <- mods.mass[, parameterEstimates(mod[[1]]), 
                       by = .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                     .(tempave_metab, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_mass = est, slope_mass.se = se)]

slopes2 <- merge(slopes.microclim, slopes.npp)
slopes2 <- merge(slopes2, slopes.seas)
slopes2 <- merge(slopes2, slopes.human)
slopes2 <- merge(slopes2, slopes.mass)


# save slopes
saveRDS(slopes2, file = here('temp', 'slopes_interactions2.rds'))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

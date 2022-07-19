## Make predictions data.frame from some of the covariate models with environmental temperature (rawT)
## From turnover_vs_temperature_GLMM.Rmd, but set up to run in the background
## so that we can also calculate SEs
#
# run as
# nohup Rscript --vanilla code/pred_GLMMmodrawTsdTTRealmCovariateAllJtu.R > logs/pred_GLMMmodrawTsdTTRealmCovariateAllJtu.Rout &
# or run by hand to, for example, start after the predictions have been made (line 93). would need to read the predictions in by hand.
# set to run on Annotate


print(paste('This is process #', Sys.getpid()))
print(Sys.time())

# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

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
modrawTsdTTRealmmicroclimAllJtu <- readRDS('temp/modrawTsdTTRealmmicroclimAllJtu.rds') # has microclimates
modrawTsdTTRealmnppAllJtu <- readRDS('temp/modrawTsdTTRealmnppAllJtu.rds') # has npp
modrawTsdTTRealmseasAllJtu <- readRDS('temp/modrawTsdTTRealmseasAllJtu.rds') # has seasonality
modrawTsdTTRealmhumanAllJtu <- readRDS('temp/modrawTsdTTRealmhumanAllJtu.rds') # has human impact
print('models loaded')

# set up prediction frame
newdat <- data.table(expand.grid(tempave = c(0, 8, 10, 13, 30), tempchange_abs = seq(0, 0.6, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial'), microclim.sc = c(-2, 2)))
newdat[, ':='(npp.sc = microclim.sc, seas.sc = microclim.sc, 
               human_bowler.sc = microclim.sc, mass.sc = microclim.sc)]
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]
newdat[, microclim := unscaleme(microclim.sc, 'microclim.sc')]
newdat[, npp := unscaleme(npp.sc, 'npp.sc')]
newdat[, seas := unscaleme(seas.sc, 'seas.sc')]
newdat[, human_bowler := unscaleme(human_bowler.sc, 'human_bowler.sc')]
newdat[, mass := unscaleme(mass.sc, 'mass.sc')]

# predict
preds.microclim <- predict(modrawTsdTTRealmmicroclimAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.microclim <- preds.microclim$fit
newdat$Jtu.sc.microclim.se <- preds.microclim$se.fit
print('finished microclim predictions')

preds.npp <- predict(modrawTsdTTRealmnppAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.npp <- preds.npp$fit
newdat$Jtu.sc.npp.se <- preds.npp$se.fit
print('finished npp predictions')

preds.seas <- predict(modrawTsdTTRealmseasAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.seas <- preds.seas$fit
newdat$Jtu.sc.seas.se <- preds.seas$se.fit
print('finished seas predictions')

preds.human <- predict(modrawTsdTTRealmhumanAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.human <- preds.human$fit
newdat$Jtu.sc.human.se <- preds.human$se.fit
print('finished human predictions')

# write out
saveRDS(newdat, file = here('temp', 'preds_rawTsdTTRealmCovariate.rds'))


# calculate slopes and SE of the slope using latent variables (since predictions have SE)
mods.microclim <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.microclim\nJtu.sc.microclim ~~ ', mean(Jtu.sc.microclim.se), '*Jtu.sc.microclim'), data.frame(Jtu.sc.microclim, duration)))), 
                         by = .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)] # takes ~30 min
print('finished microclim SEM')
mods.npp <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.npp\nJtu.sc.npp ~~ ', mean(Jtu.sc.npp.se), '*Jtu.sc.npp'), data.frame(Jtu.sc.npp, duration)))), 
                         by = .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)]
print('finished npp SEM')
mods.seas <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.seas\nJtu.sc.seas ~~ ', mean(Jtu.sc.seas.se), '*Jtu.sc.seas'), data.frame(Jtu.sc.seas, duration)))), 
                   by = .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)]
print('finished seas SEM')
mods.human <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.human\nJtu.sc.human ~~ ', mean(Jtu.sc.human.se), '*Jtu.sc.human'), data.frame(Jtu.sc.human, duration)))), 
                   by = .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)]
print('finished human SEM')


# extract slopes and SEs
slopes.microclim <- mods.microclim[, parameterEstimates(mod[[1]]), 
                         by = .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                       .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_microclim = est, slope_microclim.se = se)]
slopes.npp <- mods.npp[, parameterEstimates(mod[[1]]), 
                         by = .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                       .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_npp = est, slope_npp.se = se)]
slopes.seas <- mods.seas[, parameterEstimates(mod[[1]]), 
                       by = .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                     .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_seas = est, slope_seas.se = se)]
slopes.human <- mods.human[, parameterEstimates(mod[[1]]), 
                       by = .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                     .(tempave, tempchange_abs, microclim, npp, seas, human_bowler, mass, REALM, slope_human = est, slope_human.se = se)]

slopes2 <- merge(slopes.microclim, slopes.npp)
slopes2 <- merge(slopes2, slopes.seas)
slopes2 <- merge(slopes2, slopes.human)


# save slopes
saveRDS(slopes2, file = here('temp', 'slopes_rawinteractions2.rds'))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

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

# load functions and data ----------------

# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(here) # for relative paths
library(lavaan) # for SEM models of slopes
source(here('code', 'util.R'))



# The scaling factors
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv'))


# The model
modrawTsdTTRealmtsignmicroclimAllJtu <- readRDS('temp/modrawTsdTTRealmtsignmicroclimAllJtu.rds') # has microclimates
modrawTsdTTRealmtsignhumanAllJtu <- readRDS('temp/modrawTsdTTRealmtsignhumanAllJtu.rds') # has human impact
print('models loaded')


# Make predictions -------------------------
# set up prediction frame
newdat <- data.table(expand.grid(tempave = c(0, 8, 10, 13, 30), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial'), microclim.sc = seq(-2, 2, length.out=10)))
newdat[, ':='(human_bowler.sc = microclim.sc)]
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
newdat[, tempchange_abs := abs(tempchange)]
newdat[, tsign := signneg11(tempchange)]
newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]
newdat[, microclim := unscaleme(microclim.sc, 'microclim.sc')]
newdat[, human_bowler := unscaleme(human_bowler.sc, 'human_bowler.sc')]

# constrain human covariate to >=0 and <=10, since that is the range over which it is defined
if(newdat[, min(human_bowler)<0]){
    newdat[abs(human_bowler.sc - -2)<0.001, human_bowler := 0]
    newdat[, human_bowler.sc := scaleme(human_bowler, 'human_bowler.sc')]
}
if(newdat[, max(human_bowler)>10]){
    newdat[abs(human_bowler.sc - 2)<0.01, human_bowler := 10]
    newdat[, human_bowler.sc := scaleme(human_bowler, 'human_bowler.sc')]
}

# predict
preds.microclim <- predict(modrawTsdTTRealmtsignmicroclimAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.microclim <- preds.microclim$fit
newdat$Jtu.sc.microclim.se <- preds.microclim$se.fit
print('finished microclim predictions')

preds.human <- predict(modrawTsdTTRealmtsignhumanAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.human <- preds.human$fit
newdat$Jtu.sc.human.se <- preds.human$se.fit
print('finished human predictions')

# write out
saveRDS(newdat, file = here('temp', 'preds_rawTsdTTRealmtsignCovariate.rds'))
print(paste('Wrote preds_rawTsdTTRealmtsignCovariate.rds:', Sys.time()))


# Calculate slopes -----------------------------
# calculate slopes and SE of the slope using latent variables (since predictions have SE)
mods.microclim <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.microclim\nJtu.sc.microclim ~~ ', mean(Jtu.sc.microclim.se), '*Jtu.sc.microclim'), data.frame(Jtu.sc.microclim, duration)))), 
                         by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)] # takes ~30 min
print('finished microclim SEM')
mods.human <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.human\nJtu.sc.human ~~ ', mean(Jtu.sc.human.se), '*Jtu.sc.human'), data.frame(Jtu.sc.human, duration)))), 
                     by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)]
print('finished human SEM')


# extract slopes and SEs
slopes.microclim <- mods.microclim[, parameterEstimates(mod[[1]]), 
                                   by = .(tempave, tempchange, tempchange_abs, tsign, 
                                          microclim, human_bowler, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                                            .(tempave, tempchange, tempchange_abs, 
                                                                                              microclim, human_bowler, REALM, 
                                                                                              slope_microclim = est, slope_microclim.se = se)]
slopes.human <- mods.human[, parameterEstimates(mod[[1]]), 
                           by = .(tempave, tempchange, tempchange_abs, tsign, 
                                  microclim, human_bowler, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                                    .(tempave, tempchange, tempchange_abs, tsign, 
                                                                                      microclim, human_bowler, REALM, 
                                                                                      slope_human = est, slope_human.se = se)]

slopes2 <- merge(slopes.microclim, slopes.human)


# save slopes
saveRDS(slopes2, file = here('temp', 'slopes_rawTsdTTRealmtsignCovariate.rds'))
print(paste('Wrote slopes_rawTsdTTRealmtsignCovariate.rds:', Sys.time()))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

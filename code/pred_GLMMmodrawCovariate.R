## Make predictions data.frame from the environmental covariate models
## From turnover_vs_temperature_GLMM.Rmd, but set up to run in the background
## so that we can also calculate SEs
#
# run as
# nohup Rscript --vanilla code/pred_GLMMmodrawTsdTTRealmCovariateAllJtu.R > logs/pred_GLMMmodrawTsdTTRealmCovariateAllJtu.Rout &
# or run by hand to, for example, start after the predictions have been made (line 93). would need to read the predictions in by hand.
# set to run on Annotate


print(paste('This is process #', Sys.getpid()))
print(Sys.time())

### set arguments -----------------
n = 1000 # number of resamples to do for each timeseries


### load functions and data ----------------
# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(here) # for relative paths
source(here('code', 'util.R'))

# function for slopes and SEs from resampling
# n: number of resamples, other columns are the x, y, and se of y variables
# colnames refer to the output
slopesamp <- function(n, duration, Jtu.sc, Jtu.sc.se, colnames = c('slope', 'slope.se')){
    if(length(duration) != length(Jtu.sc)) stop('duration and Jtu.sc are not the same length')
    if(length(duration) != length(Jtu.sc.se)) stop('duration and Jtu.sc.se are not the same length')
    if(length(Jtu.sc) != length(Jtu.sc.se)) stop('Jtu.sc and Jtu.sc.se are not the same length')
    
    samp <- rep(NA, n) # will hold the slopes of the sampled data
    for(j in 1:n){
        y <- rnorm(length(Jtu.sc), mean = Jtu.sc, sd = Jtu.sc.se) # one sample
        samp[j] <- coef(lm(y ~ duration))[2] # fit line, get slope
    }
    out <- c(coef(lm(Jtu.sc ~ duration))[2], sd(samp))
    names(out) <- colnames
    return(as.list(out)) # coercing to list will allow the data.table aggregate used later to create 2 columns
}


# The scaling factors
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv'))

# The models
modmicroclim <- readRDS('temp/modrawTsdTTRealmtsignmicroclimInitAllJtu.rds') # has microclimates
modhuman <- readRDS('temp/modrawTsdTTRealmtsignhumanInitAllJtu.rds') # has human impact
print('models loaded')


### Make predictions -------------------------
# set up prediction frame
newdat <- data.table(expand.grid(tempave = c(0, 8, 10, 13, 30), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial'), microclim.sc = seq(-2, 2, length.out=10)))
newdat[, ':='(human_bowler.sc = microclim.sc)]
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat$Jtu.init <- 0.5
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
preds.microclim <- predict(modmicroclim, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.microclim <- preds.microclim$fit
newdat$Jtu.sc.microclim.se <- preds.microclim$se.fit
print('finished microclim predictions')

preds.human <- predict(modhuman, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.human <- preds.human$fit
newdat$Jtu.sc.human.se <- preds.human$se.fit
print('finished human predictions')

# write out
saveRDS(newdat, file = here('temp', 'preds_rawTsdTTRealmtsignCovariateInit.rds'))
print(paste('Wrote preds_rawTsdTTRealmtsignCovariateInit.rds:', Sys.time()))

# if reading in (eg, if running by hand)
# newdat <- readRDS(here('temp', 'preds_rawTsdTTRealmtsignCovariateInit.rds'))



### Calculate turnover rates (change in turnover per year) -----------------------------
slopes.microclim <- newdat[, slopesamp(n, duration, Jtu.sc.microclim, Jtu.sc.microclim.se, colnames = c('slope_microclim', 'slope_microclim.se')), 
                                  by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)]
slopes.human <- newdat[, slopesamp(n, duration, Jtu.sc.human, Jtu.sc.human.se, colnames = c('slope_human', 'slope_human.se')), 
                              by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)]
slopes2 <- merge(slopes.microclim, slopes.human)

# save turnover rates
saveRDS(slopes2, file = here('temp', 'slopes_rawTsdTTRealmtsignCovariateInit.rds'))
print(paste('Wrote slopes_rawTsdTTRealmtsignCovariateInit.rds:', Sys.time()))



### Calculate sensitivity of turnover rates to covariates ----------------------
# slopes2 <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsignCovariateInit.rds')) # to read in manually
slopes2 <- slopes2[tempave == 13 & tempchange > 0, ] # pick tempave and warming

sensitivity.microclim <- slopes2[, slopesamp(n, tempchange, slope_microclim, slope_microclim.se, colnames = c('sensitivity_microclim', 'sensitivity_microclim.se')), 
                                 by = .(microclim, human_bowler, REALM)]
sensitivity.human <- slopes2[, slopesamp(n, tempchange, slope_human, slope_human.se, colnames = c('sensitivity_human', 'sensitivity_human.se')), 
                                 by = .(microclim, human_bowler, REALM)]


sensitivity2 <- merge(sensitivity.microclim, sensitivity.human)

# save sensitivities
saveRDS(sensitivity2, file = here('temp', 'sensitivity_rawTsdTTRealmtsignCovariateInit.rds'))
# sensitivity2 <- readRDS(here('temp', 'sensitivity_rawTsdTTRealmtsignCovariateInit.rds')) # to read in manually
print(paste('Wrote sensitivity_rawTsdTTRealmtsignCovariateInit.rds:', Sys.time()))



print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

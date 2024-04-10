## Make predictions data.frame from the environmental covariate models
## Set up to run in the background so that we can also calculate SEs
#
# run as
# nohup Rscript --vanilla code/pred_GLMMcov.R X > logs/pred_GLMMcov_X.Rout &
# where X is the stem of the model name (see below for modname options)
# set to run on Annotate or Annotate2


print(paste('This is process #', Sys.getpid()))
print(Sys.time())

### set arguments -----------------
n = 1000 # number of resamples to do for each timeseries

### read arguments ----------------
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 1)
    stop("Have to specify a model type to predict", call. = FALSE)
if (length(args) > 1)
    stop("Have to specify only 1 model type to predict", call. = FALSE)
modname <- args[1]

MATCHMOD <- FALSE # indicator of the argument matches an accepted model name
if(modname == 'modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu'){
    micromodel <- 'modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu'
    humanmodel <- 'modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu'
    MATCHMOD <- TRUE
}
if(modname == 'modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu_marterr'){
    micromodel <- 'modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu_marterr'
    humanmodel <- 'modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu_marterr'
    MATCHMOD <- TRUE
}
if(!MATCHMOD){
    stop('Need to specify a recognized covariate model name')
}

### load functions and data ----------------
# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

library(data.table) # for handling large datasets
if(Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){ # location of the library on Annotate
    library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
} else {
    library(glmmTMB)
}
library(here) # for relative paths
source(here('code', 'util.R'))



# The scaling factors
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv'))

# The models
microinfile <- paste0('temp/', micromodel, '.rds')
humaninfile <- paste0('temp/', humanmodel, '.rds')

modmicroclim <- readRDS(microinfile) # has microclimates
modhuman <- readRDS(humaninfile) # has human impact
print('models loaded')
print(microinfile)
print(humaninfile)


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
predsfile <- paste0('preds_', modname, '.rds')
saveRDS(newdat, file = here('temp', predsfile))
print(paste('Wrote', predsfile, ':', Sys.time()))

# if reading in (eg, if running by hand)
# newdat <- readRDS(here('temp', predsfile))



### Calculate turnover rates (change in turnover per year) -----------------------------
slopes.microclim <- newdat[, slopesamp(n, duration, Jtu.sc.microclim, Jtu.sc.microclim.se, colnames = c('slope_microclim', 'slope_microclim.se')), 
                                  by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)]
slopes.human <- newdat[, slopesamp(n, duration, Jtu.sc.human, Jtu.sc.human.se, colnames = c('slope_human', 'slope_human.se')), 
                              by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)]
slopes2 <- merge(slopes.microclim, slopes.human)

# save turnover rates
slopesfile <- paste0('slopes_', modname, '.rds')
saveRDS(slopes2, file = here('temp', slopesfile))
print(paste('Wrote', slopesfile, ':', Sys.time()))


### Calculate sensitivity of turnover rates to covariates ----------------------
# slopes2 <- readRDS(here('temp', slopesfile)) # to read in manually
slopes2 <- slopes2[tempave == 13 & tempchange > 0, ] # pick tempave and warming

sensitivity.microclim <- slopes2[, slopesamp(n, tempchange, slope_microclim, slope_microclim.se, colnames = c('sensitivity_microclim', 'sensitivity_microclim.se')), 
                                 by = .(microclim, human_bowler, REALM)]
sensitivity.human <- slopes2[, slopesamp(n, tempchange, slope_human, slope_human.se, colnames = c('sensitivity_human', 'sensitivity_human.se')), 
                                 by = .(microclim, human_bowler, REALM)]

sensitivity2 <- merge(sensitivity.microclim, sensitivity.human)

# save sensitivities
sensfile <- paste0('sensitivity_', modname, '.rds')
saveRDS(sensitivity2, file = here('temp', sensfile))
# sensitivity2 <- readRDS(here('temp', sensfile)) # to read in manually
print(paste('Wrote', sensfile, ':', Sys.time()))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

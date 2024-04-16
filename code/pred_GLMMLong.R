#!/usr/bin/Rscript --vanilla


## Make predictions data.frame from Tchange and Tchange x Tave models fit to longer time series (see turnover_GLMMlong_fit.R)
#
# run as
# nohup code/pred_GLMMLong.R X > logs/pred_GLMMLong_X.Rout &
# where X is an argument (see below)
# or run by hand to, for example, start after the predictions have been made (line 93). would need to read the predictions in by hand.
# set to run on Annotate


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
predmod <- args[1]



### Functions and data ----------------------
# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(here) # for relative paths
source(here('code', 'util.R'))



# The scaling factors
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv')) # From assemble_turnover_covariates.Rmd



### Choose a model ---------------------------------
if(predmod == 'modsdTRealmtsigninitLongJtu'){ # Tchange model
    mod <- readRDS(here('temp', 'modsdTRealmtsigninitLongJtu.rds')) # From turnover_GLMMlong_fit.R
    out_preds <- 'preds_modsdTRealmtsigninitLongJtu.rds'
    out_slopes <- 'slopes_modsdTRealmtsigninitLongJtu.rds'
    doSensitivity <- FALSE # calculate sensitivity to Tave?
    print('model loaded')
} 

if(predmod == 'modrawTsdTTRealmtsigninitLongJtu'){ # Tchange x Tave model
    mod <- readRDS(here('temp', 'modrawTsdTTRealmtsigninitLongJtu.rds')) # From turnover_GLMMlong_fit.R
    out_preds <- 'preds_rawTsdTTRealmtsigninitLong.rds' # note: name is not quite of same form as previous model.
    out_slopes <- 'slopes_rawTsdTTRealmtsigninitLong.rds' # note: name is not quite of same form as previous model.
    doSensitivity <- TRUE # calculate sensitivity to Tave?
    out_sensitivity <- 'sensitivity_rawTsdTTRealmtsigninitLong.rds'
    print('model loaded')
} 

if(!exists('mod')) stop('No model loaded. Make sure argument matches one of the model names')




### Make predictions -------------------
# set up prediction frame
newdat <- data.table(expand.grid(tempave = c(0, 25), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial', 'Freshwater')))
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat$Jtu.init <- 0.5
newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
newdat[, tempchange_abs := abs(tempchange)]
newdat[, tsign := signneg11(tempchange)]
newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]

# predict with SEs
preds <- predict(mod, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc <- preds$fit
newdat$Jtu.sc.se <- preds$se.fit
print('finished  predictions')

# Write dissimilarities
saveRDS(newdat, file = here('temp', out_preds))
print(paste('Wrote', out_preds, ':', Sys.time()))


### Slope calculations -------------------
slopes <- newdat[, slopesamp(n, duration, Jtu.sc, Jtu.sc.se), 
                 by = .(tempave, tempchange, REALM)]

# Write slopes
saveRDS(slopes, file = here('temp', out_slopes))
print(paste('Wrote', out_slopes, ':', Sys.time()))
#slopes <- readRDS(file = here('temp', out_slopes)) # to read back in



### Calculate sensitivity of turnover rates to Tave ----------------------
# if we have a model with Tave
if(doSensitivity){
    sensitivity <- slopes[, slopesamp(n, abs(tempchange), slope, slope.se, colnames = c('sensitivity', 'sensitivity.se')), 
                                     by = .(tempave, tsign = sign(tempchange), REALM)]
    
    # write sensitivities
    saveRDS(sensitivity, file = here('temp', out_sensitivity))
    print(paste('Wrote', out_sensitivity, ':', Sys.time()))
}


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

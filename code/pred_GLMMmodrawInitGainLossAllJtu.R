#!/usr/bin/Rscript --vanilla


## Make predictions data.frame from sdT and TsdT models with environmental temperature (rawT), Jtu.init, and GainLoss
#
# run as
# nohup code/pred_GLMMmodrawInitGainLossAllJtu.R X > logs/pred_GLMMmodrawInitGainLossAllJtu_X.Rout &
# where X is an argument (see below)
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
library(lavaan) # for SEM models of slopes
source(here('code', 'util.R'))


# slopes and SEs from resampling
# n: number of resamples, other columns are the x, y, and se of y variables
slopesamp <- function(n, duration, Jtu.sc, Jtu.sc.se){
    if(length(duration) != length(Jtu.sc)) stop('duration and Jtu.sc are not the same length')
    if(length(duration) != length(Jtu.sc.se)) stop('duration and Jtu.sc.se are not the same length')
    if(length(Jtu.sc) != length(Jtu.sc.se)) stop('Jtu.sc and Jtu.sc.se are not the same length')
    
    samp <- rep(NA, n) # will hold the slopes of the sampled data
    for(j in 1:n){
        y <- rnorm(length(Jtu.sc), mean = Jtu.sc, sd = Jtu.sc.se) # one sample
        samp[j] <- coef(lm(y ~ duration))[2] # fit line, get slope
    }
    out <- c(coef(lm(Jtu.sc ~ duration))[2], sd(samp))
    names(out) <- c('slope', 'slope.se')
    return(as.list(out)) # coercing to list will allow the data.table aggregate used later to create 2 columns
}


# The scaling factors
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv')) # From assemble_turnover_covariates.Rmd



### Choose a model ---------------------------------
if(predmod == 'modInitGainLossAllJtu'){
    mod <- readRDS(here('temp', 'modInitGainLossAllJtu.rds')) # From turnover_vs_temperature_GLMM_fit_Jtu.init.gainloss.R
    out_preds <- 'preds_modInitGainLossAllJtu.rds'
    out_slopes <- 'slopes_modInitGainLossAllJtu.rds'
    print('model loaded')
} 

if(predmod == 'modsdTRealmtsignInitGainLossAllJtu'){
    mod <- readRDS(here('temp', 'modsdTRealmtsignInitGainLossAllJtu.rds')) # From turnover_vs_temperature_GLMM_fit_Jtu.init.gainloss.R
    out_preds <- 'preds_modsdTRealmtsignInitGainLossAllJtu.rds'
    out_slopes <- 'slopes_modsdTRealmtsignInitGainLossAllJtu.rds'
    print('model loaded')
} 

if(predmod == 'modrawTsdTTRealmtsignInitGainLossAllJtu'){
    mod <- readRDS(here('temp', 'modrawTsdTTRealmtsignInitGainLossAllJtu.rds')) # From turnover_vs_temperature_GLMM_fit_Jtu.init.gainloss.R
    out_preds <- 'preds_modrawTsdTTRealmtsignInitGainLossAllJtu.rds' 
    out_slopes <- 'slopes_modrawTsdTTRealmtsignInitGainLossAllJtu.rds' 
    print('model loaded')
} 


if(!exists('mod')) stop('No model loaded. Make sure argument matches one of the model names')




### Make predictions -------------------
# set up prediction frame
newdat <- data.table(expand.grid(tempave = seq(-20, 30, length.out = 10), 
                                 tempchange = seq(-1.5, 2, length.out=10), 
                                 duration = 1:10, 
                                 REALM = c('Marine', 'Terrestrial', 'Freshwater'),
                                 Jtu.init = seq(0, 1, length.out = 5),
                                 gainlossprop = seq(-2.5, 2.5, length.out = 5)))
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
newdat[, tempchange_abs := abs(tempchange)]
newdat[, tsign := signneg11(tempchange)]
newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]

# predict with SEs
preds <- predict(mod, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc <- preds$fit
newdat$Jtu.sc.se <- preds$se.fit
print('finished  predictions')

### Write dissimilarities ------------------
saveRDS(newdat, file = here('temp', out_preds))


### Slope calculations -------------------
slopes <- newdat[, slopesamp(n, duration, Jtu.sc, Jtu.sc.se), 
                 by = .(tempave, tempchange, Jtu.init, gainlossprop, REALM)]

### Write slopes -------------
saveRDS(slopes, file = here('temp', out_slopes))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

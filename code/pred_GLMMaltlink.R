#!/usr/bin/Rscript --vanilla


## Make predictions data.frame from Tchange and Tchange x Tave models fit with alternative link functions (see turnover_GLMMaltlink_fit.R)
#
# run as
# nohup Rscript --vanilla code/pred_GLMMaltlink.R X > logs/pred_GLMMaltlink_X.Rout &
# where X is an argument (see below)
# or as:
# pred_GLMMaltlink.sh X1 X2
# which spawns a job for each argument


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
if(Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){ # location of the library on Annotate
    library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
} else {
    library(glmmTMB)
}
library(here) # for relative paths
source(here('code', 'util.R'))


# slopes and SEs from resampling
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
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv')) # From assemble_turnover_covariates.Rmd



### Choose a model ---------------------------------
inmod <- paste0(predmod, '.rds')
mod <- readRDS(here('temp', inmod)) # From turnover_GLMMaltlink_fit.R
out_preds <- paste0('preds_', predmod, '.rds')
out_slopes <- paste0('slopes_', predmod, '.rds')
out_sensitivity <- paste0('sensitivity_', predmod, '.rds')
if(grepl('raw', predmod)){
    doSensitivity <- TRUE # calculate sensitivity to Tave if this is a Tave model
} else {
    doSensitivity <- FALSE
}

if(!exists('mod')) stop('No model loaded. Make sure argument matches one of the model names')
print('model loaded')




### Make predictions -------------------
# set up prediction frame
newdat <- data.table(expand.grid(tempave = seq(-10, 36, by = 1), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial', 'Freshwater')))
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

#!/usr/bin/Rscript --vanilla

## Resample from the uncertainty in the predictions from the sdT and TsdT models with environmental temperature (rawT)
## and calculate the SE of the slopes of the resulting timeseries
## Requires the preds_X.rds file produced by pred_GLMMmodrawXAllJtu.R script
#
# run as
# nohup code/pred_GLMMmodrawXAllJtu_resampleonly.R X > logs/pred_GLMMmodrawXAllJtu_X.Rout &
# where X is an argument (see below)
# or as:
# pred_modrawXAllJtu_resamponly.sh X1 X2
# which spawns a job for each argument
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


# Choose a model
if(predmod == 'modsdTRealmtsignAllJtu'){
    out_preds <- 'preds_modsdTRealmtsignAllJtu.rds'
    out_slopes <- 'slopes_modsdTRealmtsignAllJtu.resamp.rds'
} 

if(predmod == 'modrawTsdTTRealmtsignAllJtu'){
    out_preds <- 'preds_rawTsdTTRealmtsign.rds' # note: name is not quite of same form as previous model.
    out_slopes <- 'slopes_rawTsdTTRealmtsign.resamp.rds' # note: name is not quite of same form as previous model.
} 

if(predmod == 'modsdTRealmtsigninitAllJtu'){
    out_preds <- 'preds_modsdTRealmtsigninitAllJtu.rds'
    out_slopes <- 'slopes_modsdTRealmtsigninitAllJtu.resamp.rds'
} 

if(predmod == 'modrawTsdTTRealmtsigninitAllJtu'){
    out_preds <- 'preds_rawTsdTTRealmtsigninit.rds' # note: name is not quite of same form as previous model.
    out_slopes <- 'slopes_rawTsdTTRealmtsigninit.resamp.rds' # note: name is not quite of same form as previous model.
} 

if(!exists('out_preds')) stop('Make sure argument matches one of the model names')



### Functions and data ----------------------
# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

library(data.table) # for handling large datasets
library(here) # for relative paths
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



### Load dissimilarity predictions -------------------
newdat <- readRDS(here('temp', out_preds))


### Slope calculations -------------------
slopes <- newdat[, slopesamp(n, duration, Jtu.sc, Jtu.sc.se), 
                 by = .(tempave, tempchange, REALM)]


### Write slopes -------------
saveRDS(slopes, file = here('temp', out_slopes))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

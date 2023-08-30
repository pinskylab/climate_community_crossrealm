#!/usr/bin/Rscript --vanilla

## Resample from the uncertainty in the predictions from the covariate models with environmental temperature (rawT)
## and calculate the SE of the slopes of the resulting timeseries
## Requires the preds_rawTsdTTRealmtsignCovariate.rds file produced by pred_GLMMmodrawTsdTTRealmCovariateAllJtu.R script

print(paste('This is process #', Sys.getpid()))
print(Sys.time())

# load functions and data ----------------

# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

library(data.table) # for handling large datasets
library(here) # for relative paths
source(here('code', 'util.R'))

n = 1000 # number of resamples to do for each timeseries
newdat <- readRDS(here('temp', 'preds_rawTsdTTRealmtsignCovariate.rds'))

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

slopes.microclim.resamp <- newdat[, slopesamp(n, duration, Jtu.sc.microclim, Jtu.sc.microclim.se, colnames = c('slope_microclim', 'slope_microclim.se')), 
                              by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)]
slopes.human.resamp <- newdat[, slopesamp(n, duration, Jtu.sc.human, Jtu.sc.human.se, colnames = c('slope_human', 'slope_human.se')), 
                          by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)]
slopes2.resamp <- merge(slopes.microclim.resamp, slopes.human.resamp)

# save resamp slopes
saveRDS(slopes2.resamp, file = here('temp', 'slopes_rawTsdTTRealmtsignCovariate.resamp.rds'))
# slopes2.resamp <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsignCovariate.resamp.rds')) # to read in manually
print(paste('Wrote slopes_rawTsdTTRealmtsignCovariate.resamp.rds:', Sys.time()))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

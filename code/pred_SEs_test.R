## Sample from the predictions SEs from some of the covariate models with environmental temperature (rawT)
# fit lm
# repeat
# done to understand what the SEs on the linear predictions are
# assumes that pred_GLMMmodrawTsdTTRealmCovariateAllJtu.R has been run to produce preds_rawTsdTTRealmtsignCovariate.rds


# load functions and data ----------------

# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(here) # for relative paths
source(here('code', 'util.R'))

newdat <- readRDS(here('temp', 'preds_rawTsdTTRealmtsignCovariate.rds'))

# select one example of a time-series
i <- newdat[, tempave == 0 & tempchange == -1.5 & REALM == 'Marine' & human_bowler==0]

# sample from the SEs for microclim preds
n <- 1000
slopesamp <- rep(NA, n) # will hold the slopes of the sampled data
for(j in 1:n){
    thisdat <- newdat[i, .(duration, Jtu.sc = rnorm(10, mean = Jtu.sc.human, sd = Jtu.sc.human.se))] # one sample
    slopesamp[j] <- coef(thisdat[, lm(Jtu.sc ~ duration)])[2] # fit line, get slope
}
hist(slopesamp)

# compare to latent variable and lm slopes
slopes2 <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsignCovariate.rds')) # latent variable slopes and se
slopes2.lm <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsignCovariate.lm.rds')) # lm slopes and ses

cbind(mean(slopesamp), sd(slopesamp))
slopes2[tempave == 0 & tempchange == -1.5 & REALM == 'Marine' & human_bowler==0, .(slope_human, slope_human.se)]
slopes2.lm[tempave == 0 & tempchange == -1.5 & REALM == 'Marine' & human_bowler==0, .(slope_human, slope_human.se)]

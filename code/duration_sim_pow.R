#!/usr/bin/env Rscript

# Code to run the true-positive (power) simulations and analysis for duration_sim.Rmd
# Set up as a script to run in the background on Annotate, e.g.,
# nohup code/duration_sim_pow.R > logs/duration_sim_pow.Rout &
# (this works if code is executable, e.g., chmod u+x code/duration_sim_pow.R)

# print basic info about the job ############################
print(paste('This is process #', Sys.getpid()))
print(Sys.time())


# load libraries ############################
require(data.table)
library(here)
library(glmmTMB, lib.loc = '/usr/lib64/R/library/')


# Functions ################
# transformation for 2 categories. Eq. 1 in Douma & Weedon 2019 MEE
transform01 <- function(x) (x * (length(x) - 1) + 0.5) / (length(x))

slope <- function(x){
    mod <- lm(x ~ I(1:length(x)))
    return(as.numeric(coef(mod)[2]))
}
slopese <- function(x){
    mod <- lm(x ~ I(1:length(x)))
    return(as.numeric(sqrt(diag(vcov(mod)))[2]))
}



# Simulations #############
set.seed(5)
nreps <- 100 # number of datasets

cors <- data.table(n = 1000, minduration = rep(c(10,3), nreps*3), maxduration = rep(c(10, 102), nreps*3), effect = rep(c(0.1, 1, 10), c(nreps*2, nreps*2, nreps*2)), 
                   lm.p = NA_real_, lm.m = NA_real_,
                   glmmwgt.p = NA_real_, glmmwgt.beta = NA_real_,
                   glmmonebeta.p = NA_real_, glmmonebeta.beta = NA_real_) # holds the parameters for and summaries from each dataset. use n = 1000 timeseries per dataset.
#cors2 <- copy(cors)
#cors2[, n:=10000] # add another run with 10k timeseries per dataset
#cors <- rbind(cors, cors2) 
#rm(cors2)

cors[, name := paste0(minduration, '-', maxduration, '-', effect)]
cors[, range := maxduration - minduration]


print(paste0(nrow(cors), ' simulations to make'))
for(i in 1:nrow(cors)){
    cat(i)
    
    # make datasets with timeseries of all the same length
    if(cors[i, minduration == maxduration]){
        len <- cors[i, minduration]
        dat <- data.table(tsid = rep(1:cors[i,n], rep(len, cors[i,n])))
    }
    
    # variable length timeseries
    if(cors[i, minduration != maxduration]){
        mind <- cors[i, minduration]
        maxd <- cors[i, maxduration]
        ndur <- maxd - mind + 1
        dat <- data.table(tsid = rep(1:cors[i,n], rep(mind:maxd, cors[i,n]/ndur)))
    }
    
    # add response and explanatory variables
    dat[, ':='(obsid = 1:.N, tslope = runif(1, 0.005, 0.1), value2 = rnorm(.N), duration = .N, e = rnorm(.N, 0, 0.1)), by = tsid] # add time counter, a random temporal slope for each tsid, a random explanatory variable, and error
    dat[, slope2_abs := abs(slope(value2)), by = tsid]
    dat[, ':='(value1.logit = obsid*tslope + cors[i,effect]*tslope*obsid*slope2_abs + e)] # value1 has a positive slope like dissimilarites and is on a logit scale. Positively related to slope2 with magnitude 'effect'
    dat[, value1 := exp(value1.logit)/(exp(value1.logit)+1)] # inverse-logit for value 1
    dat[, value1.sc := transform01(value1)]
    dat[, duration := max(obsid), by = tsid]
    
    # calc slopes
    slopes <- dat[, .(slope1 = slope(value1), slope2_abs = unique(slope2_abs),
                      se1 = slopese(value1), duration = .N), by = tsid]
    
    # statistical tests
    test2 <- slopes[, summary(lm(slope1 ~ slope2_abs))$coefficients] # a simple regression
    test3 <- slopes[, summary(glmmTMB(slope1 ~ slope2_abs, disp=~se1))$coefficients$cond] # glmmTMB and dispersion

    cors[i, ':='(lm.p = test2[2,4], lm.m = test2[2,1], glmmwgt.p = test3[2,4], glmmwgt.beta = test3[2,1])]    
    
    mod5 <- dat[, glmmTMB(value1.sc ~ obsid + obsid:slope2_abs + (obsid|tsid),
                          family = beta_family(link = 'logit'))] # a one-stage beta model
    if(mod5$fit$convergence==0){ # if the one-stage beta model converged
        test5 <- summary(mod5)$coefficients$cond
        cors[i, ':='(glmmonebeta.p = test5[3,4], glmmonebeta.beta = test5[3,1])]    
    }

    # write out temp file
    write.csv(cors, gzfile(here('temp', 'simulated_ts_temp_pow.csv.gz')), row.names = FALSE)
}



# write out full file ##################
write.csv(cors, gzfile(here('output', 'simulated_ts_pow.csv.gz')), row.names = FALSE)

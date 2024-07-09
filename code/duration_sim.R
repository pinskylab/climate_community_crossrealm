#!/usr/bin/env Rscript

# Code to run the false-positive simulations and analysis for duration_sim.Rmd
# Set up as a script to run in the background on Annotate, e.g.,
# nohup code/duration_sim.R > logs/duration_sim.Rout &
# (this works if code is executable, e.g., chmod u+x code/duration_sim.R)

# print basic info about the job ############################
print(paste('This is process #', Sys.getpid()))
print(Sys.time())


# load libraries ############################
require(data.table)
library(here)
library(glmmTMB, lib.loc = '/usr/lib64/R/library/')


# Functions ################
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
nreps <- 100 # number of datasets at each duration range
ndata <- 1000 # number timeseries per dataset
nstudy <- 50 # number of studies per dataset. timeseries are nested within studies.
studytable <- data.table(tsid = 1:ndata, studyid = rep(1:nstudy, ceiling(ndata/nstudy))[1:ndata]) # table to match study IDs to timeseries IDs

cors <- data.table(n = ndata, minduration = rep(c(10,3,3,3,3), nreps), maxduration = rep(c(10, 12, 27, 52, 102), nreps), cor.p = NA_real_, cor.cor = NA_real_, lm.m = NA_real_,
                   glmmwgt.p = NA_real_, glmmwgt.beta = NA_real_, glmmonegauss.p = NA_real_, glmmonegauss.beta = NA_real_,
                   glmmonebeta.p = NA_real_, glmmonebeta.beta = NA_real_) # holds the summaries from each dataset.

cors[, name := paste0(minduration, '-', maxduration)]
cors[, range := maxduration - minduration]


print(paste0(nrow(cors), ' simulations to make'))
for(i in 1:nrow(cors)){
    cat(i)
    
    # prep to make datasets with timeseries of all the same length
    if(cors[i, minduration == maxduration]){
        len <- cors[i, minduration]
        dat <- data.table(tsid = rep(1:cors[i,n], rep(len, cors[i,n])))
        dat <- merge(dat, studytable)
    }
    
    # prep to make datasets with variable length timeseries
    if(cors[i, minduration != maxduration]){
        mind <- cors[i, minduration]
        maxd <- cors[i, maxduration]
        ndur <- maxd - mind + 1
        dat <- data.table(tsid = rep(1:cors[i,n], rep(mind:maxd, cors[i,n]/ndur)))
        dat <- merge(dat, studytable)
    }
    
    # add response and explanatory variables
    dat[, ':='(obsid = 1:.N, tslope = runif(1, 0, 0.097)), by = tsid] # time counter and a random slope for each tsid
    dat[, studyslope := runif(1, 0.01, 0.05), by = studyid] # a random slope for each studyid
    dat[, ':='(value1.logit = obsid*studyslope + obsid*tslope + rnorm(.N), value2 = rnorm(.N))] # create two random timeseries. value1 has a positive slope like dissimilarities and is on a logit scale.
    dat[, value1 := exp(value1.logit)/(exp(value1.logit)+1)] # inverse logit-transform
    dat[, duration := max(obsid), by = tsid]
    
    # calc slopes
    dat[, slope2_abs := abs(slope(value2)), by = tsid]
    slopes <- dat[, .(slope1 = slope(value1), slope2_abs = unique(slope2_abs), se1 = slopese(value1)), by = tsid]
    
    # statistical tests
    # NEED TO ADD STUDYID
    test <- slopes[, cor.test(slope1, slope2_abs)] # a pearson correlation
    test2 <- slopes[, coef(lm(slope1 ~ slope2_abs))[2]] # a simple regression
    test3 <- slopes[, summary(glmmTMB(slope1 ~ slope2_abs, disp=~se1))$coefficients$cond] # glmmTMB and dispersion: a meta-analysis model

    cors[i, ':='(cor.p = test$p.value, cor.cor = test$estimate, lm.m = test2, glmmwgt.p = test3[2,4], glmmwgt.beta = test3[2,1])]    
    
    mod4 <- dat[, glmmTMB(value1 ~ obsid + obsid:slope2_abs + (obsid|tsid))] # a one-stage Gaussian model
    if(mod4$fit$convergence==0){ # if the one-stage gaussian model converged
        test4 <- summary(mod4)$coefficients$cond
        cors[i, ':='(glmmonegauss.p = test4[3,4], glmmonegauss.beta = test4[3,1])]    
    }
    
    mod5 <- dat[, glmmTMB(value1.sc ~ obsid + obsid:slope2_abs + (obsid|tsid),
                          family = beta_family(link = 'logit'))] # a one-stage beta model
    if(mod5$fit$convergence==0){ # if the one-stage beta model converged
        test5 <- summary(mod5)$coefficients$cond
        cors[i, ':='(glmmonebeta.p = test5[3,4], glmmonebeta.beta = test5[3,1])]    
    }

    # write out temp file
    write.csv(cors, gzfile(here('temp', 'simulated_ts_temp.csv.gz')), row.names = FALSE)
}



# write out full file ##################
write.csv(cors, gzfile(here('output', 'simulated_ts.csv.gz')), row.names = FALSE)

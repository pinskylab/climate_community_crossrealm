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



# Set up simulations #############
set.seed(5)
nreps <- 100 # number of datasets
effects <- c(0.1, 10^(-2/3), 10^(-1/3), 1, 10) # effect sizes to simulate

cors <- data.table(n = 1000, minduration = rep(c(10,3), nreps*length(effects)), maxduration = rep(c(10, 102), nreps*length(effects)), effect = rep(effects, rep(nreps*2, length(effects))), 
                   lm.p = NA_real_, lm.m = NA_real_,
                   glmmwgt.p = NA_real_, glmmwgt.beta = NA_real_,
                   glmmonebeta.p = NA_real_, glmmonebeta.beta = NA_real_) # holds the parameters for and summaries from each dataset. use n = 1000 timeseries per dataset.
#cors2 <- copy(cors)
#cors2[, n:=10000] # add another run with 10k timeseries per dataset
#cors <- rbind(cors, cors2) 
#rm(cors2)

cors[, repID := 1:.N, by = .(n, minduration, maxduration, effect)] # to differentiate the rows for each replicate in a given parameter set
cors[, name := paste0(minduration, '-', maxduration, '-', signif(effect,2))]
cors[, range := maxduration - minduration]

if(file.exists(here('output', 'simulated_ts_pow.csv.gz'))){ # check for previous set of simulations
   corsold <- fread(here('temp', 'simulated_ts_temp_pow.csv.gz'))
   if(!('repID' %in% names(corsold))){
       corsold[, repID := 1:.N, by = .(n, minduration, maxduration, effect)] # add repID column if missing
   }
   corsnew <- merge(corsold, cors[, .(n, minduration, maxduration, effect, name, range, repID)], all = TRUE) # merge in the existing sims so that we add on rather than replace
   cors <- corsnew
}

cors[, todo := is.na(lm.p) | is.na(lm.m) | is.na(glmmwgt.p) | is.na(glmmwgt.beta) | is.na(glmmonebeta.p) | is.na(glmmonebeta.beta)] # mark rows to simulate (or re-sim if some model fits failed)

    
# Run simulations #######################
print(paste0(cors[, sum(todo)], ' simulations to make'))
rowstodo <- cors[, which(todo)]
for(i in rowstodo){
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
    
    print(warnings())
}



# write out full file ##################
cors[, todo := NULL] # remove the todo column
write.csv(cors, gzfile(here('output', 'simulated_ts_pow.csv.gz')), row.names = FALSE)

print(Sys.time())

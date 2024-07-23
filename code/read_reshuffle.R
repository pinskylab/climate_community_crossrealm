# Read in the model fits to data with reshuffled Tchange
# Output a summary data table

### Functions -----------
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(glmmTMB) # for beta regression
library(scales) # for defining signed sqrt axis transformation
library(here)

# Process the reshuffled model coefficient results from fit_turnover_GLMM_reshuffle.R

### Read in distribution of model coefficients from model fits ----------------
# Takes a while.
# read in previous file if desired (eg, might have a subset of the reshuffles so that not all new ones need to be read in)
# coeftable <- fread('output/coeftable_reshuffle_modTchangexYearxRealmJtu.csv') # read in previous written version?

# read in each reshuffled model fit. takes about 3 seconds per model (slow)
mods <- list.files(path = 'temp', pattern = glob2rx('modTchangexYearxRealmJtu_reshuff*.rds'), full.names=TRUE)
n <- length(mods)
n # number of models
exists('coeftable')
if(!exists('coeftable')) coeftable <- data.table(type = rep('null', 0), shuffID = rep(NA_integer_,0), terrwarm = rep(NA_real_,0), terrcool = rep(NA_real_,0), marwarm = rep(NA_real_,0), marcool = rep(NA_real_,0), freshwarm = rep(NA_real_,0), freshcool = rep(NA_real_,0)) # table to hold coefficients. only create if it doesn't yet exist (to allow adding to it)
nrow(coeftable)
for(i in 1:n){
    id <- as.numeric(gsub('temp/modTchangexYearxRealmJtu_reshuff|.rds', '', mods[i]))
    cat(paste0(' ', id))
    if(!(id %in% coeftable$shuffID)){ # only read in model file if not already in coeftable
        cat('*') # small marker for files that get read in
        temp <- readRDS(mods[i]) # slow (~ 3 sec). faster would be to read summaries from .Rout files
        if(temp$sdr$pdHess){ # checks positive definite Hessian. skip if not true
            if(is.numeric(AIC(temp))){ # check for an AIC value. skip if not true
                coefs <- fixef(temp)$cond
                
                newline <- data.table(type = 'null',
                                      shuffID = id,
                                      terrwarm =  coefs["duration:tempchange_abs.sc"] + coefs["duration:tsign1:tempchange_abs.sc"] + coefs["duration:REALMTerrestrial:tempchange_abs.sc"] + coefs["duration:REALMTerrestrial:tsign1:tempchange_abs.sc"],
                                      terrcool =  coefs["duration:tempchange_abs.sc"] + coefs["duration:REALMTerrestrial:tempchange_abs.sc"],
                                      marwarm =   coefs["duration:tempchange_abs.sc"] + coefs["duration:tsign1:tempchange_abs.sc"] + coefs["duration:REALMMarine:tempchange_abs.sc"] + coefs["duration:REALMMarine:tsign1:tempchange_abs.sc"],
                                      marcool =   coefs["duration:tempchange_abs.sc"] + coefs["duration:REALMMarine:tempchange_abs.sc"],
                                      freshwarm = coefs["duration:tempchange_abs.sc"] + coefs["duration:tsign1:tempchange_abs.sc"],
                                      freshcool = coefs["duration:tempchange_abs.sc"])
                coeftable <- rbind(coeftable, newline)
            } else {
                print(paste('No AIC for shuffID', id))
            }
        } else {
            print(paste('Non-positive definite Hessian for shuffID', id))
        }
    }
}

if(!('obs' %in% coeftable$type)){ # only read in observed model file if not already in coeftable
    mod_obs <- readRDS('temp/modTchangexYearxRealmJtu.rds')
    coefs_obs <- fixef(mod_obs)$cond
    coeftable <- rbind(coeftable, data.table(type = 'obs', shuffID = NA_integer_,
                                             terrwarm =  coefs_obs["duration:tempchange_abs.sc"] + coefs_obs["duration:tsign1:tempchange_abs.sc"] + coefs_obs["duration:REALMTerrestrial:tempchange_abs.sc"] + coefs_obs["duration:REALMTerrestrial:tsign1:tempchange_abs.sc"],
                                             terrcool =  coefs_obs["duration:tempchange_abs.sc"] + coefs_obs["duration:REALMTerrestrial:tempchange_abs.sc"],
                                             marwarm =   coefs_obs["duration:tempchange_abs.sc"] + coefs_obs["duration:tsign1:tempchange_abs.sc"] + coefs_obs["duration:REALMMarine:tempchange_abs.sc"] + coefs_obs["duration:REALMMarine:tsign1:tempchange_abs.sc"],
                                             marcool =   coefs_obs["duration:tempchange_abs.sc"] + coefs_obs["duration:REALMMarine:tempchange_abs.sc"],
                                             freshwarm = coefs_obs["duration:tempchange_abs.sc"] + coefs_obs["duration:tsign1:tempchange_abs.sc"],
                                             freshcool = coefs_obs["duration:tempchange_abs.sc"]))
}

coeftable <- coeftable[order(type, shuffID), ] # nice order
nrow(coeftable)

sum(duplicated(coeftable[, .(terrwarm, terrcool, marwarm, marcool, freshwarm, freshcool)])) # check that a model hasn't be used 2x. Should be 0

# write out
write.csv(coeftable, file = 'output/coeftable_reshuffle_modTchangexYearxRealmJtu.csv', row.names=FALSE)


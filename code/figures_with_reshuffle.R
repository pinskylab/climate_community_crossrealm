# Making tables and figures with the model fits to data with reshuffled Tchange

### Functions -----------
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(glmmTMB) # for beta regression
library(scales) # for defining signed sqrt axis transformation
library(here)





### Read in distribution of model coefficients from model fits ----------------
# Only run this if needed. Takes a while.
# read in previous file if desired (eg, might have a subset of the reshuffles so that not all new ones need to be read in)
# coeftable <- fread('output/coeftable_reshuffle_modOBsdTMERtsRealmtsigninitAllJtu.csv') # read in previous written version?

# read in each reshuffled model fit. takes about 3 seconds per model (slow)
mods <- list.files(path = 'temp', pattern = glob2rx('modOBsdTMERtsRealmtsigninitAllJtu_reshuff*.rds'), full.names=TRUE)
n <- length(mods)
n # number of models
exists('coeftable')
if(!exists('coeftable')) coeftable <- data.table(type = rep('null', 0), shuffID = rep(NA_integer_,0), terrwarm = rep(NA_real_,0), terrcool = rep(NA_real_,0), marwarm = rep(NA_real_,0), marcool = rep(NA_real_,0), freshwarm = rep(NA_real_,0), freshcool = rep(NA_real_,0)) # table to hold coefficients. only create if it doesn't yet exist (to allow adding to it)
nrow(coeftable)
for(i in 1:n){
    id <- as.numeric(gsub('temp/modOBsdTMERtsRealmtsigninitAllJtu_reshuff|.rds', '', mods[i]))
    cat(paste0(' ', id))
    if(!(id %in% coeftable$shuffID)){ # only read in model file if not already in coeftable
        cat('*') # small marker for files that get read in
        temp <- readRDS(mods[i]) # slow (~ 3 sec). faster would be to read summaries from .Rout files
        if(temp$sdr$pdHess){ # checks positive definite Hessian. skip if not true
            if(is.numeric(AIC(temp))){ # check for an AIC value. skip if not true
            coefs <- fixef(temp)$cond
            
            newline <- data.table(type = 'null',
                                  shuffID = id,
                                  terrwarm = coefs["duration:REALMTerrestrial:tsign1:tempchange_abs.sc"],
                                  terrcool = coefs["duration:REALMTerrestrial:tsign-1:tempchange_abs.sc"],
                                  marwarm = coefs["duration:REALMMarine:tsign1:tempchange_abs.sc"],
                                  marcool = coefs["duration:REALMMarine:tsign-1:tempchange_abs.sc"],
                                  freshwarm = coefs["duration:REALMFreshwater:tsign1:tempchange_abs.sc"],
                                  freshcool = coefs["duration:REALMFreshwater:tsign-1:tempchange_abs.sc"])
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
    mod_obs <- readRDS('temp/modOBsdTMERtsRealmtsigninitAllJtu.rds')
    coefs_obs <- fixef(mod_obs)$cond
    coeftable <- rbind(coeftable, data.table(type = 'obs', shuffID = NA_integer_,
                                             terrwarm = coefs_obs["duration:REALMTerrestrial:tsign1:tempchange_abs.sc"],
                                             terrcool = coefs_obs["duration:REALMTerrestrial:tsign-1:tempchange_abs.sc"],
                                             marwarm = coefs_obs["duration:REALMMarine:tsign1:tempchange_abs.sc"],
                                             marcool = coefs_obs["duration:REALMMarine:tsign-1:tempchange_abs.sc"],
                                             freshwarm = coefs_obs["duration:REALMFreshwater:tsign1:tempchange_abs.sc"],
                                             freshcool = coefs_obs["duration:REALMFreshwater:tsign-1:tempchange_abs.sc"]))
}

coeftable <- coeftable[order(type, shuffID), ] # nice order
nrow(coeftable)

sum(duplicated(coeftable[, .(terrwarm, terrcool, marwarm, marcool, freshwarm, freshcool)])) # check that a model hasn't be used 2x. Should be 0

# write out
write.csv(coeftable, file = 'output/coeftable_reshuffle_modOBsdTMERtsRealmtsigninitAllJtu.csv', row.names=FALSE)

# Analyze distribution of model coefficients --------------------------
# read in the pre-extracted coefficients (previous section)
coeftable <- fread('output/coeftable_reshuffle_modOBsdTMERtsRealmtsigninitAllJtu.csv') # read in previous written version?
nrow(coeftable) # should be >1000

# plot distributions
par(mfrow=c(2,3))
coeftable[type == 'null', ][1:1000, hist(terrwarm, xlim = range(coeftable$terrwarm))]
coeftable[type == 'obs', abline(v=terrwarm, col = 'red')]

coeftable[type == 'null', ][1:1000, hist(terrcool, xlim = range(coeftable$terrcool))]
coeftable[type == 'obs', abline(v=terrcool, col = 'red')]

coeftable[type == 'null', ][1:1000, hist(marwarm, xlim = range(coeftable$marwarm))]
coeftable[type == 'obs', abline(v=marwarm, col = 'red')]

coeftable[type == 'null', ][1:1000, hist(marcool, xlim = range(coeftable$marcool))]
coeftable[type == 'obs', abline(v=marcool, col = 'red')]

coeftable[type == 'null', ][1:1000, hist(freshwarm, xlim = range(coeftable$freshwarm))]
coeftable[type == 'obs', abline(v=freshwarm, col = 'red')]

coeftable[type == 'null', ][1:1000, hist(freshcool, xlim = range(coeftable$freshcool))]
coeftable[type == 'obs', abline(v=freshcool, col = 'red')]

# how many reshuffles >= obs?
coeftable[type=='null', sum(terrwarm >= coeftable[type=='obs', terrwarm])]
coeftable[type=='null', sum(terrcool >= coeftable[type=='obs', terrcool])]
coeftable[type=='null', sum(marwarm >= coeftable[type=='obs', marwarm])]
coeftable[type=='null', sum(marcool >= coeftable[type=='obs', marcool])]
coeftable[type=='null', sum(freshwarm >= coeftable[type=='obs', freshwarm])]
coeftable[type=='null', sum(freshcool >= coeftable[type=='obs', freshcool])]

# empirical p-values as (x+1)/(n+1)
# use first 1000
print(paste0('Terr warm p=', coeftable[type=='null',][1:1000, signif((sum(terrwarm >= coeftable[type=='obs', terrwarm])+1)/(.N+1),3)]))
print(paste0('Terr cool p=', coeftable[type=='null',][1:1000, signif((sum(terrcool >= coeftable[type=='obs', terrcool])+1)/(.N+1),3)]))
print(paste0('Mar warm p=', coeftable[type=='null',][1:1000, signif((sum(marwarm >= coeftable[type=='obs', marwarm])+1)/(.N+1),3)]))
print(paste0('Mar cool p=', coeftable[type=='null',][1:1000, signif((sum(marcool >= coeftable[type=='obs', marcool])+1)/(.N+1),3)]))
print(paste0('Fresh warm p=', coeftable[type=='null',][1:1000, signif((sum(freshwarm >= coeftable[type=='obs', freshwarm])+1)/(.N+1),3)]))
print(paste0('Fresh cool p=', coeftable[type=='null',][1:1000, signif((sum(freshcool >= coeftable[type=='obs', freshcool])+1)/(.N+1),3)]))

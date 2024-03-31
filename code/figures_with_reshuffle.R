# Making tables and figures with the bootstrap model fits

### Functions -----------
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(glmmTMB) # for beta regression
library(scales) # for defining signed sqrt axis transformation
library(here)





### Distribution of model coefficients ----------------
# read in each reshuffled model fit. takes about 3 seconds per model (slow)
mods <- list.files(path = 'temp', pattern = glob2rx('modOBsdTMERtsRealmtsigninitAllJtu_reshuff*.rds'), full.names=TRUE)
n <- length(mods)
n # number of models
exists('coeftable')
# coeftable <- fread('output/coeftable_reshuffle_modOBsdTMERtsRealmtsigninitAllJtu.csv') # read in previous written version?
if(!exists('coeftable')) coeftable <- data.table(type = rep('null', 0), shuffID = rep(NA_integer_,0), terrwarm = rep(NA_real_,0), terrcool = rep(NA_real_,0), marwarm = rep(NA_real_,0), marcool = rep(NA_real_,0), freshwarm = rep(NA_real_,0), freshcool = rep(NA_real_,0)) # table to hold coefficients. only create if it doesn't yet exist (to allow adding to it)
nrow(coeftable)
for(i in 1:n){
    id <- as.numeric(gsub('temp/modOBsdTMERtsRealmtsigninitAllJtu_reshuff|.rds', '', mods[i]))
    cat(paste0(' ', id))
    if(!(id %in% coeftable$shuffID)){ # only read in model file if not already in coeftable
        cat('*') # small marker for files that get read in
        temp <- readRDS(mods[i]) # slow (~ 3 sec). faster would be to read summaries from .Rout files
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

# write out
write.csv(coeftable, file = 'output/coeftable_reshuffle_modOBsdTMERtsRealmtsigninitAllJtu.csv', row.names=FALSE)

# plot distributions
par(mfrow=c(2,3))
coeftable[type == 'null', hist(terrwarm, xlim = range(coeftable$terrwarm))]
coeftable[type == 'obs', abline(v=terrwarm, col = 'red')]

coeftable[type == 'null', hist(terrcool, xlim = range(coeftable$terrcool))]
coeftable[type == 'obs', abline(v=terrcool, col = 'red')]

coeftable[type == 'null', hist(marwarm, xlim = range(coeftable$marwarm))]
coeftable[type == 'obs', abline(v=marwarm, col = 'red')]

coeftable[type == 'null', hist(marcool, xlim = range(coeftable$marcool))]
coeftable[type == 'obs', abline(v=marcool, col = 'red')]

coeftable[type == 'null', hist(freshwarm, xlim = range(coeftable$freshwarm))]
coeftable[type == 'obs', abline(v=freshwarm, col = 'red')]

coeftable[type == 'null', hist(freshcool, xlim = range(coeftable$freshcool))]
coeftable[type == 'obs', abline(v=freshcool, col = 'red')]


# empirical p-values
print(paste0('Terr warm p=', coeftable[, signif(sum(terrwarm >= coeftable[type=='obs', terrwarm])/.N,2)]))
print(paste0('Terr cool p=', coeftable[, signif(sum(terrcool >= coeftable[type=='obs', terrcool])/.N,2)]))
print(paste0('Mar warm p=', coeftable[, signif(sum(marwarm >= coeftable[type=='obs', marwarm])/.N,2)]))
print(paste0('Mar cool p=', coeftable[, signif(sum(marcool >= coeftable[type=='obs', marcool])/.N,2)]))
print(paste0('Fresh warm p=', coeftable[, signif(sum(freshwarm >= coeftable[type=='obs', freshwarm])/.N,2)]))
print(paste0('Fresh cool p=', coeftable[, signif(sum(freshcool >= coeftable[type=='obs', freshcool])/.N,2)]))

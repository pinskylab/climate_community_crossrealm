#!/usr/bin/Rscript --vanilla

# Script to fit Tchange glmmTMB model with ordbeta errors and reshuffled Tchange
# Set up to be run on the command line for multiple reshuffles
# Arguments are the min and max reshuffling IDs
# nohup code/fit_turnover_GLMM_reshuffle.R 1 10 > logs/fit_turnover_GLMM_reshuffle1-10.Rout &
# (this works if code is executable, e.g., chmod u+x code/fit_turnover_GLMM_reshuffle.R)
# Note: this requires a newer version of glmmTMB, e.g, 1.1.8 (installed on Annotate2 but not Annotate)

# Read command line arguments ############

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 2)
    stop("Have to specify a min and max reshuffling ID", call. = FALSE)
if (length(args) > 2)
    stop("Have to specify only a min and max reshuffling ID", call. = FALSE)

shuffIDs <- as.numeric(args[1]):as.numeric(args[2]) # used for random seeds

# print basic info about the job ############################

print(paste('This is script fit_turnover_GLMM_reshuffle.R'))
print(paste('This is process #', Sys.getpid()))
print(paste('Reshuffle IDs', paste(shuffIDs, collapse = ',')))
print(Sys.time())


# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB) # for ME ordbeta models
library(here)


# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd

trendsall[, tsign := as.factor(tsign)]

# calculate initial dissimilarities
trendsall[, minduration := min(duration), by = rarefyID]
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu)), by = .(rarefyID)] # initial dissimilarities
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])

## Choose dataset
iallJtu <-
    trendsall[, complete.cases(
        Jtu,
        tempchange_abs.sc,
        REALM,
        tempave.sc,
        duration,
        microclim.sc,
        human_bowler.sc
    )]
trendsalltrim <- trendsall[iallJtu,] # trim


print('Warnings so far:')
print(warnings())

# Reshuffling loop ###########################
tchange <- trendsalltrim[, .(tempchange_abs.sc = unique(tempchange_abs.sc), tsign = unique(tsign)), by = rarefyID]
trendsalltrim[, ':='(tempchange_abs.sc = NULL, tsign = NULL)]

for(i in 1:length(shuffIDs)){
    print(paste('Starting shuffle ID', shuffIDs[i]))
    print(Sys.time())
    
    # reshuffle Tchange across rarefyIDs
    set.seed(shuffIDs[i])
    tchange_shuff <- data.table(rarefyID = sample(tchange$rarefyID), tempchange_abs.sc = tchange$tempchange_abs.sc, tsign = tchange$tsign) # reshuffle rarefyIDs
    temp <- merge(trendsalltrim, tchange_shuff, all = TRUE, by = 'rarefyID') # add Tchange variables back onto data frame with reshuffled rarefyIDs
    
    # fit Tchange x Year x Realm model
    print(paste(nrow(temp), 'data points'))
    print(paste(temp[, length(unique(STUDY_ID))], 'studies'))
    print(paste(temp[, length(unique(rarefyID))], 'time series'))
    
    try(
        {
            mod <- glmmTMB(
                Jtu ~ duration +
                    REALM +
                    tsign +
                    Jtu.init:duration +
                    REALM:tsign:tempchange_abs.sc +
                    REALM:tsign:tempchange_abs.sc:duration +
                    (duration | STUDY_ID / rarefyID),
                data = temp,
                family = ordbeta(link = 'logit'),
                dispformula = ~ REALM
            )
            print('finished fitting')
            
            # print and save results
            print(summary(mod))
            outfile <- paste0('temp/modOBsdTMERtsRealmtsigninitAllJtu_reshuff', shuffIDs[i], '.rds')
            saveRDS(mod, file = outfile)
            print(paste0('saved ', outfile))
        })
    print(Sys.time())
    print(warnings())
}   



print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

#!/usr/bin/Rscript --vanilla

# Script to fit and predict from glmmTMB models with ordbeta errors and downsampling
# Downsampling picks t pairwise dissimilarities from each time series of length t species composition observations
# Set up to be run on the command line for multiple downsamplings at a time
# Argument is model name to run (see below for options) and min and max downsampling ID (used as random number seed), e.g.
# nohup code/fit_pred_turnover_GLMM_downsamp.R modOBInitAllJtu 1 20> logs/fit_pred_turnover_GLMM_downsamp_modOBInitAllJtu_1-20.Rout &
# Note: this requires a newer version of glmmTMB, e.g, 1.1.8 (installed on Annotate2 but not Annotate)
#
# Beware that this script will use up to all the cpus on a machine.
# Not clear why this happens, since code below tries to force data.table and glmmTMB to use one core

# Read command line arguments ############

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 3)
    stop("Have to specify a model to fit and the min and max downsampling ID", call. = FALSE)
if (length(args) > 3)
    stop("Have to specify only 1 model to fit and min and max downsampling ID", call. = FALSE)
fitmod <- args[1] # which model to fit
bootID <- as.numeric(args[2]):as.numeric(args[3]) # used for random seed

# print basic info about the job ############################

print(paste('This is script turnover_GLMM_fit_boot.R'))
print(paste('This is process #', Sys.getpid()))
print(paste('Downsampling IDs', paste(bootID, collapse=',')))
print(Sys.time())


# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB) # for ME ordbeta models
library(here)
source(here('code', 'util.R'))

setDTthreads(threads = 1) # to avoid multithreading on data.table reading

# load and prep data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv')) # The scaling factors for predictions. From assemble_turnover_covariates.Rmd

trendsall[, tsign := as.factor(tsign)]

# calculate initial dissimilarities
trendsall[, minduration := min(duration), by = rarefyID]
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu)), by = .(rarefyID)] # initial dissimilarities
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])

# Choose dataset
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


print('Warnings so far:')
print(warnings())

# Downsampling loop ###########################

for(i in 1:length(bootID)){
    MATCHMOD <- FALSE # indicator to check if the argument matched a model name

    print(paste('Starting downsampling ID', bootID[i]))
    print(Sys.time())
    
    # downsample the dataset to t dissimilarities per rarefyID
    thisiallJtu <- downsampRarefyID(dat = trendsall, inds = iallJtu, seed = bootID[i])
    
    # Tchange x Year x Realm model
    if (fitmod == 'modTchangexYearxRealmJtu') {
        if (MATCHMOD)
            stop('Model name matched more than one model!')
        print(paste(sum(thisiallJtu), 'data points'))
        print(paste(trendsall[thisiallJtu, length(unique(STUDY_ID))], 'studies'))
        print(paste(trendsall[thisiallJtu, length(unique(rarefyID))], 'time series'))
        mod <- glmmTMB(
            Jtu ~ Jtu.init*duration +
                REALM*tsign*tempchange_abs.sc*duration +
                (duration | STUDY_ID / rarefyID),
            data = trendsall[thisiallJtu, ],
            family = ordbeta(link = 'logit'),
            dispformula = ~ REALM,
            control = glmmTMBControl(parallel = 1) # prevent multithreading
        )
        MATCHMOD <- TRUE
    }
    
    # Tchange x Tave x Year x Realm model
    if (fitmod == 'modTchangexTavexYearxRealmJtu') {
        if (MATCHMOD)
            stop('Model name matched more than one model!')
        print(paste(sum(thisiallJtu), 'data points'))
        print(paste(trendsall[thisiallJtu, length(unique(STUDY_ID))], 'studies'))
        print(paste(trendsall[thisiallJtu, length(unique(rarefyID))], 'time series'))
        mod <- glmmTMB(
            Jtu ~ Jtu.init*duration +
                REALM*tsign*tempave.sc*tempchange_abs.sc*duration +
                (duration | STUDY_ID / rarefyID),
            data = trendsall[thisiallJtu, ],
            family = ordbeta(link = 'logit'),
            dispformula = ~ REALM,
            control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3), # higher iteration limit to facilitate convergence
                                     parallel = 1) # prevent multithreading
        )
        MATCHMOD <- TRUE
    }
    
    print(paste('Finished model fitting:', Sys.time()))
    
    doSensitivity <- FALSE # default is not to calculate sensitivity to Tave
    if(grepl('TchangexTave|Human|Microclim', fitmod)){ # Tchange x Tave, human, or microclim model
        doSensitivity <- TRUE # calculate sensitivity to Tave if we loaded Tave model
    } 
    doCovar <- FALSE
    if(grepl('Human|Microclim', fitmod)){
        doCovar <- TRUE # some extra steps if a covariate model
    }
 
    # print and save results, then make predictions
    if (MATCHMOD == FALSE)
        stop("Model name did not match anything", call. = FALSE)
    if (MATCHMOD) {
        
        ## print summary and save model --------------
        print(summary(mod))
        outfile <- paste0('temp/', fitmod, '_boot', bootID[i], '.rds')
        saveRDS(mod, file = outfile)
        print(paste0('saved ', outfile))
        print(Sys.time())
        print(warnings())
        
        # set up file names and paramters for predictions
        out_preds <- paste0('preds_', fitmod, '_boot', bootID[i], '.rds')
        out_slopes <- paste0('slopes_', fitmod, '_boot', bootID[i], '.rds')
        out_sensitivity <- paste0('sensitivity_', fitmod, '_boot', bootID[i], '.rds')
        n = 1000 # number of resamples to do for each timeseries
        
        # set up prediction frame
        if(!doCovar) newdat <- data.table(expand.grid(tempave = seq(-10, 36, by = 1), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial', 'Freshwater')))
        if(doCovar){ # a bit extra for the covariate models
            newdat <- data.table(expand.grid(tempave = c(0, 8, 10, 13, 30), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial'), microclim.sc = seq(-2, 2, length.out=10)))  
            newdat[, ':='(human_bowler.sc = microclim.sc)]
            newdat[, microclim := unscaleme(microclim.sc, 'microclim.sc')]
            newdat[, human_bowler := unscaleme(human_bowler.sc, 'human_bowler.sc')]
            
            # constrain human covariate to >=0 and <=10, since that is the range over which it is defined
            if(newdat[, min(human_bowler)<0]){
                newdat[abs(human_bowler.sc - -2)<0.001, human_bowler := 0]
                newdat[, human_bowler.sc := scaleme(human_bowler, 'human_bowler.sc')]
            }
            if(newdat[, max(human_bowler)>10]){
                newdat[abs(human_bowler.sc - 2)<0.01, human_bowler := 10]
                newdat[, human_bowler.sc := scaleme(human_bowler, 'human_bowler.sc')]
            }
        } 
        newdat$STUDY_ID <- 1
        newdat$rarefyID <- 1
        newdat$Jtu.init <- 0.5
        newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
        newdat[, tempchange_abs := abs(tempchange)]
        newdat[, tsign := signneg11(tempchange)]
        newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]

        # Make predictions -------------------
        # predict with SEs
        preds <- predict(mod, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
        newdat$Jtu <- preds$fit
        newdat$Jtu.se <- preds$se.fit
        print('finished  predictions')
  
        # Write dissimilarities
        saveRDS(newdat, file = here('temp', out_preds))
        print(paste('Wrote', out_preds, ':', Sys.time()))
        
        # Slope calculations -------------------
        if(!doCovar){
            slopes <- newdat[, slopesamp(n, duration, Jtu, Jtu.se), 
                         by = .(tempave, tempchange, REALM)]
        }
        if(doCovar){
            slopes <- newdat[, slopesamp(n, duration, Jtu, Jtu.se), 
                                       by = .(tempave, tempchange, microclim, human_bowler, REALM)]    
        }

        # Write slopes
        saveRDS(slopes, file = here('temp', out_slopes))
        print(paste('Wrote', out_slopes, ':', Sys.time()))
        
        # Calculate sensitivity of turnover rates to Tave ----------------------
        # if we have a model with Tave
        if(doSensitivity){
            if(!doCovar){
                sensitivity <- slopes[, slopesamp(n, abs(tempchange), slope, slope.se, colnames = c('sensitivity', 'sensitivity.se')), 
                                  by = .(tempave, tsign = sign(tempchange), REALM)]
            }
            if(doCovar){
                slopes <- slopes[tempave == 13 & tempchange > 0, ] # pick tempave and warming
                sensitivity <- slopes[, slopesamp(n, tempchange, slope, slope.se, colnames = c('sensitivity', 'sensitivity.se')), 
                                                 by = .(microclim, human_bowler, REALM)]
            }
            
            # write sensitivities
            saveRDS(sensitivity, file = here('temp', out_sensitivity))
            print(paste('Wrote', out_sensitivity, ':', Sys.time()))
        }
    }
}
print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

## Make predictions data.frame from modrawTsdTTRealmtsignmicroclimhumanAllJtu
## Adapted from modrawTsdTTRealmCovariateAllJtu.R
#
# run as
# nohup Rscript --vanilla code/pred_GLMMmodrawTsdTTRealmmicroclimhumanAllJtu.R > logs/pred_GLMMmodrawTsdTTRealmmicroclimhumanAllJtu.Rout &
# or run by hand to, for example, start after the predictions have been made (line 93). would need to read the predictions in by hand.
# set to run on Annotate


print(paste('This is process #', Sys.getpid()))
print(Sys.time())

# load functions and data ----------------

# needed to run this from the Annotate R console. Not needed in RStudio on Annotate. Not clear why.
if(Sys.getenv('RSTUDIO')=='' & Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){
    dyn.load('/usr/local/lib64/libgfortran.so.5')
}

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(here) # for relative paths
library(lavaan) # for SEM models of slopes

scaleme <- function(x, nm){
    if(!(nm %in% scalingall[,var])) stop('nm not found in scalingall')
    if(scalingall[var==nm, log]){
        x.sc <- (log(x + scalingall[var==nm, plus]) - scalingall[var == nm, center]) / scalingall[var == nm, scale]  
    } else {
        x.sc <- (x  + scalingall[var==nm, plus] - scalingall[var == nm, center]) / scalingall[var == nm, scale]
    }
    return(x.sc)
}
unscaleme <- function(x.sc, nm){
    if(!(nm %in% scalingall[,var])) stop('nm not found in scalingall')
    if(scalingall[var==nm, log]){
        x <- exp(x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center]) - scalingall[var==nm, plus]
    } else {
        x <- x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center] - scalingall[var==nm, plus]
    }
    return(x)
}

# sign of temperature change
signneg11 <- function(x){ # assign 0 a sign of 1 so that there are only 2 levels
    out <- sign(x)
    out[out == 0] <- 1
    return(out)
}


# The scaling factors
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv'))


# The model
modrawTsdTTRealmtsignmicroclimhumanAllJtu <- readRDS('temp/modrawTsdTTRealmtsignmicroclimhumanAllJtu.rds') # has microclimates and human
print('models loaded')


# Make predictions -------------------------
# set up prediction frame
newdat <- data.table(expand.grid(tempave = c(0, 8, 10, 13, 30), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial'), microclim.sc = c(-2, 2), human_bowler.sc = c(-2, 2)))
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
newdat[, tempchange_abs := abs(tempchange)]
newdat[, tsign := signneg11(tempchange)]
newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]
newdat[, microclim := unscaleme(microclim.sc, 'microclim.sc')]
newdat[, human_bowler := unscaleme(human_bowler.sc, 'human_bowler.sc')]

# predict
preds.microclimhuman <- predict(modrawTsdTTRealmtsignmicroclimhumanAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.microclimhuman <- preds.microclimhuman$fit
newdat$Jtu.sc.microclimhuman.se <- preds.microclimhuman$se.fit
print('finished microclimhuman predictions')

# write out
saveRDS(newdat, file = here('temp', 'preds_rawTsdTTRealmtsignmicroclimhuman.rds'))
print(paste('Wrote preds_rawTsdTTRealmtsignmicroclimhuman.rds:', Sys.time()))


# Calculate slopes -----------------------------
# calculate slopes and SE of the slope using latent variables (since predictions have SE)
mods.microclimhuman <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.microclimhuman\nJtu.sc.microclimhuman ~~ ', mean(Jtu.sc.microclimhuman.se), '*Jtu.sc.microclimhuman'), 
                                                    data.frame(Jtu.sc.microclimhuman, duration)))), 
                         by = .(tempave, tempchange, tempchange_abs, tsign, microclim, human_bowler, REALM)]
print('finished microclimhuman SEM')


# extract slopes and SEs
slopes.microclimhuman <- mods.microclimhuman[, parameterEstimates(mod[[1]]), 
                                   by = .(tempave, tempchange, tempchange_abs, tsign, 
                                          microclim, human_bowler, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                                                            .(tempave, tempchange, tempchange_abs, 
                                                                                              microclim, human_bowler, REALM, 
                                                                                              slope_microclimhuman = est, slope_microclimhuman.se = se)]

# save slopes
saveRDS(slopes.microclimhuman, file = here('temp', 'slopes_rawTsdTTRealmtsignmicroclimhuman.rds'))
print(paste('Wrote slopes_rawTsdTTRealmtsignmicroclimhuman.rds:', Sys.time()))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

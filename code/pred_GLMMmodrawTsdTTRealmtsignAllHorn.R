## Make predictions data.frame from the tsign models with environmental temperature (rawT)
## From turnover_vs_temperature_GLMM.Rmd, but set up to run in the background
## so that we can also calculate SEs
#
# run as
# nohup Rscript --vanilla code/pred_GLMMmodrawTsdTTRealmtsignAllHorn.R > logs/pred_GLMMmodrawTsdTTRealmtsignAllHorn.Rout &
# or run by hand to, for example, start after the predictions have been made (line 93). would need to read the predictions in by hand.
# set to run on Annotate


print(paste('This is process #', Sys.getpid()))
print(Sys.time())

### Functions and data ----------------------
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
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv')) # From assemble_turnover_covariates.Rmd

# The models
modrawTsdTTRealmtsignAllHorn <- readRDS(here('temp', 'modrawTsdTTRealmtsignAllHorn.rds')) # tempave, tempchange_abs, realm, tsign. From turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllHorn.R
print('models loaded')

### Make predictions -------------------
# set up prediction frame
newdat <- data.table(expand.grid(tempave = seq(-20, 30, length.out = 100), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial', 'Freshwater')))
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
newdat[, tempchange_abs := abs(tempchange)]
newdat[, tsign := signneg11(tempchange)]
newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]

# predict
preds.realmtsign <- predict(modrawTsdTTRealmtsignAllHorn, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Horn.sc.realmtsign <- preds.realmtsign$fit
newdat$Horn.sc.realmtsign.se <- preds.realmtsign$se.fit
print('finished realmtsign predictions')


### Write dissimilarities ------------------
saveRDS(newdat, file = here('temp', 'preds_rawTsdTTRealmtsignHorn.rds'))


### Slope calculations -------------------
# calculate slopes and SE of the slope using latent variables (since predictions have SE)
mods.realmtsign <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Horn.sc.realmtsign\nHorn.sc.realmtsign ~~ ', mean(Horn.sc.realmtsign.se), '*Horn.sc.realmtsign'), data.frame(Horn.sc.realmtsign, duration)))), 
                     by = .(tempave, tempchange, REALM)] # takes ~30 min
print('finished realmtsign SEM')


# extract slopes and SEs
slopes2 <- mods.realmtsign[, parameterEstimates(mod[[1]]), 
                           by = .(tempave, tempchange, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                               .(tempave, tempchange, REALM, slope_realmtsign = est, slope_realmtsign.se = se)]

### Write slopes -------------
saveRDS(slopes2, file = here('temp', 'slopes_rawTsdTTRealmtsignHorn.rds'))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

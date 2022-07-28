## Make predictions data.frame from the thermal_bias models with environmental temperature (rawT)
## From turnover_vs_temperature_GLMM.Rmd, but set up to run in the background
## so that we can also calculate SEs
#
# run as
# nohup Rscript --vanilla code/pred_GLMMmodrawTsdTTRealmthermal_biasAllJtu.R > logs/pred_GLMMmodrawTsdTTRealmthermal_biasAllJtu.Rout &
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
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv'))

# The models
modrawTsdTTRealmthermal_biassdTAllJtu_thermal_biasdata <- readRDS(here('temp','modrawTsdTTRealmthermal_biassdTAllJtu_thermal_biasdata.rds')) # has thermal bias:tempchange_abs:tsign
modrawTsdTTRealmthermal_biasAllJtu_thermal_biasdata <- readRDS(here('temp','modrawTsdTTRealmthermal_biasAllJtu_thermal_biasdata.rds')) # also adds thermal bias:tempave
print('models loaded')

### Make predictions -------------------
# set up prediction frame
newdat <- data.table(expand.grid(tempave = seq(-20, 30, length.out = 100), tempchange = seq(-1.5, 2, length.out=100), duration = 1:10, REALM = c('Marine', 'Terrestrial'), thermal_bias.sc = c(-2, 2)))
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
newdat[, tempchange_abs := abs(tempchange)]
newdat[, tsign := signneg11(tempchange)]
newdat[, tempchange_abs.sc := scaleme(tempchange_abs, 'tempchange_abs.sc')]
newdat[, thermal_bias := unscaleme(thermal_bias.sc, 'thermal_bias.sc')]

# predict
preds.thermal_biassdT <- predict(modrawTsdTTRealmthermal_biassdTAllJtu_thermal_biasdata, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.thermal_biassdT <- preds.thermal_biassdT$fit
newdat$Jtu.sc.thermal_biassdT.se <- preds.thermal_biassdT$se.fit
print('finished thermal_biassdT predictions')

preds.thermal_bias <- predict(modrawTsdTTRealmthermal_biasAllJtu_thermal_biasdata, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc.thermal_bias <- preds.thermal_bias$fit
newdat$Jtu.sc.thermal_bias.se <- preds.thermal_bias$se.fit
print('finished thermal_bias predictions')


### Write dissimilarities ------------------
saveRDS(newdat, file = here('temp', 'preds_rawTsdTTRealmthermal_bias.rds'))


### Slope calculations -------------------
# calculate slopes and SE of the slope using latent variables (since predictions have SE)
mods.thermal_biassdT <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.thermal_biassdT\nJtu.sc.thermal_biassdT ~~ ', mean(Jtu.sc.thermal_biassdT.se), '*Jtu.sc.thermal_biassdT'), data.frame(Jtu.sc.thermal_biassdT, duration)))), 
                               by = .(tempave, tempchange, thermal_bias, REALM)]
print('finished thermal_biassdT SEM')

mods.thermal_bias <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Jtu.sc.thermal_bias\nJtu.sc.thermal_bias ~~ ', mean(Jtu.sc.thermal_bias.se), '*Jtu.sc.thermal_bias'), data.frame(Jtu.sc.thermal_bias, duration)))), 
                            by = .(tempave, tempchange, thermal_bias, REALM)]
print('finished thermal_bias SEM')


# extract slopes and SEs
slopes.thermal_biassdT <- mods.thermal_biassdT[, parameterEstimates(mod[[1]]), 
                           by = .(tempave, tempchange, thermal_bias, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                               .(tempave, tempchange, thermal_bias, REALM, slope_thermal_biassdT = est, slope_thermal_biassdT.se = se)]
slopes.thermal_bias <- mods.thermal_bias[, parameterEstimates(mod[[1]]), 
                           by = .(tempave, tempchange, thermal_bias, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                               .(tempave, tempchange, thermal_bias, REALM, slope_thermal_bias = est, slope_thermal_bias.se = se)]

slopes2 <- merge(slopes.thermal_biassdT, slopes.thermal_bias)


### Write slopes -------------
saveRDS(slopes2, file = here('temp', 'slopes_rawTsdTTRealmthermal_bias.rds'))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

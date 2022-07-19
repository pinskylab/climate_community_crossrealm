## Make predictions data.frame from the tempchange model (environmental temperature, tempchange)
## From turnover_vs_temperature_GLMM.Rmd, but set up to run in the background
## so that we can also calculate SEs

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(here) # for relative paths

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
    if(!(nm %in% scalingall[,var])) stop(paste0(nm, ' not found in scalingall'))
    if(scalingall[var==nm, log]){
        x <- exp(x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center]) - scalingall[var==nm, plus]
    } else {
        x <- x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center] - scalingall[var==nm, plus]
    }
    return(x)
}

# The scaling factors
scalingall <- fread(here('output', 'turnover_w_covariates_scaling.csv'))

# The model
modrawTdTTRealmAllJtu <- readRDS('temp/modrawTdTTRealmAllJtu.rds') # uses duration and tempchange


# set up prediction frame.
newdat <- data.table(expand.grid(tempave = seq(-20, 30, length.out = 100), tempchange = seq(-1.5, 2.0, length.out = 100), duration = seq(1, 10, length.out=10), REALM = c('Marine', 'Terrestrial', 'Freshwater')))
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave.sc := scaleme(tempave, 'tempave.sc')]
newdat[, tempchange.sc := scaleme(tempchange, 'tempchange.sc')]

# predict
preds <- predict(modrawTdTTRealmAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc <- preds$fit
newdat$Jtu.sc.se <- preds$se.fit

saveRDS(newdat, file = here('temp', 'preds_rawTdTTRealm.rds'))

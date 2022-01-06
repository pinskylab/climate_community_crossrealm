## Make predictions data.frame from the baseline model
## From turnover_vs_temperature_GLMM.Rmd, but set up to run in the background

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models
library(here) # for relative paths

unscaleme <- function(x.sc, nm){
    if(!(nm %in% scalingall[,var])) stop('nm not found in scalingall')
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
modTsdTTRealmAllJtu <- readRDS('temp/modTsdTTRealmAllJtu.rds') # uses duration and abs(tempchange)


# set up prediction frame
newdat <- data.table(expand.grid(tempave_metab.sc = seq(-1.6, 1.5, length.out = 100), tempchange_abs.sc = seq(-0.8, 20, length.out = 100), duration = seq(1, 10, length.out=10), REALM = c('Marine', 'Terrestrial', 'Freshwater')))
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave_metab := unscaleme(tempave_metab.sc, 'tempave_metab.sc')]
newdat[, tempchange_abs := unscaleme(tempchange_abs.sc, 'tempchange_abs.sc')]

# predict
preds <- predict(modTsdTTRealmAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc <- preds$fit
newdat$Jtu.sc.se <- preds$se.fit

saveRDS(newdat, file = here('temp', 'preds_TsdTTRealm.rds'))

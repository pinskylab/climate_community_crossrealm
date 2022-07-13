## Make predictions data.frame from the baseline model (environmental temperature, tempchange_abs)
## From turnover_vs_temperature_GLMM.Rmd, but set up to run in the background
## so that we can also calculate SEs

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
modrawTsdTTRealmAllJtu <- readRDS('temp/modrawTsdTTRealmAllJtu.rds') # uses duration and abs(tempchange)


# set up prediction frame. tempchange_abs set so goes from 0 2.
newdat <- data.table(expand.grid(tempave.sc = seq(-1.6, 1.5, length.out = 100), tempchange_abs.sc = seq(-0.7953, 47, length.out = 100), duration = seq(1, 10, length.out=10), REALM = c('Marine', 'Terrestrial', 'Freshwater')))
newdat$STUDY_ID <- 1
newdat$rarefyID <- 1
newdat[, tempave := unscaleme(tempave.sc, 'tempave.sc')]
newdat[, tempchange_abs := unscaleme(tempchange_abs.sc, 'tempchange_abs.sc')]

# predict
preds <- predict(modrawTsdTTRealmAllJtu, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Jtu.sc <- preds$fit
newdat$Jtu.sc.se <- preds$se.fit

saveRDS(newdat, file = here('temp', 'preds_rawTsdTTRealm.rds'))

#!/usr/bin/Rscript --vanilla


## Make Horn predictions data.frame from sdT and TsdT models with environmental temperature (rawT)
#
# run as
# nohup Rscript --vanilla code/pred_GLMMmodrawXAllHorn.R X > logs/pred_GLMMmodrawXAllHorn_X.Rout &
# where X is an argument (see below)
# or as:
# pred_modrawXAllHorn.sh X1 X2
# which spawns a job for each argument
# or run by hand to, for example, start after the predictions have been made (line 93). would need to read the predictions in by hand.
# set to run on Annotate


print(paste('This is process #', Sys.getpid()))
print(Sys.time())

# read arguments ----------------
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 1)
    stop("Have to specify a model type to predict", call. = FALSE)
if (length(args) > 1)
    stop("Have to specify only 1 model type to predict", call. = FALSE)
predmod <- args[1]



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

# Choose a model
if(predmod == 'modsdTRealmtsignAllHorn'){
    mod <- readRDS(here('temp', 'modsdTRealmtsignAllHorn.rds')) # From turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllHorn.R
    out_preds <- 'preds_modsdTRealmtsignAllHorn.rds'
    out_slopes <- 'slopes_modsdTRealmtsignAllHorn.rds'
    print('model loaded')
} 

if(!exists('mod')) stop('No model loaded. Make sure argument matches one of the model names')




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
preds <- predict(mod, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$Horn.sc <- preds$fit
newdat$Horn.sc.se <- preds$se.fit
print('finished  predictions')

### Write dissimilarities ------------------
saveRDS(newdat, file = here('temp', out_preds))


### Slope calculations -------------------
# calculate slopes and SE of the slope using latent variables (since predictions have SE)
slopemods <- newdat[, .(mod = list(lavaan(paste0('fithat ~ duration\nfithat =~ 1*Horn.sc\nHorn.sc ~~ ', mean(Horn.sc.se), '*Horn.sc'), data.frame(Horn.sc, duration)))), 
                         by = .(tempave, tempchange, REALM)]
print('finished SEM')


# extract slopes and SEs
slopes <- slopemods[, parameterEstimates(mod[[1]]), 
                           by = .(tempave, tempchange, REALM)][lhs == 'fithat' & op == '~' & rhs == 'duration', 
                                                               .(tempave, tempchange, REALM, slope = est, slope.se = se)]

### Write slopes -------------
saveRDS(slopes, file = here('temp', out_slopes))


print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

#!/usr/bin/Rscript --vanilla

# Script to fit glmmTMB models
# Set up to be run on the command line for one model at a time
# Argument is model name to run (see below for options), e.g.
# nohup code/turnover_vs_init_GLMM_fit.R modInitAllJtu > logs/turnover_vs_init_GLMMmodInitAllJtu.Rout &
# (this works if code is executable, e.g., chmod u+x turnover_GLMM_fit.R)
# (otherwise using nohup Rscript ...)

# Read command line arguments ############

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 1)
    stop("Have to specify a model to fit", call. = FALSE)
if (length(args) > 1)
    stop("Have to specify only 1 model to fit", call. = FALSE)
fitmod <- args[1]
MATCHMOD <- FALSE # indicator to check if the argument matched a model name

# print basic info about the job ############################

print(paste('This is process #', Sys.getpid()))
print(Sys.time())



# load libraries ############################

library(data.table) # for handling large datasets
library(glmmTMB, lib.loc = "/usr/lib64/R/library") # for ME models. path is set to run on Annotate.

# load data ############################

# Turnover and covariates assembled by assemble_turnover_covariates.Rmd
trendsall <- fread('output/turnover_w_covariates.csv.gz') # From assemble_turnover_covariates.Rmd
trendsall[, minduration := min(duration), by = rarefyID]

# calc initial dissimilarities
initdiss <- trendsall[duration == minduration, .(Jtu.init = mean(Jtu.sc)), by = .(REALM, rarefyID)] # initial dissimilarities

# add initdiss back to trendsall
trendsall <- merge(trendsall, initdiss[, .(rarefyID, Jtu.init)])
trendsall[, Jtu.init.sc := scale(Jtu.init)]
nrow(trendsall)


# Models ############################

## Choose dataset
iallJtu <-
    trendsall[, complete.cases(
        Jtu.sc
    )]

# Initial dissimilarity model #################################
if (fitmod == 'modInitAllJtu') {
    if (MATCHMOD) stop('Model name matched more than one model!')
    print(paste(sum(iallJtu), 'data points'))
    mod <- glmmTMB(
        Jtu.sc ~ duration +
            Jtu.init.sc:duration +
            (duration | STUDY_ID / rarefyID),
        data = trendsall[iallJtu, ],
        family = beta_family(link = 'logit'),
        dispformula = ~ REALM
    )
    MATCHMOD <- TRUE
}

# print and save results ############################
if (!MATCHMOD) stop("Model name did not match anything", call. = FALSE)
if (MATCHMOD) {
    print(summary(mod))
    saveRDS(mod, file = paste0('temp/', fitmod, '.rds'))
    print(paste0('saved ', fitmod, '.rds'))
    print(Sys.time())
    print(warnings())
}

print(warnings())
print(paste('Ended', Sys.time(), sep = ''))

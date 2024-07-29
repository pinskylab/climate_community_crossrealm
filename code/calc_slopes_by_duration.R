#!/usr/bin/env Rscript

# Calculate slope for Gaussian white noise and temperature time-series
# prep for making figure on duration problem
print(paste('started', Sys.time()))

### Functions -----------
library(data.table) # for handling large datasets
library(here)
source(here('code', 'util.R'))

# load biotime trends
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends[, duration := year2 - year1]
trends <- trends[measure == 'Jtu' & duration>2, ] # trim to Jaccard turnover

# load temperature slopes
scalingall <- fread('output/turnover_w_covariates_scaling.csv') # covariate scaling data. From assemble_turnover_covariates.Rmd
tempchange <- fread('output/turnover_w_covariates.csv.gz') # covariate data. From assemble_turnover_covariates.Rmd
tempchange <- tempchange[, .(tempchange.sc = mean(tempchange.sc, na.rm=TRUE), duration = max(duration)), by = rarefyID] # summarize by rarefyID
tempchange[, tempchange := tempchange.sc * scalingall[var == 'tempchange.sc', scale] + scalingall[var == 'tempchange.sc', center]]

print(paste('data loaded', Sys.time()))

# make slopes of Gaussian white noise timeseries
set.seed(10)
trends[, gauss.slope := calcslopeGauss(duration), by = rarefyID]
print(paste('slopes made', Sys.time()))

# make mean predictions of turnover rate, white noise slopes, and temperature change rate
modloess <- trends[, loess(disstrend~duration)] # loess fit (slow)
print(paste('loess fit to disstrend', Sys.time()))

predsloess <- data.table(duration = 2:118)
predsloess[, c('disstrend', 'se') := predict(modloess, newdata = predsloess, se.fit = TRUE)]
print(paste('disstrend loess predictions made', Sys.time()))

modloessgauss <- trends[, loess(gauss.slope~duration)] # loess fit
print(paste('loess fit to gaussian slopes', Sys.time()))

predsloess[, c('gauss.slope', 'gauss.se') := predict(modloessgauss, newdata = predsloess, se.fit = TRUE)]
print(paste('gaussian loess predictions made', Sys.time()))


modloesstemp <- tempchange[, loess(tempchange~duration)] # loess fit
print(paste('loess fit to tempchange', Sys.time()))

predsloess[, c('tempchange', 'tempchange.se') := predict(modloesstemp, newdata = predsloess, se.fit = TRUE)]
print(paste('tempchange loss predictions made', Sys.time()))

# write out
write.csv(trends[, .(rarefyID, duration, gauss.slope)], file = gzfile('output/gaussian_slopes.csv.gz'), row.names = FALSE)
write.csv(predsloess, file = gzfile('output/slopes_vs_duration.csv.gz'), row.names = FALSE)

print(paste('finished', Sys.time()))

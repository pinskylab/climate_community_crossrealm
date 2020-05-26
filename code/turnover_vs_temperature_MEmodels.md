Drivers of variation in the community response to temperature change
across realms
================

Collaborators: Shane Blowes, Jon Chase, Helmut Hillebrand, Michael
Burrows, Amanda Bates, Uli Brose, Benoit Gauzens, Laura Antao
Assistance: Katherine Lew, Josef Hauser

# Introduction

  - Climate change is driving a widespread reorganization of ecological
    communities around the world (Parmsesan & Yohe 2003, Poloczanska et
    al. 2013),
  - but the impacts of climate change vary substantially from one
    location to another and among taxa (Molinos et al. 2016 NCC, Antao
    et al. 2020 NEE).
  - Community reorganization is substantially more common than an
    aggregate loss or gain of species (Dornelas et al. 2014 Science,
    Blowes et al. 2019 Science, Hillebrand et al. 2017 J Appl Ecol)
  - There are many hypotheses for why some communities are more
    sensitive to warming than others, including differences in
      - metabolic rates (Dillon et al. 2010 Nature),
      - thermal physiology (Deutsch et al. 2008 PNAS, Pinsky et al. 2019
        Nature),
      - microclimate availability (Burrows et al. 2019 NCC, Suggitt et
        al. 2018 NCC),
      - species mobility (Poloczanska et al. 2013 NCC, Burrows et
        al. 2011 Science, Sunday et al. 2012 NCC)
      - or generation time (Beaugrand et al. 2009 DSR II, Poloczanska et
        al. 2013 NCC),
      - consumers vs. producers (Petchey et al. 1999 Nature)
      - community composition (Stuart-Smith et al. 2015 Nature,
        Beaugrand et a. 2015 NCC, Trisos et al. 2020 Nature),
      - ecosystem productivity (Thomas et al. 2017 GCB, Brett 1971 Am
        Zoo),
      - exposure to human impacts (White & Kerr 2006 Ecography)
      - and among realms (Antao et al. 2020 NEE).
  - Scaling up from organismal effects to whole ecological communities
    is complex, and yet these scales are critical for ecosystem
    functioning and human well-being.
  - There is a need for a comprehensive test to understand where warming
    is driving and is likely to drive the most dramatic community
    turnover

# Methods

  - BioTime dataset, gridded to 96 km2 hexagons, summarized as temporal
    turnover (Blowes)
      - Temporal slope of Jaccard turnover compared to the first year
      - Same for Jaccard total
      - and Morisita-Horn turnover
  - Tested explanatory variables for differences in rate of turnover:
      - Temperature trend over the time-frame of each time-series (CRU
        TS 4.03 on land and in freshwater, ERSST v5 in the ocean)
      - Seasonality as a metric of thermal sensitivity (Deutsch et
        al. 2008 PNAS). Standard deviation of monthly temperatures.
      - Average temperature as a metric of metabolic rates (Dillon et
        al. 2010 Nature, Antao et al. 2020 Nat E\&E)
      - Mobility calculated from body mass and taxonomic group
        classifications of mobility mode (fly, run, swim, crawl,
        sessile). Fly/run/swim followed the allometric relationship in
        Hirt et al. 2017 Nat E\&E. Crawl set at 0.1 km/hr, sessile set
        to 0 km/hr. Then calculated averaged within each assemblage.
      - Net primary productivity (NPP) from the merged land/ocean
        product produced by the [Ocean
        Productivity](http://www.science.oregonstate.edu/ocean.productivity/)
        group at Oregon State using methods from Zhao et al. 2005 and
        Behrenfeld & Falkowski 1997.
  - TO DO:
      - Generation time calculated from body mass and endotherm
        vs. ectotherm classifications, following McCoy & Gillooly 2008
        ELE. Averaged across species within each assemblage.
      - Human impact calculated from Carsten Meyer’s ecosystem cube data
      - Thermal bias calculated from Species Temperature Indices (Mike
        Burrows)
      - Microclimates calculated from WorldClim and BioOracle (Laura
        Antao)
      - Consumer vs. producer classification
  - Maybe to do:
      - Species pool richness
  - Differences in temporal turnover (response variable) modeled with a
    linear mixed effects model (nlme package, lme() function). See below
    for details.

<!-- end list -->

``` r
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
#library(lme4)
library(nlme) # for ME models
library(beanplot) # for beanplots
library(maps) # for map
library(ggeffects) # marginal effect plots
```

    ## Warning: package 'ggeffects' was built under R version 3.6.2

``` r
library(gridExtra) # to combine ggplots together
library(grid) # to combine ggplots together

# tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) 
```

``` r
# Turnover and covariates assembled by turnover_vs_temperature_prep.Rmd
trends <- fread('output/turnover_w_covariates.csv.gz')

# set realm order
trends[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = TRUE)]
```

Log-transform some variables, then center and scale. Do they look ok?

``` r
trends[, temptrend_comb_allyr.sc := scale(temptrend_comb_allyr)]
trends[, temptrend_comb_allyr_abs.sc := scale(log(abs(temptrend_comb_allyr)))]
trends[, seas_comb.sc := scale(seas_comb)]
trends[, tempave_comb.sc := scale(tempave_comb)]
trends[, npp.sc := scale(npp)]
trends[, mass_geomean.sc := scale(log(mass_geomean))]
trends[, speed_geomean.sc := scale(log(speed_geomean+1))]

# histograms to examine
par(mfrow = c(3,3))
invisible(trends[, hist(temptrend_comb_allyr.sc, main = 'Temperature trend (°C/yr)')])
invisible(trends[, hist(temptrend_comb_allyr_abs.sc, main = 'log abs(Temperature trend) (°C/yr)')])
invisible(trends[, hist(seas_comb.sc, main = 'Seasonality (°C)')])
invisible(trends[, hist(tempave_comb.sc, main = 'Temperature (°C)')])
invisible(trends[, hist(npp.sc, main = 'Net primary productivity')])
invisible(trends[, hist(mass_geomean.sc, main = 'log Mass (g)')])
invisible(trends[, hist(speed_geomean.sc, main = 'log Speed (km/hr)')])
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-1.png)<!-- -->

# Results

## Where do we have data?

``` r
world <- map_data('world')
ggplot(world, aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = 'lightgray', color = 'white') +
    geom_point(data = trends, aes(rarefyID_x, rarefyID_y, group = REALM, color = REALM), size = 0.5, alpha = 0.4)  +
    scale_color_brewer(palette="Set1")
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/map-1.png)<!-- -->
Mostly northern hemisphere, but spread all over. No so much in Africa or
much of Asia.

## Plot turnover vs. temperature trends

Lines are ggplot smoother fits by
    realm.

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-1.png)<!-- -->

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-2.png)<!-- -->

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 4262 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 4262 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-3.png)<!-- -->
Trends are pretty symmetric around no trend in temperature, which
implies warming or cooling drives similar degree of community turnover.

## Compare covariates across realms

``` r
i <- trends[, !duplicated(STUDY_ID)]; sum(i)
```

    ## [1] 332

``` r
par(mfrow=c(2,3))
beanplot(rarefyID_y ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Latitude (degN)', ll = 0.05)
beanplot(tempave_comb ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature (degC)', ll = 0.05)
beanplot(seas_comb ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Seasonality (degC)', ll = 0.05)
beanplot(npp ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'NPP', ll = 0.05)
beanplot(mass_geomean ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Mass (g)', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(speed_geomean ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Speed (km/hr)', ll = 0.05, log = '')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->
Marine are in generally warmer locations (seawater doesn’t freeze)
Marine have much lower seasonality. Marine and freshwater have some very
small masses (plankton), but much of dataset is similar to terrestrial.
Marine has a lot of slow, crawling organisms, but land has plants. Land
also has birds (fast).

## Model choice: Jaccard turnover temporal trend vs. temperature trend

First examine how many data points are available

``` r
# the cases we can compare
cat('Number of data points available:\n')
```

    ## Number of data points available:

``` r
apply(trends[, .(Jtutrend, temptrend_comb_allyr, tempave_comb, REALM, seas_comb, npp, mass_geomean, speed_geomean)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##             Jtutrend temptrend_comb_allyr         tempave_comb 
    ##                53467                50335                50335 
    ##                REALM            seas_comb                  npp 
    ##                53467                50335                53314 
    ##         mass_geomean        speed_geomean 
    ##                53078                52649

``` r
i <- trends[, complete.cases(Jtutrend, temptrend_comb_allyr, REALM, seas_comb, npp, tempave_comb, mass_geomean, speed_geomean)]
cat('Overall:\n')
```

    ## Overall:

``` r
sum(i)
```

    ## [1] 49471

### Choose the variance structure

Try combinations of

  - variance scaled to a power of the number of years in the community
    time-series
  - variance scaled to a power of the abs temperature trend
  - random intercept for STUDY\_ID
  - random slope (abs temperature trend) for STUDY\_ID
  - random intercept for rarefyID (for overdispersion)

And choose the one with lowest AIC

``` r
# fit models for variance structure
fixed <- formula(Jtutrend ~ temptrend_comb_allyr_abs.sc*REALM + 
                     temptrend_comb_allyr_abs.sc*tempave_comb + 
                     temptrend_comb_allyr_abs.sc*seas_comb + 
                     temptrend_comb_allyr_abs.sc*npp + 
                     temptrend_comb_allyr_abs.sc*mass_geomean + 
                     temptrend_comb_allyr_abs.sc*speed_geomean)
i <- trends[, complete.cases(Jtutrend, temptrend_comb_allyr, REALM, seas_comb, npp, tempave_comb, mass_geomean, speed_geomean)]
mods <- vector('list', 0)
mods[[1]] <- gls(fixed, data = trends[i,])
mods[[2]] <- gls(fixed, data = trends[i,], weights = varPower(-0.5, ~nyrBT))
mods[[3]] <- gls(fixed, data = trends[i,], weights = varPower(0.5, ~ abs(temptrend_comb_allyr)))

mods[[4]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, control = lmeControl(opt = "optim"))
mods[[5]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, control = lmeControl(opt = "optim"))

mods[[6]] <- lme(fixed, data = trends[i,], random = ~temptrend_comb_allyr_abs.sc | STUDY_ID)
mods[[7]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_comb_allyr_abs.sc, rarefyID = ~1)) # includes overdispersion. new formula so that random slope is only for study level (not enough data to extend to rarefyID).

mods[[8]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, weights = varPower(-0.5, ~nyrBT))
mods[[9]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, weights = varPower(-0.5, ~nyrBT))
mods[[10]] <- lme(fixed, data = trends[i,], random = ~temptrend_comb_allyr_abs.sc|STUDY_ID, weights = varPower(-0.5, ~nyrBT))
mods[[11]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_comb_allyr_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT))

mods[[12]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, weights = varPower(-0.5, ~abs(temptrend_comb_allyr)))
mods[[13]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, weights = varPower(-0.5, ~abs(temptrend_comb_allyr)))
mods[[14]] <- lme(fixed, data = trends[i,], random = ~temptrend_comb_allyr_abs.sc|STUDY_ID, weights = varPower(-0.5, ~abs(temptrend_comb_allyr)))
mods[[15]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_comb_allyr_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~abs(temptrend_comb_allyr)))

aics <- sapply(mods, AIC)
minaics <- aics - min(aics)
minaics
```

    ##  [1] 48445.19140 18222.58517 41627.12157 38675.44804 38677.44804
    ##  [6] 37786.50724 37788.50724  1765.08077  1588.13774    20.25581
    ## [11]     0.00000 32399.98390 28678.50831 31558.81180 27907.09088

Chooses the random slopes & intercepts, overdispersion, and variance
scaled to number of years. We haven’t dealt with potential testing on
the boundary issues here yet.

### Compare fixed effects among models

Try interactions of abs temperature trend with each covariate:

  - realm
  - average
temperature
  - seasonality
  - speed
  - mass
  - NPP

<!-- end list -->

``` r
i <- trends[, complete.cases(Jtutrend, temptrend_comb_allyr, REALM, seas_comb, npp, tempave_comb, mass_geomean, speed_geomean)]

randef <- list(STUDY_ID = ~ temptrend_comb_allyr_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

mods <- vector('list', 0)
mods[[1]] <- lme(Jtutrend ~ 1, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[2]] <- lme(Jtutrend ~ temptrend_comb_allyr.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[3]] <- lme(Jtutrend ~ temptrend_comb_allyr_abs.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[4]] <- lme(Jtutrend ~ temptrend_comb_allyr_abs.sc * REALM, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[5]] <- lme(Jtutrend ~ temptrend_comb_allyr_abs.sc * tempave_comb.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[6]] <- lme(Jtutrend ~ temptrend_comb_allyr_abs.sc * seas_comb.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[7]] <- lme(Jtutrend ~ temptrend_comb_allyr_abs.sc * speed_geomean.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[8]] <- lme(Jtutrend ~ temptrend_comb_allyr_abs.sc * mass_geomean.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[9]] <- lme(Jtutrend ~ temptrend_comb_allyr_abs.sc * npp.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')

mods[[10]] <- lme(Jtutrend ~ temptrend_comb_allyr_abs.sc + speed_geomean.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')

mods[[11]] <- lme(Jtutrend ~ speed_geomean.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')

# examine models
cat('AIC:\n')
```

    ## AIC:

``` r
aics <- sapply(mods, AIC) - min(sapply(mods, AIC)); aics
```

    ##  [1] 137.2464752 139.0933902 102.0253274  82.1151257  88.3846613
    ##  [6] 100.6581440   0.1344852  57.3818067  96.4063838   0.0000000
    ## [11]  34.9180486

``` r
cat(paste0('\nChose model ', which.min(aics), ': ', paste0(mods[[which.min(aics)]]$call$fixed[c(2,1,3)], collapse = ' '), '\n'))
```

    ## 
    ## Chose model 10: Jtutrend `~` temptrend_comb_allyr_abs.sc + speed_geomean.sc

``` r
cat('\nModel terms:\n')
```

    ## 
    ## Model terms:

``` r
summary(mods[[which.min(aics)]])
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: trends[i, ] 
    ##         AIC       BIC  logLik
    ##   -153128.4 -153049.1 76573.2
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_comb_allyr_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                             StdDev     Corr  
    ## (Intercept)                 0.05418508 (Intr)
    ## temptrend_comb_allyr_abs.sc 0.01890756 0.457 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.002242457 0.3161294
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.255242 
    ## Fixed effects: Jtutrend ~ temptrend_comb_allyr_abs.sc + speed_geomean.sc 
    ##                                  Value   Std.Error    DF   t-value p-value
    ## (Intercept)                 0.04468975 0.004104368 49262 10.888340       0
    ## temptrend_comb_allyr_abs.sc 0.01378307 0.002073437 49262  6.647447       0
    ## speed_geomean.sc            0.00677595 0.000662523 49262 10.227497       0
    ##  Correlation: 
    ##                             (Intr) tm___.
    ## temptrend_comb_allyr_abs.sc  0.302       
    ## speed_geomean.sc            -0.033  0.003
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -7.8881796 -0.2609638  0.1015781  0.5795188  6.7819934 
    ## 
    ## Number of Observations: 49471
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    207                  49471

Chooses a model with temperature trend and speed, but maybe not the
interaction

## Examine the chosen model

``` r
bmod <- update(mods[[7]], method = 'REML') # re-fit with REML
summary(bmod)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC     BIC   logLik
    ##   -153083.1 -152995 76551.57
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_comb_allyr_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                             StdDev     Corr  
    ## (Intercept)                 0.05446662 (Intr)
    ## temptrend_comb_allyr_abs.sc 0.01955632 0.463 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.002229942 0.3160733
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.255114 
    ## Fixed effects: Jtutrend ~ temptrend_comb_allyr_abs.sc * speed_geomean.sc 
    ##                                                   Value   Std.Error    DF
    ## (Intercept)                                  0.04466173 0.004126729 49261
    ## temptrend_comb_allyr_abs.sc                  0.01370534 0.002129866 49261
    ## speed_geomean.sc                             0.00708839 0.000698578 49261
    ## temptrend_comb_allyr_abs.sc:speed_geomean.sc 0.00102555 0.000734488 49261
    ##                                                t-value p-value
    ## (Intercept)                                  10.822551  0.0000
    ## temptrend_comb_allyr_abs.sc                   6.434837  0.0000
    ## speed_geomean.sc                             10.146883  0.0000
    ## temptrend_comb_allyr_abs.sc:speed_geomean.sc  1.396279  0.1626
    ##  Correlation: 
    ##                                              (Intr) tm___. spd_g.
    ## temptrend_comb_allyr_abs.sc                   0.308              
    ## speed_geomean.sc                             -0.032 -0.011       
    ## temptrend_comb_allyr_abs.sc:speed_geomean.sc -0.003 -0.043  0.317
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -7.8791935 -0.2608093  0.1014711  0.5797334  6.7944542 
    ## 
    ## Number of Observations: 49471
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    207                  49471

``` r
hist(residuals(bmod))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-1.png)<!-- -->

``` r
qqnorm(bmod, ~ resid(., type = 'p'), abline = c(0,1))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-2.png)<!-- -->

``` r
qqnorm(bmod, ~ ranef(., level = 1))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-3.png)<!-- -->

``` r
qqnorm(bmod, ~ ranef(., level = 2))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-4.png)<!-- -->

``` r
plot(fitted(bmod), resid(bmod, type = 'normalized'), col = '#00000033')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-5.png)<!-- -->

``` r
plot(trends[i,nyrBT], residuals(bmod, type = 'normalized'), col = '#00000011', log = 'x')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-6.png)<!-- -->

``` r
plot(trends[i,temptrend_comb_allyr_abs.sc], residuals(bmod, type = 'normalized'), col = '#00000011')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-7.png)<!-- -->

``` r
plot(trends[i,speed_geomean.sc], residuals(bmod, type = 'normalized'), col = '#00000033')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-8.png)<!-- -->

Residuals may still be overdispersed despite the observation-level RE?
Not sure what else to try. Still some heterogeneity (cone) for nyrBT and
temperature trend, but may be created by so few data points in certain
locations on the
axis

## Plot the chosen model

``` r
ggplot(data = trends[i,], aes(x = speed_geomean.sc, y = Jtutrend, color = REALM)) +
    geom_point(shape = 16, alpha = 0.5) +
    geom_smooth(method = 'lm') +
    scale_color_brewer(palette="Set1")
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/examine%20model-1.png)<!-- -->

``` r
# fix temptrend at low or high
newdat1 <- trends[i, .(STUDY_ID, rarefyID, REALM, temptrend_comb_allyr_abs.sc = quantile(temptrend_comb_allyr_abs.sc, 0.25), speed_geomean, speed_geomean.sc, grp = 1)]
newdat2 <- trends[i, .(STUDY_ID, rarefyID, REALM, temptrend_comb_allyr_abs.sc = quantile(temptrend_comb_allyr_abs.sc, 0.75), speed_geomean, speed_geomean.sc, grp = 2)]
newdat <- rbind(newdat1, newdat2)
newdat$preds <- predict(bmod, newdata = newdat, level = 0)

ggplot(newdat, aes(speed_geomean.sc, preds, color = grp, group = grp)) +
    geom_line()
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/examine%20model-2.png)<!-- -->

``` r
ggplot(newdat, aes(speed_geomean, preds, color = grp, group = grp)) +
    geom_line()
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/examine%20model-3.png)<!-- -->

``` r
# fix speed at low or high
newdat1 <- trends[i, .(STUDY_ID, rarefyID, REALM, temptrend_comb_allyr, temptrend_comb_allyr_abs.sc, speed_geomean.sc = quantile(speed_geomean.sc, 0.25), grp = 1)]
newdat2 <- trends[i, .(STUDY_ID, rarefyID, REALM, temptrend_comb_allyr, temptrend_comb_allyr_abs.sc, speed_geomean.sc = quantile(speed_geomean.sc, 0.75), grp = 2)]
newdat <- rbind(newdat1, newdat2)
newdat$preds <- predict(bmod, newdata = newdat, level = 0)

ggplot(newdat, aes(temptrend_comb_allyr_abs.sc, preds, color = grp, group = grp)) +
    geom_line()
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/examine%20model-4.png)<!-- -->

``` r
ggplot(newdat, aes(temptrend_comb_allyr, preds, color = grp, group = grp)) +
    geom_line()
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/examine%20model-5.png)<!-- -->

# To do

  - plot residuals on top of temperature relationship
  - plot turnover vs. speed, look for outliers
  - other temporal turnover metrics (Jaccard total, M-H)

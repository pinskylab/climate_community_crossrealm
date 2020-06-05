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
trends[, tempave.sc := scale(tempave)]
trends[, tempave_metab.sc := scale(tempave_metab)]
trends[, seas.sc := scale(seas)]
trends[, microclim.sc := scale(log(microclim))]
trends[, temptrend.sc := scale(temptrend)]
trends[, temptrend_abs.sc := scale(log(abs(temptrend)))]
trends[, npp.sc := scale(log(npp))]
trends[, mass.sc := scale(log(mass_geomean))]
trends[, speed.sc := scale(log(speed_geomean+1))]
trends[, lifespan.sc := scale(log(lifespan_geomean))]
trends[, thermal_bias.sc := scale(thermal_bias)]
trends[, consumerfrac.sc := scale(consfrac)]

# histograms to examine
par(mfrow = c(3,4))
invisible(trends[, hist(tempave.sc, main = 'Environmental temperature (°C)')])
invisible(trends[, hist(tempave_metab.sc, main = 'Metabolic temperature (°C)')])
invisible(trends[, hist(seas.sc, main = 'Seasonality (°C)')])
invisible(trends[, hist(microclim.sc, main = 'log Microclimates (°C)')])
invisible(trends[, hist(temptrend.sc, main = 'Temperature trend (°C/yr)')])
invisible(trends[, hist(temptrend_abs.sc, main = 'log abs(Temperature trend) (°C/yr)')])
invisible(trends[, hist(npp.sc, main = 'log Net primary productivity')])
invisible(trends[, hist(mass.sc, main = 'log Mass (g)')])
invisible(trends[, hist(speed.sc, main = 'log Speed (km/hr)')])
invisible(trends[, hist(lifespan.sc, main = 'log Lifespan (yr)')])
invisible(trends[, hist(thermal_bias.sc, main = 'Thermal bias (°C)')])
invisible(trends[, hist(consumerfrac.sc, main = 'Consumers (fraction)')])
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-1.png)<!-- -->

First examine how many data points are available

``` r
# the cases we can compare
cat('Number of data points available:\n')
```

    ## Number of data points available:

``` r
apply(trends[, .(Jtutrend, temptrend.sc, tempave_metab.sc, REALM, seas.sc, microclim.sc, npp.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, thermal_bias.sc)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##         Jtutrend     temptrend.sc tempave_metab.sc            REALM 
    ##            53467            50335            50335            53467 
    ##          seas.sc     microclim.sc           npp.sc          mass.sc 
    ##            50335            52262            53314            53078 
    ##         speed.sc      lifespan.sc  consumerfrac.sc  thermal_bias.sc 
    ##            52649            51797            47988            49624

``` r
i <- trends[, complete.cases(Jtutrend, temptrend.sc, tempave_metab.sc, REALM, seas.sc, microclim.sc, npp.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, thermal_bias.sc)]
cat('Overall:\n')
```

    ## Overall:

``` r
sum(i)
```

    ## [1] 43559

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
fixed <- formula(Jtutrend ~ REALM + tempave_metab.sc + seas.sc + microclim.sc + npp.sc + temptrend.sc +
                     mass.sc + speed.sc + lifespan.sc + consumerfrac.sc + thermal_bias.sc)
i <- trends[, complete.cases(Jtutrend, REALM, tempave_metab.sc, seas.sc, microclim.sc, npp.sc, temptrend.sc,
                             mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, thermal_bias.sc)]
mods <- vector('list', 0)
mods[[1]] <- gls(fixed, data = trends[i,])
mods[[2]] <- gls(fixed, data = trends[i,], weights = varPower(-0.5, ~nyrBT))
mods[[3]] <- gls(fixed, data = trends[i,], weights = varPower(0.5, ~ abs(temptrend)))

mods[[4]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, control = lmeControl(opt = "optim"))
mods[[5]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, control = lmeControl(opt = "optim"))

mods[[6]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc | STUDY_ID)
mods[[7]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)) # includes overdispersion. new formula so that random slope is only for study level (not enough data to extend to rarefyID).

mods[[8]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, weights = varPower(-0.5, ~nyrBT))
mods[[9]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, weights = varPower(-0.5, ~nyrBT))
mods[[10]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc|STUDY_ID, weights = varPower(-0.5, ~nyrBT))
mods[[11]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT))

mods[[12]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, weights = varPower(-0.5, ~abs(temptrend)))
mods[[13]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, weights = varPower(-0.5, ~abs(temptrend)))
mods[[14]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc|STUDY_ID, weights = varPower(-0.5, ~abs(temptrend)))
mods[[15]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~abs(temptrend)))

aics <- sapply(mods, AIC)
minaics <- aics - min(aics)
minaics
```

    ##  [1] 40894.193939 14779.149044 33681.552359 35003.152188 35005.152188
    ##  [6] 29874.125281 29876.125282  3323.064861  3307.416479     4.542291
    ## [11]     0.000000 27707.938936 24132.563076 25160.466709 21909.263091

Chooses the random slopes & intercepts, overdispersion, and variance
scaled to number of years. We haven’t dealt with potential testing on
the boundary issues here yet.

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
beanplot(tempave ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature (degC)', ll = 0.05)
beanplot(tempave_metab ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Metabolic Temperature (degC)', ll = 0.05)
beanplot(seas ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Seasonality (degC)', ll = 0.05)
beanplot(microclim ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Microclimates (degC)', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(temptrend ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature trend (degC/yr)', ll = 0.05)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->

``` r
beanplot(npp ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'NPP', ll = 0.05)
beanplot(mass_geomean ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Mass (g)', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(speed_geomean ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Speed (km/hr)', ll = 0.05, log = '')
beanplot(lifespan_geomean ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Lifespan (yr)', ll = 0.05, log = '')
beanplot(thermal_bias ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Thermal bias (degC)', ll = 0.05, log = '')
#beanplot(consfrac ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Consumers (fraction)', ll = 0.05, log = '') # too sparse
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-2.png)<!-- -->
Marine are in generally warmer locations (seawater doesn’t freeze)
Marine have much lower seasonality. Marine and freshwater have some very
small masses (plankton), but much of dataset is similar to terrestrial.
Marine has a lot of slow, crawling organisms, but land has plants. Land
also has birds (fast).

## Model choice: Jaccard turnover temporal trend vs. covariates

  - realm
  - average metabolic temperature
  - seasonality
  - microclimates
  - NPP
  - speed
  - mass
  - lifespan
  - consumer vs. producer
  - thermal
bias

<!-- end list -->

``` r
i <- trends[, complete.cases(Jtutrend, REALM, tempave_metab.sc, seas.sc, microclim.sc, npp.sc,
                             mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, thermal_bias.sc)]

randef <- list(STUDY_ID = ~ 1, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

mods <- vector('list', 0)
mods[[1]] <- lme(Jtutrend ~ 1, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[2]] <- lme(Jtutrend ~ REALM, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[3]] <- lme(Jtutrend ~ tempave_metab.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[4]] <- lme(Jtutrend ~ seas.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[5]] <- lme(Jtutrend ~ microclim.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[6]] <- lme(Jtutrend ~ npp.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[7]] <- lme(Jtutrend ~ speed.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[8]] <- lme(Jtutrend ~ mass.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[9]] <- lme(Jtutrend ~ lifespan.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[10]] <- lme(Jtutrend ~ consumerfrac.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods[[11]] <- lme(Jtutrend ~ thermal_bias.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')

# examine models
cat('AIC:\n')
```

    ## AIC:

``` r
aics <- sapply(mods, AIC) - min(sapply(mods, AIC)); aics
```

    ##  [1] 70.044979 54.468409 69.975536 71.059149 37.894461  5.195981 19.339194
    ##  [8]  0.000000 63.511457 64.818691 69.745851

``` r
cat(paste0('\nChose model ', which.min(aics), ': ', paste0(mods[[which.min(aics)]]$call$fixed[c(2,1,3)], collapse = ' '), '\n'))
```

    ## 
    ## Chose model 8: Jtutrend `~` mass.sc

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
    ##         AIC       BIC   logLik
    ##   -125559.6 -125507.5 62785.81
    ## 
    ## Random effects:
    ##  Formula: ~1 | STUDY_ID
    ##         (Intercept)
    ## StdDev:  0.05912184
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev: 0.001206208 0.345842
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.299339 
    ## Fixed effects: Jtutrend ~ mass.sc 
    ##                   Value   Std.Error    DF   t-value p-value
    ## (Intercept)  0.04062228 0.004742349 43375  8.565855       0
    ## mass.sc     -0.00697705 0.000819578 43375 -8.512981       0
    ##  Correlation: 
    ##         (Intr)
    ## mass.sc 0.189 
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -9.33802335 -0.24004997  0.09812242  0.60777708  6.44613610 
    ## 
    ## Number of Observations: 43559
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    183                  43559

Chooses a model with mass (slower turnover for larger body size). Next
best is NPP (deltaAIC 5)

## Model choice: Jaccard turnover temporal trend vs. temperature trend

Try interactions of abs temperature trend with each covariate:

  - realm
  - average metabolic temperature
  - seasonality
  - microclimates
  - NPP
  - speed
  - mass
  - lifespan
  - consumer vs. producer
  - thermal
bias

<!-- end list -->

``` r
i <- trends[, complete.cases(Jtutrend, temptrend, REALM, tempave_metab.sc, seas.sc, microclim.sc, npp.sc, 
                             mass.sc, speed.sc, lifespan.sc, thermal_bias.sc, consumerfrac.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

mods2 <- vector('list', 0)
mods2[[1]] <- lme(Jtutrend ~ 1, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[2]] <- lme(Jtutrend ~ temptrend.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[3]] <- lme(Jtutrend ~ temptrend_abs.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[4]] <- lme(Jtutrend ~ temptrend_abs.sc * REALM, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[5]] <- lme(Jtutrend ~ temptrend_abs.sc * tempave_metab.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[6]] <- lme(Jtutrend ~ temptrend_abs.sc * seas.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[7]] <- lme(Jtutrend ~ temptrend_abs.sc * microclim.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[8]] <- lme(Jtutrend ~ temptrend_abs.sc * npp.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[9]] <- lme(Jtutrend ~ temptrend_abs.sc * speed.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[10]] <- lme(Jtutrend ~ temptrend_abs.sc * mass.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[11]] <- lme(Jtutrend ~ temptrend_abs.sc * lifespan.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[12]] <- lme(Jtutrend ~ temptrend_abs.sc * consumerfrac.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')
mods2[[13]] <- lme(Jtutrend ~ temptrend_abs.sc * thermal_bias.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')

mods2[[14]] <- lme(Jtutrend ~ temptrend_abs.sc + speed.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')

mods2[[15]] <- lme(Jtutrend ~ speed.sc, random = randef, 
                 weights = varef, data = trends[i,], method = 'ML')

# examine models
cat('AIC:\n')
```

    ## AIC:

``` r
aics <- sapply(mods2, AIC) - min(sapply(mods2, AIC)); aics
```

    ##  [1] 122.155932 123.821336  88.705187  69.506448  61.336645  88.617908
    ##  [7]   6.801944  53.238895   0.000000  28.358985  13.821183  84.418094
    ## [13]  61.351105  21.002366  53.452021

``` r
cat(paste0('\nChose model ', which.min(aics), ': ', paste0(mods2[[which.min(aics)]]$call$fixed[c(2,1,3)], collapse = ' '), '\n'))
```

    ## 
    ## Chose model 9: Jtutrend `~` temptrend_abs.sc * speed.sc

``` r
cat('\nModel terms:\n')
```

    ## 
    ## Model terms:

``` r
summary(mods2[[which.min(aics)]])
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -129026.7 -128939.9 64523.36
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.05597433 (Intr)
    ## temptrend_abs.sc 0.01853060 0.355 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 0.0009920627 0.3120502
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.255396 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc * speed.sc 
    ##                                 Value   Std.Error    DF   t-value p-value
    ## (Intercept)                0.04550969 0.004491278 43373 10.132904       0
    ## temptrend_abs.sc           0.01473795 0.002170159 43373  6.791186       0
    ## speed.sc                   0.00615244 0.000789113 43373  7.796656       0
    ## temptrend_abs.sc:speed.sc -0.00383958 0.000793227 43373 -4.840456       0
    ##  Correlation: 
    ##                           (Intr) tmpt_. spd.sc
    ## temptrend_abs.sc           0.232              
    ## speed.sc                  -0.063 -0.007       
    ## temptrend_abs.sc:speed.sc  0.003 -0.061  0.110
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -7.56673953 -0.27281439  0.08978463  0.57986869  6.83169173 
    ## 
    ## Number of Observations: 43559
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    183                  43559

Chooses a model with temperature trend and speed. More turnover at
faster speeds and faster temperature change, but less response to
temperature change at faster speeds.

Examine the chosen model

``` r
bmod <- update(mods2[[9]], method = 'REML') # re-fit with REML
summary(bmod)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -128982.4 -128895.6 64501.19
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.05621517 (Intr)
    ## temptrend_abs.sc 0.01880514 0.348 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 0.0009923436 0.3120258
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.255337 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc * speed.sc 
    ##                                 Value   Std.Error    DF   t-value p-value
    ## (Intercept)                0.04553018 0.004510611 43373 10.094016       0
    ## temptrend_abs.sc           0.01479151 0.002198582 43373  6.727751       0
    ## speed.sc                   0.00615925 0.000789288 43373  7.803553       0
    ## temptrend_abs.sc:speed.sc -0.00385260 0.000795136 43373 -4.845213       0
    ##  Correlation: 
    ##                           (Intr) tmpt_. spd.sc
    ## temptrend_abs.sc           0.228              
    ## speed.sc                  -0.063 -0.007       
    ## temptrend_abs.sc:speed.sc  0.003 -0.060  0.110
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -7.56648491 -0.27280515  0.08980588  0.57983358  6.83172649 
    ## 
    ## Number of Observations: 43559
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    183                  43559

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
plot(trends[i,temptrend_abs.sc], residuals(bmod, type = 'normalized'), col = '#00000011')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-7.png)<!-- -->

``` r
plot(trends[i,speed.sc], residuals(bmod, type = 'normalized'), col = '#00000033')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/model%20examination-8.png)<!-- -->

Residuals may still be overdispersed despite the observation-level RE?
Not sure what else to try. Still some heterogeneity (cone) for nyrBT and
temperature trend, but may be created by so few data points in certain
locations on the axis

Plot the chosen
model

``` r
ggplot(data = trends[i,], aes(x = speed.sc, y = Jtutrend, color = REALM)) +
    geom_point(shape = 16, alpha = 0.5) +
    geom_smooth(method = 'lm') +
    scale_color_brewer(palette="Set1")
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/examine%20model-1.png)<!-- -->

``` r
# fix temptrend at low or high
newdat1 <- trends[i, .(STUDY_ID, rarefyID, REALM, temptrend_abs.sc = quantile(temptrend_abs.sc, 0.25), speed_geomean, speed.sc, grp = 1)]
newdat2 <- trends[i, .(STUDY_ID, rarefyID, REALM, temptrend_abs.sc = quantile(temptrend_abs.sc, 0.75), speed_geomean, speed.sc, grp = 2)]
newdat <- rbind(newdat1, newdat2)
newdat$preds <- predict(bmod, newdata = newdat, level = 0)

ggplot(newdat, aes(speed.sc, preds, color = grp, group = grp)) +
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
newdat1 <- trends[i, .(STUDY_ID, rarefyID, REALM, temptrend, temptrend_abs.sc, speed.sc = quantile(speed.sc, 0.25), grp = 1)]
newdat2 <- trends[i, .(STUDY_ID, rarefyID, REALM, temptrend, temptrend_abs.sc, speed.sc = quantile(speed.sc, 0.75), grp = 2)]
newdat <- rbind(newdat1, newdat2)
newdat$preds <- predict(bmod, newdata = newdat, level = 0)

ggplot(newdat, aes(temptrend_abs.sc, preds, color = grp, group = grp)) +
    geom_line()
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/examine%20model-4.png)<!-- -->

``` r
ggplot(newdat, aes(temptrend, preds, color = grp, group = grp)) +
    geom_line()
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/examine%20model-5.png)<!-- -->

# To do

  - plot residuals on top of temperature relationship
  - plot turnover vs. speed, look for outliers
  - other temporal turnover metrics (Jaccard total, M-H)

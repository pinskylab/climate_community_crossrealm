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
library(nlme) # for ME models
library(maps) # for map
library(gridExtra) # to combine ggplots together
library(grid) # to combine ggplots together

options(width=500) # turn off most text wrapping

# tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) 
```

``` r
# Turnover and covariates assembled by turnover_vs_temperature_prep.Rmd
trends <- fread('output/turnover_w_covariates.csv.gz')

# set realm order
trends[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = FALSE)]

# set up sign of temperature change
trends[, tsign := factor(sign(temptrend))]

# realm that combined Terrestrial and Freshwater, for interacting with human impact
trends[, REALM2 := REALM]
levels(trends$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")

# group Marine invertebrates/plants in with All
trends[, taxa_mod2 := taxa_mod]
trends[taxa_mod == 'Marine invertebrates/plants', taxa_mod2 := 'All']
```

### Log-transform some variables, then center and scale.

``` r
trends[, tempave.sc := scale(tempave)]
trends[, tempave_metab.sc := scale(tempave_metab)]
trends[, seas.sc := scale(seas)]
trends[, microclim.sc := scale(log(microclim))]
trends[, temptrend.sc := scale(temptrend, center = FALSE)]
trends[, temptrend_abs.sc := scale(abs(temptrend), center = FALSE)] # do not center, so that 0 is still 0 temperature change
trends[, mass.sc := scale(log(mass_mean_weight))]
trends[, speed.sc := scale(log(speed_mean_weight+1))]
trends[, lifespan.sc := scale(log(lifespan_mean_weight))]
trends[, consumerfrac.sc := scale(consfrac)]
trends[, endothermfrac.sc := scale(endofrac)]
trends[, nspp.sc := scale(log(Nspp))]
trends[, thermal_bias.sc := scale(thermal_bias)]
trends[, npp.sc := scale(log(npp))]
trends[, veg.sc := scale(log(veg+1))]
trends[, human_bowler.sc := scale(log(human_bowler+1)), by = REALM2] # separate scaling by realm
trends[REALM2 == 'TerrFresh', human_footprint.sc := scale(log(human_venter+1))]
trends[REALM2 == 'Marine', human_footprint.sc := scale(log(human_halpern))]
```

### Examine how many data points are available

Just turnover

``` r
cat('Overall # time-series: ', nrow(trends), '\n')
```

    ## Overall # time-series:  53013

``` r
cat('# studies: ', trends[, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  332

``` r
cat('Data points: ', trends[, sum(nyrBT)], '\n')
```

    ## Data points:  293973

``` r
trends[, table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##        1025       48647        3341

``` r
trends[, table(taxa_mod)]
```

    ## taxa_mod
    ##                         All                  Amphibians                     Benthos                       Birds                        Fish               Invertebrates                     Mammals Marine invertebrates/plants                       Plant                    Reptiles 
    ##                        1705                         379                        4679                       13741                       28473                        2996                         525                         206                         305                           4

``` r
trends[, table(taxa_mod, REALM)]
```

    ##                              REALM
    ## taxa_mod                      Freshwater Marine Terrestrial
    ##   All                                  0   1702           3
    ##   Amphibians                           2      0         377
    ##   Benthos                              0   4679           0
    ##   Birds                                0  11099        2642
    ##   Fish                              1006  27467           0
    ##   Invertebrates                       15   2901          80
    ##   Mammals                              0    478          47
    ##   Marine invertebrates/plants          0    206           0
    ##   Plant                                1    115         189
    ##   Reptiles                             1      0           3

With all covariates (Venter/Halpern for human)

``` r
# the cases we can compare
apply(trends[, .(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_footprint.sc)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##           Jtutrend              REALM         tempave.sc   tempave_metab.sc            seas.sc       microclim.sc       temptrend.sc            mass.sc           speed.sc        lifespan.sc    consumerfrac.sc   endothermfrac.sc            nspp.sc    thermal_bias.sc             npp.sc             veg.sc human_footprint.sc 
    ##              53013              53013              49916              49916              49916              51834              49916              52820              52734              51540              47534              53013              53013              49371              52863              52890              51197

``` r
i <- trends[, complete.cases(Jtutrend, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_footprint.sc)]
cat('Overall # time-series: ', sum(i), '\n')
```

    ## Overall # time-series:  43174

``` r
cat('# studies: ', trends[i, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  227

``` r
cat('Data points: ', trends[i, sum(nyrBT)], '\n')
```

    ## Data points:  218934

``` r
trends[i, table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##         967       39491        2716

``` r
trends[i, table(taxa_mod)]
```

    ## taxa_mod
    ##           All    Amphibians       Benthos         Birds          Fish Invertebrates       Mammals         Plant      Reptiles 
    ##           521            12           590         11575         27256          2525           514           180             1

``` r
trends[i, table(taxa_mod, REALM)]
```

    ##                REALM
    ## taxa_mod        Freshwater Marine Terrestrial
    ##   All                    0    520           1
    ##   Amphibians             2      0          10
    ##   Benthos                0    590           0
    ##   Birds                  0   9100        2475
    ##   Fish                 955  26301           0
    ##   Invertebrates          9   2450          66
    ##   Mammals                0    477          37
    ##   Plant                  1     53         126
    ##   Reptiles               0      0           1

### Choose the variance structure for mixed effects modles

Try combinations of

  - variance scaled to a power of the number of years in the community
    time-series
  - variance scaled to a power of the abs temperature trend
  - random intercept for taxa\_mod
  - random intercept for STUDY\_ID
  - random slope (abs temperature trend) for taxa\_mod
  - random slope (abs temperature trend) for STUDY\_ID
  - random intercept for rarefyID (for overdispersion)

And choose the one with lowest AIC (not run: takes a long time)

``` r
# fit models for variance structure
fixed <- formula(Jtutrend ~ REALM + tempave_metab.sc + seas.sc + microclim.sc + npp.sc + temptrend_abs.sc +
                     mass.sc + speed.sc + lifespan.sc + consumerfrac.sc + thermal_bias.sc)
i <- trends[, complete.cases(Jtutrend, REALM, tempave_metab.sc, seas.sc, microclim.sc, npp.sc, temptrend_abs.sc,
                             mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, thermal_bias.sc)]
mods <- vector('list', 0)
mods[[1]] <- gls(fixed, data = trends[i,])
mods[[2]] <- gls(fixed, data = trends[i,], weights = varPower(-0.5, ~nyrBT))
mods[[3]] <- gls(fixed, data = trends[i,], weights = varPower(0.5, ~ abs(temptrend)))

mods[[4]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2, control = lmeControl(opt = "optim"))
mods[[5]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, control = lmeControl(opt = "optim"))
mods[[6]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID, control = lmeControl(opt = "optim"))
mods[[7]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, control = lmeControl(opt = "optim"))
mods[[8]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID/rarefyID, control = lmeControl(opt = "optim"))

mods[[9]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc | taxa_mod)
mods[[10]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc | STUDY_ID)
mods[[11]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc | taxa_mod2/STUDY_ID, control = lmeControl(opt = "optim"))
mods[[12]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)) # includes overdispersion. new formula so that random slope is only for study level (not enough data to extend to rarefyID).
mods[[13]] <- lme(fixed, data = trends[i,], random = list(taxa_mod2 = ~ temptrend_abs.sc, STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)) # 30+ min to fit

mods[[14]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, weights = varPower(-0.5, ~nyrBT))
mods[[15]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2, weights = varPower(-0.5, ~nyrBT))
mods[[16]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID, weights = varPower(-0.5, ~nyrBT))
mods[[17]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, weights = varPower(-0.5, ~nyrBT))
mods[[18]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID/rarefyID, weights = varPower(-0.5, ~nyrBT))
mods[[19]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc|STUDY_ID, weights = varPower(-0.5, ~nyrBT))
mods[[20]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT))
mods[[21]] <- lme(fixed, data = trends[i,], random = list(taxa_mod2 = ~ temptrend_abs.sc, STUDY_ID = ~ 1), weights = varPower(-0.5, ~nyrBT))
mods[[22]] <- lme(fixed, data = trends[i,], random = list(taxa_mod2 = ~ temptrend_abs.sc, STUDY_ID = ~ 1, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT))
mods[[23]] <- lme(fixed, data = trends[i,], random = list(taxa_mod2 = ~ temptrend_abs.sc, STUDY_ID = ~ temptrend_abs.sc), weights = varPower(-0.5, ~nyrBT)) # singular precision warning with lmeControl(opt = 'optim') and convergence error without
mods[[24]] <- lme(fixed, data = trends[i,], random = list(taxa_mod2 = ~ temptrend_abs.sc, STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT)) # singular precision warning with lmeControl(opt = 'optim') and convergence error without

mods[[25]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2, weights = varPower(-0.5, ~abs(temptrend)))
mods[[26]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, weights = varPower(-0.5, ~abs(temptrend)))
mods[[27]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, weights = varPower(-0.5, ~abs(temptrend)))
mods[[28]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID/rarefyID, weights = varPower(-0.5, ~abs(temptrend)))
mods[[29]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc|STUDY_ID, weights = varPower(-0.5, ~abs(temptrend)))
mods[[30]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc|taxa_mod2/STUDY_ID, weights = varPower(-0.5, ~abs(temptrend)), control = lmeControl(opt = "optim"))
mods[[31]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~abs(temptrend)))
mods[[32]] <- lme(fixed, data = trends[i,], random = list(taxa_mod2 = ~ temptrend_abs.sc, STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~abs(temptrend)), control = lmeControl(opt = "optim")) # singular precision warning

aics <- sapply(mods, AIC)
minaics <- aics - min(aics)
minaics
which.min(aics)
```

Chooses the random slopes (temptrend\_abs) & intercepts for STUDY\_ID,
overdispersion, and variance scaled to number of years. We haven’t dealt
with potential testing on the boundary issues here yet.

# Results

## Where do we have data?

``` r
world <- map_data('world')
ggplot(world, aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = 'lightgray', color = 'white') +
    geom_point(data = trends, aes(rarefyID_x, rarefyID_y, group = REALM, color = REALM), size = 0.5, alpha = 0.4)  +
    scale_color_brewer(palette="Set1", name = 'Realm') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=20)) +
  labs(x = 'Longitude (°)', y = 'Latitude (°)')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/map-1.png)<!-- -->

Mostly northern hemisphere, but spread all over. No so much in Africa or
much of Asia.

Average rates of turnover

``` r
trends[abs(temptrend) >= 0.5, .(mean(Jtutrend), sd(Jtutrend)/sqrt(.N))] # turnover per year for locations changing temperature
```

    ##           V1          V2
    ## 1: 0.1814753 0.004483237

``` r
trends[abs(temptrend) < 0.1, .(mean(Jtutrend), sd(Jtutrend)/sqrt(.N))] # not changing temperature
```

    ##            V1           V2
    ## 1: 0.04592881 0.0003434666

``` r
trends[temptrend >= 0.5, .(mean(Jtutrend), sd(Jtutrend)/sqrt(.N))] # warming
```

    ##           V1          V2
    ## 1: 0.1616492 0.006620151

``` r
trends[temptrend <= -0.5, .(mean(Jtutrend), sd(Jtutrend)/sqrt(.N))] # cooling
```

    ##           V1          V2
    ## 1: 0.1961448 0.005982736

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) < 35, .(mean(Jtutrend), sd(Jtutrend)/sqrt(.N))] # tropics and sub-tropics
```

    ##           V1         V2
    ## 1: 0.4165367 0.04623876

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) >= 35 & abs(rarefyID_y) < 66.56339, .(mean(Jtutrend), sd(Jtutrend)/sqrt(.N))] # temperate
```

    ##          V1          V2
    ## 1: 0.178024 0.004397925

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) >= 66.56339, .(mean(Jtutrend), sd(Jtutrend)/sqrt(.N))] # arctic
```

    ##           V1         V2
    ## 1: 0.1962083 0.05650476

## Temperature-only model (Jtutrend, Jbetatrend, Horntrend)

``` r
i <- trends[, complete.cases(Jtutrend, REALM, temptrend)]

randef <- list(STUDY_ID = ~ abs(temptrend), rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modonlyTtrend.rds')){
  modonlyTtrend <- readRDS('temp/modonlyTtrend.rds')
} else {
  modonlyTtrend <- lme(Jtutrend ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrend, file = 'temp/modonlyTtrend.rds')
}

i2 <- trends[, complete.cases(Jbetatrend, REALM, temptrend)]
if(file.exists('temp/modonlyTtrendJbeta.rds')){
  modonlyTtrendJbeta <- readRDS('temp/modonlyTtrendJbeta.rds')
} else {
  modonlyTtrendJbeta <- lme(Jbetatrend ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i2,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbeta, file = 'temp/modonlyTtrendJbeta.rds')
}

i3 <- trends[, complete.cases(Horntrend, REALM, temptrend)]
if(file.exists('temp/modonlyTtrendHorn.rds')){
  modonlyTtrendHorn <- readRDS('temp/modonlyTtrendHorn.rds')
} else {
  modonlyTtrendHorn <- lme(Horntrend ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i3,], method = 'REML')
  saveRDS(modonlyTtrendHorn, file = 'temp/modonlyTtrendHorn.rds')
}

summary(modonlyTtrend)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -158937.6 -158831.8 79480.79
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev     Corr  
    ## (Intercept)    0.04375688 (Intr)
    ## abs(temptrend) 0.29351811 -0.123
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.001792571 0.2893301
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.223838 
    ## Fixed effects: Jtutrend ~ abs(temptrend) * REALM 
    ##                                      Value  Std.Error    DF   t-value p-value
    ## (Intercept)                     -0.0010209 0.01173020 49596 -0.087032  0.9306
    ## abs(temptrend)                   0.4667858 0.12697860 49596  3.676098  0.0002
    ## REALMMarine                      0.0459010 0.01255993   314  3.654554  0.0003
    ## REALMTerrestrial                 0.0146895 0.01271420   314  1.155364  0.2488
    ## abs(temptrend):REALMMarine      -0.1025169 0.13228517 49596 -0.774969  0.4384
    ## abs(temptrend):REALMTerrestrial -0.1944451 0.13908368 49596 -1.398044  0.1621
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.424                                
    ## REALMMarine                     -0.934  0.396                         
    ## REALMTerrestrial                -0.923  0.391  0.862                  
    ## abs(temptrend):REALMMarine       0.407 -0.960 -0.407 -0.375           
    ## abs(temptrend):REALMTerrestrial  0.387 -0.913 -0.362 -0.433  0.876    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.38476416 -0.29781039  0.08738927  0.54929160  6.52538535 
    ## 
    ## Number of Observations: 49916
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    317                  49916

``` r
summary(modonlyTtrendJbeta)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##       AIC       BIC logLik
    ##   -164158 -164052.2  82091
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev     Corr  
    ## (Intercept)    0.05691594 (Intr)
    ## abs(temptrend) 0.34739161 -0.16 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 1.28009e-06 0.3525846
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.394679 
    ## Fixed effects: Jbetatrend ~ abs(temptrend) * REALM 
    ##                                      Value  Std.Error    DF   t-value p-value
    ## (Intercept)                      0.0003338 0.01476240 49596  0.022611  0.9820
    ## abs(temptrend)                   0.6078271 0.14642231 49596  4.151192  0.0000
    ## REALMMarine                      0.0588638 0.01581834   314  3.721240  0.0002
    ## REALMTerrestrial                 0.0226529 0.01594186   314  1.420971  0.1563
    ## abs(temptrend):REALMMarine      -0.1448866 0.15243618 49596 -0.950474  0.3419
    ## abs(temptrend):REALMTerrestrial -0.2667107 0.16041942 49596 -1.662584  0.0964
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.428                                
    ## REALMMarine                     -0.933  0.399                         
    ## REALMTerrestrial                -0.926  0.396  0.864                  
    ## abs(temptrend):REALMMarine       0.411 -0.961 -0.411 -0.380           
    ## abs(temptrend):REALMTerrestrial  0.390 -0.913 -0.364 -0.437  0.877    
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -7.3578768 -0.1587203  0.2033006  0.6705745  6.5270460 
    ## 
    ## Number of Observations: 49916
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    317                  49916

``` r
summary(modonlyTtrendHorn)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -146392.5 -146286.9 73208.24
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)    0.0549979 (Intr)
    ## abs(temptrend) 0.3228314 -0.095
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.007989212 0.3034899
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.217439 
    ## Fixed effects: Horntrend ~ abs(temptrend) * REALM 
    ##                                      Value  Std.Error    DF   t-value p-value
    ## (Intercept)                     -0.0032335 0.01443800 48521 -0.223955  0.8228
    ## abs(temptrend)                   0.6867279 0.14372508 48521  4.778066  0.0000
    ## REALMMarine                      0.0557766 0.01567510   273  3.558295  0.0004
    ## REALMTerrestrial                 0.0211790 0.01570463   273  1.348586  0.1786
    ## abs(temptrend):REALMMarine      -0.2257689 0.15063739 48521 -1.498757  0.1339
    ## abs(temptrend):REALMTerrestrial -0.3604184 0.15821499 48521 -2.278029  0.0227
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.399                                
    ## REALMMarine                     -0.921  0.368                         
    ## REALMTerrestrial                -0.919  0.367  0.847                  
    ## abs(temptrend):REALMMarine       0.381 -0.954 -0.378 -0.350           
    ## abs(temptrend):REALMTerrestrial  0.363 -0.908 -0.334 -0.408  0.867    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.21250285 -0.32193926  0.06521414  0.53688251  6.03752738 
    ## 
    ## Number of Observations: 48800
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    276                  48800

### Plot the temp-only coefficients

``` r
# make table of coefficients
coefs <- as.data.frame(summary(modonlyTtrend)$tTable)
coefs2 <- as.data.frame(summary(modonlyTtrendJbeta)$tTable)
coefs3 <- as.data.frame(summary(modonlyTtrendHorn)$tTable)
coefs$mod <- 'Jtu'
coefs2$mod <- 'Jbeta'
coefs3$mod <- 'Horn'
rows1 <- which(grepl('temptrend', rownames(coefs))) # extract temperature effect
cols <- c('Value', 'Std.Error', 'mod')
allcoefs <- rbind(coefs[rows1, cols], coefs2[rows1, cols], coefs3[rows1, cols])
allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to marine effects
allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to terrestrial effects

allcoefs$lCI <- allcoefs$Value - allcoefs$Std.Error # lower confidence interval
allcoefs$uCI <- allcoefs$Value + allcoefs$Std.Error
allcoefs$y <- c(3, 2, 1, 2.9, 1.9, 0.9, 2.8, 1.8, 0.8) # y-values
allcoefs$col <- c(rep('black', 3), rep('light grey', 3), rep('dark grey', 3))
allcoefs$realm <- rep(c('Freshwater', 'Marine', 'Terrestrial'), 3)

par(las = 1, mai = c(0.8, 2, 0.1, 0.1))
plot(0,0, col = 'white', xlim=c(-0.02, 0.85), ylim = c(0.7,3), 
     yaxt='n', xlab = 'Turnover per |°C/yr|', ylab ='')
axis(2, at = 3:1, labels = c('Freshwater', 'Marine', 'Terrestrial'), cex.axis = 0.7)
abline(v = 0, col = 'grey')
for(i in 1:nrow(allcoefs)){
  with(allcoefs[i, ], points(Value, y, pch = 16, col = col))
  with(allcoefs[i, ], lines(x = c(lCI, uCI), y = c(y, y), col = col))
}
legend('bottomright', col = c('black', 'dark grey', 'light grey'), lwd = 1, pch = 16, 
       legend = c('Jaccard turnover', 'Jaccard total', 'Horn-Morisita'))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/modonlyTtrendsimp%20coefs-1.png)<!-- -->
\#\#\# Nicer plots of turnover vs. temperature data and model fit
Scatterplot, violin plots, and coefficient plots all together

``` r
# on macbook: fig.width=3, fig.height=2.375, fig.retina=3, out.width=3, out.height=2.375
# on external monitor: fig.width=6, fig.height=4.5
trends[temptrend <= -0.5, temptrendtext := 'Cooling']
trends[abs(temptrend) <= 0.1, temptrendtext := 'Stable']
trends[temptrend >= 0.5, temptrendtext := 'Warming']

trends[abs(rarefyID_y) < 35, latzone := 'Subtropics']
trends[abs(rarefyID_y) >= 35 & abs(rarefyID_x) < 66.56339, latzone := 'Temperate'] 
trends[abs(rarefyID_y) >= 66.56339, latzone := 'Polar']
trends[, latzone := factor(latzone, levels = c('Subtropics', 'Temperate', 'Polar'))]

p1 <- ggplot(trends, aes(temptrend, Jtutrend, color = REALM, fill = REALM, size = nyrBT)) +
  geom_point(na.rm = TRUE, shape = 16, alpha = 0.1) + 
  geom_smooth(data=subset(trends, abs(temptrend) < 0.75), method = 'gam', formula = y ~ s(x, bs = "cs"), 
              na.rm = TRUE) +
  scale_color_brewer(palette="Set1", name = 'Realm') + 
  scale_fill_brewer(palette="Set1", name = 'Realm') + 
  labs(x = 'Temperature trend (°C/yr)', y = 'Jaccard turnover', tag = 'A') +
  scale_size_continuous(range = c(1, 8), breaks = c(2, 5, 20)) +
  guides(size = guide_legend(title = 'Years',
                             override.aes = list(linetype=0, fill = NA, alpha = 1))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        legend.key.size = unit(0.5,"line"), 
        axis.text=element_text(size=8),
        axis.title=element_text(size=10))

p2 <- ggplot(trends[!is.na(temptrendtext), ], aes(temptrendtext, Jtutrend)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = 'grey') +
  labs(x = '', y = 'Jaccard turnover', tag = 'B') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10))

p3 <- ggplot(trends[abs(temptrend) >= 0.5 & !is.na(latzone), ], aes(latzone, Jtutrend)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = 'grey') + 
  labs(x = '', y = '', tag = 'C') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=7),
        axis.title=element_text(size=10))

p4 <- ggplot(allcoefs, aes(Value, y, group = mod, color = mod)) +
  geom_errorbarh(aes(xmin = lCI, xmax = uCI, height = 0)) + 
  geom_point() + 
  labs(x = expression(atop('Temperature change effect', '(Turnover '~degree*'C'^'-1'*')')), y = '', tag = 'D') +
  scale_color_grey() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position='none',
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) + 
  coord_cartesian(xlim =c(0, 1)) +
  scale_y_continuous(name = '', breaks = c(1, 2, 3), labels = c('Terrestrial', 'Marine', 'Freshwater'))

grid.arrange(p1, p2, p3, p4, ncol = 3, layout_matrix = rbind(c(1,1,1), c(2,3,4)),
             heights=c(unit(0.66, "npc"), unit(0.34, "npc")))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/turnover%20vs.%20temperature%20big%20plot-1.png)<!-- -->

## Full models

Try static covariates plus interactions of abs temperature trend with
each covariate:

  - realm
  - speed
  - mass
  - lifespan
  - average metabolic temperature
  - consumer fraction
  - endotherm fraction
  - environmental temperature
  - seasonality
  - microclimates
  - thermal bias
  - NPP
  - vegetation
  - human footprint

Except for thermal bias: interact with temperature trend (not abs)

### Full model for Jaccard total

``` r
# using Bowler for human impact
i <- trends[, complete.cases(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modTfull1.rds')){
  modTfull1 <- readRDS('temp/modTfull1.rds')
} else {
  modTfull1 <- lme(Jtutrend ~ temptrend_abs.sc*REALM +
                     temptrend_abs.sc*tsign + 
                     temptrend_abs.sc*tempave.sc +
                     temptrend_abs.sc*tempave_metab.sc + 
                     temptrend_abs.sc*seas.sc + 
                     temptrend_abs.sc*microclim.sc + 
                     temptrend_abs.sc*mass.sc + 
                     temptrend_abs.sc*speed.sc + 
                     temptrend_abs.sc*lifespan.sc + 
                     temptrend_abs.sc*consumerfrac.sc +
                     temptrend_abs.sc*endothermfrac.sc +
                     temptrend_abs.sc*nspp.sc +
                     temptrend_abs.sc*thermal_bias.sc:tsign +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*veg.sc +
                     temptrend_abs.sc*human_bowler.sc:REALM2,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modTfull1, file = 'temp/modTfull1.rds')
}

# using Venter/Halpern for human impact
i <- trends[, complete.cases(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_footprint.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modTfullfootprint.rds')){
  modTfullfootprint <- readRDS('temp/modTfullfootprint.rds')
} else {
  modTfullfootprint <- lme(Jtutrend ~ temptrend_abs.sc*REALM + 
                     temptrend_abs.sc*tsign + 
                     temptrend_abs.sc*tempave.sc +
                     temptrend_abs.sc*tempave_metab.sc + 
                     temptrend_abs.sc*seas.sc + 
                     temptrend_abs.sc*microclim.sc + 
                     temptrend_abs.sc*mass.sc + 
                     temptrend_abs.sc*speed.sc + 
                     temptrend_abs.sc*lifespan.sc + 
                     temptrend_abs.sc*consumerfrac.sc +
                     temptrend_abs.sc*endothermfrac.sc +
                     temptrend_abs.sc*nspp.sc +
                     temptrend_abs.sc*thermal_bias.sc:tsign +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*veg.sc +
                     temptrend_abs.sc*human_footprint.sc:REALM2,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modTfullfootprint, file = 'temp/modTfullfootprint.rds')
}

summary(modTfull1)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -133022.6 -132623.7 66557.29
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.04683428 (Intr)
    ## temptrend_abs.sc 0.03852623 -0.115
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 2.152948e-06 0.2862641
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.238848 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave.sc + temptrend_abs.sc * tempave_metab.sc +      temptrend_abs.sc * seas.sc + temptrend_abs.sc * microclim.sc +      temptrend_abs.sc * mass.sc + temptrend_abs.sc * speed.sc +      temptrend_abs.sc * lifespan.sc + temptrend_abs.sc * consumerfrac.sc +      temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.00440294 0.013574325 42907   0.324358  0.7457
    ## temptrend_abs.sc                                  0.05469078 0.019107051 42907   2.862335  0.0042
    ## REALMMarine                                       0.05506689 0.014427616   222   3.816770  0.0002
    ## REALMTerrestrial                                  0.03363773 0.015797059   222   2.129366  0.0343
    ## tsign1                                           -0.01039531 0.000589710 42907 -17.627842  0.0000
    ## tempave.sc                                       -0.00430156 0.000875433 42907  -4.913633  0.0000
    ## tempave_metab.sc                                  0.00655550 0.001963041 42907   3.339462  0.0008
    ## seas.sc                                          -0.00059689 0.000446905 42907  -1.335619  0.1817
    ## microclim.sc                                      0.00054683 0.000243427 42907   2.246396  0.0247
    ## mass.sc                                           0.00003329 0.000823940 42907   0.040403  0.9678
    ## speed.sc                                         -0.00180200 0.000614968 42907  -2.930232  0.0034
    ## lifespan.sc                                      -0.00428265 0.001590956 42907  -2.691874  0.0071
    ## consumerfrac.sc                                   0.00177252 0.000965237 42907   1.836357  0.0663
    ## endothermfrac.sc                                 -0.01808343 0.004089531 42907  -4.421884  0.0000
    ## nspp.sc                                          -0.00546412 0.000450892 42907 -12.118467  0.0000
    ## npp.sc                                            0.00084605 0.000380943 42907   2.220931  0.0264
    ## veg.sc                                            0.00078347 0.000466926 42907   1.677929  0.0934
    ## temptrend_abs.sc:REALMMarine                     -0.00726761 0.019766931 42907  -0.367665  0.7131
    ## temptrend_abs.sc:REALMTerrestrial                -0.02968306 0.023018110 42907  -1.289552  0.1972
    ## temptrend_abs.sc:tsign1                           0.00884131 0.001077008 42907   8.209139  0.0000
    ## temptrend_abs.sc:tempave.sc                       0.01487416 0.002131008 42907   6.979873  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                 0.00122358 0.004767730 42907   0.256638  0.7975
    ## temptrend_abs.sc:seas.sc                          0.00375569 0.001042268 42907   3.603385  0.0003
    ## temptrend_abs.sc:microclim.sc                    -0.00284268 0.000671360 42907  -4.234210  0.0000
    ## temptrend_abs.sc:mass.sc                         -0.00508926 0.001806487 42907  -2.817211  0.0048
    ## temptrend_abs.sc:speed.sc                         0.00252254 0.000913109 42907   2.762586  0.0057
    ## temptrend_abs.sc:lifespan.sc                      0.00799909 0.003719181 42907   2.150768  0.0315
    ## temptrend_abs.sc:consumerfrac.sc                 -0.00159405 0.001421502 42907  -1.121383  0.2621
    ## temptrend_abs.sc:endothermfrac.sc                 0.01293849 0.006444364 42907   2.007722  0.0447
    ## temptrend_abs.sc:nspp.sc                         -0.00678702 0.000957142 42907  -7.090925  0.0000
    ## tsign-1:thermal_bias.sc                          -0.00032048 0.000653545 42907  -0.490367  0.6239
    ## tsign1:thermal_bias.sc                           -0.00044784 0.000539433 42907  -0.830209  0.4064
    ## temptrend_abs.sc:npp.sc                          -0.00357444 0.001040147 42907  -3.436479  0.0006
    ## temptrend_abs.sc:veg.sc                           0.00038406 0.001304025 42907   0.294519  0.7684
    ## human_bowler.sc:REALM2Marine                      0.00111969 0.000347230 42907   3.224639  0.0013
    ## human_bowler.sc:REALM2TerrFresh                   0.00009230 0.000367867 42907   0.250912  0.8019
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.00504727 0.001213634 42907   4.158809  0.0000
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00183069 0.001179299 42907   1.552354  0.1206
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.00353947 0.000943215 42907  -3.752555  0.0002
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.00177712 0.000825450 42907  -2.152905  0.0313
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2M h_.:REALM2T t_.:-1 t_.:1: t_.:_.:REALM2M
    ## temptrend_abs.sc                                 -0.414                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.925  0.387                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.864  0.353  0.797                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.036  0.021  0.008  0.012                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.008  0.007  0.002 -0.004 -0.063                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.061 -0.028 -0.051 -0.040 -0.035 -0.468                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.051  0.046  0.069 -0.014 -0.107  0.189 -0.096                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.023  0.025  0.033 -0.008  0.011 -0.022 -0.031  0.176                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                           0.011  0.010  0.002  0.002  0.023  0.018 -0.492 -0.008  0.017                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.025 -0.006 -0.041 -0.028  0.037 -0.002  0.134 -0.066 -0.055 -0.127                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.042 -0.026 -0.032 -0.045 -0.031 -0.052  0.684  0.051  0.020 -0.806  0.248                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.017  0.007  0.038  0.261  0.006 -0.012 -0.074 -0.022  0.007  0.045 -0.158 -0.121                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                  0.150 -0.053 -0.080 -0.314  0.004  0.162 -0.207  0.042  0.024 -0.002  0.108  0.014 -0.351                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                          -0.012 -0.015 -0.008 -0.025 -0.041  0.009 -0.027 -0.026 -0.087 -0.206 -0.074  0.143  0.007  0.080                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.020 -0.004 -0.025  0.011  0.016 -0.309  0.263 -0.192 -0.214 -0.050 -0.011  0.065 -0.019 -0.085 -0.178                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.117  0.149  0.121 -0.001 -0.016  0.170 -0.125  0.153  0.031  0.037 -0.055 -0.024 -0.034  0.015  0.012 -0.335                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.397 -0.958 -0.398 -0.331 -0.005 -0.006  0.022 -0.060 -0.029 -0.014  0.016  0.017 -0.013  0.034  0.022  0.020 -0.163                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.340 -0.835 -0.312 -0.417 -0.006 -0.013  0.025  0.039  0.016 -0.015  0.010  0.030 -0.089  0.138  0.040 -0.025  0.012  0.785                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.018 -0.046 -0.006 -0.005 -0.486  0.043  0.003  0.014 -0.037 -0.013 -0.011  0.013 -0.001  0.006  0.027 -0.049  0.003  0.014      0.009                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.007 -0.020 -0.009 -0.004  0.044 -0.600  0.315 -0.058  0.075  0.024  0.001  0.018 -0.020 -0.122 -0.016  0.164 -0.152  0.013      0.026     -0.068                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.018  0.036  0.015  0.017 -0.007  0.328 -0.605 -0.002 -0.024  0.253 -0.022 -0.369  0.055  0.151 -0.008 -0.156  0.071 -0.024     -0.044      0.013 -0.423                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.025 -0.067 -0.033  0.028  0.044 -0.012 -0.005 -0.627 -0.108  0.011  0.013 -0.023  0.024  0.005  0.020  0.102 -0.139  0.089     -0.079     -0.020  0.041            0.056                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.012 -0.033 -0.017  0.008 -0.040  0.068 -0.009 -0.087 -0.691 -0.008  0.015 -0.004 -0.002  0.000  0.039  0.158  0.012  0.038     -0.032      0.048 -0.136            0.112   0.146                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.007 -0.016 -0.009 -0.011 -0.015  0.025  0.271  0.008 -0.002 -0.555  0.023  0.470 -0.016  0.016  0.131  0.008 -0.004  0.020      0.029      0.024 -0.079           -0.469  -0.049            -0.017                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.007  0.009  0.016  0.008 -0.019  0.018 -0.067  0.031  0.020  0.046 -0.563 -0.108  0.093 -0.065  0.059  0.003  0.049 -0.026     -0.011     -0.003 -0.034            0.064  -0.038            -0.005            -0.021                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.016  0.041  0.013  0.019  0.022  0.020 -0.378 -0.019 -0.004  0.445 -0.051 -0.557  0.067 -0.006 -0.115 -0.021 -0.014 -0.021     -0.050     -0.020 -0.022            0.663   0.058             0.016            -0.824             0.086                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.003 -0.073 -0.010 -0.090 -0.003 -0.022  0.083  0.023  0.001 -0.031  0.094  0.100 -0.480  0.119  0.001  0.020  0.015  0.081      0.294      0.004  0.055           -0.157  -0.029            -0.005             0.059            -0.173            -0.171                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                -0.046  0.152  0.029  0.136  0.007 -0.205  0.247  0.000  0.012  0.018 -0.077  0.002  0.101 -0.395 -0.063  0.093 -0.042 -0.102     -0.365     -0.021  0.254           -0.401  -0.028            -0.062            -0.019             0.110            -0.030           -0.301                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.010  0.032  0.011  0.028  0.022 -0.026  0.021  0.010  0.034  0.133  0.051 -0.118  0.003 -0.065 -0.535  0.081 -0.037 -0.040     -0.075     -0.036 -0.002            0.033  -0.028            -0.015            -0.212            -0.093             0.205           -0.006            0.108                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.022 -0.006 -0.025 -0.003 -0.171  0.453 -0.090 -0.110 -0.077  0.030 -0.002 -0.072  0.010  0.010 -0.007 -0.087  0.041  0.012     -0.010      0.022 -0.258            0.073   0.088             0.049            -0.006             0.012             0.032           -0.019           -0.016            -0.002                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.035 -0.015 -0.046 -0.010 -0.028  0.684 -0.147 -0.237 -0.187  0.051  0.026 -0.076  0.007  0.029 -0.043 -0.125  0.050  0.016     -0.006      0.070 -0.334            0.096   0.114             0.081            -0.008            -0.008             0.028           -0.017           -0.041            -0.002             0.571                                                                                               
    ## temptrend_abs.sc:npp.sc                           0.004  0.018  0.001 -0.019 -0.013  0.117 -0.102  0.084  0.160  0.013  0.003 -0.019  0.009  0.030  0.070 -0.636  0.212 -0.041      0.044      0.039 -0.195            0.165  -0.214            -0.280             0.018             0.009             0.017           -0.035           -0.091            -0.088             0.031  0.057                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.078 -0.205 -0.082  0.008  0.009 -0.108  0.058 -0.133  0.012 -0.010  0.029 -0.001  0.010 -0.013 -0.032  0.210 -0.709  0.227     -0.027      0.017  0.188           -0.065   0.233             0.005            -0.019            -0.062             0.040           -0.011            0.044             0.035            -0.030 -0.049 -0.347                                                                                 
    ## human_bowler.sc:REALM2Marine                      0.011 -0.016 -0.013 -0.002  0.004 -0.020  0.023 -0.155 -0.090 -0.019  0.010  0.000 -0.003  0.001 -0.043 -0.118  0.023  0.022      0.004      0.014  0.023           -0.027   0.143             0.026             0.018             0.007             0.002            0.014            0.003             0.052             0.072  0.077  0.068             0.004                                                               
    ## human_bowler.sc:REALM2TerrFresh                  -0.031  0.060  0.024 -0.010  0.033 -0.257  0.096 -0.058  0.155  0.043  0.046 -0.052 -0.054 -0.071  0.016 -0.088  0.163 -0.065      0.021     -0.025  0.139           -0.078  -0.148            -0.049            -0.018            -0.021             0.015            0.019            0.063             0.035            -0.068 -0.085  0.069            -0.250            0.009                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.004 -0.005  0.004 -0.004  0.018 -0.349  0.108  0.090  0.086 -0.012  0.000  0.042 -0.019 -0.023 -0.001  0.056 -0.073 -0.012      0.026      0.103  0.678           -0.113  -0.117            -0.112             0.005            -0.007            -0.064            0.038            0.018            -0.042            -0.469 -0.367 -0.077             0.094           -0.061       0.047                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.007  0.001  0.006 -0.002  0.064 -0.384  0.105  0.096  0.116 -0.008 -0.010  0.038 -0.017 -0.021  0.009  0.066 -0.084 -0.010      0.022     -0.128  0.716           -0.091  -0.088            -0.145            -0.007             0.012            -0.054            0.034            0.010            -0.054            -0.336 -0.510 -0.094             0.116           -0.068       0.038       0.792                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.010  0.022  0.011  0.002  0.000  0.012 -0.022  0.134  0.037  0.016  0.007 -0.007  0.006  0.001  0.027  0.067  0.007 -0.025      0.001     -0.009 -0.049            0.064  -0.218            -0.064            -0.019            -0.017            -0.002           -0.016           -0.031            -0.065            -0.058 -0.068 -0.130             0.005           -0.672       0.006       0.065  0.066               
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.037 -0.117 -0.038  0.014 -0.014  0.159 -0.081 -0.121 -0.061 -0.027 -0.013  0.027  0.018  0.048  0.016  0.053 -0.304  0.128     -0.037      0.047 -0.210            0.115   0.267             0.095             0.038             0.025            -0.030           -0.028           -0.093            -0.038             0.048  0.055 -0.131             0.492            0.005      -0.631      -0.091 -0.087 -0.006        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.34720709 -0.35408008  0.05210934  0.52828824  6.60624447 
    ## 
    ## Number of Observations: 43169
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    225                  43169

``` r
summary(modTfullfootprint)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -133008.5 -132609.6 66550.27
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.04664700 (Intr)
    ## temptrend_abs.sc 0.03913752 -0.111
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 1.335773e-06 0.2861977
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.238589 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave.sc + temptrend_abs.sc * tempave_metab.sc +      temptrend_abs.sc * seas.sc + temptrend_abs.sc * microclim.sc +      temptrend_abs.sc * mass.sc + temptrend_abs.sc * speed.sc +      temptrend_abs.sc * lifespan.sc + temptrend_abs.sc * consumerfrac.sc +      temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_footprint.sc:REALM2 
    ##                                                           Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                          0.00506425 0.013546761 42907   0.373834  0.7085
    ## temptrend_abs.sc                                     0.05026013 0.019249213 42907   2.611023  0.0090
    ## REALMMarine                                          0.05419204 0.014393916   222   3.764927  0.0002
    ## REALMTerrestrial                                     0.03370578 0.015770337   222   2.137290  0.0337
    ## tsign1                                              -0.01032468 0.000590108 42907 -17.496260  0.0000
    ## tempave.sc                                          -0.00437495 0.000857197 42907  -5.103788  0.0000
    ## tempave_metab.sc                                     0.00639190 0.001958767 42907   3.263225  0.0011
    ## seas.sc                                             -0.00073859 0.000437846 42907  -1.686868  0.0916
    ## microclim.sc                                         0.00065090 0.000240500 42907   2.706453  0.0068
    ## mass.sc                                              0.00014962 0.000823894 42907   0.181605  0.8559
    ## speed.sc                                            -0.00181056 0.000615299 42907  -2.942570  0.0033
    ## lifespan.sc                                         -0.00433486 0.001591270 42907  -2.724151  0.0064
    ## consumerfrac.sc                                      0.00169485 0.000963783 42907   1.758535  0.0787
    ## endothermfrac.sc                                    -0.01809365 0.004079957 42907  -4.434766  0.0000
    ## nspp.sc                                             -0.00533365 0.000451616 42907 -11.810134  0.0000
    ## npp.sc                                               0.00094939 0.000377683 42907   2.513723  0.0120
    ## veg.sc                                               0.00045684 0.000455722 42907   1.002446  0.3161
    ## temptrend_abs.sc:REALMMarine                        -0.00244289 0.019895364 42907  -0.122787  0.9023
    ## temptrend_abs.sc:REALMTerrestrial                   -0.03136929 0.023240533 42907  -1.349766  0.1771
    ## temptrend_abs.sc:tsign1                              0.00892584 0.001076603 42907   8.290742  0.0000
    ## temptrend_abs.sc:tempave.sc                          0.01355708 0.002112646 42907   6.417109  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                    0.00326975 0.004741165 42907   0.689651  0.4904
    ## temptrend_abs.sc:seas.sc                             0.00358430 0.000997483 42907   3.593343  0.0003
    ## temptrend_abs.sc:microclim.sc                       -0.00303599 0.000673140 42907  -4.510194  0.0000
    ## temptrend_abs.sc:mass.sc                            -0.00498346 0.001806778 42907  -2.758200  0.0058
    ## temptrend_abs.sc:speed.sc                            0.00251670 0.000913958 42907   2.753629  0.0059
    ## temptrend_abs.sc:lifespan.sc                         0.00774495 0.003724768 42907   2.079309  0.0376
    ## temptrend_abs.sc:consumerfrac.sc                    -0.00177757 0.001433793 42907  -1.239765  0.2151
    ## temptrend_abs.sc:endothermfrac.sc                    0.01118765 0.006482681 42907   1.725776  0.0844
    ## temptrend_abs.sc:nspp.sc                            -0.00711565 0.000959061 42907  -7.419396  0.0000
    ## tsign-1:thermal_bias.sc                             -0.00047954 0.000650943 42907  -0.736677  0.4613
    ## tsign1:thermal_bias.sc                              -0.00058788 0.000537384 42907  -1.093963  0.2740
    ## temptrend_abs.sc:npp.sc                             -0.00457894 0.001025852 42907  -4.463549  0.0000
    ## temptrend_abs.sc:veg.sc                              0.00223365 0.001222602 42907   1.826962  0.0677
    ## human_footprint.sc:REALM2Marine                     -0.00013794 0.000245615 42907  -0.561623  0.5744
    ## human_footprint.sc:REALM2TerrFresh                   0.00032388 0.000327108 42907   0.990117  0.3221
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc             0.00506146 0.001208126 42907   4.189510  0.0000
    ## temptrend_abs.sc:tsign1:thermal_bias.sc              0.00180363 0.001177048 42907   1.532334  0.1254
    ## temptrend_abs.sc:human_footprint.sc:REALM2Marine    -0.00109822 0.000672981 42907  -1.631866  0.1027
    ## temptrend_abs.sc:human_footprint.sc:REALM2TerrFresh  0.00015767 0.000811210 42907   0.194365  0.8459
    ##  Correlation: 
    ##                                                     (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2M h_.:REALM2T t_.:-1 t_.:1: t_.:_.:REALM2M
    ## temptrend_abs.sc                                    -0.413                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                         -0.925  0.385                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                    -0.864  0.355  0.798                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                              -0.035  0.020  0.007  0.012                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                          -0.011  0.011  0.003 -0.006 -0.059                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                     0.063 -0.030 -0.053 -0.040 -0.036 -0.463                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                             -0.049  0.034  0.066 -0.013 -0.111  0.184 -0.094                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                        -0.019  0.017  0.030 -0.006  0.005  0.013 -0.046  0.180                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                              0.011  0.009  0.001  0.002  0.023  0.024 -0.495 -0.012  0.009                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                             0.024 -0.005 -0.041 -0.027  0.034  0.006  0.130 -0.059 -0.060 -0.128                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                          0.042 -0.026 -0.033 -0.046 -0.028 -0.061  0.690  0.048  0.025 -0.806  0.248                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                     -0.017  0.007  0.038  0.260  0.007 -0.020 -0.071 -0.026  0.014  0.045 -0.157 -0.122                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                     0.150 -0.052 -0.079 -0.316  0.006  0.155 -0.204  0.039  0.033 -0.001  0.110  0.013 -0.353                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                             -0.010 -0.016 -0.009 -0.026 -0.040  0.016 -0.028 -0.033 -0.100 -0.207 -0.078  0.147  0.009  0.084                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                               0.019  0.000 -0.025  0.010  0.017 -0.341  0.276 -0.214 -0.215 -0.050 -0.005  0.060 -0.023 -0.093 -0.187                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                              -0.116  0.135  0.119  0.003 -0.018  0.191 -0.137  0.126  0.017  0.037 -0.054 -0.025 -0.037  0.020  0.013 -0.330                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                         0.396 -0.958 -0.396 -0.334 -0.004 -0.011  0.024 -0.046 -0.020 -0.013  0.015  0.017 -0.013  0.032  0.024  0.016 -0.147                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                    0.342 -0.842 -0.314 -0.417 -0.007 -0.008  0.022  0.035  0.014 -0.016  0.010  0.031 -0.088  0.140  0.041 -0.023  0.001  0.793                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                              0.018 -0.043 -0.005 -0.006 -0.486  0.043  0.004  0.022 -0.033 -0.013 -0.012  0.013 -0.001  0.005  0.028 -0.050  0.009  0.011      0.011                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                          0.009 -0.027 -0.011 -0.001  0.042 -0.603  0.311 -0.068  0.059  0.023  0.002  0.019 -0.019 -0.119 -0.019  0.182 -0.185  0.022      0.019     -0.065                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                   -0.020  0.042  0.017  0.015 -0.005  0.325 -0.605  0.000 -0.017  0.256 -0.022 -0.372  0.054  0.150 -0.007 -0.169  0.093 -0.031     -0.040      0.011 -0.413                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                             0.019 -0.050 -0.026  0.026  0.050 -0.032 -0.001 -0.619 -0.086  0.019  0.016 -0.029  0.024  0.001  0.030  0.107 -0.096  0.071     -0.071     -0.030  0.060            0.054                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                        0.007 -0.017 -0.011  0.007 -0.040  0.058 -0.004 -0.059 -0.704 -0.005  0.019 -0.007 -0.003 -0.002  0.039  0.155  0.053  0.020     -0.027      0.041 -0.115            0.106   0.107                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                             0.007 -0.015 -0.009 -0.011 -0.015  0.026  0.271  0.014  0.002 -0.555  0.022  0.471 -0.016  0.016  0.132  0.009  0.000  0.018      0.029      0.023 -0.079           -0.472  -0.060            -0.026                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                           -0.007  0.008  0.016  0.008 -0.018  0.019 -0.068  0.033  0.023  0.046 -0.563 -0.108  0.094 -0.065  0.060  0.002  0.049 -0.025     -0.010     -0.003 -0.038            0.065  -0.043            -0.010            -0.020                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                        -0.016  0.042  0.013  0.019  0.022  0.019 -0.378 -0.020 -0.006  0.445 -0.050 -0.558  0.067 -0.007 -0.116 -0.019 -0.014 -0.022     -0.051     -0.020 -0.021            0.667   0.060             0.020            -0.824             0.084                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                     0.003 -0.074 -0.010 -0.090 -0.003 -0.022  0.081  0.024 -0.001 -0.030  0.094  0.098 -0.481  0.119  0.001  0.023  0.013  0.082      0.294      0.004  0.053           -0.155  -0.030            -0.001             0.058            -0.172            -0.168                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                   -0.045  0.150  0.028  0.138  0.005 -0.201  0.243 -0.004  0.004  0.017 -0.076  0.003  0.103 -0.395 -0.065  0.099 -0.055 -0.099     -0.369     -0.019  0.246           -0.394  -0.025            -0.056            -0.018             0.108            -0.029           -0.305                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                            -0.012  0.036  0.013  0.029  0.020 -0.027  0.020  0.026  0.035  0.134  0.053 -0.120  0.004 -0.067 -0.538  0.093 -0.033 -0.045     -0.076     -0.038  0.001            0.035  -0.049            -0.016            -0.214            -0.096             0.208           -0.006            0.109                                                                                                                                   
    ## tsign-1:thermal_bias.sc                              0.020 -0.002 -0.023 -0.004 -0.169  0.453 -0.086 -0.102 -0.064  0.034 -0.001 -0.074  0.008  0.006 -0.002 -0.086  0.049  0.006     -0.009      0.019 -0.256            0.071   0.069             0.042            -0.009             0.011             0.033           -0.020           -0.013            -0.006                                                                                                                 
    ## tsign1:thermal_bias.sc                               0.034 -0.012 -0.045 -0.011 -0.026  0.691 -0.144 -0.234 -0.172  0.055  0.026 -0.078  0.006  0.025 -0.036 -0.124  0.053  0.013     -0.004      0.068 -0.337            0.096   0.096             0.070            -0.010            -0.009             0.028           -0.019           -0.040            -0.007             0.567                                                                                               
    ## temptrend_abs.sc:npp.sc                              0.006  0.008 -0.001 -0.017 -0.015  0.139 -0.117  0.090  0.156  0.012  0.006 -0.018  0.011  0.036  0.076 -0.639  0.187 -0.030      0.040      0.044 -0.230            0.192  -0.223            -0.268             0.019             0.008             0.014           -0.039           -0.109            -0.101             0.029  0.053                                                                                        
    ## temptrend_abs.sc:veg.sc                              0.076 -0.191 -0.078  0.002  0.016 -0.147  0.079 -0.101  0.049 -0.005  0.026 -0.003  0.011 -0.023 -0.031  0.201 -0.689  0.209     -0.011      0.007  0.252           -0.106   0.182            -0.066            -0.025            -0.063             0.041           -0.009            0.070             0.027            -0.048 -0.060 -0.319                                                                                 
    ## human_footprint.sc:REALM2Marine                     -0.011  0.005  0.008  0.013 -0.032  0.010 -0.017  0.061  0.074 -0.011  0.025 -0.021  0.009 -0.022 -0.065  0.054 -0.025 -0.004     -0.006      0.006  0.001           -0.014  -0.017            -0.064             0.008             0.003             0.004           -0.002            0.020             0.049            -0.006  0.002 -0.023             0.014                                                               
    ## human_footprint.sc:REALM2TerrFresh                  -0.031  0.056  0.025  0.003  0.015 -0.155  0.045 -0.053  0.036  0.035  0.039 -0.043 -0.047 -0.043 -0.022 -0.019  0.172 -0.058     -0.001     -0.023  0.113           -0.048  -0.077             0.009            -0.023            -0.033             0.029            0.028            0.048             0.056            -0.034 -0.067  0.015            -0.201           -0.006                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc            -0.002 -0.013  0.001 -0.003  0.017 -0.353  0.106  0.072  0.077 -0.014  0.000  0.044 -0.019 -0.021 -0.004  0.054 -0.094 -0.003      0.022      0.107  0.680           -0.110  -0.091            -0.097             0.008            -0.006            -0.065            0.038            0.014            -0.038            -0.466 -0.363 -0.081             0.133           -0.002       0.025                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc             -0.006 -0.004  0.004 -0.001  0.063 -0.394  0.105  0.080  0.107 -0.010 -0.008  0.037 -0.018 -0.020  0.006  0.063 -0.097 -0.004      0.019     -0.126  0.723           -0.090  -0.067            -0.125            -0.006             0.011            -0.053            0.035            0.008            -0.047            -0.331 -0.509 -0.095             0.141           -0.014       0.047       0.790                      
    ## temptrend_abs.sc:human_footprint.sc:REALM2Marine     0.001 -0.008 -0.001 -0.004  0.002  0.016 -0.017  0.008 -0.063 -0.003  0.011  0.006 -0.008  0.023  0.036 -0.030  0.016  0.006      0.005     -0.004 -0.007            0.026   0.037             0.125            -0.015            -0.006            -0.009            0.007           -0.032            -0.035            -0.002 -0.022  0.071            -0.039           -0.656      -0.006       0.011  0.031               
    ## temptrend_abs.sc:human_footprint.sc:REALM2TerrFresh  0.031 -0.095 -0.030  0.001 -0.001  0.121 -0.060 -0.072  0.028 -0.021 -0.021  0.028  0.020  0.038  0.038 -0.005 -0.218  0.099     -0.005      0.036 -0.168            0.070   0.181            -0.041             0.039             0.041            -0.049           -0.032           -0.073            -0.081             0.016  0.048 -0.032             0.373            0.000      -0.649      -0.044 -0.078 -0.001        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.35903450 -0.35316683  0.05317953  0.52698821  6.60662084 
    ## 
    ## Number of Observations: 43169
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    225                  43169

#### Plot the coefficients from the full model

Bowler human footprint in black, Venter/Halpern in grey

``` r
coefs <- summary(modTfull1)$tTable
coefs2 <- summary(modTfullfootprint)$tTable
par(las = 1, mai = c(0.5, 3, 0.1, 0.1))
rows1 <- which(!grepl('Intercept', rownames(coefs)))
plot(0,0, col = 'white', xlim=c(-0.05, 0.08), ylim = c(1,length(rows1)), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(rows1):1, labels = rownames(coefs)[rows1], cex.axis = 0.7)
abline(v = 0, col = 'grey')
for(i in 1:length(rows1)){
  x = coefs[rows1[i], 1]
  se = coefs[rows1[i], 2]
  points(x, length(rows1) + 1 - i, pch = 16)
  lines(x = c(x-se, x+se), y = c(length(rows1) + 1 - i, length(rows1) + 1 - i))

  x = coefs2[rows1[i], 1]
  se = coefs2[rows1[i], 2]
  points(x, length(rows1) + 1 - i, pch = 16, col = 'grey', cex = 0.75)
  lines(x = c(x-se, x+se), y = c(length(rows1) + 1 - i, length(rows1) + 1 - i), col = 'grey')
}
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20fullTmod1-1.png)<!-- -->

#### Plot interactions

``` r
# set up the interactions to plot
ints <- data.frame(vars = c('tsign', 'tempave', 'tempave_metab', 'seas', 'microclim', 'mass', 'speed', 'lifespan', 'consumerfrac', 'endothermfrac', 'nspp', 'thermal_bias', 'npp', 'veg', 'human_bowler'),
           min = c(1, -10, 10, 0.1, 0, 0, 0, 0.3, 0, 0, 0.3, -10, 1.9, 0, 0), 
           max = c(2, 30, 40, 16, 6, 8, 2, 2, 1, 1, 2.6, 10, 3.7, 1, 9),
           log = c(F, F, F, F, F, T, T, T, F, F, T, F, T, F, F),
           len = c(2, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
           discrete = c(T, F, F, F, F, F, F, F, F, F, F, F, F, F, F),
           stringsAsFactors = FALSE)
basetab <- data.frame(REALM = 'Freshwater', REALM2 = 'TerrFresh', 
                      tempave.sc = 0, tempave_metab.sc = 0, 
                      seas.sc = 0, microclim.sc = 0, mass.sc = 0, 
                      speed.sc = 0, lifespan.sc = 0, endothermfrac.sc = 0, 
                      nspp.sc = 0, thermal_bias.sc = 0, npp.sc = 0, human_bowler.sc = 0, veg.sc = 0,
                      consumerfrac.sc = 0,
                      nyrBT = 20, STUDY_ID = 127L, rarefyID = '127_514668')

# make the data frames for each interaction to plot                
for(j in 1:nrow(ints)){
  # set up a grid of temperature trends and the interacting variable
  if(ints$log[j]) intvars <- list(temptrend = seq(-1.5, 1.5, length.out = 100), 
                                  new = 10^seq(ints$min[j], ints$max[j], length.out = ints$len[j]),
                                   var = ints$vars[j])
  if(!ints$log[j]) intvars <- list(temptrend = seq(-1.5, 1.5, length.out = 100), 
                                   new = seq(ints$min[j], ints$max[j], length.out = ints$len[j]),
                                   var = ints$vars[j])
  names(intvars) <- c('temptrend', ints$vars[j], 'var')
  thisdat <- expand.grid(intvars)
  
  # scale the interacting variable
  cent <- attr(trends[[paste0(ints$var[j], '.sc')]], 'scaled:center')
  scl <- attr(trends[[paste0(ints$var[j], '.sc')]], 'scaled:scale')
  if(!is.null(cent) & !is.null(scl)){
    if(ints$log[j]) thisdat[[paste0(ints$var[j], '.sc')]] <- (log(thisdat[[ints$var[j]]]) - cent)/scl
    if(!ints$log[j]) thisdat[[paste0(ints$var[j], '.sc')]] <- (thisdat[[ints$var[j]]] - cent)/scl
  }

  # merge with the rest of the columns
  if(ints$var[j] != 'tsign') colnamestouse <- setdiff(colnames(basetab), paste0(ints$var[j], '.sc'))
  if(ints$var[j] == 'tsign') colnamestouse <- setdiff(colnames(basetab), ints$var[j])
  thisdat <- cbind(thisdat, basetab[, colnamestouse])

  # merge with the previous iterations
  if(j == 1) newdat <- thisdat
  if(j > 1){
    colstoadd <- setdiff(colnames(thisdat), colnames(newdat))
    for(toadd in colstoadd){
      newdat[[toadd]] <- NA
    }
    
    colstoadd2 <- setdiff(colnames(newdat), colnames(thisdat))
    for(toadd in colstoadd2){
      thisdat[[toadd]] <- NA
    }
    
    newdat <- rbind(newdat, thisdat)
  } 
}

# add two extra rows so that all factor levels are represented (for predict.lme to work)
newdat <- rbind(newdat[1:2, ], newdat)
newdat$REALM <- as.character(newdat$REALM)
newdat$REALM2 <- as.character(newdat$REALM2)
newdat$REALM[1:2] <- c('Marine', 'Terrestrial')
newdat$REALM2[1:2] <- c('Marine', 'Marine')

# trim to at least some temperature change (so that tsign is -1 or 1)
newdat <- newdat[newdat$temptrend != 0,]

# scale the temperature vars
newdat$temptrend.sc <- newdat$temptrend/attr(trends$temptrend.sc, 'scaled:scale') 
newdat$temptrend_abs <- abs(newdat$temptrend)
newdat$temptrend_abs.sc <- (newdat$temptrend_abs)/attr(trends$temptrend_abs.sc, 'scaled:scale')
newdat$tsign <- factor(sign(newdat$temptrend))

# make predictions
newdat$preds <- predict(object = modTfull1, newdata = newdat, level = 0)

#remove the extra rows
newdat <- newdat[newdat$REALM == 'Freshwater', ]

# prep the plots
intplots <- vector('list', nrow(ints))
for(j in 1:nrow(ints)){
  thisplot <- ggplot(newdat[newdat$var == ints$vars[j], ], 
                     aes_string(x = 'temptrend', y = 'preds', 
                                group = ints$vars[j], 
                                color = ints$vars[j])) +
    geom_line() +
    coord_cartesian(ylim = c(0, 1)) +
    theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm'))
  if(ints$log[j] & !ints$discrete[j]){
    intplots[[j]] <- thisplot + scale_color_distiller(palette = "YlGnBu", trans = 'log')
  }
  if(!ints$log[j] & !ints$discrete[j]){
    intplots[[j]] <- thisplot + scale_color_distiller(palette = "YlGnBu", trans = 'identity')
  }
  if(ints$discrete[j]){
    intplots[[j]] <- thisplot + scale_color_brewer(palette = "Dark2")
  }
}

#grid.arrange(grobs = intplots, '+', theme(plot.margin = unit(c(0,0,0,0), 'cm'))), ncol=2)
#do.call('grid.arrange', c(intplots, ncol = 2))
grid.arrange(grobs = intplots, ncol = 3)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/interaction%20plots%20modTfull1-1.png)<!-- -->
Positive thermal bias means that the species in the assemblage have
warmer thermal niches than the environment (and vice versa).

#### Plot residuals against each predictor

``` r
resids <- resid(modTfull1)
preds <- getData(modTfull1)
col = '#00000033'
cex = 0.5
par(mfrow = c(5,4))
boxplot(resids ~ preds$REALM, cex = cex, col = col)
plot(preds$temptrend_abs.sc, resids, cex = cex, col = col)
plot(preds$tsign, resids, cex = cex, col = col)
plot(preds$tempave.sc, resids, cex = cex, col = col)
plot(preds$tempave_metab.sc, resids, cex = cex, col = col)
plot(preds$seas.sc, resids, cex = cex, col = col)
plot(preds$microclim.sc, resids, cex = cex, col = col)
plot(preds$mass.sc, resids, cex = cex, col = col)
plot(preds$speed.sc, resids, cex = cex, col = col)
plot(preds$lifespan.sc, resids, cex = cex, col = col)
plot(preds$consumerfrac.sc, resids, cex = cex, col = col)
plot(preds$endothermfrac.sc, resids, cex = cex, col = col)
plot(preds$nspp.sc, resids, cex = cex, col = col)
plot(preds$thermal_bias.sc, resids, cex = cex, col = col)
plot(preds$npp.sc, resids, cex = cex, col = col)
plot(preds$veg.sc, resids, cex = cex, col = col)
plot(preds$human_bowler.sc, resids, cex = cex, col = col)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/resids%20modTfull1-1.png)<!-- -->

#### Remove each term from the full model

``` r
if(file.exists('output/aics_from_full.csv')){
  aicsfromfull <- read.csv('output/aics_from_full.csv')

} else {
  i <- trends[, complete.cases(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                               temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                               consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
  
  randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
  varef <- varPower(-0.5, ~nyrBT)
  
  terms <- c('temptrend_abs.sc*REALM', 
             'temptrend_abs.sc*tsign',
             'temptrend_abs.sc*tempave.sc',
             'temptrend_abs.sc*tempave_metab.sc',
             'temptrend_abs.sc*seas.sc',
             'temptrend_abs.sc*microclim.sc',
             'temptrend_abs.sc*mass.sc',
             'temptrend_abs.sc*speed.sc', 
             'temptrend_abs.sc*lifespan.sc', 
             'temptrend_abs.sc*consumerfrac.sc',
             'temptrend_abs.sc*endothermfrac.sc',
             'temptrend_abs.sc*nspp.sc',
             'temptrend_abs.sc*thermal_bias.sc:tsign',
             'temptrend_abs.sc*npp.sc',
             'temptrend_abs.sc*veg.sc',
             'temptrend_abs.sc*human_bowler.sc:REALM2')
  modTdrops <- vector('list', length(terms)+2)
  names(modTdrops) <- c('full', '-temptrend_abs.sc', paste0('-', terms))
  
  # fit full model with ML for model comparison
  modTdrops[[1]] <- lme(formula(paste0('Jtutrend ~ ', paste(terms, collapse = ' + '))),
                        random = randef, weights = varef, data = trends[i,], method = 'ML')
  
  # w/out temptrend
  modTdrops[[2]] <- lme(formula(paste0('Jtutrend ~ ', paste(gsub('temptrend_abs.sc\\*', '', terms), collapse = ' + '))),
                        random = list(STUDY_ID = ~ 1, rarefyID = ~1), weights = varef, data = trends[i,], method = 'ML')
  
  for(j in 1:length(terms)){
    print(j)
    modTdrops[[j+2]] <- lme(formula(paste0('Jtutrend ~ ', paste(terms[-j], collapse = ' + '))),
                            random = randef, weights = varef, data = trends[i,], method = 'ML')
  }
  
  aics <- sapply(modTdrops, AIC)
  
  aicsfromfull <- data.frame(mod = names(aics), dAIC_Jtu = aics - aics[1])
  
  write.csv(aicsfromfull, file = 'output/aics_from_full.csv', row.names = FALSE)
}

aicsfromfull
```

    ##                                         mod      dAIC_Jtu   dAIC_Jbeta    dAIC_Horn
    ## 1                                      full    0.00000000     0.000000    0.0000000
    ## 2                         -temptrend_abs.sc 7392.04492882 12840.286245 7605.6472326
    ## 3                   -temptrend_abs.sc*REALM   14.87466466     7.098437   10.4577563
    ## 4                   -temptrend_abs.sc*tsign  307.26456844   589.870020  240.0240794
    ## 5              -temptrend_abs.sc*tempave.sc   45.08557421   109.482738   30.5043643
    ## 6        -temptrend_abs.sc*tempave_metab.sc   14.83704694   128.758735    9.1313302
    ## 7                 -temptrend_abs.sc*seas.sc   11.03831929    36.657551   22.4356323
    ## 8            -temptrend_abs.sc*microclim.sc   12.88175453    73.883018    2.6523187
    ## 9                 -temptrend_abs.sc*mass.sc    7.95995370     5.221911   -3.4124476
    ## 10               -temptrend_abs.sc*speed.sc    5.05076864    14.138661   -2.7615998
    ## 11            -temptrend_abs.sc*lifespan.sc    4.07127373    -1.204624    1.6582203
    ## 12        -temptrend_abs.sc*consumerfrac.sc   -0.07879609     1.815516   -0.9988319
    ## 13       -temptrend_abs.sc*endothermfrac.sc   13.73066680    14.590190    5.7036860
    ## 14                -temptrend_abs.sc*nspp.sc  398.21647587  1568.523765  732.0112555
    ## 15  -temptrend_abs.sc*thermal_bias.sc:tsign   22.37488755   105.350725   17.6958524
    ## 16                 -temptrend_abs.sc*npp.sc    8.71935851     3.171746    2.9179283
    ## 17                 -temptrend_abs.sc*veg.sc    3.74525115    17.061611   -2.3119391
    ## 18 -temptrend_abs.sc*human_bowler.sc:REALM2   10.86674547     1.996388   14.3718703

### Sensitivity analysis: total turnover and Horn-Morisita models

#### Fit full models for total and HM

``` r
i2 <- trends[, complete.cases(Jbetatrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
i3 <- trends[, complete.cases(Horntrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# full models
if(file.exists('temp/modTfullJbeta.rds')){
  modTfullJbeta <- readRDS('temp/modTfullJbeta.rds')
} else {
  modTfullJbeta <- lme(Jbetatrend ~ temptrend_abs.sc*REALM + 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave.sc +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*lifespan.sc + 
                         temptrend_abs.sc*consumerfrac.sc +
                         temptrend_abs.sc*endothermfrac.sc +
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*veg.sc +
                         temptrend_abs.sc*human_bowler.sc:REALM2,
                       random = randef, weights = varef, data = trends[i2,], method = 'REML')
  saveRDS(modTfullJbeta, file = 'temp/modTfullJbeta.rds')
}

if(file.exists('temp/modTfullHorn.rds')){
  modTfullHorn <- readRDS('temp/modTfullHorn.rds')
} else {
  modTfullHorn <- lme(Horntrend ~ temptrend_abs.sc*REALM + 
                        temptrend_abs.sc*tsign +
                        temptrend_abs.sc*tempave.sc +
                        temptrend_abs.sc*tempave_metab.sc + 
                        temptrend_abs.sc*seas.sc + 
                        temptrend_abs.sc*microclim.sc + 
                        temptrend_abs.sc*mass.sc + 
                        temptrend_abs.sc*speed.sc + 
                        temptrend_abs.sc*lifespan.sc + 
                        temptrend_abs.sc*consumerfrac.sc +
                        temptrend_abs.sc*endothermfrac.sc +
                        temptrend_abs.sc*nspp.sc +
                        temptrend_abs.sc*thermal_bias.sc:tsign +
                        temptrend_abs.sc*npp.sc +
                        temptrend_abs.sc*veg.sc +
                        temptrend_abs.sc*human_bowler.sc:REALM2,
                      random = randef, weights = varef, data = trends[i3,], method = 'REML')
  saveRDS(modTfullHorn, file = 'temp/modTfullHorn.rds')
}

summary(modTfullJbeta)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -139261.5 -138862.3 69676.77
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06255686 (Intr)
    ## temptrend_abs.sc 0.04777874 -0.173
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 7.124273e-05 0.3160086
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.352512 
    ## Fixed effects: Jbetatrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave.sc + temptrend_abs.sc * tempave_metab.sc +      temptrend_abs.sc * seas.sc + temptrend_abs.sc * microclim.sc +      temptrend_abs.sc * mass.sc + temptrend_abs.sc * speed.sc +      temptrend_abs.sc * lifespan.sc + temptrend_abs.sc * consumerfrac.sc +      temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.01923879 0.017457458 43218   1.102038  0.2705
    ## temptrend_abs.sc                                  0.08506960 0.022613395 43218   3.761912  0.0002
    ## REALMMarine                                       0.06411882 0.018544623   230   3.457543  0.0006
    ## REALMTerrestrial                                  0.03675475 0.020048351   230   1.833305  0.0680
    ## tsign1                                           -0.01258368 0.000516531 43218 -24.361901  0.0000
    ## tempave.sc                                       -0.00783124 0.000735146 43218 -10.652636  0.0000
    ## tempave_metab.sc                                  0.01457147 0.001659219 43218   8.782126  0.0000
    ## seas.sc                                          -0.00045585 0.000360506 43218  -1.264458  0.2061
    ## microclim.sc                                      0.00092771 0.000175205 43218   5.295006  0.0000
    ## mass.sc                                           0.00036376 0.000706034 43218   0.515221  0.6064
    ## speed.sc                                         -0.00189170 0.000526558 43218  -3.592572  0.0003
    ## lifespan.sc                                      -0.00218132 0.001354324 43218  -1.610633  0.1073
    ## consumerfrac.sc                                   0.00260054 0.001122633 43218   2.316464  0.0205
    ## endothermfrac.sc                                 -0.01827918 0.004359351 43218  -4.193096  0.0000
    ## nspp.sc                                          -0.00821385 0.000386221 43218 -21.267209  0.0000
    ## npp.sc                                            0.00075856 0.000284814 43218   2.663362  0.0077
    ## veg.sc                                            0.00165940 0.000365672 43218   4.537951  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.02928580 0.023324742 43218  -1.255568  0.2093
    ## temptrend_abs.sc:REALMTerrestrial                -0.00913657 0.026985327 43218  -0.338575  0.7349
    ## temptrend_abs.sc:tsign1                           0.01026889 0.001026980 43218   9.999116  0.0000
    ## temptrend_abs.sc:tempave.sc                       0.01438288 0.002002988 43218   7.180714  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                 0.00185922 0.004463909 43218   0.416501  0.6770
    ## temptrend_abs.sc:seas.sc                         -0.00377759 0.000955138 43218  -3.955019  0.0001
    ## temptrend_abs.sc:microclim.sc                    -0.00513173 0.000584645 43218  -8.777522  0.0000
    ## temptrend_abs.sc:mass.sc                         -0.00469213 0.001692060 43218  -2.773029  0.0056
    ## temptrend_abs.sc:speed.sc                         0.00339794 0.000872203 43218   3.895810  0.0001
    ## temptrend_abs.sc:lifespan.sc                      0.00440578 0.003496135 43218   1.260185  0.2076
    ## temptrend_abs.sc:consumerfrac.sc                 -0.00069175 0.001597973 43218  -0.432893  0.6651
    ## temptrend_abs.sc:endothermfrac.sc                 0.00524609 0.006949016 43218   0.754939  0.4503
    ## temptrend_abs.sc:nspp.sc                         -0.01505676 0.000893453 43218 -16.852330  0.0000
    ## tsign-1:thermal_bias.sc                          -0.00088046 0.000562857 43218  -1.564274  0.1178
    ## tsign1:thermal_bias.sc                           -0.00262409 0.000461084 43218  -5.691128  0.0000
    ## temptrend_abs.sc:npp.sc                          -0.00140981 0.000925201 43218  -1.523784  0.1276
    ## temptrend_abs.sc:veg.sc                          -0.00346921 0.001132709 43218  -3.062754  0.0022
    ## human_bowler.sc:REALM2TerrFresh                   0.00024360 0.000298867 43218   0.815095  0.4150
    ## human_bowler.sc:REALM2Marine                      0.00031110 0.000189206 43218   1.644260  0.1001
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.00900567 0.001157565 43218   7.779834  0.0000
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00553350 0.001128072 43218   4.905267  0.0000
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.00160427 0.000754126 43218  -2.127330  0.0334
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.00168345 0.000783751 43218  -2.147945  0.0317
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.428                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.929  0.400                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.875  0.376  0.815                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.026  0.017  0.006  0.010                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.003  0.001 -0.005 -0.007 -0.065                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.041 -0.022 -0.032 -0.031 -0.032 -0.441                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.030  0.037  0.040 -0.013 -0.104  0.144 -0.063                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.008  0.019  0.014 -0.006  0.023 -0.070 -0.005  0.119                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                           0.010  0.006  0.000 -0.002  0.020  0.008 -0.481  0.000  0.053                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.015 -0.004 -0.026 -0.016  0.035  0.003  0.136 -0.052 -0.025 -0.131                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.030 -0.023 -0.023 -0.035 -0.031 -0.043  0.685  0.045 -0.003 -0.795  0.258                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.004  0.004  0.026  0.248  0.008 -0.040 -0.049 -0.041  0.004  0.034 -0.108 -0.100                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                  0.135 -0.057 -0.069 -0.270 -0.002  0.116 -0.144  0.022  0.016  0.009  0.079  0.024 -0.340                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                          -0.008 -0.011 -0.007 -0.018 -0.040  0.000 -0.030 -0.059 -0.115 -0.200 -0.071  0.132 -0.001  0.064                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.004  0.004 -0.005  0.003  0.015 -0.248  0.216 -0.089 -0.087 -0.050 -0.011  0.063 -0.021 -0.042 -0.182                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.075  0.121  0.074 -0.004 -0.016  0.134 -0.093  0.119 -0.023  0.036 -0.058 -0.021 -0.058 -0.009  0.005 -0.269                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.412 -0.960 -0.413 -0.357 -0.004  0.003  0.015 -0.047 -0.021 -0.011  0.012  0.015 -0.012  0.036  0.019  0.008 -0.131                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.361 -0.856 -0.334 -0.432 -0.005 -0.011  0.023  0.027  0.010 -0.009  0.008  0.025 -0.093  0.136  0.030 -0.015  0.006  0.812                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.013 -0.037 -0.005 -0.004 -0.484  0.048  0.002  0.012 -0.042 -0.016 -0.011  0.016 -0.002  0.008  0.025 -0.049 -0.005  0.013      0.007                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.001 -0.008 -0.002 -0.002  0.049 -0.625  0.313 -0.050  0.091  0.033  0.002  0.013 -0.014 -0.099 -0.013  0.155 -0.134  0.000      0.017     -0.076                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.011  0.029  0.009  0.014 -0.009  0.329 -0.620 -0.013 -0.039  0.248 -0.017 -0.375  0.047  0.116 -0.009 -0.142  0.051 -0.020     -0.037      0.015 -0.397                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.016 -0.052 -0.019  0.021  0.042  0.010 -0.017 -0.648 -0.084  0.010  0.005 -0.019  0.039  0.014  0.041  0.059 -0.123  0.069     -0.053     -0.017  0.051            0.047                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.006 -0.028 -0.008  0.006 -0.046  0.096 -0.022 -0.049 -0.692 -0.027 -0.001  0.010  0.002  0.005  0.046  0.100  0.028  0.031     -0.021      0.047 -0.128            0.116   0.136                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.006 -0.012 -0.007 -0.007 -0.016  0.035  0.269  0.008 -0.016 -0.557  0.020  0.471 -0.007  0.013  0.130  0.004 -0.003  0.017      0.019      0.026 -0.084           -0.467  -0.049            -0.012                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.002  0.006  0.009  0.002 -0.018  0.019 -0.067  0.027  0.006  0.047 -0.563 -0.109  0.065 -0.049  0.059 -0.002  0.051 -0.020     -0.013     -0.003 -0.032            0.060  -0.037             0.000            -0.022                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.012  0.035  0.010  0.014  0.023  0.016 -0.388 -0.014  0.005  0.439 -0.046 -0.567  0.055 -0.013 -0.111 -0.018 -0.016 -0.019     -0.040     -0.021 -0.019            0.668   0.053             0.014            -0.815             0.086                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.000 -0.065 -0.009 -0.088 -0.004 -0.008  0.069  0.027  0.000 -0.024  0.076  0.087 -0.494  0.131  0.001  0.021  0.023  0.077      0.280      0.004  0.038           -0.131  -0.027            -0.002             0.045            -0.149            -0.144                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                -0.046  0.155  0.030  0.126  0.008 -0.179  0.216  0.008  0.019  0.012 -0.067  0.002  0.103 -0.394 -0.053  0.068 -0.022 -0.102     -0.355     -0.020  0.210           -0.343  -0.025            -0.057            -0.010             0.102            -0.032           -0.311                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.009  0.023  0.010  0.020  0.023 -0.024  0.025  0.033  0.048  0.131  0.049 -0.114  0.001 -0.056 -0.545  0.086 -0.033 -0.033     -0.054     -0.036 -0.006            0.024  -0.043            -0.022            -0.204            -0.092             0.192            0.000            0.091                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.015 -0.006 -0.017 -0.004 -0.167  0.471 -0.092 -0.109 -0.077  0.019 -0.009 -0.061  0.000  0.013  0.004 -0.081  0.038  0.012     -0.008      0.017 -0.277            0.074   0.088             0.050             0.000             0.017             0.026           -0.013           -0.017            -0.010                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.023 -0.013 -0.031 -0.010 -0.032  0.710 -0.149 -0.242 -0.170  0.036  0.018 -0.061 -0.007  0.026 -0.038 -0.125  0.048  0.017     -0.003      0.069 -0.356            0.097   0.120             0.068             0.002            -0.004             0.020           -0.009           -0.040            -0.010             0.569                                                                                               
    ## temptrend_abs.sc:npp.sc                           0.008  0.008 -0.006 -0.013 -0.010  0.081 -0.073  0.034  0.110  0.013 -0.003 -0.017  0.003  0.007  0.065 -0.654  0.183 -0.025      0.029      0.035 -0.185            0.151  -0.190            -0.250             0.019             0.012             0.019           -0.030           -0.068            -0.091             0.026  0.056                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.053 -0.155 -0.054  0.006  0.008 -0.088  0.042 -0.114  0.026 -0.010  0.032 -0.001  0.018 -0.003 -0.026  0.190 -0.761  0.170     -0.015      0.023  0.160           -0.048   0.205             0.011            -0.016            -0.060             0.037           -0.013            0.025             0.030            -0.030 -0.047 -0.316                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.026  0.050  0.017 -0.011  0.038 -0.283  0.102 -0.059  0.125  0.045  0.051 -0.057 -0.071 -0.077  0.016 -0.064  0.191 -0.054      0.015     -0.032  0.162           -0.091  -0.149            -0.037            -0.021            -0.025             0.017            0.025            0.064             0.036            -0.081 -0.098  0.053            -0.263                                                               
    ## human_bowler.sc:REALM2Marine                      0.003 -0.010 -0.002  0.000 -0.008  0.033 -0.010 -0.050 -0.074 -0.019 -0.050  0.013  0.007  0.000  0.033 -0.212  0.052  0.009      0.003      0.015  0.000           -0.018   0.087             0.039             0.022             0.036            -0.013            0.002            0.004             0.013             0.057  0.064  0.131            -0.028           -0.001                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.004  0.002  0.004 -0.001  0.015 -0.375  0.110  0.085  0.085 -0.004  0.005  0.036 -0.013 -0.023 -0.011  0.064 -0.071 -0.018      0.014      0.103  0.686           -0.092  -0.110            -0.104            -0.003            -0.007            -0.054            0.026            0.014            -0.034            -0.471 -0.379 -0.082             0.084            0.058      -0.046                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.006  0.008  0.005  0.000  0.066 -0.410  0.106  0.093  0.115  0.002 -0.005  0.030 -0.010 -0.020  0.002  0.073 -0.082 -0.018      0.011     -0.125  0.727           -0.071  -0.090            -0.136            -0.016             0.012            -0.045            0.022            0.008            -0.050            -0.342 -0.522 -0.096             0.102            0.052      -0.052       0.794                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.026 -0.086 -0.025  0.010 -0.019  0.192 -0.093 -0.110 -0.066 -0.028 -0.015  0.031  0.022  0.048  0.019  0.033 -0.312  0.093     -0.024      0.054 -0.237            0.123   0.228             0.108             0.042             0.029            -0.035           -0.027           -0.088            -0.036             0.060  0.070 -0.105             0.458           -0.684      -0.001      -0.104 -0.106               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.004  0.015  0.004  0.001  0.004 -0.017 -0.003  0.077  0.026  0.013  0.037 -0.012  0.002  0.001 -0.011  0.100 -0.006 -0.015      0.001     -0.007 -0.038            0.062  -0.180            -0.084            -0.018            -0.028             0.005           -0.007           -0.029            -0.043            -0.043 -0.057 -0.156             0.023            0.013      -0.613       0.048  0.047 -0.001        
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -7.4746814 -0.2440563  0.1302303  0.6083779  6.7313862 
    ## 
    ## Number of Observations: 43488
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    233                  43488

``` r
summary(modTfullHorn)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -125431.1 -125032.9 62761.55
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06117307 (Intr)
    ## temptrend_abs.sc 0.04027972 -0.08 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.004102246 0.2958833
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##    power 
    ## -1.22294 
    ## Fixed effects: Horntrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave.sc + temptrend_abs.sc * tempave_metab.sc +      temptrend_abs.sc * seas.sc + temptrend_abs.sc * microclim.sc +      temptrend_abs.sc * mass.sc + temptrend_abs.sc * speed.sc +      temptrend_abs.sc * lifespan.sc + temptrend_abs.sc * consumerfrac.sc +      temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.01107922 0.017093964 42251   0.648136  0.5169
    ## temptrend_abs.sc                                  0.08040481 0.020814432 42251   3.862936  0.0001
    ## REALMMarine                                       0.05795446 0.018364960   199   3.155708  0.0018
    ## REALMTerrestrial                                  0.04354349 0.019976795   199   2.179703  0.0305
    ## tsign1                                           -0.01017320 0.000652951 42251 -15.580346  0.0000
    ## tempave.sc                                       -0.00342630 0.001020285 42251  -3.358185  0.0008
    ## tempave_metab.sc                                  0.00516455 0.002240704 42251   2.304880  0.0212
    ## seas.sc                                          -0.00261568 0.000538621 42251  -4.856250  0.0000
    ## microclim.sc                                      0.00073103 0.000280597 42251   2.605278  0.0092
    ## mass.sc                                          -0.00018185 0.000909808 42251  -0.199876  0.8416
    ## speed.sc                                         -0.00013320 0.000711212 42251  -0.187291  0.8514
    ## lifespan.sc                                      -0.00374602 0.001774365 42251  -2.111191  0.0348
    ## consumerfrac.sc                                   0.00214128 0.001298915 42251   1.648515  0.0993
    ## endothermfrac.sc                                 -0.01522740 0.004924544 42251  -3.092144  0.0020
    ## nspp.sc                                          -0.00849495 0.000508028 42251 -16.721423  0.0000
    ## npp.sc                                           -0.00007975 0.000434614 42251  -0.183505  0.8544
    ## veg.sc                                            0.00069277 0.000554676 42251   1.248970  0.2117
    ## temptrend_abs.sc:REALMMarine                     -0.02246053 0.021572719 42251  -1.041154  0.2978
    ## temptrend_abs.sc:REALMTerrestrial                -0.05927379 0.025202434 42251  -2.351907  0.0187
    ## temptrend_abs.sc:tsign1                           0.00734926 0.001147683 42251   6.403564  0.0000
    ## temptrend_abs.sc:tempave.sc                       0.01344607 0.002286524 42251   5.880572  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                 0.00459335 0.005114635 42251   0.898080  0.3691
    ## temptrend_abs.sc:seas.sc                          0.00501253 0.001133346 42251   4.422771  0.0000
    ## temptrend_abs.sc:microclim.sc                    -0.00117423 0.000725047 42251  -1.619517  0.1053
    ## temptrend_abs.sc:mass.sc                         -0.00096131 0.001925826 42251  -0.499166  0.6177
    ## temptrend_abs.sc:speed.sc                         0.00096148 0.000979263 42251   0.981840  0.3262
    ## temptrend_abs.sc:lifespan.sc                      0.00105477 0.003965406 42251   0.265994  0.7902
    ## temptrend_abs.sc:consumerfrac.sc                 -0.00045406 0.001739325 42251  -0.261054  0.7941
    ## temptrend_abs.sc:endothermfrac.sc                 0.00663560 0.006986511 42251   0.949773  0.3422
    ## temptrend_abs.sc:nspp.sc                         -0.00964784 0.001018570 42251  -9.471943  0.0000
    ## tsign-1:thermal_bias.sc                          -0.00114227 0.000709046 42251  -1.610994  0.1072
    ## tsign1:thermal_bias.sc                           -0.00018020 0.000607772 42251  -0.296497  0.7669
    ## temptrend_abs.sc:npp.sc                           0.00242345 0.001102079 42251   2.198982  0.0279
    ## temptrend_abs.sc:veg.sc                          -0.00076024 0.001318386 42251  -0.576641  0.5642
    ## human_bowler.sc:REALM2TerrFresh                   0.00032460 0.000466298 42251   0.696122  0.4864
    ## human_bowler.sc:REALM2Marine                      0.00152861 0.000383633 42251   3.984576  0.0001
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.00550039 0.001292207 42251   4.256588  0.0000
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00217968 0.001260050 42251   1.729840  0.0837
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.00074368 0.000875438 42251   0.849495  0.3956
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.00109468 0.001011335 42251  -1.082414  0.2791
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.377                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.917  0.349                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.858  0.323  0.781                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.031  0.020  0.005  0.012                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.007  0.004  0.003 -0.006 -0.070                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.056 -0.025 -0.046 -0.039 -0.030 -0.495                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.046  0.038  0.062 -0.014 -0.113  0.197 -0.099                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.019  0.022  0.026 -0.006  0.006 -0.041 -0.025  0.147                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                           0.009  0.011  0.002  0.005  0.021  0.010 -0.476 -0.009  0.021                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.020 -0.005 -0.038 -0.023  0.051 -0.009  0.132 -0.079 -0.054 -0.134                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.038 -0.024 -0.030 -0.043 -0.028 -0.047  0.670  0.049  0.019 -0.805  0.248                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.007  0.007  0.008  0.270  0.009 -0.019 -0.057 -0.030  0.005  0.041 -0.129 -0.106                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                  0.146 -0.049 -0.071 -0.302  0.002  0.170 -0.197  0.044  0.022 -0.003  0.096  0.018 -0.339                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                          -0.011 -0.016 -0.007 -0.022 -0.035  0.016 -0.032 -0.013 -0.078 -0.212 -0.066  0.144  0.005  0.078                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.015  0.004 -0.018  0.008  0.013 -0.298  0.248 -0.210 -0.226 -0.041 -0.005  0.054 -0.010 -0.078 -0.194                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.111  0.125  0.113 -0.002 -0.018  0.145 -0.106  0.136  0.022  0.030 -0.053 -0.019 -0.039  0.013  0.017 -0.277                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.361 -0.957 -0.358 -0.303 -0.005 -0.005  0.019 -0.048 -0.025 -0.014  0.015  0.016 -0.011  0.033  0.024  0.011 -0.137                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.311 -0.839 -0.283 -0.384 -0.006 -0.008  0.025  0.041  0.012 -0.019  0.009  0.030 -0.085  0.133  0.041 -0.024  0.012  0.788                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.016 -0.044 -0.005 -0.005 -0.488  0.041  0.006  0.015 -0.036 -0.014 -0.015  0.012 -0.003  0.005  0.025 -0.047  0.003  0.014      0.009                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.004 -0.013 -0.007 -0.002  0.048 -0.572  0.306 -0.055  0.082  0.029  0.004  0.013 -0.014 -0.111 -0.015  0.145 -0.125  0.008      0.020     -0.068                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.015  0.034  0.013  0.017 -0.007  0.312 -0.591 -0.008 -0.022  0.253 -0.025 -0.365  0.056  0.135 -0.007 -0.144  0.052 -0.022     -0.053      0.012 -0.414                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.018 -0.058 -0.025  0.030  0.044 -0.026  0.000 -0.607 -0.100  0.010  0.017 -0.022  0.028 -0.001  0.015  0.109 -0.130  0.076     -0.080     -0.019  0.026            0.065                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.009 -0.027 -0.012  0.007 -0.039  0.068 -0.010 -0.075 -0.660 -0.008  0.013 -0.005  0.000 -0.002  0.033  0.155  0.023  0.030     -0.029      0.050 -0.152            0.113   0.136                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.008 -0.019 -0.009 -0.012 -0.013  0.026  0.263  0.005 -0.007 -0.554  0.029  0.467 -0.017  0.017  0.129  0.007  0.000  0.022      0.034      0.025 -0.090           -0.465  -0.043            -0.010                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.005  0.007  0.016  0.005 -0.028  0.019 -0.070  0.036  0.019  0.055 -0.578 -0.119  0.080 -0.062  0.053  0.003  0.046 -0.024     -0.014     -0.001 -0.037            0.063  -0.037            -0.002            -0.023                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.013  0.039  0.012  0.017  0.020  0.017 -0.368 -0.015 -0.003  0.446 -0.060 -0.555  0.070 -0.012 -0.113 -0.018 -0.016 -0.020     -0.050     -0.019 -0.013            0.658   0.052             0.010            -0.823             0.088                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.001 -0.069  0.001 -0.082 -0.005 -0.016  0.077  0.028  0.004 -0.029  0.085  0.097 -0.489  0.113 -0.001  0.004  0.014  0.077      0.262      0.008  0.057           -0.156  -0.035            -0.010             0.055            -0.159            -0.171                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                -0.040  0.155  0.028  0.128  0.005 -0.195  0.236  0.002  0.012  0.021 -0.078 -0.001  0.095 -0.365 -0.061  0.082 -0.031 -0.107     -0.371     -0.018  0.252           -0.388  -0.039            -0.064            -0.025             0.115            -0.026           -0.273                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.011  0.033  0.011  0.026  0.016 -0.026  0.021  0.002  0.032  0.133  0.042 -0.118  0.004 -0.062 -0.527  0.087 -0.037 -0.041     -0.072     -0.032 -0.005            0.033  -0.026            -0.011            -0.209            -0.086             0.205           -0.006            0.100                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.021 -0.007 -0.022 -0.005 -0.165  0.465 -0.102 -0.108 -0.089  0.028 -0.002 -0.070  0.011  0.019 -0.007 -0.096  0.040  0.011     -0.010      0.024 -0.262            0.074   0.087             0.053            -0.004             0.013             0.029           -0.020           -0.016             0.001                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.034 -0.015 -0.042 -0.013 -0.029  0.650 -0.130 -0.230 -0.201  0.047  0.028 -0.072  0.008  0.028 -0.041 -0.124  0.039  0.015     -0.002      0.071 -0.313            0.084   0.105             0.081            -0.007            -0.010             0.026           -0.016           -0.034             0.002             0.588                                                                                               
    ## temptrend_abs.sc:npp.sc                           0.009  0.002 -0.005 -0.018 -0.012  0.110 -0.094  0.080  0.155  0.008  0.004 -0.014  0.001  0.028  0.075 -0.614  0.164 -0.023      0.040      0.038 -0.177            0.159  -0.203            -0.280             0.018             0.007             0.017           -0.012           -0.080            -0.090             0.034  0.057                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.066 -0.189 -0.068  0.010  0.009 -0.099  0.048 -0.128  0.021 -0.004  0.030 -0.007  0.013 -0.013 -0.032  0.182 -0.640  0.209     -0.027      0.021  0.177           -0.051   0.232            -0.022            -0.023            -0.062             0.045           -0.009            0.037             0.038            -0.030 -0.043 -0.304                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.025  0.044  0.019 -0.013  0.035 -0.255  0.114 -0.061  0.151  0.040  0.044 -0.044 -0.057 -0.070  0.009 -0.081  0.131 -0.048      0.018     -0.021  0.134           -0.077  -0.122            -0.034            -0.016            -0.022             0.012            0.019            0.059             0.033            -0.071 -0.073  0.043            -0.205                                                               
    ## human_bowler.sc:REALM2Marine                      0.011 -0.017 -0.014 -0.003  0.007 -0.023  0.024 -0.146 -0.101 -0.018  0.008  0.000 -0.001  0.003 -0.033 -0.137  0.023  0.023      0.006      0.015  0.026           -0.028   0.138             0.026             0.019             0.006            -0.002            0.011           -0.002             0.047             0.085  0.079  0.078             0.004            0.011                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.004 -0.001  0.003 -0.002  0.020 -0.328  0.098  0.087  0.088 -0.007  0.002  0.037 -0.016 -0.019  0.002  0.047 -0.064 -0.013      0.016      0.101  0.682           -0.102  -0.123            -0.123            -0.006            -0.010            -0.052            0.045            0.021            -0.047            -0.466 -0.354 -0.069             0.093            0.046      -0.060                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.008  0.006  0.006  0.001  0.067 -0.355  0.088  0.091  0.117 -0.004 -0.009  0.033 -0.014 -0.015  0.012  0.058 -0.072 -0.011      0.010     -0.134  0.710           -0.075  -0.086            -0.153            -0.018             0.011            -0.044            0.039            0.010            -0.060            -0.338 -0.496 -0.083             0.112            0.033      -0.066       0.795                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.028 -0.104 -0.029  0.015 -0.017  0.158 -0.089 -0.121 -0.054 -0.026 -0.013  0.024  0.018  0.046  0.018  0.042 -0.263  0.114     -0.041      0.050 -0.229            0.132   0.283             0.079             0.038             0.026            -0.027           -0.026           -0.099            -0.037             0.052  0.050 -0.100             0.480           -0.562       0.005      -0.103 -0.092               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.009  0.025  0.012  0.003  0.004  0.014 -0.022  0.121  0.036  0.015  0.011 -0.007  0.005 -0.001  0.023  0.072  0.008 -0.031     -0.004     -0.010 -0.051            0.066  -0.220            -0.063            -0.021            -0.018             0.001           -0.017           -0.019            -0.066            -0.062 -0.064 -0.135            -0.002            0.006      -0.658       0.067  0.065 -0.013        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.76846562 -0.38338664  0.03973782  0.54230508  6.19339158 
    ## 
    ## Number of Observations: 42490
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    202                  42490

#### plot coefs for full Total and HM models

``` r
coefs2 <- summary(modTfullJbeta)$tTable
coefs3 <- summary(modTfullHorn)$tTable
varstoplot <- unique(c(rownames(coefs2), rownames(coefs3)))

rows1 <- which(!grepl('Intercept', varstoplot) | grepl(':', varstoplot)) # vars to plot in first graph
rows1_2 <- which(rownames(coefs2) %in% varstoplot[rows1]) # rows in coefs2
rows1_3 <- which(rownames(coefs3) %in% varstoplot[rows1]) # rows in coefs3
xlims1 <- range(c(coefs2[rows1_2,1] - coefs2[rows1_2,2], 
                  coefs2[rows1_2,1] + coefs2[rows1_2,2], 
                  coefs3[rows1_3,1] - coefs3[rows1_3,2], 
                  coefs3[rows1_3,1] + coefs3[rows1_3,2]))

cols <- c('black', 'grey') # for Jbeta and Horn models, respectively
offs1 <- 0.1 # offset vertically for the two models
offs2 <- 0.01 # offset vertically for the two models (plot 2)

par(las = 1, mai = c(0.5, 3, 0.1, 0.1))

plot(0,0, col = 'white', xlim=xlims1, ylim = c(1,length(rows1)), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(rows1):1, labels = varstoplot[rows1], cex.axis = 0.7)
abline(v = 0, col = 'grey')
for(i in 1:length(rows1)){
  if(varstoplot[rows1[i]] %in% rownames(coefs2)){
    x = coefs2[rownames(coefs2) == varstoplot[rows1[i]], 1]
    se = coefs2[rownames(coefs2) == varstoplot[rows1[i]], 2]
    points(x, length(rows1) + 1 - i + offs1, pch = 16, col = cols[1])
    lines(x = c(x-se, x+se), y = c(length(rows1) + 1 - i + offs1, length(rows1) + 1 - i + offs1), col = cols[1])
  }
  if(varstoplot[rows1[i]] %in% rownames(coefs3)){
    x = coefs3[rownames(coefs3) == varstoplot[rows1[i]], 1]
    se = coefs3[rownames(coefs3) == varstoplot[rows1[i]], 2]
    points(x, length(rows1) + 1 - i - offs1, pch = 16, col = cols[2])
    lines(x = c(x-se, x+se), y = c(length(rows1) + 1 - i - offs1, length(rows1) + 1 - i - offs1), col = cols[2])
  }
}
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20total%20and%20HM%20mod%20coefs-1.png)<!-- -->

Black is for Jaccard total turnover (pres/abs), grey is for
Morisita-Horn turnover (considers abundance)

#### Delete terms from the full Jbeta and Horn models, plot dAICs

``` r
if(file.exists('output/aics_from_full.csv')){
  aicsfromfull <- read.csv('output/aics_from_full.csv')
  
  if(all(c('dAIC_Jbeta', 'dAIC_Horn') %in% colnames(aicsfromfull))){
    runJbetaHorn <- FALSE
  } else {
    runJbetaHorn <- TRUE
  }
} else {
  runJbetaHorn <- TRUE
}

if(runJbetaHorn){
  i <- trends[, complete.cases(Jbetatrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                               temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                               consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
  i2 <- trends[, complete.cases(Horntrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                               temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                               consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
  
  randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
  varef <- varPower(-0.5, ~nyrBT)
  
  terms <- c('temptrend_abs.sc*REALM', 
             'temptrend_abs.sc*tsign',
             'temptrend_abs.sc*tempave.sc',
             'temptrend_abs.sc*tempave_metab.sc',
             'temptrend_abs.sc*seas.sc',
             'temptrend_abs.sc*microclim.sc',
             'temptrend_abs.sc*mass.sc',
             'temptrend_abs.sc*speed.sc', 
             'temptrend_abs.sc*lifespan.sc', 
             'temptrend_abs.sc*consumerfrac.sc',
             'temptrend_abs.sc*endothermfrac.sc',
             'temptrend_abs.sc*nspp.sc',
             'temptrend_abs.sc*thermal_bias.sc:tsign',
             'temptrend_abs.sc*npp.sc',
             'temptrend_abs.sc*veg.sc',
             'temptrend_abs.sc*human_bowler.sc:REALM2')
  
  modTJbetadrops <- vector('list', length(terms)+2)
  modTHorndrops <- vector('list', length(terms)+2)

  names(modTJbetadrops) <- c('full', '-temptrend_abs.sc', paste0('-', terms))
  names(modTHorndrops) <- c('full', '-temptrend_abs.sc', paste0('-', terms))
  
  # fit full model with ML for model comparison
  modTJbetadrops[[1]] <- lme(formula(paste0('Jbetatrend ~ ', paste(terms, collapse = ' + '))),
                        random = randef, weights = varef, data = trends[i,], method = 'ML')
  modTHorndrops[[1]] <- lme(formula(paste0('Horntrend ~ ', paste(terms, collapse = ' + '))),
                        random = randef, weights = varef, data = trends[i2,], method = 'ML')
  
  # w/out temptrend
  modTJbetadrops[[2]] <- lme(formula(paste0('Jbetatrend ~ ', paste(gsub('temptrend_abs.sc\\*', '', terms), collapse = ' + '))),
                        random = list(STUDY_ID = ~ 1, rarefyID = ~1), weights = varef, data = trends[i,], method = 'ML')
  modTHorndrops[[2]] <- lme(formula(paste0('Horntrend ~ ', paste(gsub('temptrend_abs.sc\\*', '', terms), collapse = ' + '))),
                        random = list(STUDY_ID = ~ 1, rarefyID = ~1), weights = varef, data = trends[i2,], method = 'ML')
  
  for(j in 1:length(terms)){
    print(j)
    modTJbetadrops[[j+2]] <- lme(formula(paste0('Jbetatrend ~ ', paste(terms[-j], collapse = ' + '))),
                                 random = randef, weights = varef, data = trends[i,], method = 'ML', 
                                 control = lmeContrl(returnObject = TRUE)) # return fitted object even if a convergence error
    modTHorndrops[[j+2]] <- lme(formula(paste0('Horntrend ~ ', paste(terms[-j], collapse = ' + '))),
                            random = randef, weights = varef, data = trends[i2,], method = 'ML')
  }
  
  aicsJbeta <- sapply(modTJbetadrops, AIC)
  aicsHorn <- sapply(modTHorndrops, AIC)

  if(exists('aicsfromfull')){
    aicsfromfull$dAIC_Jbeta <- aicsJbeta - aicsJbeta[1]
    aicsfromfull$dAIC_Horn <- aicsHorn - aicsHorn[1]
  } else {
    aicsfromfull <- data.frame(mod = names(aics), 
                               dAIC_Jbeta = aicsJbeta - aicsJbeta[1],
                               dAIC_Horn <- aicsHorn - aicsHorn[1])
  }
  
  write.csv(aicsfromfull, file = 'output/aics_from_full.csv', row.names = FALSE)
}

aicsfromfull
```

    ##                                         mod      dAIC_Jtu   dAIC_Jbeta    dAIC_Horn
    ## 1                                      full    0.00000000     0.000000    0.0000000
    ## 2                         -temptrend_abs.sc 7392.04492882 12840.286245 7605.6472326
    ## 3                   -temptrend_abs.sc*REALM   14.87466466     7.098437   10.4577563
    ## 4                   -temptrend_abs.sc*tsign  307.26456844   589.870020  240.0240794
    ## 5              -temptrend_abs.sc*tempave.sc   45.08557421   109.482738   30.5043643
    ## 6        -temptrend_abs.sc*tempave_metab.sc   14.83704694   128.758735    9.1313302
    ## 7                 -temptrend_abs.sc*seas.sc   11.03831929    36.657551   22.4356323
    ## 8            -temptrend_abs.sc*microclim.sc   12.88175453    73.883018    2.6523187
    ## 9                 -temptrend_abs.sc*mass.sc    7.95995370     5.221911   -3.4124476
    ## 10               -temptrend_abs.sc*speed.sc    5.05076864    14.138661   -2.7615998
    ## 11            -temptrend_abs.sc*lifespan.sc    4.07127373    -1.204624    1.6582203
    ## 12        -temptrend_abs.sc*consumerfrac.sc   -0.07879609     1.815516   -0.9988319
    ## 13       -temptrend_abs.sc*endothermfrac.sc   13.73066680    14.590190    5.7036860
    ## 14                -temptrend_abs.sc*nspp.sc  398.21647587  1568.523765  732.0112555
    ## 15  -temptrend_abs.sc*thermal_bias.sc:tsign   22.37488755   105.350725   17.6958524
    ## 16                 -temptrend_abs.sc*npp.sc    8.71935851     3.171746    2.9179283
    ## 17                 -temptrend_abs.sc*veg.sc    3.74525115    17.061611   -2.3119391
    ## 18 -temptrend_abs.sc*human_bowler.sc:REALM2   10.86674547     1.996388   14.3718703

``` r
# transform for a plot
aicsfromfulllong <- reshape(aicsfromfull, direction = 'long',
                            varying = c('dAIC_Jtu', 'dAIC_Jbeta', 'dAIC_Horn'),
                            v.names = 'dAIC',
                            idvar = 'mod',
                            timevar = 'type',
                            times = c('Jtu', 'Jbeta', 'Horn'))

trans = function(x) sign(x)*sqrt(abs(x))
aicsfromfulllong$dAIC_tr <- trans(aicsfromfulllong$dAIC)

# plot
xlims <- range(aicsfromfulllong$dAIC_tr)
xticks <- c(-10, 0, 10, 100, 1000, 10000)
par(mai = c(0.5, 3, 0.1, 0.1))
with(aicsfromfulllong[aicsfromfulllong$type == 'Jtu',], plot(dAIC_tr, nrow(aicsfromfull):1, 
                                                           col = 'light grey', xlim = xlims, yaxt = 'n', ylab = '', xaxt = 'n'))
with(aicsfromfulllong[aicsfromfulllong$type == 'Jbeta',], points(dAIC_tr, nrow(aicsfromfull):1 - 0.1, col = 'dark grey'))
with(aicsfromfulllong[aicsfromfulllong$type == 'Horn',], points(dAIC_tr, nrow(aicsfromfull):1 - 0.2, col = 'black'))
axis(2, at = nrow(aicsfromfull):1, labels = aicsfromfull$mod, las = 1, cex.axis = 0.7)
axis(1, at = trans(xticks), labels = xticks, cex.axis = 0.5)
abline(v = 0, lty =2, col = 'grey')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/term%20deletion%20from%20Jbeta%20and%20Horn%20models-1.png)<!-- -->
Light grey is for Jaccard turnover, dark grey is for Jaccard total,
black is for Morisita-Horn

### Plot interaction coefficients from all full models

``` r
# fig.width = 3, fig.height = 5, out.height=2.5, out.width=3, fig.retina =3 for macbook screen
# double that for external monitor
coefs <- as.data.table(summary(modTfull1)$tTable)
coefs2 <- as.data.table(summary(modTfullJbeta)$tTable)
coefs3 <- as.data.table(summary(modTfullHorn)$tTable)

coefs$mod <- 'Jtu'
coefs2$mod <- 'Jbeta'
coefs3$mod <- 'Horn'

coefs$var <- rownames(summary(modTfull1)$tTable)
coefs2$var <- rownames(summary(modTfullJbeta)$tTable)
coefs3$var <- rownames(summary(modTfullHorn)$tTable)

# extract temperature effects and bind model coefs together
cols <- c('var', 'Value', 'Std.Error', 'mod')

allcoefsfull <- rbind(coefs[grep('temptrend|REALM', var), ..cols], 
                      coefs2[grep('temptrend|REALM', var), ..cols], 
                      coefs3[grep('temptrend|REALM', var), ..cols])
allcoefsfull$var[allcoefsfull$var == 'temptrend_abs.sc'] <- 'temptrend_abs.sc:REALMFreshwater'

# add average temperature effect (across realms) to realm-specific temperature effects
meantempeffect <- allcoefsfull[var == 'temptrend_abs.sc:REALMFreshwater', mean(Value), by = mod]$V1
allcoefsfull[var == 'temptrend_abs.sc:REALMMarine', Value := Value + meantempeffect]
allcoefsfull[var == 'temptrend_abs.sc:REALMTerrestrial', Value := Value + meantempeffect]

# remove non-temperature effects
allcoefsfull <- allcoefsfull[grepl(':', allcoefsfull$var) & grepl('temptrend', allcoefsfull$var), ]

# add info for plotting
allcoefsfull$lCI <- allcoefsfull$Value - allcoefsfull$Std.Error # lower confidence interval
allcoefsfull$uCI <- allcoefsfull$Value + allcoefsfull$Std.Error
nvar <- nrow(allcoefsfull)/3
allcoefsfull$y <- 1:nvar + rep(c(0, 0.1, 0.2), c(nvar, nvar, nvar)) # y-values

# clean up some variable names for nicer plotting
allcoefsfull$varname <- gsub('temptrend_abs.sc:|temptrend.sc:', '', allcoefsfull$var)
allcoefsfull$varname <- gsub('REALM|REALM2', '', allcoefsfull$varname)
allcoefsfull$varname <- gsub('.sc', '', allcoefsfull$varname)
allcoefsfull$varname <- gsub('_bowler', '', allcoefsfull$varname)
allcoefsfull$varname <- gsub('tsign1', 'warming', allcoefsfull$varname)
allcoefsfull$varname <- gsub('tsign-1', 'cooling', allcoefsfull$varname)

# set x-axes
xlims1 <- c(-0.01, 0.11) # for realms
xlims2 <- c(-0.007, 0.017) # for traits
xlims3 <- c(-0.01, 0.02) # for environment
xlims4 <- c(-0.016, 0.01) # for community
xlims5 <- c(-0.005, 0.0025) # for human

ddg <- 0.5 # vertical dodge for each model

# choose which variables are in which graph
set1 <- c('Terrestrial', 'Marine', 'Freshwater')
set2 <- c('mass', 'speed', 'lifespan', 'consumerfrac', 'endothermfrac', 'tempave_metab')
set3 <- c('seas', 'microclim', 'tempave')
set4 <- c('npp', 'nspp', 'tsign-1:thermal_bias', 'tsign1:thermal_bias')
set5 <- c('human_bowler:TerrFresh', 'human_bowler:Marine')

p1 <- ggplot(subset(allcoefsfull, varname %in% set1), 
                    aes(varname, Value, group = mod, color = mod)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'light grey') +
  geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0, position = position_dodge(ddg)) + 
  geom_point(position = position_dodge(ddg)) + 
  labs(y = 'Temperature change effect', x = '', tag = 'A') +
  scale_color_grey() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position='none',
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) + 
  coord_flip(ylim = xlims1)

p2 <- ggplot(subset(allcoefsfull, varname %in% set2), 
                    aes(varname, Value, group = mod, color = mod)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'light grey') +
  geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0, position = position_dodge(ddg)) + 
  geom_point(position = position_dodge(ddg)) + 
  labs(y = 'Interaction with temperature change effect', x = '', tag = 'B') +
  scale_color_grey() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position='none',
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) + 
  coord_flip(ylim = xlims2)

p3 <- ggplot(subset(allcoefsfull, varname %in% set3), 
                    aes(varname, Value, group = mod, color = mod)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'light grey') +
  geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0, position = position_dodge(ddg)) + 
  geom_point(position = position_dodge(ddg)) + 
  labs(y = 'Interaction with temperature change effect', x = '', tag = 'C') +
  scale_color_grey() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position='none',
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) + 
  coord_flip(ylim = xlims3)

p4 <- ggplot(subset(allcoefsfull, varname %in% set4), 
                    aes(varname, Value, group = mod, color = mod)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'light grey') +
  geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0, position = position_dodge(ddg)) + 
  geom_point(position = position_dodge(ddg)) + 
  labs(y = 'Interaction with temperature change effect', x = '', tag = 'D') +
  scale_color_grey() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position='none',
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) + 
  coord_flip(ylim = xlims4)

p5 <- ggplot(subset(allcoefsfull, varname %in% set5), 
                    aes(varname, Value, group = mod, color = mod)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'light grey') +
  geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0, position = position_dodge(ddg)) + 
  geom_point(position = position_dodge(ddg)) + 
  labs(y = 'Interaction with temperature change effect', x = '', tag = 'E') +
  scale_color_grey() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position='none',
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) + 
  coord_flip(ylim = xlims5)

grid.arrange(p1, p2, p3, p4, p5, ncol = 2, layout_matrix = rbind(c(1,2), c(3,4), c(5, NA)))
```

<img src="turnover_vs_temperature_MEmodels_files/figure-gfm/plot coefs for all models-1.png" width="6" height="5" />

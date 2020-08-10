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
  - Explanatory variables considered for differences in rate of
    turnover:
      - Temperature trend over the time-frame of each time-series (CRU
        TS 4.03 on land and in freshwater, ERSST v5 in the ocean)
      - Seasonality as a metric of thermal sensitivity (Deutsch et
        al. 2008 PNAS). Standard deviation of monthly temperatures.
      - Microclimates calculated from WorldClim and BioOracle (Laura
        Antao)
      - Average temperature
      - Body mass, collated from databases and literature searches
      - Metabolic temperature, from average temperature if ectotherms
        (Dillon et al. 2010 Nature, Antao et al. 2020 Nat E\&E)
      - Mobility calculated from body mass and taxonomic group
        classifications of mobility mode (fly, run, swim, crawl,
        sessile). Fly/run/swim followed the allometric relationship in
        Hirt et al. 2017 Nat E\&E. Crawl set at 0.1 km/hr, sessile set
        to 0 km/hr. Then calculated averaged within each assemblage.
      - Generation time calculated from body mass and endotherm
        vs. ectotherm classifications, following McCoy & Gillooly 2008
        ELE. Averaged across species within each assemblage.
      - Consumer vs. producer classification by species
      - Endotherm vs. ectotherm classification by species
      - Species richness, calculated as the number of species in the
        assemblage
      - Net primary productivity (NPP) from the merged land/ocean
        product produced by the [Ocean
        Productivity](http://www.science.oregonstate.edu/ocean.productivity/)
        group at Oregon State using methods from Zhao et al. 2005 and
        Behrenfeld & Falkowski 1997.
      - Human impact calculated from Bowler et al. 2020 (also try data
        from Venter et al. 2016 and Halpern et al. 2008)
      - Thermal bias calculated from Species Temperature Indices (Mike
        Burrows)
      - Vegetation cover index, calculated from %tree cover and
        %non-tree veg cover (latter counted as 1/2), from vegetation
        continuous fields (Ruben Remelgado)
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
library(RColorBrewer)

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

With all covariates (Bowler for human)

``` r
# the cases we can compare
apply(trends[, .(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##         Jtutrend            REALM       tempave.sc tempave_metab.sc          seas.sc     microclim.sc     temptrend.sc          mass.sc         speed.sc      lifespan.sc  consumerfrac.sc endothermfrac.sc          nspp.sc  thermal_bias.sc           npp.sc           veg.sc  human_bowler.sc 
    ##            53013            53013            49916            49916            49916            51834            49916            52820            52734            51540            47534            53013            53013            49371            52863            52890            53013

``` r
i <- trends[, complete.cases(Jtutrend, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
cat('Overall # time-series: ', sum(i), '\n')
```

    ## Overall # time-series:  43493

``` r
cat('# studies: ', trends[i, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  235

``` r
cat('Data points: ', trends[i, sum(nyrBT)], '\n')
```

    ## Data points:  221956

``` r
trends[i, table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##         978       39738        2777

``` r
trends[i, table(taxa_mod)]
```

    ## taxa_mod
    ##           All    Amphibians       Benthos         Birds          Fish Invertebrates       Mammals         Plant      Reptiles 
    ##           521            12           590         11753         27345          2559           515           196             2

``` r
trends[i, table(taxa_mod, REALM)]
```

    ##                REALM
    ## taxa_mod        Freshwater Marine Terrestrial
    ##   All                    0    520           1
    ##   Amphibians             2      0          10
    ##   Benthos                0    590           0
    ##   Birds                  0   9221        2532
    ##   Fish                 966  26379           0
    ##   Invertebrates          9   2484          66
    ##   Mammals                0    477          38
    ##   Plant                  1     67         128
    ##   Reptiles               0      0           2

### Choose the variance structure for mixed effects models

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

Average rates of turnover (with year 1)

``` r
trends[abs(temptrend) >= 0.5, .(ave = mean(Jtutrend), sd = sd(Jtutrend)/sqrt(.N))] # turnover per year for locations changing temperature
```

    ##          ave          sd
    ## 1: 0.1814753 0.004483237

``` r
trends[abs(temptrend) < 0.1, .(ave = mean(Jtutrend), sd = sd(Jtutrend)/sqrt(.N))] # not changing temperature
```

    ##           ave           sd
    ## 1: 0.04592881 0.0003434666

``` r
trends[temptrend >= 0.5, .(ave = mean(Jtutrend), sd = sd(Jtutrend)/sqrt(.N))] # warming
```

    ##          ave          sd
    ## 1: 0.1616492 0.006620151

``` r
trends[temptrend <= -0.5, .(ave = mean(Jtutrend), sd = sd(Jtutrend)/sqrt(.N))] # cooling
```

    ##          ave          sd
    ## 1: 0.1961448 0.005982736

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) < 35, .(ave = mean(Jtutrend), sd = sd(Jtutrend)/sqrt(.N))] # tropics and sub-tropics
```

    ##          ave         sd
    ## 1: 0.4165367 0.04623876

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) >= 35 & abs(rarefyID_y) < 66.56339, .(ave = mean(Jtutrend), sd = sd(Jtutrend)/sqrt(.N))] # temperate
```

    ##         ave          sd
    ## 1: 0.178024 0.004397925

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) >= 66.56339, .(ave = mean(Jtutrend), sd = sd(Jtutrend)/sqrt(.N))] # arctic
```

    ##          ave         sd
    ## 1: 0.1962083 0.05650476

Average rates of turnover (without year 1)

``` r
trends[abs(temptrend) >= 0.5, .(ave = mean(Jtutrendrem0, na.rm=TRUE), 
                                sd = sd(Jtutrendrem0, na.rm=TRUE)/sqrt(.N))] # turnover per year for locations changing temperature
```

    ##           ave          sd
    ## 1: 0.01886862 0.008306934

``` r
trends[abs(temptrend) < 0.1, .(ave = mean(Jtutrendrem0, na.rm=TRUE), 
                               sd = sd(Jtutrendrem0, na.rm=TRUE)/sqrt(.N))] # not changing temperature
```

    ##            ave          sd
    ## 1: 0.006235547 0.000464845

``` r
trends[temptrend >= 0.5, .(ave = mean(Jtutrendrem0, na.rm=TRUE), 
                           sd = sd(Jtutrendrem0, na.rm=TRUE)/sqrt(.N))] # warming
```

    ##           ave         sd
    ## 1: 0.01569674 0.01435269

``` r
trends[temptrend <= -0.5, .(ave = mean(Jtutrendrem0, na.rm=TRUE), 
                            sd = sd(Jtutrendrem0, na.rm=TRUE)/sqrt(.N))] # cooling
```

    ##           ave         sd
    ## 1: 0.02036127 0.01026449

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) < 35, 
       .(ave = mean(Jtutrendrem0, na.rm=TRUE), 
         sd = sd(Jtutrendrem0, na.rm=TRUE)/sqrt(.N))] # tropics and sub-tropics
```

    ##    ave sd
    ## 1: 0.2 NA

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) >= 35 & abs(rarefyID_y) < 66.56339, 
       .(ave = mean(Jtutrendrem0, na.rm=TRUE), 
         sd = sd(Jtutrendrem0, na.rm=TRUE)/sqrt(.N))] # temperate
```

    ##           ave          sd
    ## 1: 0.01838431 0.008387229

``` r
trends[abs(temptrend) >= 0.5 & abs(rarefyID_y) >= 66.56339, 
       .(ave = mean(Jtutrendrem0, na.rm=TRUE), 
         sd = sd(Jtutrendrem0, na.rm=TRUE)/sqrt(.N))] # arctic
```

    ##    ave sd
    ## 1: NaN NA

## Temperature-only model (Jtutrend, Jbetatrend, Horntrend)

### With year 1

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

### Without year 1

``` r
i4 <- trends[, complete.cases(Jtutrendrem0, REALM, temptrend)]

randef <- list(STUDY_ID = ~ abs(temptrend), rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modonlyTtrendrem0.rds')){
  modonlyTtrendrem0 <- readRDS('temp/modonlyTtrendrem0.rds')
} else {
  modonlyTtrendrem0 <- lme(Jtutrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i4,], method = 'REML')
  saveRDS(modonlyTtrendrem0, file = 'temp/modonlyTtrendrem0.rds')
}

i5 <- trends[, complete.cases(Jbetatrendrem0, REALM, temptrend)]
if(file.exists('temp/modonlyTtrendJbetarem0.rds')){
  modonlyTtrendJbetarem0 <- readRDS('temp/modonlyTtrendJbetarem0.rds')
} else {
  modonlyTtrendJbetarem0 <- lme(Jbetatrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i5,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetarem0, file = 'temp/modonlyTtrendJbetarem0.rds')
}

i6 <- trends[, complete.cases(Horntrendrem0, REALM, temptrend)]
if(file.exists('temp/modonlyTtrendHornrem0.rds')){
  modonlyTtrendHornrem0 <- readRDS('temp/modonlyTtrendHornrem0.rds')
} else {
  modonlyTtrendHornrem0 <- lme(Horntrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i6,], method = 'REML')
  saveRDS(modonlyTtrendHornrem0, file = 'temp/modonlyTtrendHornrem0.rds')
}

summary(modonlyTtrendrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i4, ] 
    ##         AIC       BIC   logLik
    ##   -104458.6 -104356.5 52241.32
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev      Corr  
    ## (Intercept)    0.007657053 (Intr)
    ## abs(temptrend) 0.197650298 -0.935
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01097849 2.000161
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.121022 
    ## Fixed effects: Jtutrendrem0 ~ abs(temptrend) * REALM 
    ##                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                      0.00206219 0.00357904 36452  0.5761868  0.5645
    ## abs(temptrend)                   0.00391191 0.09140164 36452  0.0427991  0.9659
    ## REALMMarine                      0.00405551 0.00378073   289  1.0726812  0.2843
    ## REALMTerrestrial                -0.00020593 0.00397072   289 -0.0518610  0.9587
    ## abs(temptrend):REALMMarine       0.01183742 0.09631002 36452  0.1229095  0.9022
    ## abs(temptrend):REALMTerrestrial  0.08779171 0.10065203 36452  0.8722299  0.3831
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.807                                
    ## REALMMarine                     -0.947  0.764                         
    ## REALMTerrestrial                -0.901  0.727  0.853                  
    ## abs(temptrend):REALMMarine       0.766 -0.949 -0.814 -0.690           
    ## abs(temptrend):REALMTerrestrial  0.733 -0.908 -0.693 -0.810  0.862    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.57101792 -0.22892007 -0.01904447  0.27434430  5.64177564 
    ## 
    ## Number of Observations: 36747
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    292                  36747

``` r
summary(modonlyTtrendJbetarem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i5, ] 
    ##         AIC       BIC   logLik
    ##   -134862.5 -134760.4 67443.27
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev      Corr  
    ## (Intercept)    0.005965236 (Intr)
    ## abs(temptrend) 0.160013131 -0.038
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.003861919 0.9688156
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.887978 
    ## Fixed effects: Jbetatrendrem0 ~ abs(temptrend) * REALM 
    ##                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                      0.00415073 0.00262295 36452  1.5824664  0.1136
    ## abs(temptrend)                   0.10667907 0.06952719 36452  1.5343504  0.1250
    ## REALMMarine                      0.00275950 0.00278740   289  0.9899899  0.3230
    ## REALMTerrestrial                 0.00086624 0.00289642   289  0.2990741  0.7651
    ## abs(temptrend):REALMMarine      -0.08361477 0.07330769 36452 -1.1406003  0.2540
    ## abs(temptrend):REALMTerrestrial -0.01111424 0.07655288 36452 -0.1451839  0.8846
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.543                                
    ## REALMMarine                     -0.941  0.511                         
    ## REALMTerrestrial                -0.906  0.492  0.852                  
    ## abs(temptrend):REALMMarine       0.515 -0.948 -0.524 -0.467           
    ## abs(temptrend):REALMTerrestrial  0.494 -0.908 -0.464 -0.542  0.861    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.57430409 -0.30896466 -0.02462766  0.32539894  8.23224466 
    ## 
    ## Number of Observations: 36747
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    292                  36747

``` r
summary(modonlyTtrendHornrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i6, ] 
    ##      AIC      BIC logLik
    ##   -99770 -99668.1  49897
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev     Corr  
    ## (Intercept)    0.01201968 (Intr)
    ## abs(temptrend) 0.25314398 0.016 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01882813 2.441788
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.318769 
    ## Fixed effects: Horntrendrem0 ~ abs(temptrend) * REALM 
    ##                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                      0.00450157 0.00551034 35745  0.8169315  0.4140
    ## abs(temptrend)                   0.16293554 0.11769902 35745  1.3843406  0.1663
    ## REALMMarine                      0.00475924 0.00583586   254  0.8155169  0.4155
    ## REALMTerrestrial                 0.00211647 0.00613279   254  0.3451070  0.7303
    ## abs(temptrend):REALMMarine      -0.04174830 0.12448593 35745 -0.3353656  0.7374
    ## abs(temptrend):REALMTerrestrial -0.03620605 0.12997739 35745 -0.2785565  0.7806
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.472                                
    ## REALMMarine                     -0.944  0.446                         
    ## REALMTerrestrial                -0.899  0.424  0.848                  
    ## abs(temptrend):REALMMarine       0.447 -0.945 -0.455 -0.401           
    ## abs(temptrend):REALMTerrestrial  0.428 -0.906 -0.404 -0.473  0.856    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.75457482 -0.22483445 -0.02254226  0.23615730  5.73488343 
    ## 
    ## Number of Observations: 36005
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    257                  36005

### Plot the temp-only coefficients

``` r
colors <- brewer.pal(6, 'Dark2')

# make table of coefficients
coefs <- as.data.frame(summary(modonlyTtrend)$tTable)
coefs2 <- as.data.frame(summary(modonlyTtrendJbeta)$tTable)
coefs3 <- as.data.frame(summary(modonlyTtrendHorn)$tTable)
coefs4 <- as.data.frame(summary(modonlyTtrendrem0)$tTable)
coefs5 <- as.data.frame(summary(modonlyTtrendJbetarem0)$tTable)
coefs6 <- as.data.frame(summary(modonlyTtrendHornrem0)$tTable)
coefs$mod <- 'Jtu'
coefs2$mod <- 'Jbeta'
coefs3$mod <- 'Horn'
coefs4$mod <- 'Jturem0'
coefs5$mod <- 'Jbetarem0'
coefs6$mod <- 'Hornrem0'
rows1 <- which(grepl('temptrend', rownames(coefs))) # extract temperature effect
cols <- c('Value', 'Std.Error', 'mod')
allcoefs <- rbind(coefs[rows1, cols], coefs2[rows1, cols], coefs3[rows1, cols], 
                  coefs4[rows1, cols], coefs5[rows1, cols], coefs6[rows1, cols])
allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to marine effects
allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to terrestrial effects

allcoefs$lCI <- allcoefs$Value - allcoefs$Std.Error # lower confidence interval
allcoefs$uCI <- allcoefs$Value + allcoefs$Std.Error
allcoefs$y <- c(3, 2, 1) + rep(c(0, -0.1, -0.2, -0.3, -0.4, -0.5), c(3, 3, 3, 3, 3, 3)) # y-values
allcoefs$col <- c(rep(colors[1], 3), rep(colors[2], 3), rep(colors[3], 3), 
                  rep(colors[4], 3), rep(colors[5], 3), rep(colors[6], 3))
allcoefs$realm <- rep(c('Freshwater', 'Marine', 'Terrestrial'), 6)

par(las = 1, mai = c(0.8, 2, 0.1, 0.1))
plot(0,0, col = 'white', xlim=c(-0.1, 0.85), ylim = c(0.5,3), 
     yaxt='n', xlab = 'Turnover per |°C/yr|', ylab ='')
axis(2, at = 3:1, labels = c('Freshwater', 'Marine', 'Terrestrial'), cex.axis = 0.7)
abline(v = 0, col = 'grey')
for(i in 1:nrow(allcoefs)){
  with(allcoefs[i, ], points(Value, y, pch = 16, col = col))
  with(allcoefs[i, ], lines(x = c(lCI, uCI), y = c(y, y), col = col))
}
legend('bottomright', col = colors, lwd = 1, pch = 16, 
       legend = c('Jaccard turnover', 'Jaccard total', 'Horn-Morisita',
                  'Jaccard turnover rem0', 'Jaccard total rem0', 'Horn-Morisita rem0'))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/modonlyTtrendsimp%20coefs-1.png)<!-- -->

### Nicer plots of turnover vs. temperature data and model fit

Scatterplot, violin plots, and coefficient plots all together Off for
now

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

### Fit full models

#### Full model for Jaccard turnover

Try Bowler or Venter/Halpern human impact

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
i1.2 <- trends[, complete.cases(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
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
                   random = randef, weights = varef, data = trends[i1.2,], method = 'REML')
  saveRDS(modTfullfootprint, file = 'temp/modTfullfootprint.rds')
}

summary(modTfull1)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -134114.8 -133715.6 67103.42
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.04686628 (Intr)
    ## temptrend_abs.sc 0.03743575 -0.109
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.000833882 0.2868584
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.239456 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave.sc + temptrend_abs.sc * tempave_metab.sc +      temptrend_abs.sc * seas.sc + temptrend_abs.sc * microclim.sc +      temptrend_abs.sc * mass.sc + temptrend_abs.sc * speed.sc +      temptrend_abs.sc * lifespan.sc + temptrend_abs.sc * consumerfrac.sc +      temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.00624996 0.013522842 43221   0.462178  0.6440
    ## temptrend_abs.sc                                  0.05256641 0.018691270 43221   2.812351  0.0049
    ## REALMMarine                                       0.05079819 0.014313195   232   3.549046  0.0005
    ## REALMTerrestrial                                  0.03116727 0.015711282   232   1.983751  0.0485
    ## tsign1                                           -0.01029608 0.000585304 43221 -17.591001  0.0000
    ## tempave.sc                                       -0.00390212 0.000869560 43221  -4.487461  0.0000
    ## tempave_metab.sc                                  0.00699218 0.002000215 43221   3.495715  0.0005
    ## seas.sc                                          -0.00057830 0.000440756 43221  -1.312074  0.1895
    ## microclim.sc                                      0.00049592 0.000225298 43221   2.201167  0.0277
    ## mass.sc                                          -0.00005477 0.000819793 43221  -0.066806  0.9467
    ## speed.sc                                          0.00073186 0.000804142 43221   0.910118  0.3628
    ## lifespan.sc                                      -0.00338539 0.001719600 43221  -1.968706  0.0490
    ## consumerfrac.sc                                   0.00132224 0.000926806 43221   1.426659  0.1537
    ## endothermfrac.sc                                 -0.01635500 0.004085713 43221  -4.002973  0.0001
    ## nspp.sc                                          -0.00572662 0.000444697 43221 -12.877559  0.0000
    ## npp.sc                                            0.00073215 0.000354391 43221   2.065941  0.0388
    ## veg.sc                                            0.00084978 0.000440302 43221   1.929995  0.0536
    ## temptrend_abs.sc:REALMMarine                     -0.00704472 0.019290557 43221  -0.365190  0.7150
    ## temptrend_abs.sc:REALMTerrestrial                -0.03046043 0.022303675 43221  -1.365714  0.1720
    ## temptrend_abs.sc:tsign1                           0.00884845 0.001074661 43221   8.233710  0.0000
    ## temptrend_abs.sc:tempave.sc                       0.01313199 0.002126877 43221   6.174310  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                 0.00609489 0.004817364 43221   1.265192  0.2058
    ## temptrend_abs.sc:seas.sc                          0.00348024 0.001030909 43221   3.375891  0.0007
    ## temptrend_abs.sc:microclim.sc                    -0.00273933 0.000651675 43221  -4.203515  0.0000
    ## temptrend_abs.sc:mass.sc                         -0.00590075 0.001803839 43221  -3.271220  0.0011
    ## temptrend_abs.sc:speed.sc                        -0.00734962 0.001611263 43221  -4.561403  0.0000
    ## temptrend_abs.sc:lifespan.sc                      0.01263258 0.003853999 43221   3.277785  0.0010
    ## temptrend_abs.sc:consumerfrac.sc                 -0.00043569 0.001332391 43221  -0.327001  0.7437
    ## temptrend_abs.sc:endothermfrac.sc                 0.01528921 0.006304006 43221   2.425316  0.0153
    ## temptrend_abs.sc:nspp.sc                         -0.00608065 0.000938280 43221  -6.480630  0.0000
    ## tsign-1:thermal_bias.sc                          -0.00044560 0.000642677 43221  -0.693354  0.4881
    ## tsign1:thermal_bias.sc                           -0.00030504 0.000537044 43221  -0.567992  0.5700
    ## temptrend_abs.sc:npp.sc                          -0.00334774 0.000999382 43221  -3.349807  0.0008
    ## temptrend_abs.sc:veg.sc                          -0.00017089 0.001218044 43221  -0.140299  0.8884
    ## human_bowler.sc:REALM2TerrFresh                   0.00014064 0.000365563 43221   0.384721  0.7004
    ## human_bowler.sc:REALM2Marine                      0.00064985 0.000273742 43221   2.373966  0.0176
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.00446804 0.001203202 43221   3.713462  0.0002
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00131604 0.001166689 43221   1.128010  0.2593
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.00183998 0.000807812 43221  -2.277735  0.0227
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.00291421 0.000882591 43221  -3.301879  0.0010
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.408                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.928  0.383                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.865  0.350  0.801                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.037  0.022  0.009  0.012                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.003  0.001 -0.002 -0.002 -0.071                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.052 -0.023 -0.045 -0.041 -0.023 -0.469                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.045  0.042  0.063 -0.014 -0.109  0.183 -0.099                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.016  0.022  0.025 -0.008  0.012 -0.042 -0.026  0.142                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                           0.014  0.009 -0.001  0.001  0.020  0.026 -0.498  0.000  0.034                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.002 -0.005  0.013  0.021 -0.047  0.103 -0.279  0.100  0.068  0.140                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.032 -0.022 -0.027 -0.044 -0.015 -0.087  0.707  0.013 -0.009 -0.774 -0.454                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.017  0.009  0.030  0.272  0.017 -0.019 -0.038 -0.043 -0.007  0.015 -0.070 -0.047                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                  0.146 -0.050 -0.074 -0.308  0.008  0.142 -0.165  0.027  0.011 -0.014 -0.169  0.067 -0.331                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                          -0.009 -0.017 -0.013 -0.027 -0.039  0.004 -0.015 -0.045 -0.099 -0.214  0.007  0.145 -0.002  0.085                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.010  0.006 -0.014  0.014  0.014 -0.276  0.213 -0.147 -0.157 -0.041  0.047  0.033 -0.028 -0.083 -0.180                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.109  0.141  0.113  0.000 -0.016  0.145 -0.104  0.130  0.001  0.033  0.046 -0.026 -0.049  0.007  0.004 -0.283                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.392 -0.960 -0.394 -0.330 -0.005 -0.001  0.019 -0.055 -0.026 -0.014 -0.007  0.017 -0.010  0.033  0.026  0.011 -0.155                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.340 -0.843 -0.314 -0.415 -0.005 -0.012  0.022  0.038  0.015 -0.014 -0.007  0.027 -0.101  0.134  0.042 -0.028  0.009  0.796                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.019 -0.046 -0.007 -0.005 -0.486  0.044  0.000  0.014 -0.037 -0.013  0.017  0.007 -0.006  0.004  0.026 -0.049 -0.001  0.015      0.007                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.003 -0.009 -0.006 -0.004  0.048 -0.595  0.309 -0.054  0.078  0.024 -0.059  0.029 -0.013 -0.109 -0.016  0.153 -0.137  0.005      0.024     -0.069                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.014  0.027  0.013  0.014 -0.011  0.326 -0.582 -0.010 -0.029  0.250  0.078 -0.346  0.058  0.134 -0.009 -0.139  0.064 -0.022     -0.037      0.016 -0.441                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.021 -0.061 -0.029  0.029  0.043 -0.010 -0.005 -0.618 -0.092  0.008 -0.035 -0.008  0.031  0.011  0.030  0.081 -0.131  0.084     -0.078     -0.019  0.048            0.047                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.009 -0.030 -0.012  0.010 -0.040  0.075 -0.012 -0.065 -0.673 -0.015 -0.027  0.008  0.001  0.007  0.039  0.127  0.026  0.034     -0.032      0.048 -0.132            0.108   0.136                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.007 -0.015 -0.009 -0.010 -0.014  0.024  0.265  0.007 -0.010 -0.552 -0.030  0.436 -0.011  0.019  0.130  0.006 -0.004  0.021      0.027      0.024 -0.075           -0.471  -0.044            -0.011                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.005  0.020 -0.004 -0.002  0.024 -0.087  0.129 -0.032 -0.036 -0.045 -0.536  0.187  0.031  0.092 -0.021  0.003 -0.052  0.003      0.001     -0.005  0.157           -0.223   0.069             0.041             0.059                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.013  0.032  0.013  0.016  0.015  0.038 -0.377 -0.010  0.007  0.433  0.148 -0.529  0.055 -0.027 -0.108 -0.018 -0.003 -0.020     -0.042     -0.017 -0.059            0.679   0.034             0.002            -0.807            -0.303                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.004 -0.074 -0.006 -0.100 -0.010 -0.008  0.067  0.036  0.004 -0.018  0.060  0.061 -0.477  0.102  0.010  0.032  0.035  0.074      0.304      0.006  0.028           -0.152  -0.043            -0.003             0.052            -0.094            -0.140                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                -0.043  0.142  0.028  0.131  0.004 -0.189  0.218  0.011  0.026  0.023  0.132 -0.033  0.086 -0.397 -0.064  0.087 -0.028 -0.096     -0.344     -0.019  0.234           -0.365  -0.044            -0.075            -0.027            -0.192             0.018           -0.256                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.012  0.033  0.014  0.029  0.022 -0.025  0.018  0.020  0.040  0.132 -0.024 -0.106  0.008 -0.065 -0.530  0.080 -0.032 -0.045     -0.076     -0.037 -0.001            0.038  -0.039            -0.016            -0.211             0.019             0.199           -0.019            0.110                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.022 -0.008 -0.024 -0.001 -0.168  0.474 -0.113 -0.098 -0.073  0.036  0.077 -0.094  0.004  0.000 -0.002 -0.089  0.042  0.013     -0.011      0.019 -0.267            0.085   0.081             0.046            -0.007            -0.063             0.044           -0.009           -0.007            -0.005                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.035 -0.017 -0.045 -0.009 -0.033  0.690 -0.156 -0.228 -0.178  0.054  0.058 -0.093  0.005  0.016 -0.038 -0.126  0.047  0.018     -0.005      0.070 -0.331            0.102   0.105             0.069            -0.008            -0.055             0.040           -0.008           -0.030            -0.006             0.578                                                                                               
    ## temptrend_abs.sc:npp.sc                           0.011  0.004 -0.006 -0.022 -0.013  0.098 -0.076  0.056  0.133  0.009  0.000 -0.010  0.017  0.024  0.066 -0.626  0.179 -0.029      0.049      0.040 -0.190            0.153  -0.195            -0.265             0.018            -0.037             0.024           -0.049           -0.081            -0.085             0.032  0.058                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.073 -0.195 -0.076  0.007  0.009 -0.097  0.052 -0.121  0.023 -0.008 -0.042  0.007  0.022 -0.006 -0.027  0.185 -0.704  0.216     -0.022      0.021  0.177           -0.069   0.222            -0.005            -0.016             0.081             0.018           -0.038            0.027             0.031            -0.032 -0.048 -0.315                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.030  0.056  0.023 -0.009  0.036 -0.274  0.102 -0.067  0.144  0.049 -0.023 -0.050 -0.046 -0.076  0.018 -0.075  0.161 -0.061      0.019     -0.026  0.149           -0.083  -0.139            -0.034            -0.019             0.013             0.013            0.018            0.068             0.033            -0.078 -0.089  0.049            -0.237                                                               
    ## human_bowler.sc:REALM2Marine                      0.008 -0.015 -0.009 -0.003 -0.004  0.007  0.008 -0.094 -0.077 -0.022 -0.001  0.011 -0.001  0.006  0.000 -0.172  0.037  0.016      0.005      0.015  0.012           -0.020   0.109             0.026             0.020            -0.010            -0.004            0.011            0.002             0.030             0.066  0.068  0.102            -0.011            0.007                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.005  0.002  0.004 -0.005  0.018 -0.356  0.117  0.081  0.078 -0.011 -0.053  0.050 -0.010 -0.015 -0.007  0.061 -0.073 -0.016      0.024      0.106  0.684           -0.138  -0.102            -0.104             0.008             0.125            -0.091            0.013            0.002            -0.034            -0.469 -0.363 -0.086             0.096            0.054      -0.053                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.008  0.008  0.007 -0.003  0.068 -0.387  0.111  0.087  0.109 -0.007 -0.042  0.045 -0.009 -0.013  0.003  0.071 -0.081 -0.015      0.021     -0.129  0.717           -0.113  -0.073            -0.137            -0.005             0.104            -0.079            0.012           -0.006            -0.047            -0.338 -0.507 -0.101             0.113            0.044      -0.059       0.790                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.034 -0.111 -0.035  0.013 -0.017  0.173 -0.086 -0.116 -0.057 -0.027  0.010  0.025  0.018  0.051  0.017  0.037 -0.292  0.120     -0.034      0.051 -0.232            0.128   0.259             0.089             0.039            -0.026            -0.025           -0.028           -0.102            -0.035             0.055  0.060 -0.103             0.471           -0.623       0.001      -0.105 -0.100               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.008  0.021  0.008  0.002  0.004 -0.005 -0.010  0.094  0.023  0.016 -0.013 -0.009  0.008 -0.001  0.003  0.086  0.002 -0.023      0.000     -0.009 -0.045            0.062  -0.201            -0.069            -0.019             0.011             0.002           -0.016           -0.031            -0.052            -0.051 -0.060 -0.147             0.010            0.010      -0.623       0.057  0.057 -0.007        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.34716400 -0.35364009  0.05228794  0.52742556  6.59480022 
    ## 
    ## Number of Observations: 43493
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    235                  43493

``` r
summary(modTfullfootprint)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i1.2, ] 
    ##         AIC       BIC   logLik
    ##   -132997.1 -132598.2 66544.57
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.04731816 (Intr)
    ## temptrend_abs.sc 0.03945538 -0.122
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 2.288017e-06 0.2861346
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.238289 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave.sc + temptrend_abs.sc * tempave_metab.sc +      temptrend_abs.sc * seas.sc + temptrend_abs.sc * microclim.sc +      temptrend_abs.sc * mass.sc + temptrend_abs.sc * speed.sc +      temptrend_abs.sc * lifespan.sc + temptrend_abs.sc * consumerfrac.sc +      temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_footprint.sc:REALM2 
    ##                                                           Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                          0.00614488 0.013718684 42910   0.447920  0.6542
    ## temptrend_abs.sc                                     0.04784413 0.019396502 42910   2.466637  0.0136
    ## REALMMarine                                          0.05315323 0.014571491   224   3.647755  0.0003
    ## REALMTerrestrial                                     0.03211621 0.015952565   224   2.013231  0.0453
    ## tsign1                                              -0.01029118 0.000590575 42910 -17.425700  0.0000
    ## tempave.sc                                          -0.00417666 0.000862401 42910  -4.843062  0.0000
    ## tempave_metab.sc                                     0.00684280 0.002016805 42910   3.392893  0.0007
    ## seas.sc                                             -0.00079182 0.000439357 42910  -1.802236  0.0715
    ## microclim.sc                                         0.00064317 0.000240691 42910   2.672195  0.0075
    ## mass.sc                                             -0.00020556 0.000825405 42910  -0.249036  0.8033
    ## speed.sc                                             0.00094404 0.000808816 42910   1.167187  0.2431
    ## lifespan.sc                                         -0.00347648 0.001729701 42910  -2.009875  0.0445
    ## consumerfrac.sc                                      0.00110799 0.000956573 42910   1.158286  0.2468
    ## endothermfrac.sc                                    -0.01763629 0.004146569 42910  -4.253224  0.0000
    ## nspp.sc                                             -0.00543570 0.000450500 42910 -12.065933  0.0000
    ## npp.sc                                               0.00093928 0.000377888 42910   2.485596  0.0129
    ## veg.sc                                               0.00045953 0.000455763 42910   1.008275  0.3133
    ## temptrend_abs.sc:REALMMarine                        -0.00131806 0.020036337 42910  -0.065783  0.9476
    ## temptrend_abs.sc:REALMTerrestrial                   -0.03006858 0.023344785 42910  -1.288021  0.1977
    ## temptrend_abs.sc:tsign1                              0.00888970 0.001076800 42910   8.255666  0.0000
    ## temptrend_abs.sc:tempave.sc                          0.01199939 0.002140077 42910   5.606992  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                    0.00827095 0.004835214 42910   1.710566  0.0872
    ## temptrend_abs.sc:seas.sc                             0.00334333 0.000999735 42910   3.344214  0.0008
    ## temptrend_abs.sc:microclim.sc                       -0.00312228 0.000673529 42910  -4.635699  0.0000
    ## temptrend_abs.sc:mass.sc                            -0.00535362 0.001809402 42910  -2.958779  0.0031
    ## temptrend_abs.sc:speed.sc                           -0.00762283 0.001623110 42910  -4.696436  0.0000
    ## temptrend_abs.sc:lifespan.sc                         0.01255989 0.003880775 42910   3.236438  0.0012
    ## temptrend_abs.sc:consumerfrac.sc                    -0.00040636 0.001399341 42910  -0.290392  0.7715
    ## temptrend_abs.sc:endothermfrac.sc                    0.01419937 0.006587460 42910   2.155515  0.0311
    ## temptrend_abs.sc:nspp.sc                            -0.00699961 0.000955097 42910  -7.328689  0.0000
    ## tsign-1:thermal_bias.sc                             -0.00039022 0.000653171 42910  -0.597427  0.5502
    ## tsign1:thermal_bias.sc                              -0.00050174 0.000538447 42910  -0.931821  0.3514
    ## temptrend_abs.sc:npp.sc                             -0.00433863 0.001023871 42910  -4.237478  0.0000
    ## temptrend_abs.sc:veg.sc                              0.00193537 0.001224079 42910   1.581079  0.1139
    ## human_footprint.sc:REALM2TerrFresh                   0.00032637 0.000327150 42910   0.997614  0.3185
    ## human_footprint.sc:REALM2Marine                     -0.00010411 0.000245884 42910  -0.423410  0.6720
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc             0.00437477 0.001216408 42910   3.596466  0.0003
    ## temptrend_abs.sc:tsign1:thermal_bias.sc              0.00119593 0.001182277 42910   1.011548  0.3118
    ## temptrend_abs.sc:human_footprint.sc:REALM2TerrFresh  0.00025085 0.000810935 42910   0.309339  0.7571
    ## temptrend_abs.sc:human_footprint.sc:REALM2Marine    -0.00102449 0.000673358 42910  -1.521470  0.1281
    ##  Correlation: 
    ##                                                     (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                    -0.417                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                         -0.925  0.390                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                    -0.865  0.359  0.799                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                              -0.036  0.020  0.008  0.012                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                          -0.011  0.009  0.005 -0.005 -0.063                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                     0.057 -0.028 -0.049 -0.039 -0.027 -0.474                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                             -0.047  0.034  0.064 -0.013 -0.113  0.192 -0.109                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                        -0.017  0.016  0.028 -0.007  0.004  0.021 -0.055  0.183                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                              0.015  0.008 -0.003  0.001  0.021  0.037 -0.501 -0.006  0.011                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                             0.003 -0.006  0.013  0.020 -0.045  0.099 -0.272  0.093  0.071  0.137                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                          0.032 -0.022 -0.026 -0.044 -0.014 -0.099  0.710  0.014  0.006 -0.774 -0.449                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                     -0.015  0.007  0.031  0.263  0.016 -0.026 -0.036 -0.041  0.003  0.016 -0.067 -0.049                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                     0.145 -0.051 -0.076 -0.312  0.010  0.134 -0.163  0.029  0.027 -0.010 -0.168  0.063 -0.329                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                             -0.008 -0.017 -0.012 -0.028 -0.038  0.018 -0.019 -0.037 -0.104 -0.217  0.010  0.150 -0.003  0.090                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                               0.018  0.001 -0.024  0.013  0.015 -0.336  0.252 -0.208 -0.210 -0.042  0.048  0.031 -0.033 -0.099 -0.188                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                              -0.114  0.133  0.116  0.003 -0.018  0.195 -0.139  0.126  0.018  0.036  0.042 -0.028 -0.052  0.018  0.010 -0.331                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                         0.400 -0.958 -0.401 -0.339 -0.004 -0.011  0.024 -0.046 -0.020 -0.012 -0.008  0.016 -0.009  0.032  0.025  0.016 -0.145                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                    0.346 -0.845 -0.318 -0.420 -0.006 -0.008  0.018  0.034  0.016 -0.015 -0.006  0.028 -0.095  0.139  0.042 -0.027  0.000  0.797                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                              0.018 -0.043 -0.005 -0.005 -0.486  0.044  0.000  0.022 -0.033 -0.012  0.017  0.006 -0.004  0.003  0.027 -0.049  0.009  0.011      0.010                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                          0.008 -0.023 -0.012 -0.002  0.044 -0.604  0.312 -0.070  0.053  0.019 -0.059  0.035 -0.013 -0.106 -0.022  0.183 -0.189  0.021      0.019     -0.065                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                   -0.018  0.035  0.017  0.011 -0.008  0.332 -0.585  0.002 -0.014  0.251  0.076 -0.350  0.055  0.135 -0.005 -0.164  0.101 -0.030     -0.032      0.012 -0.435                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                             0.018 -0.048 -0.026  0.025  0.051 -0.037  0.006 -0.617 -0.089  0.017 -0.039 -0.015  0.029  0.006  0.030  0.107 -0.098  0.069     -0.071     -0.031  0.069            0.040                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                        0.006 -0.016 -0.011  0.008 -0.039  0.053  0.000 -0.060 -0.703 -0.006 -0.034  0.001 -0.002  0.001  0.040  0.153  0.051  0.021     -0.030      0.041 -0.106            0.099   0.110                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                             0.006 -0.013 -0.008 -0.010 -0.014  0.021  0.266  0.014  0.001 -0.553 -0.026  0.436 -0.011  0.017  0.133  0.009 -0.002  0.018      0.028      0.023 -0.070           -0.473  -0.056            -0.024                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                           -0.005  0.020 -0.006 -0.002  0.025 -0.090  0.127 -0.037 -0.040 -0.041 -0.537  0.182  0.030  0.091 -0.024  0.002 -0.052  0.004      0.000     -0.005  0.160           -0.220   0.075             0.045             0.056                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                        -0.013  0.033  0.013  0.016  0.016  0.042 -0.378 -0.013  0.001  0.432  0.143 -0.529  0.053 -0.025 -0.110 -0.018 -0.001 -0.020     -0.044     -0.018 -0.065            0.684   0.039             0.008            -0.804            -0.295                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                     0.004 -0.078 -0.008 -0.099 -0.010 -0.010  0.065  0.034  0.002 -0.019  0.058  0.060 -0.476  0.105  0.012  0.035  0.032  0.077      0.317      0.006  0.030           -0.143  -0.045             0.000             0.052            -0.092            -0.135                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                   -0.043  0.144  0.027  0.137  0.002 -0.180  0.211  0.002  0.010  0.021  0.127 -0.030  0.086 -0.403 -0.067  0.100 -0.049 -0.096     -0.364     -0.017  0.214           -0.350  -0.034            -0.062            -0.025            -0.183             0.016           -0.280                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                            -0.013  0.037  0.014  0.030  0.019 -0.028  0.020  0.026  0.035  0.135 -0.029 -0.108  0.011 -0.067 -0.537  0.092 -0.030 -0.047     -0.077     -0.039  0.001            0.037  -0.051            -0.016            -0.215             0.026             0.201           -0.023            0.113                                                                                                                                   
    ## tsign-1:thermal_bias.sc                              0.020 -0.003 -0.022 -0.002 -0.172  0.458 -0.104 -0.095 -0.058  0.043  0.075 -0.101  0.002 -0.006 -0.001 -0.083  0.052  0.006     -0.010      0.021 -0.260            0.080   0.065             0.038            -0.011            -0.059             0.048           -0.010           -0.001            -0.007                                                                                                                 
    ## tsign1:thermal_bias.sc                               0.033 -0.013 -0.043 -0.009 -0.029  0.692 -0.158 -0.226 -0.166  0.064  0.055 -0.101  0.004  0.013 -0.033 -0.123  0.057  0.013     -0.005      0.069 -0.340            0.106   0.092             0.066            -0.013            -0.053             0.043           -0.012           -0.027            -0.010             0.569                                                                                               
    ## temptrend_abs.sc:npp.sc                              0.007  0.006 -0.001 -0.020 -0.016  0.140 -0.107  0.089  0.155  0.010 -0.006 -0.010  0.021  0.037  0.077 -0.637  0.191 -0.031      0.046      0.044 -0.237            0.191  -0.226            -0.267             0.017            -0.041             0.022           -0.051           -0.103            -0.101             0.031  0.055                                                                                        
    ## temptrend_abs.sc:veg.sc                              0.074 -0.187 -0.077  0.001  0.017 -0.152  0.084 -0.102  0.046 -0.006 -0.039  0.006  0.022 -0.019 -0.031  0.203 -0.689  0.206     -0.009      0.006  0.259           -0.123   0.185            -0.062            -0.021             0.086             0.017           -0.036            0.059             0.024            -0.051 -0.064 -0.326                                                                                 
    ## human_footprint.sc:REALM2TerrFresh                  -0.032  0.057  0.026  0.004  0.015 -0.157  0.045 -0.053  0.037  0.036 -0.030 -0.036 -0.041 -0.042 -0.020 -0.021  0.172 -0.058     -0.003     -0.023  0.115           -0.050  -0.075             0.009            -0.022             0.031             0.022            0.024            0.045             0.054            -0.037 -0.071  0.016            -0.199                                                               
    ## human_footprint.sc:REALM2Marine                     -0.012  0.005  0.008  0.013 -0.031  0.006 -0.008  0.058  0.072 -0.013 -0.043 -0.005  0.017 -0.018 -0.064  0.052 -0.025 -0.003     -0.006      0.006  0.002           -0.013  -0.017            -0.063             0.008             0.008             0.002           -0.003            0.018             0.050            -0.009 -0.001 -0.022             0.015           -0.006                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc            -0.002 -0.010  0.000 -0.005  0.019 -0.358  0.116  0.069  0.072 -0.017 -0.050  0.056 -0.011 -0.011 -0.007  0.057 -0.099 -0.003      0.025      0.107  0.686           -0.138  -0.081            -0.089             0.015             0.120            -0.099            0.018           -0.008            -0.034            -0.468 -0.366 -0.090             0.141            0.029      -0.002                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc             -0.005 -0.002  0.004 -0.003  0.065 -0.399  0.114  0.076  0.101 -0.013 -0.041  0.049 -0.011 -0.011  0.002  0.066 -0.101 -0.004      0.023     -0.126  0.727           -0.116  -0.058            -0.118             0.001             0.101            -0.084            0.018           -0.012            -0.042            -0.334 -0.510 -0.104             0.148            0.051      -0.014       0.792                      
    ## temptrend_abs.sc:human_footprint.sc:REALM2TerrFresh  0.031 -0.096 -0.030 -0.001 -0.002  0.123 -0.059 -0.071  0.028 -0.021  0.023  0.023  0.019  0.036  0.037 -0.004 -0.216  0.099     -0.001      0.037 -0.173            0.072   0.179            -0.041             0.038            -0.047            -0.038           -0.029           -0.069            -0.079             0.019  0.051 -0.034             0.370           -0.648       0.000      -0.051 -0.085               
    ## temptrend_abs.sc:human_footprint.sc:REALM2Marine     0.001 -0.008  0.000 -0.005  0.001  0.017 -0.018  0.009 -0.062 -0.003  0.004  0.004 -0.006  0.022  0.037 -0.031  0.018  0.005      0.005     -0.004 -0.011            0.031   0.035             0.124            -0.016            -0.023            -0.001            0.008           -0.027            -0.036            -0.001 -0.021  0.072            -0.041           -0.007      -0.655       0.007  0.029  0.000        
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -8.3520451 -0.3535263  0.0534660  0.5265006  6.6044778 
    ## 
    ## Number of Observations: 43174
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    227                  43174

#### Full model for Jaccard turnover, removing year 1

Try Bowler or Venter/Halpern human impact

``` r
# using Bowler for human impact
i1 <- trends[, complete.cases(Jtutrendrem0, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, human_bowler.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modTfullrem0.rds')){
  modTfullrem0 <- readRDS('temp/modTfullrem0.rds')
} else {
  modTfullrem0 <- lme(Jtutrendrem0 ~ temptrend_abs.sc*REALM +
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
                   random = randef, weights = varef, data = trends[i1,], method = 'REML')
  saveRDS(modTfullrem0, file = 'temp/modTfullrem0.rds')
}

# using Venter/Halpern for human impact
i2 <- trends[, complete.cases(Jtutrendrem0, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, human_footprint.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modTfullfootprintrem0.rds')){
  modTfullfootprintrem0 <- readRDS('temp/modTfullfootprintrem0.rds')
} else {
  modTfullfootprintrem0 <- lme(Jtutrendrem0 ~ temptrend_abs.sc*REALM + 
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
                   random = randef, weights = varef, data = trends[i2,], method = 'REML')
  saveRDS(modTfullfootprintrem0, file = 'temp/modTfullfootprintrem0.rds')
}

summary(modTfullrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i1, ] 
    ##         AIC       BIC   logLik
    ##   -82585.38 -82201.62 41338.69
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.007404753 (Intr)
    ## temptrend_abs.sc 0.012100255 0.189 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev: 0.003350251 1.507894
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.902506 
    ## Fixed effects: Jtutrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave.sc + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * lifespan.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.001162278 0.004106608 30817  0.283026  0.7772
    ## temptrend_abs.sc                                  0.019813085 0.012634930 30817  1.568120  0.1169
    ## REALMMarine                                       0.006986296 0.004310491   215  1.620766  0.1065
    ## REALMTerrestrial                                  0.000819649 0.004332310   215  0.189195  0.8501
    ## tsign1                                           -0.002596164 0.000800320 30817 -3.243909  0.0012
    ## tempave.sc                                        0.001543837 0.001079243 30817  1.430481  0.1526
    ## tempave_metab.sc                                 -0.007257093 0.002445200 30817 -2.967894  0.0030
    ## seas.sc                                           0.000683894 0.000519392 30817  1.316720  0.1879
    ## microclim.sc                                      0.000394706 0.000315293 30817  1.251870  0.2106
    ## mass.sc                                           0.001578620 0.001014467 30817  1.556107  0.1197
    ## speed.sc                                         -0.000307793 0.000868426 30817 -0.354427  0.7230
    ## lifespan.sc                                      -0.004356582 0.001964540 30817 -2.217609  0.0266
    ## consumerfrac.sc                                   0.000024176 0.000381498 30817  0.063372  0.9495
    ## endothermfrac.sc                                  0.003115145 0.002163318 30817  1.439985  0.1499
    ## nspp.sc                                          -0.001086322 0.000535655 30817 -2.028024  0.0426
    ## npp.sc                                           -0.000882242 0.000464527 30817 -1.899229  0.0575
    ## veg.sc                                            0.001418114 0.000605285 30817  2.342886  0.0191
    ## temptrend_abs.sc:REALMMarine                     -0.014626080 0.013152425 30817 -1.112044  0.2661
    ## temptrend_abs.sc:REALMTerrestrial                 0.003502901 0.012853270 30817  0.272530  0.7852
    ## temptrend_abs.sc:tsign1                          -0.001525316 0.002371375 30817 -0.643220  0.5201
    ## temptrend_abs.sc:tempave.sc                      -0.004419464 0.003569618 30817 -1.238078  0.2157
    ## temptrend_abs.sc:tempave_metab.sc                 0.020270453 0.008383109 30817  2.418011  0.0156
    ## temptrend_abs.sc:seas.sc                         -0.002238972 0.001721619 30817 -1.300503  0.1934
    ## temptrend_abs.sc:microclim.sc                    -0.002902935 0.001135893 30817 -2.555643  0.0106
    ## temptrend_abs.sc:mass.sc                         -0.003355291 0.003463984 30817 -0.968622  0.3327
    ## temptrend_abs.sc:speed.sc                         0.000036211 0.002806841 30817  0.012901  0.9897
    ## temptrend_abs.sc:lifespan.sc                      0.005281893 0.006832518 30817  0.773052  0.4395
    ## temptrend_abs.sc:consumerfrac.sc                 -0.001395529 0.001114227 30817 -1.252464  0.2104
    ## temptrend_abs.sc:endothermfrac.sc                -0.010479837 0.006901477 30817 -1.518492  0.1289
    ## temptrend_abs.sc:nspp.sc                          0.000583290 0.001699162 30817  0.343281  0.7314
    ## tsign-1:thermal_bias.sc                          -0.000427445 0.000813340 30817 -0.525543  0.5992
    ## tsign1:thermal_bias.sc                            0.000961820 0.000622615 30817  1.544807  0.1224
    ## temptrend_abs.sc:npp.sc                           0.004454219 0.001709018 30817  2.606302  0.0092
    ## temptrend_abs.sc:veg.sc                          -0.005593917 0.002227841 30817 -2.510914  0.0120
    ## human_bowler.sc:REALM2TerrFresh                   0.000919504 0.000479466 30817  1.917766  0.0551
    ## human_bowler.sc:REALM2Marine                      0.000433065 0.000414458 30817  1.044894  0.2961
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.000816471 0.002086341 30817 -0.391341  0.6955
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.000378278 0.001870718 30817  0.202210  0.8398
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.003071154 0.001540405 30817 -1.993731  0.0462
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.001910919 0.001547697 30817 -1.234686  0.2170
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.613                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.940  0.600                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.692  0.344  0.639                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.180  0.089  0.044  0.033                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.027  0.010  0.037  0.010 -0.061                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.121 -0.064 -0.131 -0.086 -0.013 -0.573                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.130  0.091  0.183 -0.087 -0.105  0.243 -0.078                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.098  0.105  0.112 -0.015  0.047  0.022 -0.072  0.146                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                          -0.015  0.056  0.026  0.072  0.021  0.001 -0.509  0.008  0.006                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.080 -0.035 -0.018 -0.071 -0.041  0.055 -0.224  0.052  0.049  0.115                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.079 -0.073 -0.090 -0.088 -0.041 -0.038  0.631  0.045  0.003 -0.848 -0.377                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.081  0.052  0.042  0.345  0.046  0.027 -0.042 -0.111 -0.016  0.104 -0.336 -0.026                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                 -0.026  0.006  0.031 -0.191  0.020  0.506 -0.542  0.106  0.063 -0.096 -0.289  0.117 -0.078                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                           0.033 -0.078 -0.101 -0.150 -0.043  0.017  0.010 -0.089 -0.113 -0.164  0.018  0.181 -0.052  0.168                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.037 -0.012 -0.063  0.044  0.004 -0.244  0.243 -0.084 -0.187 -0.036  0.043  0.033 -0.031 -0.217 -0.172                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.475  0.482  0.506  0.007 -0.002  0.129 -0.132  0.080  0.031  0.046  0.040 -0.055 -0.010  0.082  0.008 -0.273                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.599 -0.964 -0.613 -0.321 -0.029 -0.031  0.069 -0.112 -0.107 -0.054  0.010  0.070 -0.034 -0.006  0.104  0.055 -0.516                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.353 -0.648 -0.329 -0.540 -0.017  0.007  0.044  0.134  0.021 -0.080  0.058  0.072 -0.187  0.091  0.148 -0.068  0.021  0.600                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.107 -0.178 -0.050 -0.019 -0.519  0.079 -0.008  0.030 -0.046 -0.021  0.019  0.034 -0.027  0.013  0.031 -0.035 -0.023  0.069      0.019                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.047 -0.021 -0.061 -0.005  0.062 -0.730  0.460 -0.202  0.011  0.034 -0.020 -0.005 -0.010 -0.454 -0.029  0.176 -0.120  0.057     -0.035     -0.106                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.055  0.071  0.060  0.053 -0.018  0.446 -0.685  0.021 -0.006  0.307  0.127 -0.382  0.066  0.406 -0.013 -0.177  0.102 -0.077     -0.080      0.021 -0.498                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.077 -0.124 -0.105  0.124  0.067 -0.146  0.038 -0.732 -0.105  0.006 -0.017 -0.038  0.103 -0.073  0.069  0.065 -0.074  0.164     -0.214     -0.051  0.285           -0.036                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.088 -0.145 -0.094  0.017 -0.055  0.018  0.040 -0.081 -0.812  0.006 -0.025 -0.005  0.000 -0.039  0.064  0.182 -0.051  0.152     -0.038      0.048 -0.028            0.056   0.146                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.039 -0.051 -0.038 -0.062 -0.017  0.031  0.297  0.016  0.024 -0.603 -0.058  0.524 -0.064  0.068  0.109 -0.012 -0.017  0.054      0.122      0.028 -0.051           -0.539  -0.053            -0.049                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.058  0.057  0.034  0.080  0.034 -0.046  0.161  0.000 -0.011 -0.067 -0.630  0.222  0.224  0.160 -0.065 -0.009 -0.036 -0.017     -0.059     -0.017  0.064           -0.292   0.050             0.018             0.140                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.055  0.081  0.059  0.053  0.029 -0.001 -0.370 -0.044 -0.018  0.517  0.200 -0.598  0.037 -0.078 -0.129 -0.005  0.031 -0.077     -0.123     -0.036  0.008            0.661   0.064             0.034            -0.877            -0.394                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.048 -0.102 -0.029 -0.190 -0.033 -0.003  0.048  0.101 -0.006 -0.080  0.239  0.042 -0.622 -0.001  0.095  0.030 -0.002  0.077      0.304      0.021  0.003           -0.157  -0.130             0.009             0.140            -0.253            -0.152                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                 0.019  0.018 -0.014  0.055  0.008 -0.489  0.451 -0.071  0.002  0.068  0.203 -0.084 -0.045 -0.656 -0.090  0.193 -0.080 -0.010     -0.180     -0.023  0.539           -0.555   0.087            -0.052            -0.091            -0.288             0.120            0.016                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.080  0.097  0.095  0.135  0.040 -0.051  0.019  0.056  0.057  0.107 -0.077 -0.128  0.083 -0.119 -0.631  0.110 -0.011 -0.136     -0.215     -0.046  0.031            0.022  -0.065            -0.044            -0.153             0.060             0.202           -0.126            0.166                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.048 -0.007 -0.034 -0.001 -0.152  0.398 -0.164 -0.070 -0.028  0.018  0.044 -0.066  0.022  0.061  0.005 -0.035  0.028  0.002      0.005      0.012 -0.258            0.130   0.049             0.024            -0.011            -0.037             0.038            0.001           -0.084            -0.010                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.070 -0.027 -0.080 -0.009 -0.024  0.645 -0.271 -0.183 -0.117  0.022  0.042 -0.075  0.037  0.148 -0.008 -0.082  0.054  0.010      0.014      0.069 -0.362            0.175   0.099             0.066             0.000            -0.035             0.027           -0.003           -0.148            -0.035             0.502                                                                                               
    ## temptrend_abs.sc:npp.sc                          -0.015  0.071  0.038 -0.056  0.011  0.132 -0.141  0.048  0.185  0.008  0.000 -0.016  0.021  0.115  0.095 -0.760  0.239 -0.127      0.109      0.013 -0.202            0.177  -0.167            -0.275             0.027             0.005             0.000           -0.048           -0.194            -0.120             0.013  0.056                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.409 -0.551 -0.430  0.016 -0.012 -0.092  0.094 -0.068 -0.040 -0.029 -0.038  0.043  0.002 -0.056 -0.019  0.221 -0.868  0.593     -0.035      0.043  0.129           -0.097   0.116             0.093             0.003             0.039            -0.020           -0.004            0.078             0.002            -0.023 -0.055 -0.336                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.114  0.165  0.114 -0.052  0.037 -0.276  0.118 -0.088  0.166  0.059 -0.011 -0.061 -0.058 -0.143  0.004 -0.084  0.221 -0.174      0.044     -0.041  0.190           -0.110  -0.026            -0.122            -0.024            -0.002             0.029            0.031            0.147             0.036            -0.044 -0.079  0.082            -0.270                                                               
    ## human_bowler.sc:REALM2Marine                      0.046 -0.042 -0.056 -0.034  0.003  0.004  0.017 -0.177 -0.057 -0.006  0.001 -0.009 -0.025  0.008 -0.021 -0.152  0.037  0.049      0.026      0.029 -0.005           -0.021   0.154             0.027            -0.004            -0.014             0.023            0.032            0.013             0.055             0.087  0.120  0.124            -0.031            0.030                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.012  0.006 -0.017  0.010 -0.011 -0.339  0.190  0.043  0.044 -0.009 -0.017  0.035  0.009 -0.141 -0.032  0.023 -0.034 -0.005     -0.033      0.127  0.516           -0.161  -0.033            -0.066             0.025             0.037            -0.038           -0.017            0.111             0.011            -0.480 -0.355 -0.032             0.036            0.012      -0.078                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.006  0.039 -0.015  0.001  0.056 -0.397  0.195  0.058  0.086  0.005 -0.007  0.028  0.002 -0.143 -0.015  0.042 -0.059 -0.009     -0.026     -0.097  0.603           -0.137  -0.044            -0.126             0.003             0.010            -0.027            0.003            0.095            -0.020            -0.323 -0.612 -0.064             0.065            0.036      -0.119       0.615                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.149 -0.242 -0.156  0.040 -0.029  0.219 -0.112  0.017 -0.133 -0.047 -0.003  0.048  0.036  0.140  0.020  0.060 -0.299  0.251     -0.031      0.062 -0.252            0.131   0.000             0.177             0.051            -0.012            -0.048           -0.028           -0.180            -0.041             0.026  0.052 -0.117             0.378           -0.816      -0.025      -0.033 -0.073               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.036  0.049  0.043  0.026 -0.007 -0.004 -0.015  0.160  0.035 -0.003 -0.018  0.018  0.031  0.008  0.019  0.129 -0.030 -0.056     -0.025     -0.036 -0.017            0.049  -0.188            -0.066             0.020             0.020            -0.038           -0.041           -0.051            -0.071            -0.063 -0.101 -0.190             0.055           -0.025      -0.742       0.065  0.097  0.034        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -7.01782639 -0.26254472 -0.02396528  0.31237827  5.40866015 
    ## 
    ## Number of Observations: 31072
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    218                  31072

``` r
summary(modTfullfootprintrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -81632.18 -81248.82 40862.09
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.007176109 (Intr)
    ## temptrend_abs.sc 0.015663961 0.105 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev: 0.003466799 1.508598
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.903419 
    ## Fixed effects: Jtutrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave.sc + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * lifespan.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_footprint.sc:REALM2 
    ##                                                            Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                          0.002188386 0.004203789 30554  0.520575  0.6027
    ## temptrend_abs.sc                                     0.014349083 0.013515177 30554  1.061701  0.2884
    ## REALMMarine                                          0.005988774 0.004433031   207  1.350943  0.1782
    ## REALMTerrestrial                                     0.001134408 0.004404415   207  0.257562  0.7970
    ## tsign1                                              -0.002542309 0.000815299 30554 -3.118255  0.0018
    ## tempave.sc                                           0.001895972 0.001090728 30554  1.738263  0.0822
    ## tempave_metab.sc                                    -0.008159882 0.002490382 30554 -3.276559  0.0011
    ## seas.sc                                              0.000668306 0.000527638 30554  1.266598  0.2053
    ## microclim.sc                                         0.000341125 0.000325501 30554  1.047999  0.2946
    ## mass.sc                                              0.001333621 0.001029497 30554  1.295410  0.1952
    ## speed.sc                                            -0.000139640 0.000891453 30554 -0.156643  0.8755
    ## lifespan.sc                                         -0.004068677 0.001996438 30554 -2.037968  0.0416
    ## consumerfrac.sc                                      0.000164062 0.000408409 30554  0.401711  0.6879
    ## endothermfrac.sc                                     0.003749363 0.002203426 30554  1.701606  0.0888
    ## nspp.sc                                             -0.001061973 0.000552874 30554 -1.920824  0.0548
    ## npp.sc                                              -0.000732529 0.000485644 30554 -1.508366  0.1315
    ## veg.sc                                               0.001009293 0.000651152 30554  1.550012  0.1211
    ## temptrend_abs.sc:REALMMarine                        -0.007899392 0.014134205 30554 -0.558885  0.5762
    ## temptrend_abs.sc:REALMTerrestrial                    0.002162586 0.014115226 30554  0.153209  0.8782
    ## temptrend_abs.sc:tsign1                             -0.000748845 0.002413939 30554 -0.310217  0.7564
    ## temptrend_abs.sc:tempave.sc                         -0.006941302 0.003703389 30554 -1.874311  0.0609
    ## temptrend_abs.sc:tempave_metab.sc                    0.022711210 0.008555783 30554  2.654487  0.0079
    ## temptrend_abs.sc:seas.sc                            -0.002939493 0.001777935 30554 -1.653319  0.0983
    ## temptrend_abs.sc:microclim.sc                       -0.002966419 0.001164908 30554 -2.546483  0.0109
    ## temptrend_abs.sc:mass.sc                            -0.002233566 0.003514929 30554 -0.635451  0.5251
    ## temptrend_abs.sc:speed.sc                            0.001020249 0.002933502 30554  0.347792  0.7280
    ## temptrend_abs.sc:lifespan.sc                         0.002470059 0.006997851 30554  0.352974  0.7241
    ## temptrend_abs.sc:consumerfrac.sc                    -0.001703091 0.001213276 30554 -1.403713  0.1604
    ## temptrend_abs.sc:endothermfrac.sc                   -0.013909493 0.007255637 30554 -1.917060  0.0552
    ## temptrend_abs.sc:nspp.sc                             0.000131555 0.001786340 30554  0.073645  0.9413
    ## tsign-1:thermal_bias.sc                             -0.000481573 0.000845080 30554 -0.569855  0.5688
    ## tsign1:thermal_bias.sc                               0.000916580 0.000632451 30554  1.449251  0.1473
    ## temptrend_abs.sc:npp.sc                              0.003257921 0.001773047 30554  1.837470  0.0662
    ## temptrend_abs.sc:veg.sc                             -0.003200610 0.002358494 30554 -1.357056  0.1748
    ## human_footprint.sc:REALM2TerrFresh                  -0.000027522 0.000470363 30554 -0.058512  0.9533
    ## human_footprint.sc:REALM2Marine                     -0.000252934 0.000319225 30554 -0.792338  0.4282
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc            -0.000569406 0.002196856 30554 -0.259191  0.7955
    ## temptrend_abs.sc:tsign1:thermal_bias.sc              0.000166531 0.001984830 30554  0.083902  0.9331
    ## temptrend_abs.sc:human_footprint.sc:REALM2TerrFresh  0.001365182 0.001639706 30554  0.832577  0.4051
    ## temptrend_abs.sc:human_footprint.sc:REALM2Marine    -0.001766829 0.001118003 30554 -1.580344  0.1140
    ##  Correlation: 
    ##                                                     (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                    -0.624                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                         -0.940  0.606                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                    -0.673  0.357  0.615                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                              -0.179  0.082  0.041  0.041                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                          -0.054  0.035  0.068 -0.001 -0.055                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                     0.142 -0.078 -0.152 -0.076 -0.025 -0.572                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                             -0.130  0.085  0.182 -0.098 -0.107  0.247 -0.073                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                        -0.091  0.083  0.105 -0.014  0.042  0.063 -0.091  0.156                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                             -0.015  0.053  0.023  0.082  0.022  0.014 -0.510  0.010 -0.002                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                             0.083 -0.041 -0.016 -0.091 -0.040  0.057 -0.216  0.051  0.054  0.100                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                          0.083 -0.071 -0.096 -0.085 -0.043 -0.052  0.637  0.045  0.011 -0.837 -0.370                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                     -0.076  0.047  0.032  0.359  0.047  0.031 -0.037 -0.111 -0.011  0.121 -0.358 -0.010                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                    -0.046  0.018  0.050 -0.203  0.027  0.494 -0.546  0.104  0.088 -0.090 -0.281  0.103 -0.079                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                              0.034 -0.071 -0.096 -0.150 -0.042  0.021  0.006 -0.100 -0.124 -0.168  0.024  0.179 -0.054  0.169                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                               0.074 -0.038 -0.104  0.026  0.005 -0.280  0.264 -0.132 -0.206 -0.043  0.045  0.031 -0.058 -0.224 -0.176                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                              -0.505  0.473  0.536  0.013 -0.007  0.174 -0.164  0.096  0.025  0.046  0.035 -0.053 -0.001  0.112  0.018 -0.310                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                         0.607 -0.959 -0.620 -0.327 -0.023 -0.060  0.084 -0.111 -0.088 -0.051  0.009  0.067 -0.025 -0.020  0.098  0.081 -0.506                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                    0.350 -0.651 -0.319 -0.559 -0.023  0.018  0.032  0.138  0.024 -0.082  0.069  0.067 -0.207  0.112  0.138 -0.049  0.015  0.593                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                              0.101 -0.161 -0.041 -0.024 -0.518  0.078 -0.006  0.040 -0.032 -0.019  0.017  0.033 -0.021  0.010  0.026 -0.036 -0.011  0.054      0.026                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                          0.073 -0.058 -0.093  0.001  0.057 -0.736  0.448 -0.210 -0.016  0.025 -0.034  0.007 -0.008 -0.428 -0.018  0.192 -0.170  0.097     -0.037     -0.102                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                   -0.074  0.090  0.080  0.045 -0.009  0.440 -0.686  0.015  0.006  0.308  0.120 -0.391  0.051  0.410 -0.008 -0.190  0.129 -0.094     -0.067      0.018 -0.473                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                             0.080 -0.131 -0.112  0.127  0.070 -0.163  0.034 -0.736 -0.103  0.007 -0.031 -0.040  0.106 -0.068  0.084  0.092 -0.092  0.177     -0.206     -0.064  0.292           -0.025                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                        0.067 -0.107 -0.073  0.021 -0.051 -0.007  0.051 -0.081 -0.822  0.010 -0.034 -0.009  0.006 -0.051  0.072  0.194 -0.016  0.120     -0.046      0.024  0.014            0.048   0.161                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                             0.036 -0.043 -0.036 -0.065 -0.018  0.028  0.299  0.015  0.028 -0.604 -0.047  0.522 -0.068  0.064  0.113 -0.009 -0.012  0.045      0.119      0.026 -0.042           -0.542  -0.049            -0.057                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                           -0.060  0.061  0.028  0.087  0.035 -0.053  0.150 -0.014 -0.017 -0.053 -0.636  0.211  0.243  0.162 -0.062 -0.010 -0.033 -0.012     -0.065     -0.013  0.075           -0.279   0.060             0.026             0.130                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                        -0.056  0.077  0.061  0.052  0.031  0.005 -0.375 -0.041 -0.022  0.507  0.190 -0.600  0.016 -0.066 -0.125 -0.002  0.022 -0.070     -0.116     -0.034 -0.006            0.672   0.060             0.041            -0.868            -0.381                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                     0.042 -0.100 -0.015 -0.218 -0.033 -0.004  0.042  0.102 -0.001 -0.090  0.259  0.029 -0.619  0.008  0.091  0.046  0.006  0.069      0.336      0.018 -0.001           -0.148  -0.130             0.002             0.145            -0.271            -0.138                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                    0.033  0.007 -0.031  0.082  0.000 -0.458  0.441 -0.062 -0.017  0.060  0.199 -0.070 -0.029 -0.661 -0.095  0.189 -0.109  0.004     -0.221     -0.023  0.478           -0.539   0.072            -0.034            -0.083            -0.288             0.102           -0.013                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                            -0.070  0.086  0.085  0.132  0.034 -0.040  0.022  0.075  0.062  0.113 -0.075 -0.124  0.080 -0.123 -0.642  0.119 -0.023 -0.128     -0.202     -0.040  0.009            0.021  -0.084            -0.054            -0.156             0.057             0.194           -0.120            0.168                                                                                                                                   
    ## tsign-1:thermal_bias.sc                              0.040  0.004 -0.018 -0.003 -0.169  0.387 -0.150 -0.059 -0.018  0.017  0.050 -0.071  0.016  0.040  0.014 -0.023  0.033 -0.011      0.005      0.017 -0.260            0.117   0.030             0.017            -0.007            -0.040             0.037            0.005           -0.062            -0.019                                                                                                                 
    ## tsign1:thermal_bias.sc                               0.057 -0.014 -0.062 -0.002 -0.026  0.656 -0.270 -0.159 -0.102  0.027  0.045 -0.081  0.046  0.132 -0.001 -0.074  0.058 -0.005      0.007      0.066 -0.386            0.173   0.067             0.057             0.001            -0.037             0.027           -0.009           -0.130            -0.036             0.497                                                                                               
    ## temptrend_abs.sc:npp.sc                             -0.039  0.093  0.063 -0.040  0.007  0.156 -0.156  0.093  0.194  0.012 -0.007 -0.012  0.033  0.121  0.102 -0.756  0.245 -0.152      0.089      0.013 -0.234            0.201  -0.213            -0.283             0.031            -0.004            -0.005           -0.060           -0.198            -0.134             0.004  0.044                                                                                        
    ## temptrend_abs.sc:veg.sc                              0.432 -0.547 -0.453  0.007 -0.002 -0.139  0.123 -0.089 -0.017 -0.027 -0.033  0.038  0.003 -0.084 -0.027  0.241 -0.860  0.590     -0.026      0.031  0.192           -0.133   0.149             0.050            -0.005             0.042            -0.011           -0.012            0.112             0.013            -0.032 -0.062 -0.348                                                                                 
    ## human_footprint.sc:REALM2TerrFresh                  -0.127  0.174  0.130 -0.001  0.013 -0.173  0.064 -0.084  0.075  0.041 -0.021 -0.041 -0.063 -0.075 -0.007 -0.048  0.217 -0.177     -0.009     -0.034  0.140           -0.064   0.012            -0.062            -0.025             0.020             0.024            0.046            0.088             0.020            -0.008 -0.069  0.063            -0.251                                                               
    ## human_footprint.sc:REALM2Marine                     -0.007  0.000  0.004  0.000 -0.022  0.033  0.020  0.036  0.052  0.000 -0.043  0.005  0.008 -0.002 -0.070  0.039 -0.012 -0.004      0.003      0.017 -0.030           -0.025  -0.022            -0.050             0.005             0.018            -0.011            0.005            0.010             0.060             0.010  0.052 -0.035             0.010           -0.006                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc             0.028 -0.018 -0.037  0.002 -0.004 -0.351  0.183  0.028  0.039 -0.011 -0.027  0.042  0.010 -0.120 -0.027  0.021 -0.058  0.021     -0.015      0.125  0.540           -0.148  -0.016            -0.052             0.028             0.045            -0.047           -0.017            0.073             0.001            -0.481 -0.368 -0.039             0.067           -0.018      -0.027                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc              0.020  0.016 -0.033 -0.008  0.056 -0.418  0.189  0.029  0.075  0.002 -0.021  0.035 -0.002 -0.120 -0.006  0.033 -0.072  0.014     -0.010     -0.088  0.636           -0.124  -0.015            -0.101             0.007             0.022            -0.034            0.004            0.055            -0.031            -0.325 -0.614 -0.064             0.085            0.053      -0.058       0.641                      
    ## temptrend_abs.sc:human_footprint.sc:REALM2TerrFresh  0.145 -0.225 -0.150  0.006 -0.006  0.147 -0.069  0.024 -0.047 -0.032  0.014  0.033  0.053  0.082  0.011  0.037 -0.255  0.229      0.010      0.044 -0.166            0.069   0.007             0.055             0.040            -0.026            -0.036           -0.047           -0.101            -0.021            -0.005  0.049 -0.084             0.332           -0.854       0.005       0.015 -0.064               
    ## temptrend_abs.sc:human_footprint.sc:REALM2Marine    -0.007 -0.010  0.006  0.013  0.011 -0.021 -0.023 -0.009 -0.050  0.001  0.012 -0.008  0.009  0.008  0.051 -0.046  0.013  0.010     -0.011     -0.035  0.051            0.027   0.078             0.095            -0.012            -0.052             0.022            0.005            0.008            -0.059            -0.018 -0.058  0.102            -0.035            0.000      -0.766       0.032  0.070 -0.013        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -7.02311163 -0.26176466 -0.02363408  0.31203019  5.40892698 
    ## 
    ## Number of Observations: 30801
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    210                  30801

#### Full models for total and HM

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
    ##       AIC       BIC   logLik
    ##   -139247 -138847.8 69669.52
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06277468 (Intr)
    ## temptrend_abs.sc 0.04847031 -0.175
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 0.0002039828 0.3158091
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.351927 
    ## Fixed effects: Jbetatrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave.sc + temptrend_abs.sc * tempave_metab.sc +      temptrend_abs.sc * seas.sc + temptrend_abs.sc * microclim.sc +      temptrend_abs.sc * mass.sc + temptrend_abs.sc * speed.sc +      temptrend_abs.sc * lifespan.sc + temptrend_abs.sc * consumerfrac.sc +      temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.02055179 0.017533485 43221   1.172145  0.2411
    ## temptrend_abs.sc                                  0.08174432 0.022862703 43221   3.575444  0.0004
    ## REALMMarine                                       0.06280016 0.018621020   232   3.372541  0.0009
    ## REALMTerrestrial                                  0.03541199 0.020123807   232   1.759707  0.0798
    ## tsign1                                           -0.01258024 0.000517175 43221 -24.324918  0.0000
    ## tempave.sc                                       -0.00754879 0.000741422 43221 -10.181497  0.0000
    ## tempave_metab.sc                                  0.01461943 0.001713148 43221   8.533662  0.0000
    ## seas.sc                                          -0.00048603 0.000363170 43221  -1.338304  0.1808
    ## microclim.sc                                      0.00096011 0.000177505 43221   5.408898  0.0000
    ## mass.sc                                           0.00013820 0.000706775 43221   0.195539  0.8450
    ## speed.sc                                          0.00164995 0.000703435 43221   2.345566  0.0190
    ## lifespan.sc                                      -0.00195295 0.001475711 43221  -1.323394  0.1857
    ## consumerfrac.sc                                   0.00202958 0.001114234 43221   1.821500  0.0685
    ## endothermfrac.sc                                 -0.01853424 0.004401711 43221  -4.210690  0.0000
    ## nspp.sc                                          -0.00834030 0.000385732 43221 -21.621989  0.0000
    ## npp.sc                                            0.00073920 0.000287350 43221   2.572462  0.0101
    ## veg.sc                                            0.00167302 0.000366368 43221   4.566490  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.02720034 0.023575037 43221  -1.153777  0.2486
    ## temptrend_abs.sc:REALMTerrestrial                -0.00741544 0.027208680 43221  -0.272539  0.7852
    ## temptrend_abs.sc:tsign1                           0.01027236 0.001027266 43221   9.999704  0.0000
    ## temptrend_abs.sc:tempave.sc                       0.01258974 0.002030624 43221   6.199935  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                 0.00757527 0.004568951 43221   1.657989  0.0973
    ## temptrend_abs.sc:seas.sc                         -0.00407666 0.000958144 43221  -4.254745  0.0000
    ## temptrend_abs.sc:microclim.sc                    -0.00530237 0.000587355 43221  -9.027535  0.0000
    ## temptrend_abs.sc:mass.sc                         -0.00522109 0.001695750 43221  -3.078925  0.0021
    ## temptrend_abs.sc:speed.sc                        -0.00926713 0.001557226 43221  -5.951049  0.0000
    ## temptrend_abs.sc:lifespan.sc                      0.01041188 0.003667429 43221   2.839012  0.0045
    ## temptrend_abs.sc:consumerfrac.sc                  0.00101552 0.001569687 43221   0.646959  0.5177
    ## temptrend_abs.sc:endothermfrac.sc                 0.00892540 0.007069432 43221   1.262535  0.2068
    ## temptrend_abs.sc:nspp.sc                         -0.01483113 0.000890388 43221 -16.656921  0.0000
    ## tsign-1:thermal_bias.sc                          -0.00074413 0.000565287 43221  -1.316372  0.1881
    ## tsign1:thermal_bias.sc                           -0.00247830 0.000462769 43221  -5.355387  0.0000
    ## temptrend_abs.sc:npp.sc                          -0.00121354 0.000926209 43221  -1.310218  0.1901
    ## temptrend_abs.sc:veg.sc                          -0.00372755 0.001134713 43221  -3.285015  0.0010
    ## human_bowler.sc:REALM2TerrFresh                   0.00027144 0.000299147 43221   0.907388  0.3642
    ## human_bowler.sc:REALM2Marine                      0.00031037 0.000193119 43221   1.607166  0.1080
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.00813092 0.001167607 43221   6.963746  0.0000
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00473688 0.001134858 43221   4.173983  0.0000
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.00156197 0.000754427 43221  -2.070401  0.0384
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.00167583 0.000786454 43221  -2.130864  0.0331
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.430                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.929  0.402                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.876  0.377  0.815                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.026  0.017  0.006  0.009                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.003  0.000 -0.004 -0.005 -0.070                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.038 -0.020 -0.030 -0.031 -0.022 -0.454                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.030  0.037  0.040 -0.012 -0.107  0.154 -0.084                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.008  0.018  0.014 -0.006  0.020 -0.060 -0.022  0.127                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                           0.012  0.006 -0.002 -0.002  0.019  0.021 -0.486  0.007  0.056                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.000 -0.005  0.009  0.017 -0.049  0.106 -0.277  0.102  0.072  0.129                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.024 -0.019 -0.019 -0.035 -0.015 -0.086  0.707  0.005 -0.027 -0.760 -0.455                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.003  0.005  0.024  0.252  0.015 -0.044 -0.025 -0.051 -0.001  0.014 -0.049 -0.047                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                  0.133 -0.056 -0.068 -0.269  0.003  0.099 -0.111  0.012  0.009  0.001 -0.142  0.068 -0.328                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                          -0.007 -0.011 -0.009 -0.019 -0.038  0.002 -0.023 -0.060 -0.115 -0.208  0.014  0.134 -0.009  0.067                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.004  0.004 -0.005  0.006  0.013 -0.245  0.196 -0.088 -0.090 -0.044  0.048  0.035 -0.029 -0.049 -0.182                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.075  0.119  0.073 -0.004 -0.016  0.141 -0.098  0.121 -0.019  0.034  0.052 -0.029 -0.069 -0.012  0.002 -0.271                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.413 -0.960 -0.415 -0.359 -0.004  0.003  0.015 -0.046 -0.022 -0.010 -0.005  0.014 -0.010  0.036  0.020  0.008 -0.129                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.363 -0.858 -0.336 -0.433 -0.005 -0.011  0.020  0.026  0.010 -0.008 -0.006  0.023 -0.100  0.135  0.031 -0.018  0.005  0.814                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.013 -0.037 -0.005 -0.004 -0.484  0.050 -0.002  0.013 -0.041 -0.015  0.019  0.008 -0.004  0.006  0.025 -0.048 -0.004  0.012      0.007                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.000 -0.005 -0.002 -0.002  0.052 -0.625  0.315 -0.052  0.083  0.028 -0.064  0.031 -0.009 -0.087 -0.016  0.157 -0.140 -0.001      0.017     -0.076                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.009  0.023  0.008  0.010 -0.012  0.335 -0.597 -0.010 -0.034  0.244  0.069 -0.350  0.049  0.105 -0.006 -0.141  0.060 -0.019     -0.030      0.015 -0.420                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.016 -0.050 -0.020  0.021  0.044  0.003 -0.008 -0.647 -0.088  0.007 -0.037 -0.005  0.041  0.019  0.040  0.061 -0.127  0.068     -0.052     -0.018  0.060            0.031                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.006 -0.026 -0.009  0.007 -0.044  0.089 -0.014 -0.053 -0.693 -0.028 -0.029  0.018  0.002  0.009  0.045  0.103  0.024  0.030     -0.023      0.047 -0.119            0.104   0.139                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.006 -0.010 -0.007 -0.006 -0.015  0.030  0.262  0.008 -0.017 -0.555 -0.024  0.434 -0.005  0.014  0.131  0.005 -0.005  0.016      0.018      0.026 -0.073           -0.469  -0.045            -0.010                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.004  0.018 -0.003 -0.001  0.026 -0.096  0.128 -0.035 -0.041 -0.044 -0.547  0.187  0.022  0.079 -0.024  0.006 -0.057  0.001      0.000     -0.006  0.161           -0.223   0.074             0.046             0.063                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.010  0.026  0.009  0.012  0.017  0.042 -0.383 -0.008  0.014  0.425  0.140 -0.531  0.045 -0.028 -0.104 -0.019 -0.002 -0.017     -0.033     -0.018 -0.067            0.687   0.031             0.000            -0.794            -0.311                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.001 -0.070 -0.008 -0.096 -0.010  0.003  0.055  0.035  0.002 -0.014  0.053  0.053 -0.490  0.123  0.010  0.029  0.040  0.075      0.301      0.006  0.018           -0.120  -0.038            -0.003             0.038            -0.080            -0.111                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                -0.045  0.150  0.029  0.125  0.005 -0.159  0.185  0.014  0.026  0.017  0.119 -0.030  0.094 -0.401 -0.055  0.069 -0.017 -0.100     -0.349     -0.017  0.181           -0.303  -0.033            -0.063            -0.017            -0.164             0.010           -0.296                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.009  0.024  0.011  0.020  0.023 -0.025  0.026  0.032  0.046  0.132 -0.029 -0.101  0.006 -0.057 -0.544  0.084 -0.031 -0.034     -0.056     -0.037 -0.006            0.027  -0.045            -0.021            -0.205             0.020             0.186           -0.014            0.095                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.015 -0.007 -0.017 -0.003 -0.170  0.475 -0.110 -0.101 -0.071  0.028  0.080 -0.089 -0.005  0.003  0.005 -0.079  0.042  0.012     -0.008      0.019 -0.281            0.083   0.084             0.047            -0.003            -0.066             0.043           -0.003           -0.005            -0.011                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.022 -0.014 -0.030 -0.008 -0.036  0.711 -0.163 -0.233 -0.165  0.045  0.061 -0.087 -0.010  0.015 -0.035 -0.123  0.053  0.017     -0.004      0.070 -0.360            0.107   0.115             0.066            -0.002            -0.063             0.039           -0.002           -0.027            -0.012             0.572                                                                                               
    ## temptrend_abs.sc:npp.sc                           0.008  0.006 -0.006 -0.015 -0.011  0.084 -0.069  0.035  0.112  0.012  0.001 -0.012  0.009  0.009  0.066 -0.654  0.186 -0.025      0.033      0.036 -0.190            0.152  -0.193            -0.252             0.017            -0.034             0.026           -0.038           -0.065            -0.090             0.028  0.059                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.053 -0.152 -0.054  0.005  0.009 -0.094  0.050 -0.116  0.021 -0.012 -0.048  0.012  0.026  0.002 -0.025  0.191 -0.761  0.167     -0.014      0.022  0.169           -0.063   0.208             0.015            -0.012             0.080             0.014           -0.034            0.017             0.027            -0.034 -0.052 -0.320                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.027  0.050  0.018 -0.011  0.037 -0.283  0.098 -0.058  0.126  0.049 -0.024 -0.054 -0.066 -0.077  0.019 -0.065  0.193 -0.054      0.014     -0.032  0.163           -0.090  -0.148            -0.037            -0.021             0.016             0.014            0.023            0.063             0.035            -0.083 -0.100  0.054            -0.263                                                               
    ## human_bowler.sc:REALM2Marine                      0.004 -0.010 -0.004 -0.001 -0.006  0.032 -0.003 -0.056 -0.074 -0.025  0.003  0.022  0.002  0.004  0.027 -0.209  0.048  0.010      0.004      0.015 -0.001           -0.017   0.088             0.037             0.023            -0.012            -0.011            0.007            0.002             0.017             0.058  0.066  0.129            -0.026            0.002                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.005  0.005  0.004 -0.003  0.018 -0.381  0.119  0.082  0.079 -0.007 -0.058  0.050 -0.006 -0.014 -0.014  0.067 -0.077 -0.018      0.016      0.102  0.692           -0.122  -0.099            -0.096             0.006             0.129            -0.093            0.008           -0.007            -0.032            -0.473 -0.382 -0.088             0.092            0.060      -0.048                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.006  0.010  0.005 -0.001  0.067 -0.414  0.115  0.089  0.110 -0.002 -0.047  0.043 -0.004 -0.011 -0.001  0.075 -0.087 -0.018      0.014     -0.124  0.732           -0.098  -0.081            -0.130            -0.008             0.109            -0.080            0.008           -0.012            -0.047            -0.346 -0.524 -0.103             0.109            0.054      -0.054       0.797                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.026 -0.086 -0.025  0.009 -0.019  0.192 -0.090 -0.111 -0.066 -0.029  0.010  0.029  0.022  0.048  0.018  0.034 -0.312  0.093     -0.022      0.054 -0.238            0.123   0.228             0.107             0.041            -0.027            -0.028           -0.025           -0.086            -0.034             0.061  0.072 -0.106             0.457           -0.683      -0.002      -0.108 -0.110               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.005  0.015  0.004  0.002  0.003 -0.017 -0.005  0.078  0.025  0.016 -0.012 -0.014  0.006 -0.001 -0.007  0.099 -0.004 -0.015      0.000     -0.006 -0.037            0.060  -0.180            -0.082            -0.018             0.010             0.005           -0.011           -0.027            -0.046            -0.044 -0.059 -0.155             0.022            0.011      -0.615       0.049  0.049 -0.001        
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -7.3991424 -0.2444604  0.1291713  0.6077538  6.7430991 
    ## 
    ## Number of Observations: 43493
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    235                  43493

``` r
summary(modTfullHorn)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -125404.7 -125006.5 62748.35
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06157166 (Intr)
    ## temptrend_abs.sc 0.04045651 -0.087
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.004171557 0.2960534
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.223338 
    ## Fixed effects: Horntrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave.sc + temptrend_abs.sc * tempave_metab.sc +      temptrend_abs.sc * seas.sc + temptrend_abs.sc * microclim.sc +      temptrend_abs.sc * mass.sc + temptrend_abs.sc * speed.sc +      temptrend_abs.sc * lifespan.sc + temptrend_abs.sc * consumerfrac.sc +      temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.01104912 0.017196483 42254   0.642522  0.5205
    ## temptrend_abs.sc                                  0.07927630 0.020905866 42254   3.792060  0.0001
    ## REALMMarine                                       0.05780589 0.018467290   199   3.130177  0.0020
    ## REALMTerrestrial                                  0.04259399 0.020090957   199   2.120058  0.0352
    ## tsign1                                           -0.01011127 0.000653242 42254 -15.478591  0.0000
    ## tempave.sc                                       -0.00361185 0.001028626 42254  -3.511336  0.0004
    ## tempave_metab.sc                                  0.00655184 0.002310628 42254   2.835523  0.0046
    ## seas.sc                                          -0.00274154 0.000541223 42254  -5.065448  0.0000
    ## microclim.sc                                      0.00072139 0.000281288 42254   2.564581  0.0103
    ## mass.sc                                          -0.00054617 0.000912681 42254  -0.598425  0.5496
    ## speed.sc                                         -0.00140467 0.000895023 42254  -1.569424  0.1166
    ## lifespan.sc                                      -0.00188846 0.001928273 42254  -0.979352  0.3274
    ## consumerfrac.sc                                   0.00219813 0.001297007 42254   1.694769  0.0901
    ## endothermfrac.sc                                 -0.01412500 0.004972945 42254  -2.840369  0.0045
    ## nspp.sc                                          -0.00856226 0.000507427 42254 -16.873893  0.0000
    ## npp.sc                                           -0.00010123 0.000436009 42254  -0.232167  0.8164
    ## veg.sc                                            0.00067701 0.000556803 42254   1.215887  0.2240
    ## temptrend_abs.sc:REALMMarine                     -0.02179149 0.021656945 42254  -1.006213  0.3143
    ## temptrend_abs.sc:REALMTerrestrial                -0.05819644 0.025301200 42254  -2.300146  0.0214
    ## temptrend_abs.sc:tsign1                           0.00730497 0.001147993 42254   6.363251  0.0000
    ## temptrend_abs.sc:tempave.sc                       0.01278689 0.002317696 42254   5.517071  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                 0.00683279 0.005235099 42254   1.305189  0.1918
    ## temptrend_abs.sc:seas.sc                          0.00497334 0.001135970 42254   4.378054  0.0000
    ## temptrend_abs.sc:microclim.sc                    -0.00120060 0.000726011 42254  -1.653692  0.0982
    ## temptrend_abs.sc:mass.sc                         -0.00111234 0.001928691 42254  -0.576732  0.5641
    ## temptrend_abs.sc:speed.sc                        -0.00264311 0.001721359 42254  -1.535477  0.1247
    ## temptrend_abs.sc:lifespan.sc                      0.00289752 0.004139239 42254   0.700012  0.4839
    ## temptrend_abs.sc:consumerfrac.sc                  0.00003823 0.001726789 42254   0.022141  0.9823
    ## temptrend_abs.sc:endothermfrac.sc                 0.00706745 0.007077812 42254   0.998536  0.3180
    ## temptrend_abs.sc:nspp.sc                         -0.00955301 0.001015550 42254  -9.406728  0.0000
    ## tsign-1:thermal_bias.sc                          -0.00122716 0.000712091 42254  -1.723316  0.0848
    ## tsign1:thermal_bias.sc                           -0.00026068 0.000609565 42254  -0.427647  0.6689
    ## temptrend_abs.sc:npp.sc                           0.00255364 0.001103631 42254   2.313855  0.0207
    ## temptrend_abs.sc:veg.sc                          -0.00085135 0.001321141 42254  -0.644408  0.5193
    ## human_bowler.sc:REALM2TerrFresh                   0.00033710 0.000468005 42254   0.720299  0.4713
    ## human_bowler.sc:REALM2Marine                      0.00152806 0.000384373 42254   3.975450  0.0001
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.00518868 0.001304502 42254   3.977522  0.0001
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00190649 0.001268519 42254   1.502923  0.1329
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.00077440 0.000876091 42254   0.883922  0.3767
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.00110828 0.001011989 42254  -1.095151  0.2735
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.380                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.917  0.352                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.858  0.326  0.781                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.032  0.020  0.007  0.013                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.007  0.003  0.004 -0.005 -0.074                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.050 -0.023 -0.042 -0.039 -0.022 -0.506                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.043  0.037  0.059 -0.014 -0.114  0.205 -0.112                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.017  0.022  0.025 -0.007  0.006 -0.034 -0.034  0.148                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                           0.012  0.010 -0.001  0.003  0.021  0.023 -0.485 -0.004  0.023                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.005 -0.005  0.011  0.013 -0.049  0.105 -0.269  0.097  0.063  0.148                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.029 -0.020 -0.024 -0.040 -0.015 -0.087  0.691  0.018  0.003 -0.778 -0.448                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.004  0.007  0.002  0.268  0.019 -0.026 -0.021 -0.047 -0.006  0.014 -0.062 -0.040                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                  0.142 -0.048 -0.069 -0.299  0.004  0.153 -0.163  0.037  0.018 -0.011 -0.145  0.060 -0.317                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                          -0.010 -0.017 -0.009 -0.023 -0.031  0.016 -0.021 -0.019 -0.081 -0.222  0.000  0.149 -0.004  0.083                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.015  0.004 -0.018  0.009  0.011 -0.292  0.228 -0.205 -0.224 -0.034  0.046  0.028 -0.014 -0.083 -0.195                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.110  0.123  0.112 -0.003 -0.017  0.148 -0.106  0.135  0.022  0.028  0.040 -0.022 -0.049  0.013  0.014 -0.276                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.364 -0.957 -0.362 -0.306 -0.006 -0.005  0.019 -0.048 -0.025 -0.014 -0.007  0.015 -0.008  0.032  0.025  0.011 -0.136                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.313 -0.840 -0.285 -0.387 -0.006 -0.008  0.025  0.041  0.012 -0.018 -0.004  0.027 -0.085  0.133  0.041 -0.024  0.013  0.789                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.016 -0.044 -0.006 -0.005 -0.489  0.042  0.003  0.015 -0.036 -0.013  0.018  0.007 -0.006  0.004  0.024 -0.046  0.003  0.014      0.009                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.004 -0.009 -0.008 -0.002  0.050 -0.571  0.306 -0.056  0.076  0.025 -0.061  0.030 -0.011 -0.099 -0.017  0.144 -0.130  0.008      0.018     -0.068                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.014  0.028  0.013  0.016 -0.010  0.316 -0.574 -0.007 -0.017  0.248  0.078 -0.347  0.050  0.120 -0.005 -0.142  0.057 -0.021     -0.050      0.012 -0.434                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.018 -0.056 -0.025  0.029  0.045 -0.031  0.006 -0.606 -0.101  0.009 -0.035 -0.011  0.032  0.002  0.016  0.109 -0.131  0.076     -0.080     -0.020  0.034            0.052                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.009 -0.026 -0.012  0.007 -0.038  0.064 -0.005 -0.075 -0.660 -0.008 -0.026  0.001  0.003  0.001  0.033  0.155  0.021  0.030     -0.029      0.050 -0.144            0.103   0.138                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.007 -0.017 -0.008 -0.011 -0.013  0.021  0.260  0.005 -0.008 -0.552 -0.034  0.435 -0.012  0.019  0.131  0.007 -0.001  0.021      0.033      0.025 -0.080           -0.465  -0.040            -0.008                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.004  0.020 -0.005 -0.001  0.025 -0.081  0.123 -0.028 -0.034 -0.047 -0.533  0.186  0.020  0.084 -0.017 -0.001 -0.047  0.004     -0.002     -0.004  0.162           -0.216   0.064             0.039             0.058                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.011  0.031  0.011  0.016  0.015  0.038 -0.370 -0.011  0.005  0.434  0.151 -0.531  0.054 -0.029 -0.108 -0.018 -0.006 -0.018     -0.046     -0.017 -0.058            0.676   0.034            -0.001            -0.803            -0.298                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.001 -0.070  0.003 -0.082 -0.012 -0.006  0.054  0.036  0.009 -0.017  0.046  0.057 -0.488  0.098  0.009  0.005  0.025  0.073      0.262      0.008  0.039           -0.130  -0.046            -0.013             0.048            -0.069            -0.131                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                -0.039  0.149  0.026  0.127  0.003 -0.176  0.206  0.005  0.017  0.025  0.123 -0.032  0.079 -0.372 -0.064  0.082 -0.028 -0.104     -0.365     -0.017  0.221           -0.344  -0.046            -0.070            -0.032            -0.180             0.019           -0.242                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.011  0.034  0.012  0.026  0.015 -0.027  0.021  0.003  0.032  0.134 -0.021 -0.108  0.011 -0.063 -0.527  0.087 -0.034 -0.042     -0.073     -0.033 -0.006            0.035  -0.028            -0.010            -0.210             0.016             0.200           -0.021            0.106                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.021 -0.008 -0.021 -0.005 -0.169  0.469 -0.119 -0.101 -0.084  0.039  0.081 -0.098  0.006  0.007 -0.007 -0.092  0.043  0.011     -0.009      0.025 -0.266            0.083   0.083             0.050            -0.008            -0.065             0.046           -0.013           -0.004             0.001                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.034 -0.016 -0.040 -0.012 -0.034  0.651 -0.145 -0.222 -0.195  0.059  0.061 -0.099  0.008  0.016 -0.039 -0.121  0.043  0.015     -0.003      0.072 -0.316            0.092   0.100             0.078            -0.011            -0.055             0.042           -0.013           -0.021             0.000             0.590                                                                                               
    ## temptrend_abs.sc:npp.sc                           0.009  0.001 -0.005 -0.018 -0.012  0.110 -0.090  0.079  0.155  0.007 -0.004 -0.009  0.002  0.028  0.076 -0.614  0.165 -0.023      0.040      0.038 -0.180            0.163  -0.204            -0.281             0.017            -0.035             0.027           -0.009           -0.074            -0.089             0.034  0.058                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.065 -0.186 -0.067  0.010  0.010 -0.103  0.052 -0.128  0.020 -0.005 -0.039  0.001  0.019 -0.010 -0.031  0.181 -0.639  0.207     -0.028      0.020  0.186           -0.064   0.235            -0.019            -0.020             0.082             0.024           -0.025            0.029             0.034            -0.034 -0.048 -0.306                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.026  0.045  0.021 -0.012  0.033 -0.255  0.111 -0.059  0.152  0.043 -0.021 -0.041 -0.051 -0.071  0.012 -0.081  0.131 -0.049      0.018     -0.021  0.133           -0.076  -0.122            -0.034            -0.017             0.012             0.010            0.015            0.059             0.032            -0.072 -0.075  0.043            -0.204                                                               
    ## human_bowler.sc:REALM2Marine                      0.011 -0.018 -0.013 -0.003  0.007 -0.024  0.025 -0.146 -0.101 -0.018 -0.006  0.002  0.001  0.003 -0.033 -0.137  0.023  0.023      0.006      0.015  0.025           -0.025   0.138             0.026             0.019            -0.010             0.001            0.013           -0.001             0.048             0.084  0.078  0.079             0.004            0.011                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.005  0.002  0.003 -0.002  0.022 -0.332  0.106  0.084  0.083 -0.011 -0.055  0.051 -0.013 -0.010  0.000  0.047 -0.069 -0.012      0.015      0.099  0.689           -0.127  -0.114            -0.116             0.001             0.132            -0.088            0.034           -0.001            -0.046            -0.468 -0.357 -0.073             0.102            0.047      -0.061                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.008  0.008  0.005  0.001  0.069 -0.358  0.096  0.088  0.113 -0.008 -0.045  0.046 -0.013 -0.006  0.009  0.058 -0.076 -0.011      0.009     -0.134  0.715           -0.097  -0.078            -0.148            -0.011             0.110            -0.076            0.033           -0.011            -0.058            -0.341 -0.497 -0.087             0.121            0.034      -0.067       0.798                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.028 -0.105 -0.029  0.015 -0.016  0.158 -0.087 -0.122 -0.053 -0.026  0.011  0.021  0.016  0.045  0.017  0.042 -0.262  0.114     -0.040      0.050 -0.229            0.133   0.282             0.078             0.037            -0.027            -0.020           -0.020           -0.096            -0.035             0.053  0.052 -0.099             0.479           -0.561       0.006      -0.106 -0.095               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.010  0.025  0.012  0.003  0.003  0.013 -0.020  0.120  0.035  0.015 -0.012 -0.004  0.008 -0.001  0.024  0.072  0.008 -0.032     -0.004     -0.010 -0.050            0.063  -0.219            -0.062            -0.020             0.013            -0.001           -0.021           -0.020            -0.068            -0.062 -0.065 -0.135            -0.002            0.006      -0.658       0.067  0.066 -0.013        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.73851930 -0.38297739  0.03926292  0.54007605  6.19209050 
    ## 
    ## Number of Observations: 42493
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    202                  42493

#### Full models for total and HM (remove year 1)

``` r
i2 <- trends[, complete.cases(Jbetatrendrem0, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, human_bowler.sc)]
i3 <- trends[, complete.cases(Horntrendrem0, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, human_bowler.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# full models
if(file.exists('temp/modTfullJbetarem0.rds')){
  modTfullJbetarem0 <- readRDS('temp/modTfullJbetarem0.rds')
} else {
  modTfullJbetarem0 <- lme(Jbetatrendrem0 ~ temptrend_abs.sc*REALM + 
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
  saveRDS(modTfullJbetarem0, file = 'temp/modTfullJbetarem0.rds')
}

if(file.exists('temp/modTfullHornrem0.rds')){
  modTfullHornrem0 <- readRDS('temp/modTfullHornrem0.rds')
} else {
  modTfullHornrem0 <- lme(Horntrendrem0 ~ temptrend_abs.sc*REALM + 
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
  saveRDS(modTfullHornrem0, file = 'temp/modTfullHornrem0.rds')
}

summary(modTfullJbetarem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -107160.4 -106776.6 53626.19
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.007723766 (Intr)
    ## temptrend_abs.sc 0.017251054 0.308 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.002065633 0.8886333
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.821858 
    ## Fixed effects: Jbetatrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave.sc + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * lifespan.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.000118264 0.003494525 30817  0.033843  0.9730
    ## temptrend_abs.sc                                  0.031047178 0.010641227 30817  2.917631  0.0035
    ## REALMMarine                                       0.008686766 0.003676059   215  2.363065  0.0190
    ## REALMTerrestrial                                  0.004213563 0.003972879   215  1.060582  0.2901
    ## tsign1                                           -0.001650156 0.000563878 30817 -2.926443  0.0034
    ## tempave.sc                                        0.000026540 0.000774039 30817  0.034288  0.9726
    ## tempave_metab.sc                                 -0.006656000 0.001743847 30817 -3.816849  0.0001
    ## seas.sc                                           0.000547462 0.000373125 30817  1.467234  0.1423
    ## microclim.sc                                      0.000512306 0.000220708 30817  2.321190  0.0203
    ## mass.sc                                          -0.000060974 0.000720367 30817 -0.084644  0.9325
    ## speed.sc                                         -0.000271284 0.000655030 30817 -0.414156  0.6788
    ## lifespan.sc                                      -0.002371749 0.001412585 30817 -1.679013  0.0932
    ## consumerfrac.sc                                  -0.000094853 0.000327142 30817 -0.289945  0.7719
    ## endothermfrac.sc                                  0.004347454 0.001715705 30817  2.533918  0.0113
    ## nspp.sc                                          -0.001072200 0.000393780 30817 -2.722841  0.0065
    ## npp.sc                                           -0.000782973 0.000332241 30817 -2.356640  0.0184
    ## veg.sc                                            0.000955358 0.000414190 30817  2.306568  0.0211
    ## temptrend_abs.sc:REALMMarine                     -0.029609940 0.011102993 30817 -2.666843  0.0077
    ## temptrend_abs.sc:REALMTerrestrial                -0.007754553 0.011861403 30817 -0.653764  0.5133
    ## temptrend_abs.sc:tsign1                          -0.000889873 0.001630157 30817 -0.545882  0.5852
    ## temptrend_abs.sc:tempave.sc                      -0.001592587 0.002633115 30817 -0.604830  0.5453
    ## temptrend_abs.sc:tempave_metab.sc                 0.014037203 0.005955256 30817  2.357112  0.0184
    ## temptrend_abs.sc:seas.sc                         -0.002460440 0.001245871 30817 -1.974876  0.0483
    ## temptrend_abs.sc:microclim.sc                    -0.002606973 0.000797090 30817 -3.270614  0.0011
    ## temptrend_abs.sc:mass.sc                          0.001706338 0.002389242 30817  0.714175  0.4751
    ## temptrend_abs.sc:speed.sc                        -0.000871498 0.002097323 30817 -0.415529  0.6778
    ## temptrend_abs.sc:lifespan.sc                     -0.000446948 0.004837709 30817 -0.092388  0.9264
    ## temptrend_abs.sc:consumerfrac.sc                 -0.000263869 0.000926039 30817 -0.284944  0.7757
    ## temptrend_abs.sc:endothermfrac.sc                -0.006695729 0.005346416 30817 -1.252377  0.2104
    ## temptrend_abs.sc:nspp.sc                          0.002299914 0.001248942 30817  1.841489  0.0656
    ## tsign-1:thermal_bias.sc                          -0.000828663 0.000583994 30817 -1.418958  0.1559
    ## tsign1:thermal_bias.sc                           -0.000442575 0.000455332 30817 -0.971982  0.3311
    ## temptrend_abs.sc:npp.sc                           0.004466331 0.001232879 30817  3.622684  0.0003
    ## temptrend_abs.sc:veg.sc                          -0.004986820 0.001522925 30817 -3.274501  0.0011
    ## human_bowler.sc:REALM2TerrFresh                   0.000708740 0.000326658 30817  2.169669  0.0300
    ## human_bowler.sc:REALM2Marine                     -0.000153042 0.000289829 30817 -0.528042  0.5975
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.000059901 0.001562561 30817 -0.038335  0.9694
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.002330171 0.001449163 30817  1.607942  0.1079
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.001933907 0.001040743 30817 -1.858198  0.0632
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine     0.001272455 0.001093840 30817  1.163292  0.2447
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.478                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.934  0.465                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.741  0.288  0.680                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.146  0.068  0.034  0.026                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.017  0.007  0.026  0.004 -0.065                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.112 -0.049 -0.117 -0.071 -0.010 -0.547                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.110  0.070  0.160 -0.073 -0.106  0.226 -0.089                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.080  0.087  0.098 -0.015  0.037  0.017 -0.069  0.169                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                          -0.007  0.052  0.022  0.050  0.018  0.006 -0.496  0.013  0.015                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.055 -0.033  0.000 -0.033 -0.049  0.067 -0.217  0.066  0.062  0.112                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.077 -0.060 -0.084 -0.079 -0.032 -0.048  0.639  0.030 -0.003 -0.825 -0.385                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.071  0.042  0.040  0.330  0.042  0.014 -0.042 -0.100 -0.023  0.069 -0.267 -0.032                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                  0.029 -0.005 -0.009 -0.257  0.023  0.415 -0.481  0.093  0.049 -0.088 -0.300  0.111 -0.140                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                           0.023 -0.070 -0.086 -0.114 -0.036  0.016 -0.007 -0.085 -0.122 -0.171  0.015  0.167 -0.033  0.165                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.025 -0.006 -0.048  0.037  0.009 -0.246  0.231 -0.094 -0.189 -0.032  0.046  0.029 -0.034 -0.184 -0.168                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.381  0.387  0.407  0.005 -0.005  0.139 -0.130  0.096  0.032  0.046  0.048 -0.054 -0.017  0.065  0.008 -0.286                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.464 -0.954 -0.469 -0.267 -0.020 -0.026  0.050 -0.091 -0.092 -0.053  0.004  0.058 -0.025  0.011  0.096  0.047 -0.417                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.283 -0.720 -0.262 -0.428 -0.011  0.000  0.026  0.119  0.026 -0.060  0.040  0.051 -0.152  0.098  0.114 -0.055  0.019  0.665                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.082 -0.138 -0.035 -0.012 -0.510  0.078 -0.015  0.029 -0.033 -0.017  0.018  0.024 -0.021  0.011  0.022 -0.043 -0.021  0.049      0.014                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.038 -0.019 -0.054 -0.007  0.062 -0.717  0.425 -0.187  0.013  0.028 -0.033  0.006 -0.010 -0.354 -0.023  0.183 -0.133  0.045     -0.012     -0.101                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.045  0.059  0.048  0.037 -0.019  0.425 -0.673  0.028 -0.006  0.293  0.089 -0.371  0.066  0.370 -0.006 -0.172  0.095 -0.060     -0.056      0.029 -0.454                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.068 -0.098 -0.097  0.099  0.064 -0.129  0.041 -0.738 -0.125  0.002 -0.026 -0.027  0.089 -0.059  0.069  0.083 -0.095  0.138     -0.175     -0.045  0.263           -0.042                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.072 -0.116 -0.083  0.019 -0.046  0.023  0.035 -0.104 -0.813 -0.004 -0.036  0.002  0.009 -0.025  0.074  0.192 -0.052  0.126     -0.043      0.031 -0.023            0.058   0.170                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.033 -0.041 -0.035 -0.047 -0.016  0.030  0.289  0.010  0.017 -0.596 -0.044  0.507 -0.051  0.058  0.114 -0.015 -0.014  0.049      0.086      0.021 -0.045           -0.521  -0.043            -0.039                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.051  0.053  0.023  0.056  0.037 -0.060  0.131 -0.016 -0.030 -0.061 -0.608  0.192  0.182  0.164 -0.045  0.001 -0.042 -0.007     -0.034     -0.014  0.090           -0.261   0.070             0.040             0.125                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.045  0.071  0.051  0.041  0.026  0.010 -0.369 -0.031 -0.008  0.499  0.162 -0.579  0.047 -0.054 -0.125 -0.011  0.025 -0.064     -0.092     -0.025 -0.018            0.670   0.048             0.021            -0.853            -0.367                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.035 -0.105 -0.016 -0.156 -0.029 -0.001  0.054  0.089  0.008 -0.056  0.208  0.049 -0.550  0.004  0.066  0.036  0.009  0.082      0.318      0.016  0.004           -0.154  -0.107            -0.008             0.112            -0.235            -0.147                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                 0.012  0.062 -0.003  0.072  0.001 -0.407  0.415 -0.060  0.009  0.058  0.208 -0.065 -0.021 -0.574 -0.093  0.157 -0.065 -0.046     -0.249     -0.026  0.407           -0.517   0.065            -0.059            -0.072            -0.297             0.086           -0.053                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.059  0.073  0.076  0.095  0.030 -0.043  0.029  0.058  0.062  0.112 -0.053 -0.120  0.054 -0.121 -0.625  0.099 -0.018 -0.117     -0.155     -0.035  0.002            0.009  -0.067            -0.054            -0.157             0.033             0.187           -0.078            0.166                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.042 -0.004 -0.033  0.005 -0.156  0.419 -0.149 -0.077 -0.039  0.017  0.055 -0.069  0.021  0.025  0.010 -0.043  0.034  0.002     -0.011      0.010 -0.274            0.111   0.055             0.032            -0.006            -0.043             0.034           -0.005           -0.043            -0.016                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.064 -0.016 -0.073  0.000 -0.025  0.666 -0.245 -0.201 -0.129  0.025  0.046 -0.073  0.033  0.100 -0.014 -0.091  0.063  0.002     -0.011      0.062 -0.382            0.156   0.111             0.076             0.002            -0.044             0.027           -0.010           -0.100            -0.030             0.523                                                                                               
    ## temptrend_abs.sc:npp.sc                          -0.009  0.057  0.031 -0.044  0.005  0.130 -0.128  0.064  0.193  0.007 -0.002 -0.015  0.022  0.084  0.087 -0.766  0.251 -0.107      0.085      0.023 -0.212            0.175  -0.192            -0.286             0.029            -0.008             0.013           -0.054           -0.154            -0.108             0.017  0.059                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.328 -0.443 -0.347  0.011 -0.008 -0.100  0.089 -0.088 -0.042 -0.029 -0.045  0.041  0.007 -0.042 -0.020  0.236 -0.866  0.481     -0.032      0.039  0.143           -0.091   0.145             0.096             0.002             0.051            -0.017           -0.014            0.059             0.004            -0.027 -0.063 -0.351                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.097  0.136  0.096 -0.035  0.040 -0.283  0.117 -0.075  0.162  0.056 -0.016 -0.059 -0.041 -0.128  0.000 -0.087  0.221 -0.145      0.036     -0.044  0.182           -0.110  -0.049            -0.121            -0.022             0.002             0.025            0.018            0.133             0.043            -0.059 -0.093  0.089            -0.273                                                               
    ## human_bowler.sc:REALM2Marine                      0.036 -0.030 -0.044 -0.026  0.005  0.003  0.017 -0.154 -0.060 -0.015 -0.002 -0.001 -0.021  0.015 -0.022 -0.138  0.033  0.037      0.018      0.013  0.002           -0.024   0.133             0.032             0.007            -0.010             0.012            0.030            0.007             0.052             0.083  0.104  0.106            -0.026            0.024                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.009  0.000 -0.017 -0.005 -0.008 -0.356  0.164  0.044  0.052 -0.008 -0.029  0.035 -0.006 -0.084 -0.028  0.044 -0.048 -0.004      0.005      0.124  0.578           -0.126  -0.033            -0.063             0.022             0.060            -0.052           -0.002            0.030            -0.013            -0.476 -0.378 -0.059             0.050            0.023      -0.065                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.006  0.020 -0.018 -0.013  0.051 -0.410  0.167  0.064  0.091  0.003 -0.020  0.030 -0.010 -0.082 -0.007  0.061 -0.074 -0.003      0.014     -0.085  0.662           -0.103  -0.047            -0.114             0.003             0.040            -0.042            0.011            0.014            -0.048            -0.336 -0.595 -0.088             0.081            0.043      -0.090       0.690                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.124 -0.202 -0.130  0.029 -0.031  0.224 -0.112 -0.001 -0.131 -0.043  0.002  0.045  0.023  0.121  0.023  0.064 -0.304  0.209     -0.031      0.063 -0.247            0.126   0.029             0.176             0.047            -0.013            -0.046           -0.020           -0.157            -0.038             0.040  0.065 -0.124             0.387           -0.813      -0.020      -0.052 -0.086               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.030  0.035  0.036  0.021 -0.008 -0.001 -0.018  0.139  0.041  0.007 -0.011  0.007  0.026  0.000  0.020  0.108 -0.022 -0.042     -0.015     -0.013 -0.032            0.054  -0.167            -0.077             0.005             0.012            -0.023           -0.036           -0.041            -0.066            -0.061 -0.086 -0.159             0.045           -0.019      -0.740       0.046  0.062  0.027        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.38808782 -0.33122885 -0.03014344  0.33937884  8.28289240 
    ## 
    ## Number of Observations: 31072
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    218                  31072

``` r
summary(modTfullHornrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -81320.49 -80937.74 40706.24
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.01569909 (Intr)
    ## temptrend_abs.sc 0.02631579 0.237 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01343706 2.039905
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.179888 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave.sc + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * lifespan.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * endothermfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.00307249 0.007259783 30165  0.423221  0.6721
    ## temptrend_abs.sc                                  0.04642807 0.017215537 30165  2.696870  0.0070
    ## REALMMarine                                       0.00738080 0.007719060   187  0.956179  0.3402
    ## REALMTerrestrial                                  0.00740448 0.007847874   187  0.943502  0.3466
    ## tsign1                                           -0.00137424 0.000875682 30165 -1.569342  0.1166
    ## tempave.sc                                       -0.00079975 0.001619533 30165 -0.493817  0.6214
    ## tempave_metab.sc                                 -0.01353582 0.003370856 30165 -4.015543  0.0001
    ## seas.sc                                          -0.00256570 0.000855156 30165 -3.000266  0.0027
    ## microclim.sc                                      0.00039749 0.000419454 30165  0.947626  0.3433
    ## mass.sc                                          -0.00036278 0.001179479 30165 -0.307577  0.7584
    ## speed.sc                                         -0.00097035 0.001102377 30165 -0.880236  0.3787
    ## lifespan.sc                                      -0.00274683 0.002401726 30165 -1.143688  0.2528
    ## consumerfrac.sc                                   0.00065512 0.000642918 30165  1.018985  0.3082
    ## endothermfrac.sc                                  0.00940778 0.003406765 30165  2.761499  0.0058
    ## nspp.sc                                          -0.00093338 0.000701328 30165 -1.330879  0.1832
    ## npp.sc                                           -0.00259143 0.000622604 30165 -4.162250  0.0000
    ## veg.sc                                            0.00151622 0.001040297 30165  1.457487  0.1450
    ## temptrend_abs.sc:REALMMarine                     -0.03014951 0.018058683 30165 -1.669530  0.0950
    ## temptrend_abs.sc:REALMTerrestrial                 0.00025785 0.019516013 30165  0.013212  0.9895
    ## temptrend_abs.sc:tsign1                          -0.00686914 0.002392671 30165 -2.870907  0.0041
    ## temptrend_abs.sc:tempave.sc                      -0.00599890 0.004447399 30165 -1.348855  0.1774
    ## temptrend_abs.sc:tempave_metab.sc                 0.01915242 0.009567514 30165  2.001818  0.0453
    ## temptrend_abs.sc:seas.sc                         -0.00417756 0.002212697 30165 -1.887994  0.0590
    ## temptrend_abs.sc:microclim.sc                    -0.00190254 0.001343852 30165 -1.415738  0.1569
    ## temptrend_abs.sc:mass.sc                          0.00480460 0.003535036 30165  1.359137  0.1741
    ## temptrend_abs.sc:speed.sc                         0.00706982 0.003206633 30165  2.204749  0.0275
    ## temptrend_abs.sc:lifespan.sc                     -0.00526800 0.007264787 30165 -0.725142  0.4684
    ## temptrend_abs.sc:consumerfrac.sc                 -0.00159195 0.001651676 30165 -0.963838  0.3351
    ## temptrend_abs.sc:endothermfrac.sc                -0.02456350 0.008757694 30165 -2.804791  0.0050
    ## temptrend_abs.sc:nspp.sc                          0.00409981 0.001981499 30165  2.069044  0.0386
    ## tsign-1:thermal_bias.sc                          -0.00246793 0.000928815 30165 -2.657078  0.0079
    ## tsign1:thermal_bias.sc                           -0.00054521 0.000796312 30165 -0.684665  0.4936
    ## temptrend_abs.sc:npp.sc                           0.00951433 0.002009022 30165  4.735800  0.0000
    ## temptrend_abs.sc:veg.sc                          -0.00712066 0.002593485 30165 -2.745594  0.0060
    ## human_bowler.sc:REALM2TerrFresh                   0.00274467 0.000865756 30165  3.170263  0.0015
    ## human_bowler.sc:REALM2Marine                      0.00042854 0.000539496 30165  0.794327  0.4270
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.00093583 0.002422917 30165 -0.386241  0.6993
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00181564 0.002288292 30165  0.793448  0.4275
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.00414011 0.001834509 30165 -2.256793  0.0240
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine     0.00157953 0.001751877 30165  0.901623  0.3673
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s npp.sc veg.sc t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:t. tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.438                                                                                                                                                                                                                                                                                                                                                                                                                          
    ## REALMMarine                                      -0.942  0.425                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## REALMTerrestrial                                 -0.717  0.272  0.656                                                                                                                                                                                                                                                                                                                                                                                                            
    ## tsign1                                           -0.095  0.060  0.014  0.023                                                                                                                                                                                                                                                                                                                                                                                                     
    ## tempave.sc                                       -0.036  0.014  0.042  0.012 -0.061                                                                                                                                                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.103 -0.038 -0.100 -0.082 -0.013 -0.662                                                                                                                                                                                                                                                                                                                                                                                       
    ## seas.sc                                          -0.121  0.053  0.165 -0.099 -0.110  0.213 -0.074                                                                                                                                                                                                                                                                                                                                                                                
    ## microclim.sc                                     -0.059  0.048  0.071 -0.017  0.031 -0.034 -0.041  0.086                                                                                                                                                                                                                                                                                                                                                                         
    ## mass.sc                                          -0.011  0.037  0.026  0.054  0.018  0.012 -0.444  0.018  0.013                                                                                                                                                                                                                                                                                                                                                                  
    ## speed.sc                                          0.043 -0.028  0.007 -0.004 -0.052  0.097 -0.214  0.063  0.058  0.153                                                                                                                                                                                                                                                                                                                                                           
    ## lifespan.sc                                       0.066 -0.045 -0.071 -0.091 -0.028 -0.068  0.589  0.034  0.007 -0.826 -0.394                                                                                                                                                                                                                                                                                                                                                    
    ## consumerfrac.sc                                  -0.039  0.034  0.014  0.296  0.046  0.028 -0.063 -0.117 -0.021  0.040 -0.192 -0.052                                                                                                                                                                                                                                                                                                                                             
    ## endothermfrac.sc                                  0.023 -0.004 -0.009 -0.236  0.021  0.518 -0.551  0.084  0.030 -0.080 -0.254  0.082 -0.140                                                                                                                                                                                                                                                                                                                                      
    ## nspp.sc                                           0.011 -0.066 -0.061 -0.107 -0.007  0.042 -0.014 -0.020 -0.069 -0.191  0.010  0.181 -0.017  0.162                                                                                                                                                                                                                                                                                                                               
    ## npp.sc                                            0.032  0.013 -0.047  0.061 -0.005 -0.240  0.209 -0.205 -0.244 -0.032  0.036  0.026 -0.001 -0.170 -0.200                                                                                                                                                                                                                                                                                                                        
    ## veg.sc                                           -0.477  0.315  0.495  0.019 -0.012  0.128 -0.102  0.111  0.008  0.020  0.024 -0.022 -0.008  0.077  0.019 -0.205                                                                                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:REALMMarine                      0.427 -0.953 -0.430 -0.254 -0.020 -0.026  0.035 -0.066 -0.054 -0.038 -0.001  0.042 -0.024  0.014  0.084  0.016 -0.335                                                                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMTerrestrial                 0.248 -0.695 -0.230 -0.416 -0.007 -0.007  0.038  0.147  0.027 -0.055  0.024  0.053 -0.133  0.087  0.111 -0.076  0.027  0.639                                                                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:tsign1                           0.055 -0.120 -0.018 -0.002 -0.517  0.056 -0.013  0.012 -0.029 -0.014  0.021  0.018 -0.014  0.004  0.011 -0.034 -0.007  0.041      0.005                                                                                                                                                                                                                                                                                        
    ## temptrend_abs.sc:tempave.sc                       0.041 -0.032 -0.053 -0.020  0.055 -0.610  0.405 -0.103  0.058  0.035 -0.050  0.014 -0.031 -0.323 -0.030  0.139 -0.111  0.050      0.024     -0.107                                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.040  0.042  0.037  0.044 -0.016  0.396 -0.626 -0.013 -0.014  0.272  0.106 -0.359  0.069  0.354  0.007 -0.142  0.066 -0.035     -0.081      0.036 -0.503                                                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:seas.sc                          0.060 -0.096 -0.082  0.125  0.057 -0.096  0.017 -0.661 -0.074 -0.002 -0.026 -0.028  0.096 -0.041  0.035  0.126 -0.102  0.130     -0.222     -0.034  0.156            0.012                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:microclim.sc                     0.026 -0.062 -0.033  0.015 -0.052  0.052  0.003 -0.032 -0.723  0.006 -0.030 -0.008  0.005 -0.001  0.044  0.181  0.020  0.073     -0.048      0.027 -0.108            0.092   0.098                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                          0.026 -0.038 -0.028 -0.039 -0.011  0.028  0.250 -0.003  0.006 -0.592 -0.064  0.495 -0.032  0.051  0.108 -0.011 -0.003  0.043      0.083      0.023 -0.073           -0.471  -0.025            -0.024                                                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                        -0.026  0.046  0.004  0.032  0.034 -0.054  0.118 -0.004 -0.029 -0.071 -0.584  0.197  0.124  0.134 -0.048 -0.011 -0.032  0.003     -0.024     -0.010  0.106           -0.253   0.062             0.033             0.125                                                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                     -0.036  0.064  0.041  0.037  0.023  0.013 -0.331 -0.026 -0.005  0.495  0.181 -0.580  0.055 -0.045 -0.119 -0.006  0.010 -0.052     -0.089     -0.023 -0.017            0.630   0.042             0.010            -0.841            -0.358                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.023 -0.098 -0.013 -0.121 -0.026 -0.019  0.063  0.096  0.013 -0.037  0.141  0.060 -0.537  0.009  0.040 -0.006 -0.009  0.076      0.289      0.011  0.041           -0.154  -0.124            -0.012             0.089            -0.179            -0.151                                                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc                 0.009  0.075  0.006  0.070 -0.001 -0.370  0.401 -0.021  0.021  0.059  0.177 -0.051 -0.004 -0.524 -0.094  0.136 -0.055 -0.068     -0.230     -0.031  0.444           -0.555   0.005            -0.097            -0.083            -0.267             0.079           -0.057                                                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                         -0.046  0.074  0.053  0.087  0.012 -0.047  0.030  0.013  0.037  0.106 -0.051 -0.118  0.038 -0.107 -0.602  0.109 -0.018 -0.110     -0.151     -0.019 -0.002           -0.006  -0.048            -0.022            -0.141             0.036             0.178           -0.056            0.159                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.040 -0.003 -0.033  0.006 -0.168  0.410 -0.164 -0.098 -0.070  0.021  0.082 -0.091  0.036  0.035  0.003 -0.061  0.038  0.000     -0.018      0.025 -0.276            0.116   0.064             0.048            -0.002            -0.057             0.041           -0.022           -0.038            -0.009                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.066 -0.010 -0.071 -0.008 -0.009  0.520 -0.175 -0.171 -0.158  0.041  0.075 -0.100  0.047  0.030 -0.015 -0.115  0.044 -0.003     -0.006      0.048 -0.305            0.116   0.070             0.069            -0.001            -0.051             0.039           -0.025           -0.049            -0.020             0.573                                                                                               
    ## temptrend_abs.sc:npp.sc                           0.017  0.016 -0.006 -0.055  0.017  0.112 -0.101  0.085  0.177 -0.004 -0.009 -0.006  0.005  0.073  0.090 -0.665  0.126 -0.056      0.100      0.011 -0.163            0.154  -0.208            -0.311             0.035            -0.006             0.011           -0.011           -0.129            -0.095             0.019  0.073                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.304 -0.460 -0.315  0.020 -0.012 -0.096  0.062 -0.111  0.011 -0.005 -0.026  0.008  0.003 -0.047 -0.020  0.175 -0.665  0.492     -0.054      0.045  0.156           -0.075   0.185             0.001            -0.012             0.048            -0.001           -0.005            0.060             0.008            -0.028 -0.056 -0.267                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.074  0.090  0.075 -0.045  0.025 -0.225  0.154 -0.015  0.097  0.038 -0.015 -0.022 -0.072 -0.132 -0.013 -0.061  0.141 -0.098      0.047     -0.024  0.116           -0.086  -0.091            -0.024            -0.016            -0.003             0.016            0.013            0.102             0.043            -0.038 -0.025  0.033            -0.209                                                               
    ## human_bowler.sc:REALM2Marine                      0.038 -0.030 -0.053 -0.021  0.010 -0.030  0.030 -0.163 -0.077 -0.011 -0.009 -0.016 -0.017  0.000 -0.067 -0.145  0.020  0.044      0.016      0.005  0.013           -0.017   0.139             0.020             0.006            -0.006             0.024            0.027           -0.005             0.086             0.104  0.110  0.104            -0.013            0.021                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.013 -0.019 -0.023 -0.012  0.005 -0.297  0.137  0.063  0.072  0.001 -0.041  0.039 -0.028 -0.058 -0.015  0.022 -0.052  0.016      0.029      0.101  0.603           -0.141  -0.067            -0.117            -0.007             0.079            -0.051            0.037            0.042            -0.033            -0.485 -0.351 -0.028             0.074            0.011      -0.066                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.014 -0.011 -0.027 -0.012  0.039 -0.299  0.106  0.052  0.088  0.004 -0.029  0.038 -0.032 -0.023  0.005  0.056 -0.067  0.026      0.024     -0.079  0.626           -0.076  -0.037            -0.154            -0.020             0.052            -0.044            0.041           -0.016            -0.075            -0.346 -0.543 -0.065             0.104            0.002      -0.086       0.715                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.116 -0.226 -0.123  0.047 -0.022  0.166 -0.127 -0.095 -0.058 -0.033  0.009  0.020  0.032  0.114  0.040  0.048 -0.256  0.239     -0.076      0.068 -0.229            0.168   0.175             0.087             0.035            -0.017            -0.023           -0.032           -0.174            -0.058             0.030  0.023 -0.098             0.445           -0.591      -0.011      -0.046 -0.052               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.032  0.037  0.046  0.018 -0.008  0.022 -0.026  0.135  0.033  0.005 -0.006  0.017  0.022  0.006  0.060  0.097 -0.005 -0.055     -0.016      0.002 -0.024            0.037  -0.170            -0.061             0.004             0.003            -0.033           -0.034           -0.009            -0.104            -0.072 -0.084 -0.168             0.028           -0.013      -0.737       0.055  0.068  0.012        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.79037325 -0.27332859 -0.02733247  0.28919907  5.91687491 
    ## 
    ## Number of Observations: 30392
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    190                  30392

### Plots from the full models

#### Plot the coefficients

``` r
coefs <- summary(modTfull1)$tTable
coefs2 <- summary(modTfullfootprint)$tTable
coefs3 <- summary(modTfullJbeta)$tTable
coefs4 <- summary(modTfullHorn)$tTable
coefs5 <- summary(modTfullrem0)$tTable
coefs6 <- summary(modTfullJbetarem0)$tTable
coefs7 <- summary(modTfullHornrem0)$tTable

varstoplot <- unique(c(rownames(coefs), rownames(coefs2), rownames(coefs3), rownames(coefs4),
                       rownames(coefs5), rownames(coefs6), rownames(coefs7)))
varstoplot <- varstoplot[which(!grepl('Intercept', varstoplot) | grepl(':', varstoplot))] # vars to plot

rows1_1 <- which(rownames(coefs) %in% varstoplot) # rows in coefs
rows1_2 <- which(rownames(coefs2) %in% varstoplot)
rows1_3 <- which(rownames(coefs3) %in% varstoplot)
rows1_4 <- which(rownames(coefs4) %in% varstoplot)
rows1_5 <- which(rownames(coefs5) %in% varstoplot)
rows1_6 <- which(rownames(coefs6) %in% varstoplot)
rows1_7 <- which(rownames(coefs7) %in% varstoplot)
xlims <- range(c(coefs[rows1_1,1] - coefs[rows1_1,2], coefs[rows1_1,1] + coefs[rows1_1,2], 
                  coefs2[rows1_2,1] - coefs2[rows1_2,2], coefs2[rows1_2,1] + coefs2[rows1_2,2], 
                  coefs3[rows1_3,1] - coefs3[rows1_3,2], coefs3[rows1_3,1] + coefs3[rows1_3,2],
                  coefs4[rows1_4,1] - coefs4[rows1_4,2], coefs4[rows1_4,1] + coefs4[rows1_4,2],
                  coefs5[rows1_5,1] - coefs5[rows1_5,2], coefs5[rows1_5,1] + coefs5[rows1_5,2], 
                  coefs6[rows1_6,1] - coefs6[rows1_6,2], coefs6[rows1_6,1] + coefs6[rows1_6,2],
                  coefs7[rows1_7,1] - coefs7[rows1_7,2], coefs7[rows1_7,1] + coefs7[rows1_7,2]))


cols <- brewer.pal(7, 'Dark2') # for Jtu Bowler, Jtu footprint, Jbeta and Horn models, Jtu rem0, Jbeta rem0, Horn rem0
pchs <- c(16, 15, 16, 16, 16, 16, 16)
offs <- c(0.3, 0.2, 0.1, 0, -0.1, -0.2, -0.3) # offset vertically for each model



par(las = 1, mai = c(0.5, 4, 0.1, 0.1))

plot(0,0, col = 'white', xlim = xlims, ylim = c(1,length(varstoplot)), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(varstoplot):1, labels = varstoplot, cex.axis = 0.7)
abline(v = 0, col = 'grey', lty = 2)
abline(h = 1:length(varstoplot), col = 'grey', lty = 3)
for(i in 1:length(varstoplot)){
  if(varstoplot[i] %in% rownames(coefs)){
    x = coefs[rownames(coefs) == varstoplot[i], 1]
    se = coefs[rownames(coefs) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[1], pch = pchs[1], col = cols[1])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[1], length(varstoplot) + 1 - i + offs[1]), col = cols[1])
  }
  if(varstoplot[i] %in% rownames(coefs2)){
    x = coefs2[rownames(coefs2) == varstoplot[i], 1]
    se = coefs2[rownames(coefs2) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[2], pch = pchs[2], col = cols[2])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[2], length(varstoplot) + 1 - i + offs[2]), col = cols[2])
  }
  if(varstoplot[i] %in% rownames(coefs3)){
    x = coefs3[rownames(coefs3) == varstoplot[i], 1]
    se = coefs3[rownames(coefs3) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[3], pch = pchs[3], col = cols[3])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[3], length(varstoplot) + 1 - i + offs[3]), col = cols[3])
  }
  if(varstoplot[i] %in% rownames(coefs4)){
    x = coefs4[rownames(coefs4) == varstoplot[i], 1]
    se = coefs4[rownames(coefs4) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[4], pch = pchs[4], col = cols[4])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[4], length(varstoplot) + 1 - i + offs[4]), col = cols[4])
  }
  if(varstoplot[i] %in% rownames(coefs5)){
    x = coefs5[rownames(coefs5) == varstoplot[i], 1]
    se = coefs5[rownames(coefs5) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[5], pch = pchs[5], col = cols[5])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[5], length(varstoplot) + 1 - i + offs[5]), col = cols[5])
  }
  if(varstoplot[i] %in% rownames(coefs6)){
    x = coefs6[rownames(coefs6) == varstoplot[i], 1]
    se = coefs6[rownames(coefs6) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[6], pch = pchs[6], col = cols[6])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[6], length(varstoplot) + 1 - i + offs[6]), col = cols[6])
  }
  if(varstoplot[i] %in% rownames(coefs7)){
    x = coefs7[rownames(coefs7) == varstoplot[i], 1]
    se = coefs7[rownames(coefs7) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[7], pch = pchs[7], col = cols[7])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[7], length(varstoplot) + 1 - i + offs[7]), col = cols[7])
  }
}
legend('topleft', col = cols, pch = 16, lwd = 1, legend = c('Jtu', 'Jtu footprint', 'Jbeta', 'Horn', 'Jtu rem0', 'Jbeta rem0', 'Horn rem0'))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20fullTmod1-1.png)<!-- -->

#### Plot interactions (Jaccard turnover)

``` r
# set up the interactions to plot
ints <- data.frame(vars = c('tsign', 'tempave', 'tempave_metab', 'seas', 'microclim', 'mass', 'speed', 'lifespan', 'consumerfrac', 'endothermfrac', 'nspp', 'thermal_bias', 'npp', 'veg', 'human_bowler', 'human_bowler'),
           min = c(1, -10, 10, 0.1, 0, 0, 0, 0.3, 0, 0, 0.3, -10, 1.9, 0, 0, 0), 
           max = c(2, 30, 40, 16, 6, 8, 2, 2, 1, 1, 2.6, 10, 3.7, 1, 9, 9),
           log = c(F, F, F, F, F, T, T, T, F, F, T, F, T, F, F, F),
           len = c(2, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
           discrete = c(T, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F),
           REALM = c(rep('Freshwater', 15), 'Marine'),
           REALM2 = c(rep('TerrFresh', 15), 'Marine'),
           stringsAsFactors = FALSE)
basetab <- data.frame(tempave.sc = 0, tempave_metab.sc = 0, 
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

  # add realm
  thisdat$REALM <- ints$REALM[j]
  thisdat$REALM2 <- ints$REALM2[j]
  
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

# character so that new levels can be added
newdat$REALM <- as.character(newdat$REALM)
newdat$REALM2 <- as.character(newdat$REALM2)

# add extra rows so that all factor levels are represented (for predict.lme to work)
newdat <- rbind(newdat[1:4, ], newdat)
newdat$REALM[1:4] <- c('Marine', 'Marine', 'Terrestrial', 'Terrestrial')
newdat$REALM2[1:4] <- c('Marine', 'Marine', 'TerrFresh', 'TerrFresh')
newdat$temptrend[1:4] <- c(-1, 1, -1, 1)

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
newdat <- newdat[5:nrow(newdat), ]

# prep the plots
intplots <- vector('list', nrow(ints))
for(j in 1:length(intplots)){
  subs <- newdat$var == ints$vars[j] & newdat$temptrend > 0 # select warming side
  xvar <- 'temptrend_abs'
  title <- ints$vars[j]
  if(ints$vars[j] %in% c('tsign')){
    subs <- newdat$var == ints$vars[j]
  } 
  if(ints$vars[j] %in% c('thermal_bias')){
    subs <- newdat$var == ints$vars[j]
    xvar <- 'temptrend'
  } 
  if(ints$vars[j] %in% c('human_bowler')){
    subs <- newdat$var == ints$vars[j] & newdat$temptrend > 0 & newdat$REALM2 == ints$REALM2[j]
    title <- paste0('human:', ints$REALM2[j])
  } 

  thisplot <- ggplot(newdat[subs, ], 
                     aes_string(x = xvar, y = 'preds', 
                                group = ints$vars[j], 
                                color = ints$vars[j])) +
    geom_line() +
    coord_cartesian(ylim = c(0, 1)) +
    theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm')) +
    labs(title = title)
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

#### Plot interactions (Jaccard turnover) without year 1

``` r
# set up the interactions to plot
ints <- data.frame(vars = c('tsign', 'tempave', 'tempave_metab', 'seas', 'microclim', 'mass', 'speed', 'lifespan', 'consumerfrac', 'endothermfrac', 'nspp', 'thermal_bias', 'npp', 'veg', 'human_bowler', 'human_bowler'),
           min = c(1, -10, 10, 0.1, 0, 0, 0, 0.3, 0, 0, 0.3, -10, 1.9, 0, 0, 0), 
           max = c(2, 30, 40, 16, 6, 8, 2, 2, 1, 1, 2.6, 10, 3.7, 1, 9, 9),
           log = c(F, F, F, F, F, T, T, T, F, F, T, F, T, F, F, F),
           len = c(2, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
           discrete = c(T, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F),
           REALM = c(rep('Freshwater', 15), 'Marine'),
           REALM2 = c(rep('TerrFresh', 15), 'Marine'),
           stringsAsFactors = FALSE)
basetab <- data.frame(tempave.sc = 0, tempave_metab.sc = 0, 
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

  # add realm
  thisdat$REALM <- ints$REALM[j]
  thisdat$REALM2 <- ints$REALM2[j]
  
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

# character so that new levels can be added
newdat$REALM <- as.character(newdat$REALM)
newdat$REALM2 <- as.character(newdat$REALM2)

# add extra rows so that all factor levels are represented (for predict.lme to work)
newdat <- rbind(newdat[1:4, ], newdat)
newdat$REALM[1:4] <- c('Marine', 'Marine', 'Terrestrial', 'Terrestrial')
newdat$REALM2[1:4] <- c('Marine', 'Marine', 'TerrFresh', 'TerrFresh')
newdat$temptrend[1:4] <- c(-1, 1, -1, 1)

# trim to at least some temperature change (so that tsign is -1 or 1)
newdat <- newdat[newdat$temptrend != 0,]

# scale the temperature vars
newdat$temptrend.sc <- newdat$temptrend/attr(trends$temptrend.sc, 'scaled:scale') 
newdat$temptrend_abs <- abs(newdat$temptrend)
newdat$temptrend_abs.sc <- (newdat$temptrend_abs)/attr(trends$temptrend_abs.sc, 'scaled:scale')
newdat$tsign <- factor(sign(newdat$temptrend))

# make predictions
newdat$preds <- predict(object = modTfullrem0, newdata = newdat, level = 0)

#remove the extra rows
newdat <- newdat[5:nrow(newdat), ]

# prep the plots
intplots <- vector('list', nrow(ints))
for(j in 1:length(intplots)){
  subs <- newdat$var == ints$vars[j] & newdat$temptrend > 0 # select warming side
  xvar <- 'temptrend_abs'
  title <- ints$vars[j]
  if(ints$vars[j] %in% c('tsign')){
    subs <- newdat$var == ints$vars[j]
  } 
  if(ints$vars[j] %in% c('thermal_bias')){
    subs <- newdat$var == ints$vars[j]
    xvar <- 'temptrend'
  } 
  if(ints$vars[j] %in% c('human_bowler')){
    subs <- newdat$var == ints$vars[j] & newdat$temptrend > 0 & newdat$REALM2 == ints$REALM2[j]
    title <- paste0('human:', ints$REALM2[j])
  } 

  thisplot <- ggplot(newdat[subs, ], 
                     aes_string(x = xvar, y = 'preds', 
                                group = ints$vars[j], 
                                color = ints$vars[j])) +
    geom_line() +
    coord_cartesian(ylim = c(-0.6, 0.6)) +
    theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm')) +
    labs(title = title)
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

![](turnover_vs_temperature_MEmodels_files/figure-gfm/interaction%20plots%20modTfullrem0-1.png)<!-- -->

#### Plot residuals against each predictor (Jaccard turnover)

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

### Remove each term from the full model

#### Jaccard turnover

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

#### Jbeta and Horn models

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

#### Plot deltaAICs for all 3 models

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

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20dAICs-1.png)<!-- -->
Light grey is for Jaccard turnover, dark grey is for Jaccard total,
black is for Morisita-Horn

### Plot interaction coefficients from all full models

``` r
# fig.width = 3, fig.height = 5, out.height=2.5, out.width=3, fig.retina =3 for macbook screen
# double that for external monitor
coefs <- as.data.table(summary(modTfull1)$tTable)
coefs2 <- as.data.table(summary(modTfullJbeta)$tTable)
coefs3 <- as.data.table(summary(modTfullHorn)$tTable)
coefs4 <- as.data.table(summary(modTfullrem0)$tTable)
coefs5 <- as.data.table(summary(modTfullJbetarem0)$tTable)
coefs6 <- as.data.table(summary(modTfullHornrem0)$tTable)

coefs$mod <- 'Jtu'
coefs2$mod <- 'Jbeta'
coefs3$mod <- 'Horn'
coefs4$mod <- 'Jtu rem0'
coefs5$mod <- 'Jbeta rem0'
coefs6$mod <- 'Horn rem0'

coefs$var <- rownames(summary(modTfull1)$tTable)
coefs2$var <- rownames(summary(modTfullJbeta)$tTable)
coefs3$var <- rownames(summary(modTfullHorn)$tTable)
coefs4$var <- rownames(summary(modTfullrem0)$tTable)
coefs5$var <- rownames(summary(modTfullJbetarem0)$tTable)
coefs6$var <- rownames(summary(modTfullHornrem0)$tTable)

# extract temperature effects and bind model coefs together
cols <- c('var', 'Value', 'Std.Error', 'mod')

allcoefsfull <- rbind(coefs[grep('temptrend|REALM', var), ..cols], 
                      coefs2[grep('temptrend|REALM', var), ..cols], 
                      coefs3[grep('temptrend|REALM', var), ..cols],
                      coefs4[grep('temptrend|REALM', var), ..cols],
                      coefs5[grep('temptrend|REALM', var), ..cols],
                      coefs6[grep('temptrend|REALM', var), ..cols])
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
xlims2 <- c(-0.05, 0.05) # for traits
xlims3 <- c(-0.01, 0.02) # for environment
xlims4 <- c(-0.016, 0.01) # for community
xlims5 <- c(-0.006, 0.005) # for human

ddg <- 0.5 # vertical dodge for each model

# choose which variables are in which graph
set1 <- c('Terrestrial', 'Marine', 'Freshwater')
set2 <- c('mass', 'speed', 'lifespan', 'consumerfrac', 'endothermfrac', 'tempave_metab')
set3 <- c('seas', 'microclim', 'tempave')
set4 <- c('npp', 'nspp', 'tsign-1:thermal_bias', 'tsign1:thermal_bias')
set5 <- c('human:TerrFresh', 'human:Marine')

p1 <- ggplot(subset(allcoefsfull, varname %in% set1), 
                    aes(varname, Value, group = mod, color = mod)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'light grey') +
  geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0, position = position_dodge(ddg)) + 
  geom_point(position = position_dodge(ddg)) + 
  labs(y = 'Temperature change effect', x = '', tag = 'A') +
  scale_color_brewer(palette = 'Dark2') + 
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
  scale_color_brewer(palette = 'Dark2') +
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
  scale_color_brewer(palette = 'Dark2') +
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
  scale_color_brewer(palette = 'Dark2') +
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
  scale_color_brewer(palette = 'Dark2') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position='bottom',
        axis.text=element_text(size=7),
        axis.title=element_text(size=7)) + 
  coord_flip(ylim = xlims5)

grid.arrange(p1, p2, p3, p4, p5, ncol = 2, layout_matrix = rbind(c(1,2), c(3,4), c(5, NA)))
```

<img src="turnover_vs_temperature_MEmodels_files/figure-gfm/plot coefs for all models-1.png" width="6" height="5" />

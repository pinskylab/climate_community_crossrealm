Drivers of variation in the community response to temperature change
across realms
================

Collaborators: Shane Blowes, Jon Chase, Helmut Hillebrand, Michael
Burrows, Amanda Bates, Uli Brose, Benoit Gauzens, Laura Antao, Ruben
Remelgado, Carsten Meyer, Myriam Hirt, maybe others Assistance:
Katherine Lew, Josef Hauser

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
        (NOT including the first year compared to itself)
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
library(MASS) # for stepAIC
library(piecewiseSEM) # for rsquared() for nlme models
```

    ## Registered S3 methods overwritten by 'lme4':
    ##   method                          from
    ##   cooks.distance.influence.merMod car 
    ##   influence.merMod                car 
    ##   dfbeta.influence.merMod         car 
    ##   dfbetas.influence.merMod        car

    ## 
    ##   This is piecewiseSEM version 2.1.0.
    ## 
    ## 
    ##   Questions or bugs can be addressed to <LefcheckJ@si.edu>.

``` r
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

# calculate duration
trends[, duration := maxyrBT - minyrBT + 1]

# trim to data with >= 3 yrs
trends <- trends[nyrBT >= 3, ]
```

### Log-transform some variables, then center and scale.

``` r
trends[, tempave.sc := scale(tempave)]
trends[, tempave_metab.sc := scale(tempave_metab)]
trends[, seas.sc := scale(seas)]
trends[, microclim.sc := scale(log(microclim))]
trends[, temptrend.sc := scale(temptrend, center = FALSE)] # do not center
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
trends[, duration.sc := scale(log(duration))]
trends[, human_bowler.sc := scale(log(human_bowler+1)), by = REALM2] # separate scaling by realm
trends[REALM2 == 'TerrFresh', human_footprint.sc := scale(log(human_venter+1))]
trends[REALM2 == 'Marine', human_footprint.sc := scale(log(human_halpern))]
```

### Examine how many data points are available

Just turnover

``` r
cat('Overall # time-series: ', nrow(trends), '\n')
```

    ## Overall # time-series:  38188

``` r
cat('# studies: ', trends[, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  268

``` r
cat('Data points: ', trends[, sum(nyrBT)], '\n')
```

    ## Data points:  258879

``` r
trends[, table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##         628       34752        2808

``` r
trends[, table(taxa_mod)]
```

    ## taxa_mod
    ##                         All                  Amphibians                     Benthos                       Birds                        Fish               Invertebrates                     Mammals Marine invertebrates/plants                       Plant                    Reptiles 
    ##                        1442                         352                        4305                        8674                       21198                        1391                         504                         104                         214                           4

``` r
trends[, table(taxa_mod, REALM)]
```

    ##                              REALM
    ## taxa_mod                      Freshwater Marine Terrestrial
    ##   All                                  0   1439           3
    ##   Amphibians                           2      0         350
    ##   Benthos                              0   4305           0
    ##   Birds                                0   6499        2175
    ##   Fish                               610  20588           0
    ##   Invertebrates                       14   1298          79
    ##   Mammals                              0    459          45
    ##   Marine invertebrates/plants          0    104           0
    ##   Plant                                1     60         153
    ##   Reptiles                             1      0           3

With all covariates (Bowler for human)

``` r
# the cases we can compare
apply(trends[, .(Jtutrendrem0, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##     Jtutrendrem0            REALM       tempave.sc tempave_metab.sc          seas.sc     microclim.sc     temptrend.sc          mass.sc         speed.sc      lifespan.sc  consumerfrac.sc endothermfrac.sc          nspp.sc  thermal_bias.sc           npp.sc           veg.sc  human_bowler.sc 
    ##            38188            38188            36005            36005            36005            37289            36005            38113            38127            37299            38188            38188            38188            35580            38090            38095            38188

``` r
i <- trends[, complete.cases(Jtutrendrem0, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
cat('Overall # time-series: ', sum(i), '\n')
```

    ## Overall # time-series:  35327

``` r
cat('# studies: ', trends[i, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  200

``` r
cat('Data points: ', trends[i, sum(nyrBT)], '\n')
```

    ## Data points:  238085

``` r
trends[i, table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##         608       32420        2299

``` r
trends[i, table(taxa_mod)]
```

    ## taxa_mod
    ##                         All                  Amphibians                     Benthos                       Birds                        Fish               Invertebrates                     Mammals Marine invertebrates/plants                       Plant                    Reptiles 
    ##                        1417                          12                        4283                        7164                       20413                        1297                         495                         104                         140                           2

``` r
trends[i, table(taxa_mod, REALM)]
```

    ##                              REALM
    ## taxa_mod                      Freshwater Marine Terrestrial
    ##   All                                  0   1415           2
    ##   Amphibians                           2      0          10
    ##   Benthos                              0   4283           0
    ##   Birds                                0   5084        2080
    ##   Fish                               597  19816           0
    ##   Invertebrates                        8   1224          65
    ##   Mammals                              0    459          36
    ##   Marine invertebrates/plants          0    104           0
    ##   Plant                                1     35         104
    ##   Reptiles                             0      0           2

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
fixed <- formula(Jtutrendrem0 ~ temptrend_abs.sc*REALM +
                     temptrend_abs.sc*tsign + 
                     temptrend_abs.sc*tempave_metab.sc + 
                     temptrend_abs.sc*seas.sc + 
                     temptrend_abs.sc*microclim.sc + 
                     temptrend_abs.sc*mass.sc + 
                     temptrend_abs.sc*speed.sc + 
                     temptrend_abs.sc*consumerfrac.sc +
                     temptrend_abs.sc*nspp.sc +
                     temptrend_abs.sc*thermal_bias.sc:tsign +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*veg.sc +
                     temptrend_abs.sc*duration.sc +
                     temptrend_abs.sc*human_bowler.sc:REALM2)
i <- trends[, complete.cases(Jtutrendrem0, temptrend_abs.sc, REALM, tsign, tempave_metab.sc, seas.sc, 
                             microclim.sc, mass.sc, speed.sc, consumerfrac.sc, nspp.sc,
                             thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
mods <- vector('list', 0)
mods[[1]] <- gls(fixed, data = trends[i,])
mods[[2]] <- gls(fixed, data = trends[i,], weights = varPower(-0.5, ~nyrBT))
mods[[3]] <- gls(fixed, data = trends[i,], weights = varPower(0.5, ~temptrend_abs.sc))

mods[[4]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2, control = lmeControl(opt = "optim"))
mods[[5]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, control = lmeControl(opt = "optim"))
mods[[6]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID, control = lmeControl(opt = "optim"))
mods[[7]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, control = lmeControl(opt = "optim"))
mods[[8]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID/rarefyID, control = lmeControl(opt = "optim"))

mods[[9]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc | taxa_mod)
mods[[10]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)) # includes overdispersion. new formula so that random slope is only for study level (not enough data to extend to rarefyID).

mods[[11]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT))
mods[[12]] <- lme(fixed, data = trends[i,], random = list(taxa_mod2 = ~ temptrend_abs.sc, STUDY_ID = ~ 1, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT))

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

Mostly northern hemisphere, but spread all over. Not so much in Africa
or much of Asia.

## Average rates of turnover

``` r
trends[abs(temptrend) >= 0.5, temptrendtext1 := 'Changing >=0.5']
trends[abs(temptrend) <= 0.05, temptrendtext1 := 'Stable <=0.05']
trends[temptrend <= -0.5, temptrendtext2 := 'Cooling <= -0.5']
trends[temptrend >= 0.5, temptrendtext2 := 'Warming >= 0.5']
trends[abs(temptrend) >= 0.5 & abs(temptrend) < 1, temptrendtext3 := 'Slow Changing >=0.5 & <1']
trends[abs(temptrend) >= 1, temptrendtext3 := 'Fast Changing >=1']

# reshape to long format
measurenms <- c('Jtutrendrem0', 'Jtutrendrem0_se', 'Jbetatrendrem0', 'Jbetatrendrem0_se', 'Horntrendrem0', 'Horntrendrem0_se', 'Jtutrendz', 'Jbetatrendz', 'Horntrendz', 'Jtulast', 'Jbetalast', 'Hornlast', 'Jtuexp', 'Jbetaexp', 'Hornexp', 'Jtumm', 'Jbetamm', 'Hornmm')
idnms <- setdiff(names(trends), measurenms)
trends2 <- melt(trends, id.vars = idnms, measure.vars = measurenms)

# changing vs. stable
trendsum1 <- trends2[!is.na(temptrendtext1), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = temptrendtext1, type = variable)] # turnover per year for locations changing temperature


# warming vs. cooling
trendsum2 <- trends2[!is.na(temptrendtext2), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = temptrendtext2, type = variable)] # inc. direction

# slow vs. fast changing
trendsum3 <- trends2[!is.na(temptrendtext3), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = temptrendtext3, type = variable)] # inc. rate

# combine
trendsum4 <- rbind(trendsum1, trendsum2, trendsum3)
setorder(trendsum4, type, text)

write.csv(trendsum4, file = 'output/trendsummary.csv', row.names = FALSE)

trendsum4
```

    ##                          text         type          ave           se     n
    ##   1:           Changing >=0.5 Jtutrendrem0  0.019110224 0.0120121948   373
    ##   2:          Cooling <= -0.5 Jtutrendrem0  0.020361270 0.0135748363   255
    ##   3:        Fast Changing >=1 Jtutrendrem0  0.045305535 0.0295638380    39
    ##   4: Slow Changing >=0.5 & <1 Jtutrendrem0  0.016051490 0.0129653010   334
    ##   5:            Stable <=0.05 Jtutrendrem0  0.005197229 0.0005253826 25094
    ##  ---                                                                      
    ## 104:          Cooling <= -0.5       Hornmm  1.194348117 0.1256875210   254
    ## 105:        Fast Changing >=1       Hornmm  1.398230131 0.3094748400    39
    ## 106: Slow Changing >=0.5 & <1       Hornmm  1.138034377 0.1346812870   333
    ## 107:            Stable <=0.05       Hornmm 10.247314407 0.4161829010 24994
    ## 108:           Warming >= 0.5       Hornmm  1.102813567 0.2865771292   118

### Plots of turnover rates

``` r
p1 <- ggplot(trends[!is.na(text), ], aes(temptrendtext1, Horntrendrem0)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = 'grey') +
  labs(x = '', y = 'Slope (Horn)', tag = 'A', title = 'Rate of temperature change') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10)) +
  geom_abline(intercept = 0, slope = 0)

p2 <- ggplot(trends[!is.na(text), ], aes(temptrendtext2, Horntrendrem0)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = 'grey') +
  labs(x = '', y = 'Slope (Horn)', tag = 'B', title = 'Rate of temperature change') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10)) +
  geom_abline(intercept = 0, slope = 0)

p3 <- ggplot(trendsum4[type %in% c('Jtutrendrem0', 'Jbetatrendrem0', 'Horntrendrem0')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Slope', title = 'Slope of dissimilarity vs. time') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  coord_cartesian(ylim = c(-0.01,0.08))

p4 <- ggplot(trendsum4[type %in% c('Jtutrendz', 'Jbetatrendz', 'Horntrendz')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Standardized slope', title = 'Null model standardized slope') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p5 <- ggplot(trendsum4[type %in% c('Jtulast', 'Jbetalast', 'Hornlast')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Dissimilarity', title = 'First to last year') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p6 <- ggplot(trendsum4[type %in% c('Jtuexp', 'Jbetaexp', 'Hornexp')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Years', title = 'Exponential half-saturation') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p7 <- ggplot(trendsum4[type %in% c('Jtumm', 'Jbetamm', 'Hornmm')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Years', title = 'MM half-saturation') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 2, 
             layout_matrix = rbind(c(1,2), c(3,3), c(4,4), c(5,5), c(6,6), c(7,7)),
             heights=c(unit(0.1, "npc"), unit(0.18, "npc"), 
                       unit(0.18, "npc"), unit(0.18, "npc"), 
                       unit(0.18, "npc"), unit(0.18, "npc")))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-1.png)<!-- -->

## Temperature-only models

### Fit the models

``` r
randef <- list(STUDY_ID = ~ abs(temptrend), rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# *trendrem0
if(file.exists('temp/modonlyTtrendrem0.rds')){
  modonlyTtrendrem0 <- readRDS('temp/modonlyTtrendrem0.rds')
} else {
  i <- trends[, complete.cases(Jtutrendrem0, REALM, temptrend)]
  modonlyTtrendrem0 <- lme(Jtutrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrendrem0, file = 'temp/modonlyTtrendrem0.rds')
}

if(file.exists('temp/modonlyTtrendJbetarem0.rds')){
  modonlyTtrendJbetarem0 <- readRDS('temp/modonlyTtrendJbetarem0.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendrem0, REALM, temptrend)]
  modonlyTtrendJbetarem0 <- lme(Jbetatrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetarem0, file = 'temp/modonlyTtrendJbetarem0.rds')
}

if(file.exists('temp/modonlyTtrendHornrem0.rds')){
  modonlyTtrendHornrem0 <- readRDS('temp/modonlyTtrendHornrem0.rds')
} else {
  i <- trends[, complete.cases(Horntrendrem0, REALM, temptrend)]
  modonlyTtrendHornrem0 <- lme(Horntrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrendHornrem0, file = 'temp/modonlyTtrendHornrem0.rds')
}

# *trendz
# doesn't converge
# if(file.exists('temp/modonlyTtrendJtuz.rds')){
#   modonlyTtrendJtuz <- readRDS('temp/modonlyTtrendJtuz.rds')
# } else {
#   i <- trends[, complete.cases(Jtutrendz, REALM, temptrend)]
#   modonlyTtrendJtuz <- lme(Jtutrendz ~ abs(temptrend)*REALM,
#                    random = randef, weights = varef, data = trends[i,], method = 'REML', 
#                    control=lmeControl(msMaxIter = 100, maxIter = 100))
#   saveRDS(modonlyTtrendJtuz, file = 'temp/modonlyTtrendJtuz.rds')
# }

if(file.exists('temp/modonlyTtrendJbetaz.rds')){
  modonlyTtrendJbetaz <- readRDS('temp/modonlyTtrendJbetaz.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendz, REALM, temptrend)]
  modonlyTtrendJbetaz <- lme(Jbetatrendz ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetaz, file = 'temp/modonlyTtrendJbetaz.rds')
}

if(file.exists('temp/modonlyTtrendHornz.rds')){
  modonlyTtrendHornz <- readRDS('temp/modonlyTtrendHornz.rds')
} else {
  i <- trends[, complete.cases(Horntrendz, REALM, temptrend)]
  modonlyTtrendHornz <- lme(Horntrendz ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrendHornz, file = 'temp/modonlyTtrendHornz.rds')
}

# *exp
if(file.exists('temp/modonlyTtrendJtuexp.rds')){
  modonlyTtrendJtuexp <- readRDS('temp/modonlyTtrendJtuexp.rds')
} else {
  i <- trends[, complete.cases(Jtuexp, REALM, temptrend)]
  modonlyTtrendJtuexp <- lme(log(Jtuexp) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML',
                   control=lmeControl(msMaxIter = 100, maxIter = 100, opt = 'optim'))
  saveRDS(modonlyTtrendJtuexp, file = 'temp/modonlyTtrendJtuexp.rds')
}

if(file.exists('temp/modonlyTtrendJbetaexp.rds')){
  modonlyTtrendJbetaexp <- readRDS('temp/modonlyTtrendJbetaexp.rds')
} else {
  i <- trends[, complete.cases(Jbetaexp, REALM, temptrend)]
  modonlyTtrendJbetaexp <- lme(log(Jbetaexp) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetaexp, file = 'temp/modonlyTtrendJbetaexp.rds')
}

if(file.exists('temp/modonlyTtrendHornexp.rds')){
  modonlyTtrendHornexp <- readRDS('temp/modonlyTtrendHornexp.rds')
} else {
  i <- trends[, complete.cases(Hornexp, REALM, temptrend)]
  modonlyTtrendHornexp <- lme(log(Hornexp) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML',
                   control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modonlyTtrendHornexp, file = 'temp/modonlyTtrendHornexp.rds')
}


# *mm
if(file.exists('temp/modonlyTtrendJtumm.rds')){
  modonlyTtrendJtumm <- readRDS('temp/modonlyTtrendJtumm.rds')
} else {
  i <- trends[, complete.cases(Jtumm, REALM, temptrend)]
  modonlyTtrendJtumm <- lme(log(Jtumm+1) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrendJtumm, file = 'temp/modonlyTtrendJtumm.rds')
}

if(file.exists('temp/modonlyTtrendJbetamm.rds')){
  modonlyTtrendJbetamm <- readRDS('temp/modonlyTtrendJbetamm.rds')
} else {
  i <- trends[, complete.cases(Jbetamm, REALM, temptrend)]
  modonlyTtrendJbetamm <- lme(log(Jbetamm+1) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetamm, file = 'temp/modonlyTtrendJbetamm.rds')
}

if(file.exists('temp/modonlyTtrendHornmm.rds')){
  modonlyTtrendHornmm <- readRDS('temp/modonlyTtrendHornmm.rds')
} else {
  i <- trends[, complete.cases(Hornmm, REALM, temptrend)]
  modonlyTtrendHornmm <- lme(log(Hornmm+1) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendHornmm, file = 'temp/modonlyTtrendHornmm.rds')
}
```

### Summary: Slope of dissimilarity

``` r
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

### Summary: Standardized slope of dissimilarity

``` r
# summary(modonlyTtrendJtuz)
summary(modonlyTtrendJbetaz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##      AIC      BIC    logLik
    ##   105975 106075.7 -52975.52
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)     3.191278 (Intr)
    ## abs(temptrend) 17.052110 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:    1.141586 785627979
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -18.77909 
    ## Fixed effects: Jbetatrendz ~ abs(temptrend) * REALM 
    ##                                     Value Std.Error    DF    t-value p-value
    ## (Intercept)                      0.932458  1.159832 32428  0.8039600  0.4214
    ## abs(temptrend)                   0.526978  7.927470 32428  0.0664749  0.9470
    ## REALMMarine                     -0.180835  1.219576   180 -0.1482770  0.8823
    ## REALMTerrestrial                 0.845357  1.224345   180  0.6904560  0.4908
    ## abs(temptrend):REALMMarine      -5.110251  8.291061 32428 -0.6163566  0.5377
    ## abs(temptrend):REALMTerrestrial -4.893105  8.692704 32428 -0.5628979  0.5735
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.199                                
    ## REALMMarine                     -0.951  0.190                         
    ## REALMTerrestrial                -0.947  0.189  0.901                  
    ## abs(temptrend):REALMMarine       0.191 -0.956 -0.191 -0.181           
    ## abs(temptrend):REALMTerrestrial  0.182 -0.912 -0.173 -0.222  0.872    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -6.878401e+05 -5.275427e-04 -1.431198e-08  2.544332e-06  1.855959e+04 
    ## 
    ## Number of Observations: 32614
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    183                  32614

``` r
summary(modonlyTtrendHornz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##      AIC      BIC    logLik
    ##   167291 167391.6 -83633.48
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)     5.322111 (Intr)
    ## abs(temptrend) 37.053865 -0.39 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    3.129784 1.439185
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -4.601958 
    ## Fixed effects: Horntrendz ~ abs(temptrend) * REALM 
    ##                                      Value Std.Error    DF    t-value p-value
    ## (Intercept)                       1.536020  2.198935 32339  0.6985294  0.4849
    ## abs(temptrend)                   10.313924 17.463360 32339  0.5906036  0.5548
    ## REALMMarine                       0.393828  2.294248   171  0.1716589  0.8639
    ## REALMTerrestrial                  0.449837  2.323632   171  0.1935921  0.8467
    ## abs(temptrend):REALMMarine      -16.091723 18.273936 32339 -0.8805833  0.3785
    ## abs(temptrend):REALMTerrestrial  -4.990263 19.002943 32339 -0.2626047  0.7929
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.516                                
    ## REALMMarine                     -0.958  0.495                         
    ## REALMTerrestrial                -0.946  0.488  0.907                  
    ## abs(temptrend):REALMMarine       0.493 -0.956 -0.510 -0.467           
    ## abs(temptrend):REALMTerrestrial  0.474 -0.919 -0.455 -0.521  0.878    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -5.155624e-02 -7.990445e-05 -2.495224e-06  1.888799e-05  9.571167e-02 
    ## 
    ## Number of Observations: 32516
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    174                  32516

### Summary: Half-saturation exponential fits

``` r
#summary(modonlyTtrendJtuexp)
summary(modonlyTtrendJbetaexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   86640.16 86737.37 -43308.08
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)    0.9162896 (Intr)
    ## abs(temptrend) 2.7519536 -0.78 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.147013 2.019918
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.5356668 
    ## Fixed effects: log(Jbetaexp) ~ abs(temptrend) * REALM 
    ##                                     Value Std.Error    DF   t-value p-value
    ## (Intercept)                      0.514941 0.3178907 24122  1.619867  0.1053
    ## abs(temptrend)                  -5.857213 1.6869493 24122 -3.472074  0.0005
    ## REALMMarine                     -1.184732 0.3381153   230 -3.503929  0.0006
    ## REALMTerrestrial                -0.361446 0.3493468   230 -1.034635  0.3019
    ## abs(temptrend):REALMMarine       3.126224 1.7609791 24122  1.775276  0.0759
    ## abs(temptrend):REALMTerrestrial  2.857828 1.8409701 24122  1.552349  0.1206
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.596                                
    ## REALMMarine                     -0.940  0.560                         
    ## REALMTerrestrial                -0.910  0.542  0.856                  
    ## abs(temptrend):REALMMarine       0.571 -0.958 -0.598 -0.519           
    ## abs(temptrend):REALMTerrestrial  0.546 -0.916 -0.513 -0.599  0.878    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.21781048 -0.31408890  0.02948519  0.37833591  2.08374807 
    ## 
    ## Number of Observations: 24358
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    233                  24358

``` r
summary(modonlyTtrendHornexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   91586.04 91682.69 -45781.02
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev   Corr  
    ## (Intercept)    1.148612 (Intr)
    ## abs(temptrend) 1.557478 -0.752
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.659098   241.26
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -5.069527 
    ## Fixed effects: log(Hornexp) ~ abs(temptrend) * REALM 
    ##                                     Value Std.Error    DF   t-value p-value
    ## (Intercept)                      0.806719 0.4137573 23040  1.949740  0.0512
    ## abs(temptrend)                  -6.802335 1.6448215 23040 -4.135607  0.0000
    ## REALMMarine                     -0.769217 0.4377119   217 -1.757359  0.0803
    ## REALMTerrestrial                 0.284607 0.4495483   217  0.633094  0.5273
    ## abs(temptrend):REALMMarine       3.918334 1.6868767 23040  2.322834  0.0202
    ## abs(temptrend):REALMTerrestrial  4.191470 1.7269449 23040  2.427101  0.0152
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.446                                
    ## REALMMarine                     -0.945  0.421                         
    ## REALMTerrestrial                -0.920  0.410  0.870                  
    ## abs(temptrend):REALMMarine       0.435 -0.975 -0.447 -0.400           
    ## abs(temptrend):REALMTerrestrial  0.425 -0.952 -0.401 -0.451  0.929    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.256505e+00 -1.265922e-02  1.461626e-06  2.322678e-02  1.576644e+00 
    ## 
    ## Number of Observations: 23263
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    220                  23263

### Summary: Half-saturation MM fits

``` r
summary(modonlyTtrendJtumm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC   logLik
    ##   106495.8 106597.7 -53235.9
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)    0.8077487 (Intr)
    ## abs(temptrend) 2.0197857 -0.95 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:    1.060674 0.2603565
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -6.504057 
    ## Fixed effects: log(Jtumm + 1) ~ abs(temptrend) * REALM 
    ##                                      Value Std.Error    DF   t-value p-value
    ## (Intercept)                      1.2469474 0.2642367 35624  4.719055  0.0000
    ## abs(temptrend)                  -2.1831560 1.0045077 35624 -2.173359  0.0298
    ## REALMMarine                     -0.5286037 0.2810840   254 -1.880590  0.0612
    ## REALMTerrestrial                 0.2427881 0.2859494   254  0.849060  0.3966
    ## abs(temptrend):REALMMarine       1.4043803 1.0474772 35624  1.340726  0.1800
    ## abs(temptrend):REALMTerrestrial -0.1519933 1.0773398 35624 -0.141082  0.8878
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.713                                
    ## REALMMarine                     -0.940  0.670                         
    ## REALMTerrestrial                -0.924  0.659  0.869                  
    ## abs(temptrend):REALMMarine       0.683 -0.959 -0.721 -0.632           
    ## abs(temptrend):REALMTerrestrial  0.665 -0.932 -0.625 -0.718  0.894    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -2.020779e-02 -6.743656e-06 -6.577972e-08  1.102194e-06  7.561373e-04 
    ## 
    ## Number of Observations: 35884
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    257                  35884

``` r
summary(modonlyTtrendJbetamm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC     BIC    logLik
    ##   75562.64 75664.5 -37769.32
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)    0.8884553 (Intr)
    ## abs(temptrend) 1.7297848 -0.972
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:   0.4492306 0.7432453
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.2098961 
    ## Fixed effects: log(Jbetamm + 1) ~ abs(temptrend) * REALM 
    ##                                      Value Std.Error    DF   t-value p-value
    ## (Intercept)                      1.2546199 0.2336348 35656  5.370006  0.0000
    ## abs(temptrend)                  -2.3026681 0.6987494 35656 -3.295413  0.0010
    ## REALMMarine                     -0.9064721 0.2533266   254 -3.578274  0.0004
    ## REALMTerrestrial                -0.1271445 0.2519655   254 -0.504611  0.6143
    ## abs(temptrend):REALMMarine       1.9840897 0.7339632 35656  2.703255  0.0069
    ## abs(temptrend):REALMTerrestrial  0.7690663 0.7462662 35656  1.030552  0.3028
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.715                                
    ## REALMMarine                     -0.922  0.660                         
    ## REALMTerrestrial                -0.927  0.663  0.855                  
    ## abs(temptrend):REALMMarine       0.681 -0.952 -0.729 -0.631           
    ## abs(temptrend):REALMTerrestrial  0.670 -0.936 -0.618 -0.721  0.891    
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -4.3629146 -0.4699047 -0.2537973  0.2518276  5.0759496 
    ## 
    ## Number of Observations: 35916
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    257                  35916

``` r
summary(modonlyTtrendHornmm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##      AIC      BIC    logLik
    ##   112157 112258.8 -56066.48
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev   Corr  
    ## (Intercept)    1.195801 (Intr)
    ## abs(temptrend) 1.395683 -0.945
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.099728 13.50166
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.794175 
    ## Fixed effects: log(Hornmm + 1) ~ abs(temptrend) * REALM 
    ##                                      Value Std.Error    DF   t-value p-value
    ## (Intercept)                      1.4549714 0.3356369 35617  4.334956  0.0000
    ## abs(temptrend)                  -2.7052076 0.9172299 35617 -2.949323  0.0032
    ## REALMMarine                     -0.7159379 0.3614531   254 -1.980721  0.0487
    ## REALMTerrestrial                 0.4560898 0.3618769   254  1.260345  0.2087
    ## abs(temptrend):REALMMarine       2.3639551 0.9458474 35617  2.499299  0.0124
    ## abs(temptrend):REALMTerrestrial  1.1245810 0.9605334 35617  1.170788  0.2417
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.519                                
    ## REALMMarine                     -0.929  0.482                         
    ## REALMTerrestrial                -0.927  0.482  0.861                  
    ## abs(temptrend):REALMMarine       0.504 -0.970 -0.531 -0.467           
    ## abs(temptrend):REALMTerrestrial  0.496 -0.955 -0.460 -0.529  0.926    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -1.04398348 -0.09256434 -0.01407641  0.02435333  2.54821772 
    ## 
    ## Number of Observations: 35877
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    257                  35877

### Plot the temp-only coefficients

``` r
colors <- brewer.pal(3, 'Dark2')

# make table of coefficients
coefs1 <- as.data.frame(summary(modonlyTtrendrem0)$tTable)
coefs2 <- as.data.frame(summary(modonlyTtrendJbetarem0)$tTable)
coefs3 <- as.data.frame(summary(modonlyTtrendHornrem0)$tTable)
coefs1$mod <- 'Jtu'
coefs2$mod <- 'Jbeta'
coefs3$mod <- 'Horn'
rows1 <- which(grepl('temptrend', rownames(coefs1))) # extract temperature effect
cols <- c('Value', 'Std.Error', 'mod')
allcoefs <- rbind(coefs1[rows1, cols], coefs2[rows1, cols], coefs3[rows1, cols])
allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to marine effects
allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to terrestrial effects

allcoefs$lCI <- allcoefs$Value - allcoefs$Std.Error # lower confidence interval
allcoefs$uCI <- allcoefs$Value + allcoefs$Std.Error
allcoefs$y <- c(3, 2, 1) + rep(c(0, -0.1, -0.2), c(3, 3, 3)) # y-values
allcoefs$col <- c(rep(colors[1], 3), rep(colors[2], 3), rep(colors[3], 3))
allcoefs$realm <- rep(c('Freshwater', 'Marine', 'Terrestrial'), 3)

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
       legend = c('Jaccard turnover rem0', 'Jaccard total rem0', 'Horn-Morisita rem0'))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/modonlyTtrendsimp%20coefs-1.png)<!-- -->

## Full models

Try static covariates plus interactions of abs temperature trend with
each covariate:

  - realm
  - speed
  - mass
  - average metabolic temperature
  - consumer fraction
  - environmental temperature
  - seasonality
  - microclimates
  - thermal bias
  - NPP
  - vegetation
  - duration
  - human footprint

Except for thermal bias: interact with temperature trend (not abs)

### Fit/load full models

#### Bowler vs Venter/Halpern human impact

Bowler has lower AIC.

``` r
# using Bowler for human impact
i1 <- trends[, complete.cases(Jtutrendrem0, REALM, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc, human_footprint.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modTfullbowlerrem0.rds')){
  modTfullbowlerrem0 <- readRDS('temp/modTfullbowlerrem0.rds')
} else {

modTfullbowlerrem0 <- lme(Jtutrendrem0 ~ temptrend_abs.sc*REALM +
                     temptrend_abs.sc*tsign + 
                     temptrend_abs.sc*tempave_metab.sc + 
                     temptrend_abs.sc*seas.sc + 
                     temptrend_abs.sc*microclim.sc + 
                     temptrend_abs.sc*mass.sc + 
                     temptrend_abs.sc*speed.sc + 
                     temptrend_abs.sc*consumerfrac.sc +
                     temptrend_abs.sc*nspp.sc +
                     temptrend_abs.sc*thermal_bias.sc:tsign +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*veg.sc +
                     temptrend_abs.sc*duration.sc +
                     temptrend_abs.sc*human_bowler.sc:REALM2,
                   random = randef, weights = varef, data = trends[i1,], method = 'REML')
  saveRDS(modTfullbowlerrem0, file = 'temp/modTfullbowlerrem0.rds')
}

# using Venter/Halpern for human impact
if(file.exists('temp/modTfullfootprintrem0.rds')){
  modTfullfootprintrem0 <- readRDS('temp/modTfullfootprintrem0.rds')
} else {
modTfullfootprintrem0 <- lme(Jtutrendrem0 ~ temptrend_abs.sc*REALM + 
                     temptrend_abs.sc*tsign + 
                     temptrend_abs.sc*tempave_metab.sc + 
                     temptrend_abs.sc*seas.sc + 
                     temptrend_abs.sc*microclim.sc + 
                     temptrend_abs.sc*mass.sc + 
                     temptrend_abs.sc*speed.sc + 
                     temptrend_abs.sc*consumerfrac.sc +
                     temptrend_abs.sc*nspp.sc +
                     temptrend_abs.sc*thermal_bias.sc:tsign +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*veg.sc +
                     temptrend_abs.sc*duration.sc +
                     temptrend_abs.sc*human_footprint.sc:REALM2,
                   random = randef, weights = varef, data = trends[i1,], method = 'REML',
                   control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullfootprintrem0, file = 'temp/modTfullfootprintrem0.rds')

}
AIC(modTfullbowlerrem0, modTfullfootprintrem0)
```

    ##                       df       AIC
    ## modTfullbowlerrem0    42 -101579.8
    ## modTfullfootprintrem0 42 -101568.7

#### Fit/load full models

``` r
randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# rem0
if(file.exists('temp/modTfullJturem0.rds')){
  modTfullJturem0 <- readRDS('temp/modTfullJturem0.rds')
} else {
  i <- trends[, complete.cases(Jtutrendrem0, REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]
  modTfullJturem0 <- lme(Jtutrendrem0 ~ temptrend_abs.sc*REALM + 
                           temptrend_abs.sc*tsign +
                           temptrend_abs.sc*tempave_metab.sc + 
                           temptrend_abs.sc*seas.sc + 
                           temptrend_abs.sc*microclim.sc + 
                           temptrend_abs.sc*mass.sc + 
                           temptrend_abs.sc*speed.sc + 
                           temptrend_abs.sc*consumerfrac.sc +
                           temptrend_abs.sc*nspp.sc +
                           temptrend_abs.sc*thermal_bias.sc:tsign +
                           temptrend_abs.sc*npp.sc +
                           temptrend_abs.sc*veg.sc +
                           temptrend_abs.sc*duration.sc +
                           temptrend_abs.sc*human_bowler.sc:REALM2,
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML')
  saveRDS(modTfullJturem0, file = 'temp/modTfullJturem0.rds')
}

if(file.exists('temp/modTfullJbetarem0.rds')){
  modTfullJbetarem0 <- readRDS('temp/modTfullJbetarem0.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendrem0, REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]
  modTfullJbetarem0 <- lme(Jbetatrendrem0 ~ temptrend_abs.sc*REALM + 
                             temptrend_abs.sc*tsign +
                             temptrend_abs.sc*tempave_metab.sc + 
                             temptrend_abs.sc*seas.sc + 
                             temptrend_abs.sc*microclim.sc + 
                             temptrend_abs.sc*mass.sc + 
                             temptrend_abs.sc*speed.sc + 
                             temptrend_abs.sc*consumerfrac.sc +
                             temptrend_abs.sc*nspp.sc +
                             temptrend_abs.sc*thermal_bias.sc:tsign +
                             temptrend_abs.sc*npp.sc +
                             temptrend_abs.sc*veg.sc +
                             temptrend_abs.sc*duration.sc +
                             temptrend_abs.sc*human_bowler.sc:REALM2,
                           random = randef, weights = varef, data = trends[i,], 
                           method = 'REML')
  saveRDS(modTfullJbetarem0, file = 'temp/modTfullJbetarem0.rds')
}

if(file.exists('temp/modTfullHornrem0.rds')){
  modTfullHornrem0 <- readRDS('temp/modTfullHornrem0.rds')
} else {
  i <- trends[, complete.cases(Horntrendrem0, REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]
  modTfullHornrem0 <- lme(Horntrendrem0 ~ temptrend_abs.sc*REALM + 
                        temptrend_abs.sc*tsign +
                        temptrend_abs.sc*tempave_metab.sc + 
                        temptrend_abs.sc*seas.sc + 
                        temptrend_abs.sc*microclim.sc + 
                        temptrend_abs.sc*mass.sc + 
                        temptrend_abs.sc*speed.sc + 
                        temptrend_abs.sc*consumerfrac.sc +
                        temptrend_abs.sc*nspp.sc +
                        temptrend_abs.sc*thermal_bias.sc:tsign +
                        temptrend_abs.sc*npp.sc +
                        temptrend_abs.sc*veg.sc +
                        temptrend_abs.sc*duration.sc +
                        temptrend_abs.sc*human_bowler.sc:REALM2,
                      random = randef, weights = varef, data = trends[i,], 
                      method = 'REML')
  saveRDS(modTfullHornrem0, file = 'temp/modTfullHornrem0.rds')
}

# exp
if(file.exists('temp/modTfullJtuexp.rds')){
  modTfullJtuexp <- readRDS('temp/modTfullJtuexp.rds')
} else {
  i <- trends[, complete.cases(Jtuexp, REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]
  modTfullJtuexp <- lme(log(Jtuexp) ~ temptrend_abs.sc*REALM + 
                           temptrend_abs.sc*tsign +
                           temptrend_abs.sc*tempave_metab.sc + 
                           temptrend_abs.sc*seas.sc + 
                           temptrend_abs.sc*microclim.sc + 
                           temptrend_abs.sc*mass.sc + 
                           temptrend_abs.sc*speed.sc + 
                           temptrend_abs.sc*consumerfrac.sc +
                           temptrend_abs.sc*nspp.sc +
                           temptrend_abs.sc*thermal_bias.sc:tsign +
                           temptrend_abs.sc*npp.sc +
                           temptrend_abs.sc*veg.sc +
                           temptrend_abs.sc*duration.sc +
                           temptrend_abs.sc*human_bowler.sc:REALM2,
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                   control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))
  saveRDS(modTfullJtuexp, file = 'temp/modTfullJtuexp.rds')
}

if(file.exists('temp/modTfullJbetaexp.rds')){
  modTfullJbetaexp <- readRDS('temp/modTfullJbetaexp.rds')
} else {
  i <- trends[, complete.cases(Jbetaexp, REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]
  modTfullJbetaexp <- lme(log(Jbetaexp) ~ temptrend_abs.sc*REALM + 
                             temptrend_abs.sc*tsign +
                             temptrend_abs.sc*tempave_metab.sc + 
                             temptrend_abs.sc*seas.sc + 
                             temptrend_abs.sc*microclim.sc + 
                             temptrend_abs.sc*mass.sc + 
                             temptrend_abs.sc*speed.sc + 
                             temptrend_abs.sc*consumerfrac.sc +
                             temptrend_abs.sc*nspp.sc +
                             temptrend_abs.sc*thermal_bias.sc:tsign +
                             temptrend_abs.sc*npp.sc +
                             temptrend_abs.sc*veg.sc +
                             temptrend_abs.sc*duration.sc +
                             temptrend_abs.sc*human_bowler.sc:REALM2,
                           random = randef, weights = varef, data = trends[i,], 
                           method = 'REML')
  saveRDS(modTfullJbetaexp, file = 'temp/modTfullJbetaexp.rds')
}

if(file.exists('temp/modTfullHornexp.rds')){
  modTfullHornexp <- readRDS('temp/modTfullHornexp.rds')
} else {
  i <- trends[, complete.cases(Hornexp, REALM, tempave_metab.sc, 
                               seas.sc, microclim.sc, 
                               temptrend_abs.sc, mass.sc, speed.sc, 
                               consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                               veg.sc, duration.sc, human_bowler.sc)]
  modTfullHornexp <- lme(log(Hornexp) ~ temptrend_abs.sc*REALM + 
                           temptrend_abs.sc*tsign +
                           temptrend_abs.sc*tempave_metab.sc + 
                           temptrend_abs.sc*seas.sc + 
                           temptrend_abs.sc*microclim.sc + 
                           temptrend_abs.sc*mass.sc + 
                           temptrend_abs.sc*speed.sc + 
                           temptrend_abs.sc*consumerfrac.sc +
                           temptrend_abs.sc*nspp.sc +
                           temptrend_abs.sc*thermal_bias.sc:tsign +
                           temptrend_abs.sc*npp.sc +
                           temptrend_abs.sc*veg.sc +
                           temptrend_abs.sc*duration.sc +
                           temptrend_abs.sc*human_bowler.sc:REALM2,
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                         control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullHornexp, file = 'temp/modTfullHornexp.rds')
}

# mm
if(file.exists('temp/modTfullJtumm.rds')){
  modTfullJtumm <- readRDS('temp/modTfullJtumm.rds')
} else {
  i <- trends[, complete.cases(Jtumm, REALM, tempave_metab.sc, 
                               seas.sc, microclim.sc, 
                               temptrend_abs.sc, mass.sc, speed.sc, 
                               consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                               veg.sc, duration.sc, human_bowler.sc)]
  modTfullJtumm <- lme(log(Jtumm+1) ~ temptrend_abs.sc*REALM + 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*consumerfrac.sc +
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*veg.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc:REALM2,
                       random = randef, weights = varef, data = trends[i,], 
                       method = 'REML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullJtumm, file = 'temp/modTfullJtumm.rds')
}

if(file.exists('temp/modTfullJbetamm.rds')){
  modTfullJbetamm <- readRDS('temp/modTfullJbetamm.rds')
} else {
  i <- trends[, complete.cases(Jbetamm, REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                               temptrend_abs.sc, mass.sc, speed.sc, 
                               consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                               veg.sc, duration.sc, human_bowler.sc)]
  modTfullJbetamm <- lme(log(Jbetamm+1) ~ temptrend_abs.sc*REALM + 
                           temptrend_abs.sc*tsign +
                           temptrend_abs.sc*tempave_metab.sc + 
                           temptrend_abs.sc*seas.sc + 
                           temptrend_abs.sc*microclim.sc + 
                           temptrend_abs.sc*mass.sc + 
                           temptrend_abs.sc*speed.sc + 
                           temptrend_abs.sc*consumerfrac.sc +
                           temptrend_abs.sc*nspp.sc +
                           temptrend_abs.sc*thermal_bias.sc:tsign +
                           temptrend_abs.sc*npp.sc +
                           temptrend_abs.sc*veg.sc +
                           temptrend_abs.sc*duration.sc +
                           temptrend_abs.sc*human_bowler.sc:REALM2,
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                         control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullJbetamm, file = 'temp/modTfullJbetamm.rds')
}

if(file.exists('temp/modTfullHornmm.rds')){
  modTfullHornmm <- readRDS('temp/modTfullHornmm.rds')
} else {
  i <- trends[, complete.cases(Hornmm, REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]
  modTfullHornmm <- lme(log(Hornmm+1) ~ temptrend_abs.sc*REALM + 
                        temptrend_abs.sc*tsign +
                        temptrend_abs.sc*tempave_metab.sc + 
                        temptrend_abs.sc*seas.sc + 
                        temptrend_abs.sc*microclim.sc + 
                        temptrend_abs.sc*mass.sc + 
                        temptrend_abs.sc*speed.sc + 
                        temptrend_abs.sc*consumerfrac.sc +
                        temptrend_abs.sc*nspp.sc +
                        temptrend_abs.sc*thermal_bias.sc:tsign +
                        temptrend_abs.sc*npp.sc +
                        temptrend_abs.sc*veg.sc +
                        temptrend_abs.sc*duration.sc +
                        temptrend_abs.sc*human_bowler.sc:REALM2,
                        random = randef, weights = varef, data = trends[i,], 
                        method = 'REML',
                        control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullHornmm, file = 'temp/modTfullHornmm.rds')
}
```

#### Summary rem0

``` r
summary(modTfullJturem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -102597.3 -102240.7 51340.65
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.009402743 (Intr)
    ## temptrend_abs.sc 0.024575808 -0.966
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01091239 2.006389
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.129145 
    ## Fixed effects: Jtutrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.001950838 0.005305143 35753  0.367726  0.7131
    ## temptrend_abs.sc                                  0.007653289 0.013585719 35753  0.563333  0.5732
    ## REALMMarine                                       0.007417093 0.005640906   228  1.314876  0.1899
    ## REALMTerrestrial                                  0.002223840 0.004973589   228  0.447130  0.6552
    ## tsign1                                           -0.002302113 0.000641868 35753 -3.586585  0.0003
    ## tempave_metab.sc                                 -0.003624603 0.000910719 35753 -3.979936  0.0001
    ## seas.sc                                           0.000489656 0.000656766 35753  0.745555  0.4559
    ## microclim.sc                                      0.000626108 0.000335958 35753  1.863648  0.0624
    ## mass.sc                                          -0.001096185 0.000535195 35753 -2.048195  0.0405
    ## speed.sc                                          0.001324237 0.000615458 35753  2.151628  0.0314
    ## consumerfrac.sc                                  -0.000531934 0.000321813 35753 -1.652930  0.0984
    ## nspp.sc                                          -0.000190341 0.000540015 35753 -0.352473  0.7245
    ## npp.sc                                           -0.000270938 0.000464632 35753 -0.583123  0.5598
    ## veg.sc                                            0.000271663 0.000922497 35753  0.294487  0.7684
    ## duration.sc                                      -0.000523613 0.000496888 35753 -1.053784  0.2920
    ## temptrend_abs.sc:REALMMarine                     -0.010742575 0.014318824 35753 -0.750241  0.4531
    ## temptrend_abs.sc:REALMTerrestrial                 0.008041095 0.013532164 35753  0.594221  0.5524
    ## temptrend_abs.sc:tsign1                          -0.001257613 0.001853754 35753 -0.678414  0.4975
    ## temptrend_abs.sc:tempave_metab.sc                 0.006417096 0.002628940 35753  2.440944  0.0147
    ## temptrend_abs.sc:seas.sc                         -0.000026729 0.001641484 35753 -0.016283  0.9870
    ## temptrend_abs.sc:microclim.sc                    -0.003297789 0.000961165 35753 -3.431032  0.0006
    ## temptrend_abs.sc:mass.sc                         -0.000629142 0.001398167 35753 -0.449976  0.6527
    ## temptrend_abs.sc:speed.sc                        -0.001808939 0.001786365 35753 -1.012637  0.3112
    ## temptrend_abs.sc:consumerfrac.sc                  0.004202315 0.000760913 35753  5.522729  0.0000
    ## temptrend_abs.sc:nspp.sc                          0.000362491 0.001460092 35753  0.248266  0.8039
    ## tsign-1:thermal_bias.sc                          -0.000913472 0.000659212 35753 -1.385703  0.1658
    ## tsign1:thermal_bias.sc                           -0.000315733 0.000463084 35753 -0.681805  0.4954
    ## temptrend_abs.sc:npp.sc                           0.004476882 0.001370542 35753  3.266506  0.0011
    ## temptrend_abs.sc:veg.sc                          -0.003792672 0.002175542 35753 -1.743323  0.0813
    ## temptrend_abs.sc:duration.sc                     -0.001195179 0.001018796 35753 -1.173129  0.2408
    ## human_bowler.sc:REALM2TerrFresh                   0.000326064 0.000726739 35753  0.448667  0.6537
    ## human_bowler.sc:REALM2Marine                      0.000422750 0.000400649 35753  1.055164  0.2914
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.000597994 0.001407691 35753  0.424805  0.6710
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.002208910 0.001182941 35753  1.867303  0.0619
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.002730704 0.001484981 35753 -1.838882  0.0659
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.003910491 0.001102998 35753 -3.545329  0.0004
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.797                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.955  0.757                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.683  0.574  0.634                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.095  0.063  0.007  0.012                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.061 -0.037 -0.053 -0.242  0.074                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.150  0.143  0.202 -0.097 -0.070  0.138                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.068  0.047  0.077  0.006  0.015 -0.157  0.097                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.061 -0.042 -0.036 -0.038 -0.019  0.100  0.127 -0.001                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.107 -0.093 -0.057 -0.084 -0.073 -0.194 -0.055  0.070 -0.437                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.028  0.022  0.017  0.133  0.014 -0.135 -0.055  0.017  0.016 -0.144                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.020  0.000 -0.074 -0.066  0.049 -0.154 -0.037 -0.113 -0.010  0.184  0.094                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.041 -0.044 -0.061  0.054  0.041  0.037 -0.136 -0.187 -0.027  0.125 -0.045 -0.188                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.578  0.387  0.599  0.024 -0.006  0.007  0.080 -0.010  0.021 -0.019 -0.010  0.024 -0.178                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.037  0.008  0.032  0.016 -0.160  0.085 -0.006 -0.016 -0.057 -0.009  0.005 -0.252  0.062 -0.012                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.764 -0.951 -0.804 -0.534 -0.014  0.030 -0.184 -0.056  0.020  0.047 -0.015  0.046  0.064 -0.406 -0.018                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.542 -0.693 -0.499 -0.820 -0.018  0.182  0.078  0.024  0.033  0.068 -0.104  0.025 -0.047  0.019 -0.005  0.639                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.067 -0.140 -0.022 -0.014 -0.507 -0.056  0.012 -0.006 -0.005  0.041 -0.007 -0.029 -0.038 -0.018  0.031  0.044      0.011                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.037  0.042  0.028  0.171 -0.058 -0.788 -0.106  0.091 -0.066  0.134  0.116  0.154 -0.031 -0.012 -0.031 -0.027     -0.228      0.052                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.145 -0.158 -0.183  0.095  0.053 -0.100 -0.777 -0.099 -0.081  0.032  0.029  0.038  0.123 -0.109 -0.027  0.207     -0.162     -0.021  0.109                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.038 -0.062 -0.047  0.020 -0.014  0.088 -0.068 -0.760  0.028 -0.055  0.002  0.074  0.183  0.006  0.011  0.077     -0.051      0.020 -0.040   0.113                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.040  0.048  0.010  0.020 -0.002 -0.086 -0.066  0.019 -0.701  0.315  0.008  0.031 -0.007 -0.007  0.044 -0.014     -0.004      0.001  0.168   0.052            -0.043                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.084  0.104  0.040  0.073  0.069  0.167  0.031 -0.041  0.266 -0.717  0.085 -0.165 -0.064 -0.006  0.001 -0.055     -0.094     -0.035 -0.331  -0.003             0.043            -0.475                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.030 -0.037 -0.018 -0.117 -0.007  0.105  0.038  0.009  0.006  0.120 -0.773 -0.069  0.022  0.005 -0.035  0.021      0.127      0.032 -0.124  -0.063            -0.027            -0.043            -0.120                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                          0.009  0.008  0.038  0.037 -0.033  0.143  0.035  0.069  0.033 -0.159 -0.046 -0.749  0.118 -0.037  0.195 -0.070     -0.068      0.029 -0.138  -0.040            -0.043             0.003             0.183             0.085                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.061 -0.041 -0.057 -0.046 -0.078  0.132 -0.217 -0.015 -0.043 -0.006 -0.001 -0.035  0.001 -0.012  0.083  0.050      0.034     -0.042 -0.114   0.159             0.007             0.021             0.018            -0.005            0.018                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.095 -0.074 -0.103 -0.097  0.063  0.223 -0.330 -0.100 -0.043 -0.009 -0.024 -0.044 -0.057  0.000  0.084  0.073      0.084      0.001 -0.216   0.205             0.042             0.007             0.030             0.020            0.032             0.357                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.039  0.065  0.059 -0.051 -0.028 -0.007  0.116  0.168 -0.004 -0.064  0.004  0.113 -0.741  0.143 -0.043 -0.103      0.085      0.003  0.030  -0.182            -0.312             0.039             0.037            -0.017           -0.114            -0.002  0.055                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.411 -0.507 -0.431  0.022 -0.008 -0.018 -0.115  0.003 -0.018  0.010  0.011 -0.024  0.173 -0.733  0.003  0.537     -0.033      0.052  0.004   0.149             0.028             0.004             0.000            -0.011            0.017             0.020 -0.012 -0.250                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.011  0.042 -0.036 -0.049  0.036  0.028 -0.034 -0.013  0.064 -0.027 -0.008  0.267 -0.043  0.019 -0.519  0.031      0.039     -0.130  0.016   0.012             0.034            -0.073            -0.002             0.026           -0.359            -0.036 -0.086  0.020             0.001                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.126  0.159  0.132 -0.048  0.020  0.028  0.050  0.092  0.037 -0.039 -0.045 -0.014 -0.112  0.230 -0.022 -0.166      0.031     -0.026 -0.025  -0.115            -0.012             0.001             0.011             0.030            0.028             0.052  0.097  0.061            -0.289            0.012                                                               
    ## human_bowler.sc:REALM2Marine                      0.053 -0.043 -0.050 -0.038 -0.009  0.079 -0.197  0.019 -0.029  0.062 -0.076 -0.020 -0.193  0.025  0.003  0.041      0.027      0.000 -0.051   0.152            -0.017             0.029            -0.054             0.077            0.029             0.111  0.145  0.135            -0.015            0.024            0.027                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.044  0.041  0.056  0.036 -0.072 -0.146  0.203  0.057  0.002  0.021  0.008  0.015 -0.021  0.011 -0.060 -0.068     -0.053      0.238  0.202  -0.240            -0.089             0.046            -0.087             0.012           -0.002            -0.503 -0.329  0.040            -0.027            0.043           -0.051      -0.085                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.059  0.074  0.061  0.063 -0.009 -0.232  0.208  0.076  0.008  0.030  0.008  0.055  0.050 -0.011 -0.054 -0.072     -0.088     -0.010  0.340  -0.251            -0.170             0.048            -0.159             0.014           -0.056            -0.279 -0.720 -0.030             0.010            0.057           -0.064      -0.117       0.435                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.208 -0.278 -0.219  0.032 -0.011 -0.060 -0.138 -0.049 -0.019  0.036  0.028  0.021  0.108 -0.341  0.020  0.293     -0.044      0.051  0.071   0.184             0.064             0.013            -0.050            -0.031           -0.035            -0.037 -0.060 -0.151             0.495           -0.017           -0.663      -0.004       0.081  0.077               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.041  0.048  0.039  0.034 -0.003 -0.077  0.141 -0.005  0.033 -0.060  0.134  0.026  0.143 -0.013  0.005 -0.045     -0.033     -0.013  0.094  -0.171            -0.034            -0.020             0.042            -0.142           -0.039            -0.075 -0.117 -0.181             0.028           -0.029           -0.012      -0.775       0.083  0.124  0.005        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.64180587 -0.23617979 -0.02235224  0.26661608  5.42260805 
    ## 
    ## Number of Observations: 36017
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    231                  36017

``` r
summary(modTfullJbetarem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##         AIC     BIC   logLik
    ##   -132984.6 -132628 66534.29
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.00722145 (Intr)
    ## temptrend_abs.sc 0.01729128 -0.058
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.003744118 0.9597427
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.892596 
    ## Fixed effects: Jbetatrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.001213002 0.003525875 35753  0.344029  0.7308
    ## temptrend_abs.sc                                  0.024516024 0.009529501 35753  2.572645  0.0101
    ## REALMMarine                                       0.007673212 0.003735223   228  2.054285  0.0411
    ## REALMTerrestrial                                  0.007744382 0.003640043   228  2.127552  0.0344
    ## tsign1                                           -0.001148029 0.000395440 35753 -2.903171  0.0037
    ## tempave_metab.sc                                 -0.003215071 0.000595344 35753 -5.400356  0.0000
    ## seas.sc                                           0.000212133 0.000380001 35753  0.558244  0.5767
    ## microclim.sc                                      0.000441077 0.000203585 35753  2.166543  0.0303
    ## mass.sc                                          -0.001423132 0.000366399 35753 -3.884106  0.0001
    ## speed.sc                                          0.000645357 0.000391480 35753  1.648506  0.0993
    ## consumerfrac.sc                                  -0.000038161 0.000199882 35753 -0.190917  0.8486
    ## nspp.sc                                          -0.000444728 0.000347046 35753 -1.281467  0.2000
    ## npp.sc                                           -0.000443237 0.000279941 35753 -1.583321  0.1134
    ## veg.sc                                            0.000532127 0.000482030 35753  1.103929  0.2696
    ## duration.sc                                      -0.001405260 0.000323951 35753 -4.337878  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.026279282 0.009987226 35753 -2.631289  0.0085
    ## temptrend_abs.sc:REALMTerrestrial                -0.010579909 0.009882467 35753 -1.070574  0.2844
    ## temptrend_abs.sc:tsign1                          -0.001015758 0.001245579 35753 -0.815491  0.4148
    ## temptrend_abs.sc:tempave_metab.sc                 0.007723972 0.001866659 35753  4.137860  0.0000
    ## temptrend_abs.sc:seas.sc                         -0.000384825 0.001024830 35753 -0.375501  0.7073
    ## temptrend_abs.sc:microclim.sc                    -0.001965408 0.000612832 35753 -3.207090  0.0013
    ## temptrend_abs.sc:mass.sc                          0.001166200 0.000972457 35753  1.199230  0.2304
    ## temptrend_abs.sc:speed.sc                        -0.001269077 0.001204881 35753 -1.053280  0.2922
    ## temptrend_abs.sc:consumerfrac.sc                  0.001020976 0.000499736 35753  2.043029  0.0411
    ## temptrend_abs.sc:nspp.sc                          0.001421707 0.000978638 35753  1.452741  0.1463
    ## tsign-1:thermal_bias.sc                          -0.000496281 0.000425242 35753 -1.167055  0.2432
    ## tsign1:thermal_bias.sc                           -0.000216836 0.000288028 35753 -0.752832  0.4516
    ## temptrend_abs.sc:npp.sc                           0.004035336 0.000892803 35753  4.519850  0.0000
    ## temptrend_abs.sc:veg.sc                          -0.003801762 0.001356247 35753 -2.803148  0.0051
    ## temptrend_abs.sc:duration.sc                     -0.000608213 0.000663935 35753 -0.916072  0.3596
    ## human_bowler.sc:REALM2TerrFresh                   0.000589613 0.000370040 35753  1.593377  0.1111
    ## human_bowler.sc:REALM2Marine                     -0.000185665 0.000243705 35753 -0.761845  0.4462
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.000148615 0.000943342 35753 -0.157541  0.8748
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.001952104 0.000768829 35753  2.539062  0.0111
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.001846694 0.000906396 35753 -2.037403  0.0416
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.000500133 0.000712112 35753 -0.702324  0.4825
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.573                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.947  0.551                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.752  0.385  0.701                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.084  0.047  0.004  0.001                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.100 -0.026 -0.084 -0.244  0.098                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.131  0.077  0.178 -0.069 -0.074  0.088                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.078  0.060  0.090  0.022 -0.001 -0.183  0.138                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.095 -0.010 -0.071 -0.043  0.001  0.089  0.091  0.007                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.081 -0.063 -0.041 -0.088 -0.086 -0.110 -0.014  0.039 -0.493                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.012  0.017  0.010  0.096 -0.022 -0.105 -0.031  0.020  0.012 -0.089                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.007 -0.032 -0.062 -0.038  0.049 -0.200 -0.033 -0.135 -0.054  0.212  0.100                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.025 -0.013 -0.040  0.032  0.066  0.073 -0.077 -0.140 -0.028  0.163 -0.052 -0.198                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.442  0.365  0.462  0.011 -0.014 -0.003  0.073 -0.003  0.018 -0.024 -0.018  0.023 -0.204                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.041  0.009  0.031 -0.003 -0.144  0.129 -0.003 -0.024 -0.042 -0.033  0.012 -0.240  0.077 -0.012                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.554 -0.951 -0.570 -0.361 -0.009  0.023 -0.099 -0.066  0.000  0.039 -0.010  0.063  0.038 -0.386 -0.008                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.377 -0.740 -0.350 -0.552  0.001  0.148  0.116  0.024 -0.011  0.076 -0.073  0.049 -0.055  0.021  0.007  0.690                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.055 -0.127 -0.022  0.002 -0.442 -0.055  0.003 -0.006 -0.004  0.033  0.004 -0.025 -0.056 -0.025  0.020  0.040      0.006                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.001  0.056 -0.009  0.132 -0.104 -0.539 -0.130  0.018 -0.066  0.142  0.085  0.123 -0.020 -0.009 -0.054 -0.041     -0.239      0.047                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.079 -0.108 -0.106  0.104  0.050 -0.090 -0.726 -0.122 -0.063  0.014  0.028  0.045  0.090 -0.086 -0.036  0.149     -0.173     -0.018  0.080                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.053 -0.078 -0.061  0.002 -0.009  0.112 -0.091 -0.780  0.012 -0.029 -0.008  0.083  0.167 -0.009  0.011  0.090     -0.058      0.009  0.029   0.157                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.018  0.040  0.012 -0.009 -0.004 -0.077 -0.056  0.024 -0.577  0.282  0.015  0.014 -0.015 -0.010  0.055 -0.003      0.016      0.009  0.140   0.042            -0.046                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.070  0.093  0.057  0.080  0.074  0.129  0.055 -0.019  0.240 -0.601  0.071 -0.148 -0.068  0.004  0.014 -0.042     -0.098     -0.031 -0.308   0.027             0.032            -0.446                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.010 -0.032 -0.009 -0.076  0.016  0.077  0.031  0.012  0.003  0.080 -0.717 -0.065  0.026  0.009 -0.023  0.019      0.109      0.018 -0.106  -0.070            -0.031            -0.050            -0.119                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.023  0.016  0.048  0.045 -0.026  0.119  0.035  0.070  0.011 -0.153 -0.045 -0.640  0.108 -0.028  0.167 -0.078     -0.068      0.024 -0.134  -0.035            -0.050            -0.002             0.166             0.084                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.066 -0.022 -0.062 -0.049 -0.068  0.135 -0.240 -0.022 -0.040  0.027  0.006 -0.048  0.004 -0.012  0.085  0.029      0.020     -0.046 -0.045   0.153             0.009             0.019            -0.008            -0.009            0.017                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.122 -0.044 -0.126 -0.104  0.066  0.252 -0.437 -0.133 -0.017  0.026 -0.030 -0.081 -0.042  0.000  0.105  0.035      0.057      0.016 -0.124   0.210             0.051             0.013            -0.009             0.019            0.025             0.340                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.002  0.057  0.015 -0.054 -0.040 -0.016  0.065  0.155  0.001 -0.085  0.016  0.109 -0.733  0.163 -0.056 -0.093      0.085      0.025  0.053  -0.193            -0.305             0.047             0.029            -0.021           -0.114             0.006  0.081                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.351 -0.443 -0.367  0.020 -0.001 -0.016 -0.088 -0.009 -0.010  0.006  0.013 -0.026  0.181 -0.813  0.012  0.474     -0.031      0.047 -0.014   0.141             0.056             0.005             0.012            -0.016            0.013             0.012 -0.026 -0.281                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.002  0.038 -0.035 -0.045  0.039  0.035 -0.020 -0.014  0.063 -0.007 -0.003  0.245 -0.049  0.013 -0.446  0.035      0.042     -0.120  0.027  -0.003             0.012            -0.074            -0.015             0.033           -0.374            -0.031 -0.071  0.035            -0.001                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.106  0.136  0.110 -0.035  0.018  0.015  0.029  0.119  0.019 -0.034 -0.036 -0.005 -0.141  0.260 -0.024 -0.148      0.041     -0.024 -0.012  -0.119            -0.049             0.003             0.005             0.028            0.025             0.047  0.119  0.092            -0.313            0.021                                                               
    ## human_bowler.sc:REALM2Marine                      0.054 -0.032 -0.050 -0.048 -0.012  0.102 -0.187  0.026 -0.026  0.091 -0.087 -0.009 -0.195  0.032  0.005  0.029      0.025      0.007 -0.038   0.142            -0.021             0.024            -0.076             0.085            0.018             0.114  0.139  0.145            -0.024            0.031            0.033                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.029  0.034  0.033  0.033 -0.081 -0.090  0.181  0.035 -0.001  0.021  0.001  0.015 -0.013  0.018 -0.072 -0.057     -0.050      0.260  0.221  -0.254            -0.069             0.032            -0.097             0.020           -0.014            -0.429 -0.282  0.041            -0.037            0.061           -0.056      -0.072                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.037  0.069  0.035  0.050 -0.027 -0.145  0.211  0.054  0.002  0.021  0.007  0.048  0.053 -0.016 -0.068 -0.065     -0.084     -0.014  0.364  -0.303            -0.137             0.027            -0.151             0.017           -0.064            -0.223 -0.658 -0.019             0.001            0.069           -0.058      -0.099       0.446                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.154 -0.222 -0.164  0.035 -0.013 -0.037 -0.103 -0.082 -0.010  0.025  0.023  0.014  0.116 -0.361  0.024  0.236     -0.047      0.044  0.060   0.152             0.111             0.001            -0.037            -0.028           -0.023            -0.030 -0.066 -0.171             0.474           -0.022           -0.752      -0.011       0.082  0.074               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.039  0.039  0.037  0.037  0.002 -0.074  0.150 -0.003  0.025 -0.072  0.133  0.006  0.147 -0.023  0.006 -0.035     -0.031     -0.012  0.083  -0.160            -0.042            -0.015             0.056            -0.132           -0.028            -0.080 -0.110 -0.171             0.035           -0.027           -0.022      -0.766       0.075  0.110  0.010        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.38434700 -0.32179127 -0.03059932  0.31554016  8.33570233 
    ## 
    ## Number of Observations: 36017
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    231                  36017

``` r
summary(modTfullHornrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -98048.08 -97692.28 49066.04
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.01556817 (Intr)
    ## temptrend_abs.sc 0.02706589 -0.065
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01880474 2.445174
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.327824 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.00442119 0.008034858 35094  0.550252  0.5822
    ## temptrend_abs.sc                                  0.03434598 0.015794007 35094  2.174621  0.0297
    ## REALMMarine                                       0.00621500 0.008542686   197  0.727523  0.4678
    ## REALMTerrestrial                                  0.01359399 0.007940609   197  1.711958  0.0885
    ## tsign1                                           -0.00198702 0.000770042 35094 -2.580409  0.0099
    ## tempave_metab.sc                                 -0.00808278 0.001213522 35094 -6.660593  0.0000
    ## seas.sc                                          -0.00227420 0.000865427 35094 -2.627835  0.0086
    ## microclim.sc                                      0.00024108 0.000406888 35094  0.592507  0.5535
    ## mass.sc                                          -0.00122157 0.000671530 35094 -1.819079  0.0689
    ## speed.sc                                          0.00105412 0.000806486 35094  1.307056  0.1912
    ## consumerfrac.sc                                   0.00066311 0.000411324 35094  1.612132  0.1069
    ## nspp.sc                                           0.00043209 0.000680545 35094  0.634918  0.5255
    ## npp.sc                                           -0.00144005 0.000580318 35094 -2.481479  0.0131
    ## veg.sc                                            0.00085712 0.001262366 35094  0.678978  0.4972
    ## duration.sc                                      -0.00231804 0.000624845 35094 -3.709790  0.0002
    ## temptrend_abs.sc:REALMMarine                     -0.02276829 0.016595284 35094 -1.371973  0.1701
    ## temptrend_abs.sc:REALMTerrestrial                -0.00212318 0.017002259 35094 -0.124876  0.9006
    ## temptrend_abs.sc:tsign1                          -0.00540794 0.001924297 35094 -2.810344  0.0050
    ## temptrend_abs.sc:tempave_metab.sc                 0.00327068 0.003256155 35094  1.004461  0.3152
    ## temptrend_abs.sc:seas.sc                         -0.00068724 0.001816124 35094 -0.378409  0.7051
    ## temptrend_abs.sc:microclim.sc                    -0.00053626 0.001063892 35094 -0.504051  0.6142
    ## temptrend_abs.sc:mass.sc                          0.00094160 0.001537309 35094  0.612498  0.5402
    ## temptrend_abs.sc:speed.sc                         0.00055379 0.002051901 35094  0.269890  0.7872
    ## temptrend_abs.sc:consumerfrac.sc                  0.00540880 0.000862301 35094  6.272519  0.0000
    ## temptrend_abs.sc:nspp.sc                          0.00151827 0.001614751 35094  0.940250  0.3471
    ## tsign-1:thermal_bias.sc                          -0.00179711 0.000779247 35094 -2.306215  0.0211
    ## tsign1:thermal_bias.sc                           -0.00071119 0.000595760 35094 -1.193761  0.2326
    ## temptrend_abs.sc:npp.sc                           0.00730939 0.001550168 35094  4.715222  0.0000
    ## temptrend_abs.sc:veg.sc                          -0.00536541 0.002235699 35094 -2.399880  0.0164
    ## temptrend_abs.sc:duration.sc                      0.00073139 0.001122652 35094  0.651481  0.5147
    ## human_bowler.sc:REALM2TerrFresh                   0.00220308 0.001035351 35094  2.127859  0.0334
    ## human_bowler.sc:REALM2Marine                      0.00030223 0.000492898 35094  0.613173  0.5398
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.00017833 0.001469446 35094 -0.121361  0.9034
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00278199 0.001284046 35094  2.166583  0.0303
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.00369256 0.001531019 35094 -2.411828  0.0159
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.00306650 0.001238429 35094 -2.476119  0.0133
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.501                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.954  0.479                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.737  0.361  0.685                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.065  0.046 -0.001  0.011                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.101 -0.030 -0.080 -0.245  0.030                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.109  0.050  0.157 -0.095 -0.088  0.103                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.067  0.042  0.076  0.017  0.014 -0.185  0.101                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.067 -0.007 -0.044 -0.030 -0.007  0.089  0.096  0.007                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.085 -0.062 -0.039 -0.090 -0.069 -0.078 -0.011  0.059 -0.421                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.017  0.024  0.010  0.088  0.033 -0.120 -0.052  0.016 -0.014 -0.086                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.008 -0.037 -0.051 -0.053  0.052 -0.146  0.006 -0.090 -0.060  0.179  0.094                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.026  0.001 -0.038  0.058  0.012  0.024 -0.189 -0.215 -0.036  0.120 -0.037 -0.198                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.518  0.276  0.536  0.028 -0.012  0.007  0.086  0.005  0.013 -0.005 -0.003  0.021 -0.163                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.021  0.012  0.029  0.021 -0.149  0.077 -0.030 -0.022 -0.027  0.011  0.008 -0.241  0.053 -0.006                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.485 -0.949 -0.497 -0.337 -0.014  0.029 -0.069 -0.049 -0.004  0.035 -0.018  0.063  0.018 -0.292 -0.016                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.322 -0.727 -0.296 -0.510 -0.003  0.161  0.123  0.019 -0.017  0.082 -0.062  0.062 -0.057  0.030 -0.014  0.674                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.041 -0.112 -0.010  0.003 -0.490 -0.038  0.007 -0.009  0.005  0.031 -0.012 -0.021 -0.035 -0.008  0.033  0.035      0.006                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.008  0.041 -0.001  0.142 -0.063 -0.500 -0.094  0.047 -0.052  0.118  0.079  0.088 -0.044 -0.013 -0.016 -0.027     -0.273      0.036                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.059 -0.087 -0.081  0.111  0.052 -0.064 -0.641 -0.093 -0.060  0.013  0.019  0.025  0.125 -0.089 -0.016  0.123     -0.209     -0.027  0.073                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.029 -0.045 -0.034  0.002 -0.024  0.103 -0.051 -0.705  0.016 -0.040 -0.007  0.051  0.172  0.014  0.011  0.057     -0.060      0.005  0.021   0.115                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.010  0.023  0.005 -0.014  0.012 -0.062 -0.054  0.020 -0.580  0.245  0.026  0.011 -0.011 -0.009  0.047  0.012      0.032     -0.003  0.103   0.050            -0.051                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.047  0.100  0.033  0.075  0.048  0.099  0.041 -0.036  0.222 -0.581  0.071 -0.127 -0.042 -0.015 -0.008 -0.043     -0.109     -0.022 -0.285   0.012             0.040            -0.394                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.013 -0.035 -0.012 -0.064 -0.016  0.087  0.031  0.014  0.019  0.074 -0.719 -0.063  0.013 -0.001 -0.021  0.023      0.094      0.032 -0.099  -0.060            -0.037            -0.060            -0.116                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.025  0.029  0.040  0.054 -0.019  0.077  0.011  0.047  0.009 -0.134 -0.044 -0.615  0.100 -0.021  0.150 -0.081     -0.080      0.012 -0.121  -0.028            -0.020             0.007             0.159             0.088                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.060 -0.021 -0.056 -0.059 -0.088  0.204 -0.213 -0.040 -0.045  0.004 -0.007 -0.044  0.008 -0.005  0.073  0.027      0.027     -0.029 -0.061   0.132             0.022             0.020             0.003             0.002            0.014                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.103 -0.035 -0.105 -0.109  0.070  0.342 -0.321 -0.140 -0.016 -0.022 -0.017 -0.065 -0.057 -0.007  0.085  0.028      0.053      0.009 -0.115   0.148             0.048             0.019             0.011             0.014            0.022             0.393                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.010  0.023 -0.006 -0.055 -0.005 -0.012  0.086  0.158 -0.004 -0.058  0.005  0.098 -0.652  0.098 -0.036 -0.052      0.088      0.004  0.077  -0.208            -0.321             0.052             0.013             0.000           -0.100            -0.016  0.045                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.306 -0.426 -0.317  0.025 -0.011 -0.024 -0.102  0.009 -0.013  0.000  0.003 -0.021  0.158 -0.620  0.002  0.456     -0.060      0.054  0.006   0.173             0.006             0.007             0.023            -0.006            0.012             0.015 -0.016 -0.242                                                                                                  
    ## temptrend_abs.sc:duration.sc                      0.000  0.040 -0.036 -0.046  0.023  0.076 -0.039 -0.012  0.052 -0.023  0.004  0.220 -0.031  0.008 -0.448  0.035      0.050     -0.087 -0.003   0.003             0.015            -0.073            -0.008             0.041           -0.332            -0.021 -0.068  0.045            -0.001                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.069  0.066  0.071 -0.031  0.013  0.029  0.034  0.071  0.030 -0.029 -0.028 -0.005 -0.101  0.145 -0.014 -0.074      0.046     -0.013 -0.013  -0.097            -0.003            -0.001             0.007             0.016            0.023             0.041  0.083  0.046            -0.203            0.013                                                               
    ## human_bowler.sc:REALM2Marine                      0.042 -0.026 -0.048 -0.032 -0.004  0.102 -0.191 -0.005 -0.030  0.043 -0.062 -0.040 -0.187  0.022  0.014  0.032      0.016      0.003 -0.048   0.146            -0.015             0.033            -0.049             0.079            0.038             0.133  0.162  0.121            -0.014            0.036            0.019                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.022  0.025  0.024  0.037 -0.065 -0.109  0.161  0.042  0.008  0.032  0.001  0.020 -0.034  0.008 -0.056 -0.044     -0.062      0.252  0.264  -0.241            -0.086             0.013            -0.106             0.025           -0.026            -0.464 -0.268  0.070            -0.027            0.062           -0.035      -0.085                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.025  0.044  0.023  0.057 -0.018 -0.166  0.155  0.054  0.006  0.057 -0.003  0.045  0.032 -0.011 -0.053 -0.041     -0.104     -0.006  0.406  -0.220            -0.144             0.010            -0.176             0.023           -0.076            -0.239 -0.590  0.005             0.015            0.060           -0.046      -0.117       0.479                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.126 -0.217 -0.135  0.045 -0.010 -0.051 -0.127 -0.042 -0.014  0.025  0.011  0.022  0.097 -0.272  0.015  0.234     -0.086      0.046  0.088   0.224             0.048            -0.001            -0.036            -0.016           -0.039            -0.035 -0.057 -0.146             0.507           -0.011           -0.541      -0.001       0.091  0.092               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.028  0.034  0.034  0.028 -0.006 -0.076  0.137  0.000  0.029 -0.045  0.136  0.031  0.125 -0.011 -0.001 -0.039     -0.025      0.003  0.090  -0.166            -0.037            -0.028             0.046            -0.156           -0.050            -0.088 -0.118 -0.166             0.027           -0.025           -0.011      -0.737       0.086  0.124  0.001        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.67542034 -0.23386259 -0.02454167  0.23164547  5.75602618 
    ## 
    ## Number of Observations: 35327
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    200                  35327

``` r
rsquared(modTfullJturem0)
```

    ##       Response   family     link method     Marginal  Conditional
    ## 1 Jtutrendrem0 gaussian identity   none 3.800694e-05 0.0002287881

``` r
rsquared(modTfullJbetarem0)
```

    ##         Response   family     link method     Marginal  Conditional
    ## 1 Jbetatrendrem0 gaussian identity   none 0.0001307538 0.0005037384

``` r
rsquared(modTfullHornrem0)
```

    ##        Response   family     link method    Marginal Conditional
    ## 1 Horntrendrem0 gaussian identity   none 5.11082e-05 0.000265083

#### Summary exp

``` r
summary(modTfullJtuexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   93829.02 94169.38 -46872.51
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.8595795 (Intr)
    ## temptrend_abs.sc 0.1287652 -0.75 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.244636  1.31512
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.1319939 
    ## Fixed effects: log(Jtuexp) ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.3864821 0.4175784 24264   0.925532  0.3547
    ## temptrend_abs.sc                                  0.3788824 0.2116710 24264   1.789959  0.0735
    ## REALMMarine                                      -0.3361403 0.4424627   174  -0.759703  0.4485
    ## REALMTerrestrial                                  1.1247181 0.3993170   174   2.816604  0.0054
    ## tsign1                                           -0.0800175 0.0336012 24264  -2.381390  0.0173
    ## tempave_metab.sc                                 -0.3937535 0.0546561 24264  -7.204204  0.0000
    ## seas.sc                                          -0.1574622 0.0353950 24264  -4.448718  0.0000
    ## microclim.sc                                     -0.1175621 0.0184129 24264  -6.384756  0.0000
    ## mass.sc                                           0.0004837 0.0272330 24264   0.017761  0.9858
    ## speed.sc                                          0.0572922 0.0356703 24264   1.606158  0.1083
    ## consumerfrac.sc                                   0.0538704 0.0214783 24264   2.508137  0.0121
    ## nspp.sc                                          -0.3347941 0.0271760 24264 -12.319496  0.0000
    ## npp.sc                                           -0.0226195 0.0256037 24264  -0.883448  0.3770
    ## veg.sc                                            0.0283953 0.0685632 24264   0.414148  0.6788
    ## duration.sc                                       0.4862212 0.0209236 24264  23.237962  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.3212715 0.2140331 24264  -1.501037  0.1334
    ## temptrend_abs.sc:REALMTerrestrial                -0.2928269 0.2228085 24264  -1.314254  0.1888
    ## temptrend_abs.sc:tsign1                           0.0434855 0.0359460 24264   1.209743  0.2264
    ## temptrend_abs.sc:tempave_metab.sc                 0.0447311 0.0395536 24264   1.130899  0.2581
    ## temptrend_abs.sc:seas.sc                          0.0213296 0.0287205 24264   0.742663  0.4577
    ## temptrend_abs.sc:microclim.sc                     0.0436269 0.0195268 24264   2.234209  0.0255
    ## temptrend_abs.sc:mass.sc                          0.0164231 0.0256233 24264   0.640944  0.5216
    ## temptrend_abs.sc:speed.sc                        -0.0047824 0.0331303 24264  -0.144353  0.8852
    ## temptrend_abs.sc:consumerfrac.sc                 -0.0521381 0.0226744 24264  -2.299429  0.0215
    ## temptrend_abs.sc:nspp.sc                          0.0896528 0.0254677 24264   3.520252  0.0004
    ## tsign-1:thermal_bias.sc                          -0.0950776 0.0317883 24264  -2.990965  0.0028
    ## tsign1:thermal_bias.sc                           -0.0378730 0.0255918 24264  -1.479891  0.1389
    ## temptrend_abs.sc:npp.sc                           0.0464496 0.0272289 24264   1.705896  0.0880
    ## temptrend_abs.sc:veg.sc                          -0.0920762 0.0393699 24264  -2.338743  0.0194
    ## temptrend_abs.sc:duration.sc                      0.0798909 0.0263741 24264   3.029149  0.0025
    ## human_bowler.sc:REALM2TerrFresh                   0.1237597 0.0583907 24264   2.119511  0.0341
    ## human_bowler.sc:REALM2Marine                      0.0365949 0.0222430 24264   1.645238  0.0999
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0532641 0.0242091 24264   2.200170  0.0278
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0289804 0.0225008 24264   1.287971  0.1978
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0536190 0.0258276 24264  -2.076037  0.0379
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0401650 0.0275629 24264  -1.457213  0.1451
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.454                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.959  0.433                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.746  0.313  0.696                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.056  0.055 -0.002  0.027                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.094 -0.018 -0.074 -0.189 -0.044                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.084  0.034  0.134 -0.095 -0.075  0.078                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.059  0.057  0.068 -0.006 -0.008 -0.101  0.084                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.062  0.006 -0.033 -0.035  0.010  0.137  0.063  0.006                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.079 -0.058 -0.028 -0.051 -0.055 -0.129 -0.004  0.067 -0.444                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.025  0.030  0.015  0.103  0.060 -0.125 -0.085  0.008 -0.049 -0.075                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.016 -0.059 -0.047 -0.072  0.052 -0.131  0.012 -0.024 -0.075  0.178  0.082                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.050  0.001 -0.069  0.068 -0.009 -0.002 -0.306 -0.301 -0.040  0.078 -0.008 -0.168                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.532  0.295  0.554  0.012 -0.015  0.001  0.107  0.031  0.005 -0.005  0.001  0.015 -0.173                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.001  0.019  0.028  0.049 -0.170  0.014 -0.063 -0.016 -0.043  0.025 -0.001 -0.282  0.040  0.000                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.456 -0.958 -0.467 -0.306 -0.017  0.040 -0.078 -0.062 -0.014  0.032 -0.019  0.072  0.022 -0.326 -0.037                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.268 -0.618 -0.246 -0.432 -0.016  0.107  0.168  0.018 -0.044  0.081 -0.046  0.148 -0.106  0.037 -0.054  0.590                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.029 -0.148 -0.001 -0.010 -0.545  0.013 -0.006 -0.019  0.008  0.039 -0.032 -0.007 -0.004  0.017  0.098  0.058      0.021                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.022  0.026  0.014  0.089 -0.050 -0.398 -0.091  0.052 -0.169  0.206  0.098  0.037 -0.031  0.001  0.024 -0.022     -0.254      0.051                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.043 -0.086 -0.064  0.103  0.034 -0.033 -0.520 -0.080 -0.028  0.009  0.040 -0.001  0.140 -0.094 -0.024  0.128     -0.390      0.020  0.124                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.032 -0.065 -0.036  0.002 -0.015  0.083 -0.061 -0.545 -0.002 -0.043 -0.013  0.006  0.122 -0.005 -0.033  0.065     -0.075      0.063 -0.044   0.156                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                          0.000 -0.036 -0.004 -0.017  0.005 -0.116 -0.038  0.004 -0.557  0.319  0.099  0.034  0.003 -0.013  0.075  0.051      0.114     -0.031  0.361   0.007            -0.040                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.026  0.100  0.020  0.048  0.048  0.136  0.026 -0.044  0.341 -0.520 -0.008 -0.138 -0.014 -0.005 -0.019 -0.087     -0.199     -0.078 -0.505  -0.001             0.076            -0.636                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.001 -0.023  0.001 -0.025 -0.025  0.078  0.051 -0.009  0.088 -0.020 -0.573 -0.058 -0.010  0.009 -0.001  0.026      0.045      0.046 -0.187  -0.048             0.001            -0.228             0.055                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.029  0.151  0.026  0.082  0.001  0.020 -0.013  0.007 -0.001 -0.130 -0.039 -0.497  0.084 -0.023  0.070 -0.140     -0.311     -0.054  0.030  -0.016             0.035             0.033             0.251             0.092                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.034 -0.012 -0.033 -0.055 -0.039  0.285 -0.125 -0.055  0.000 -0.066 -0.010 -0.043 -0.020  0.007  0.043  0.022      0.021     -0.036 -0.043   0.115             0.024             0.003             0.031            -0.002            0.032                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.074 -0.064 -0.074 -0.093  0.020  0.439 -0.188 -0.180  0.029 -0.105 -0.023 -0.031 -0.054 -0.008  0.044  0.051      0.038      0.085 -0.086   0.142             0.084             0.002             0.042             0.006            0.000             0.374                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.003 -0.026  0.001 -0.064  0.009 -0.011  0.110  0.116  0.002 -0.015 -0.011  0.070 -0.524  0.093 -0.016 -0.006      0.195     -0.034  0.110  -0.192            -0.213             0.060            -0.042             0.033           -0.061             0.003  0.018                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.267 -0.509 -0.278  0.018 -0.005 -0.010 -0.125  0.001 -0.003 -0.001  0.010 -0.027  0.178 -0.529 -0.011  0.566     -0.132      0.048 -0.004   0.278            -0.029            -0.001             0.028            -0.026            0.019             0.013  0.012 -0.328                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.021  0.189 -0.014 -0.028 -0.012  0.115 -0.061  0.007  0.062 -0.046  0.023  0.104  0.011  0.005 -0.339  0.014      0.043     -0.035 -0.052   0.048             0.007            -0.123             0.003             0.040           -0.094             0.002 -0.058 -0.026             0.013                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.063  0.032  0.067 -0.050  0.009  0.023  0.026  0.074  0.024 -0.018 -0.027 -0.004 -0.088  0.136 -0.006 -0.039      0.055     -0.007 -0.013  -0.068            -0.019            -0.005             0.017             0.013            0.015             0.025  0.044  0.053            -0.123            0.000                                                               
    ## human_bowler.sc:REALM2Marine                      0.031 -0.042 -0.035 -0.012 -0.001  0.069 -0.198 -0.084 -0.026  0.013 -0.014 -0.032 -0.138  0.010 -0.014  0.053      0.028      0.024 -0.075   0.182             0.034             0.046            -0.045             0.047            0.016             0.131  0.149  0.124            -0.004            0.035            0.005                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.004  0.008  0.008  0.019 -0.061 -0.068  0.117  0.010 -0.006  0.061 -0.007  0.022 -0.005  0.000 -0.037 -0.032     -0.030      0.246  0.088  -0.232            -0.034             0.037            -0.126             0.014           -0.040            -0.559 -0.249  0.022            -0.045            0.026           -0.003      -0.137                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.039  0.082  0.037  0.041  0.074 -0.116  0.136  0.065 -0.002  0.058  0.004  0.003  0.015  0.022 -0.072 -0.062     -0.082     -0.222  0.147  -0.207            -0.108             0.015            -0.105            -0.003            0.017            -0.234 -0.578  0.015            -0.028            0.058           -0.032      -0.160       0.450                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.076 -0.244 -0.083  0.046 -0.009 -0.026 -0.110 -0.049 -0.006  0.018  0.014  0.007  0.105 -0.180  0.000  0.273     -0.158      0.048  0.048   0.306             0.045            -0.020            -0.025            -0.014           -0.043            -0.027 -0.032 -0.182             0.543            0.002           -0.550       0.010       0.049  0.054               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.028  0.077  0.029  0.018  0.020 -0.051  0.120  0.045  0.019 -0.032  0.058  0.003  0.106  0.002  0.005 -0.084     -0.052     -0.028  0.136  -0.284            -0.061            -0.065             0.078            -0.149           -0.002            -0.095 -0.106 -0.237             0.019           -0.020            0.005      -0.608       0.193  0.223 -0.025        
    ## 
    ## Standardized Within-Group Residuals:
    ##          Min           Q1          Med           Q3          Max 
    ## -2.885226292 -0.380500449 -0.001856242  0.466334961  2.171027700 
    ## 
    ## Number of Observations: 24474
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    177                  24474

``` r
summary(modTfullJbetaexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   83728.68 84068.01 -41822.34
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.7664634 (Intr)
    ## temptrend_abs.sc 0.2104152 -0.605
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.128957 1.976458
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.5576954 
    ## Fixed effects: log(Jbetaexp) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.0511188 0.3544972 23661   0.144201  0.8853
    ## temptrend_abs.sc                                  0.2423180 0.2159690 23661   1.122004  0.2619
    ## REALMMarine                                      -0.4545082 0.3765213   183  -1.207125  0.2289
    ## REALMTerrestrial                                  1.0638992 0.3423216   183   3.107894  0.0022
    ## tsign1                                           -0.1117519 0.0292229 23661  -3.824121  0.0001
    ## tempave_metab.sc                                 -0.3670294 0.0481278 23661  -7.626140  0.0000
    ## seas.sc                                          -0.1939765 0.0308779 23661  -6.282047  0.0000
    ## microclim.sc                                     -0.0749242 0.0160586 23661  -4.665681  0.0000
    ## mass.sc                                          -0.0175319 0.0238479 23661  -0.735153  0.4623
    ## speed.sc                                         -0.0076774 0.0310290 23661  -0.247426  0.8046
    ## consumerfrac.sc                                  -0.0189446 0.0182447 23661  -1.038363  0.2991
    ## nspp.sc                                          -0.4976847 0.0240866 23661 -20.662299  0.0000
    ## npp.sc                                           -0.0148788 0.0224377 23661  -0.663114  0.5073
    ## veg.sc                                            0.0130916 0.0568179 23661   0.230412  0.8178
    ## duration.sc                                       0.4566283 0.0190243 23661  24.002350  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.1706954 0.2206148 23661  -0.773726  0.4391
    ## temptrend_abs.sc:REALMTerrestrial                -0.5288397 0.2380377 23661  -2.221663  0.0263
    ## temptrend_abs.sc:tsign1                           0.0308175 0.0333701 23661   0.923506  0.3558
    ## temptrend_abs.sc:tempave_metab.sc                 0.1430030 0.0476167 23661   3.003215  0.0027
    ## temptrend_abs.sc:seas.sc                          0.0707029 0.0274229 23661   2.578246  0.0099
    ## temptrend_abs.sc:microclim.sc                     0.0163360 0.0188590 23661   0.866218  0.3864
    ## temptrend_abs.sc:mass.sc                         -0.0034863 0.0240136 23661  -0.145179  0.8846
    ## temptrend_abs.sc:speed.sc                         0.0074005 0.0322485 23661   0.229485  0.8185
    ## temptrend_abs.sc:consumerfrac.sc                  0.0308452 0.0215601 23661   1.430665  0.1525
    ## temptrend_abs.sc:nspp.sc                          0.1893196 0.0248692 23661   7.612622  0.0000
    ## tsign-1:thermal_bias.sc                          -0.0576764 0.0279701 23661  -2.062070  0.0392
    ## tsign1:thermal_bias.sc                            0.0123483 0.0226538 23661   0.545086  0.5857
    ## temptrend_abs.sc:npp.sc                           0.0570277 0.0270578 23661   2.107623  0.0351
    ## temptrend_abs.sc:veg.sc                          -0.0495616 0.0358338 23661  -1.383097  0.1666
    ## temptrend_abs.sc:duration.sc                      0.0939020 0.0239036 23661   3.928353  0.0001
    ## human_bowler.sc:REALM2TerrFresh                   0.1380299 0.0473741 23661   2.913612  0.0036
    ## human_bowler.sc:REALM2Marine                      0.0263521 0.0194029 23661   1.358151  0.1744
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0281297 0.0236133 23661   1.191268  0.2336
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0229682 0.0223267 23661   1.028735  0.3036
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0384541 0.0225661 23661  -1.704067  0.0884
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0454927 0.0258495 23661  -1.759906  0.0784
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.488                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.955  0.462                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.750  0.365  0.699                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.061  0.061  0.001  0.027                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.099 -0.035 -0.077 -0.194 -0.040                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.106  0.044  0.158 -0.079 -0.085  0.059                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.071  0.065  0.079 -0.001 -0.004 -0.120  0.086                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.060  0.005 -0.030 -0.031  0.013  0.120  0.058  0.002                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.080 -0.063 -0.027 -0.053 -0.061 -0.099  0.009  0.067 -0.436                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.019  0.027  0.010  0.095  0.051 -0.111 -0.081  0.010 -0.045 -0.073                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.001 -0.037 -0.038 -0.058  0.049 -0.133  0.019 -0.042 -0.079  0.168  0.078                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.066 -0.014 -0.083  0.057 -0.003  0.006 -0.283 -0.295 -0.035  0.081 -0.020 -0.177                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.523  0.261  0.545  0.024 -0.011  0.006  0.119  0.047  0.011 -0.006 -0.005  0.021 -0.198                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.002  0.016  0.026  0.045 -0.168  0.027 -0.060 -0.010 -0.032  0.020  0.005 -0.285  0.048 -0.005                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.482 -0.953 -0.496 -0.350 -0.027  0.046 -0.085 -0.071 -0.018  0.031 -0.016  0.060  0.039 -0.286 -0.030                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.311 -0.661 -0.286 -0.486 -0.019  0.146  0.130  0.009 -0.038  0.077 -0.056  0.099 -0.069  0.022 -0.038  0.624                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.032 -0.153 -0.004 -0.013 -0.539  0.008  0.002 -0.017 -0.008  0.044 -0.018 -0.011 -0.017  0.019  0.089  0.072      0.036                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.034  0.054  0.024  0.115 -0.044 -0.447 -0.061  0.063 -0.104  0.136  0.079  0.055 -0.066  0.010  0.012 -0.045     -0.325      0.050                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.051 -0.081 -0.072  0.097  0.040 -0.025 -0.536 -0.093 -0.024 -0.007  0.044 -0.002  0.133 -0.097 -0.024  0.120     -0.327      0.000  0.088                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.044 -0.081 -0.047  0.001 -0.010  0.077 -0.064 -0.559  0.009 -0.035 -0.023  0.021  0.145 -0.026 -0.037  0.085     -0.059      0.046 -0.044   0.170                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.002 -0.037 -0.004 -0.023 -0.009 -0.090 -0.026  0.012 -0.575  0.304  0.088  0.034  0.000 -0.014  0.058  0.062      0.101      0.003  0.215   0.002            -0.059                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.028  0.102  0.016  0.055  0.050  0.109  0.006 -0.041  0.305 -0.539  0.019 -0.119 -0.009 -0.011 -0.016 -0.070     -0.163     -0.080 -0.379   0.020             0.063            -0.547                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.001 -0.027 -0.002 -0.028 -0.013  0.064  0.049 -0.012  0.068  0.007 -0.623 -0.040  0.004  0.010 -0.018  0.025      0.078      0.024 -0.135  -0.065             0.013            -0.171            -0.026                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.016  0.098  0.020  0.065 -0.003  0.035 -0.021  0.028  0.017 -0.119 -0.033 -0.521  0.079 -0.030  0.102 -0.110     -0.212     -0.033 -0.027  -0.001            -0.003             0.004             0.224             0.070                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.037 -0.009 -0.036 -0.051 -0.036  0.283 -0.152 -0.046  0.000 -0.061 -0.006 -0.041 -0.023  0.010  0.049  0.018      0.017     -0.035 -0.066   0.126             0.013            -0.006             0.036            -0.003            0.028                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.078 -0.056 -0.079 -0.091  0.026  0.428 -0.224 -0.170  0.029 -0.095 -0.015 -0.046 -0.039 -0.003  0.051  0.043      0.047      0.077 -0.123   0.147             0.071            -0.011             0.051             0.005            0.024             0.387                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.007  0.001  0.011 -0.049 -0.002 -0.026  0.095  0.136 -0.004 -0.019  0.004  0.071 -0.546  0.110 -0.013 -0.038      0.123     -0.007  0.174  -0.177            -0.261             0.066            -0.056             0.009           -0.056             0.012  0.010                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.256 -0.452 -0.266  0.016 -0.012 -0.007 -0.129 -0.017 -0.014 -0.002  0.017 -0.032  0.208 -0.517 -0.014  0.502     -0.103      0.057 -0.031   0.272             0.007             0.016             0.041            -0.037            0.037             0.002  0.000 -0.365                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.023  0.158 -0.015 -0.026 -0.008  0.110 -0.062  0.002  0.047 -0.038  0.023  0.132  0.008  0.011 -0.358  0.015      0.045     -0.029 -0.041   0.040             0.016            -0.088            -0.014             0.035           -0.147             0.003 -0.057 -0.019             0.009                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.078  0.040  0.081 -0.025  0.013  0.018  0.020  0.083  0.023 -0.020 -0.029  0.000 -0.095  0.144 -0.010 -0.046      0.038     -0.006 -0.008  -0.065            -0.007            -0.006             0.015             0.016            0.011             0.029  0.052  0.044            -0.126            0.002                                                               
    ## human_bowler.sc:REALM2Marine                      0.031 -0.033 -0.035 -0.015 -0.001  0.074 -0.189 -0.068 -0.024  0.014 -0.004 -0.040 -0.133  0.014 -0.016  0.043      0.029      0.027 -0.072   0.154             0.031             0.047            -0.040             0.035            0.025             0.130  0.147  0.083            -0.005            0.041            0.009                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.005 -0.006  0.009  0.020 -0.070 -0.105  0.145  0.004 -0.021  0.071 -0.003  0.030 -0.001 -0.008 -0.040 -0.012     -0.029      0.251  0.163  -0.247            -0.009             0.072            -0.158             0.006           -0.064            -0.548 -0.283  0.018            -0.022            0.033           -0.011      -0.133                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.041  0.058  0.039  0.043  0.060 -0.155  0.149  0.058 -0.014  0.067  0.002  0.029  0.008  0.019 -0.064 -0.037     -0.086     -0.186  0.224  -0.196            -0.082             0.052            -0.149            -0.004           -0.054            -0.255 -0.589  0.035            -0.009            0.054           -0.038      -0.149       0.498                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.084 -0.223 -0.089  0.040 -0.010 -0.031 -0.109 -0.049 -0.006  0.019  0.019  0.003  0.110 -0.185  0.003  0.246     -0.137      0.046  0.050   0.294             0.036            -0.013            -0.026            -0.024           -0.025            -0.037 -0.045 -0.172             0.538           -0.007           -0.532       0.002       0.075  0.079               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.026  0.060  0.027  0.021  0.020 -0.059  0.107  0.042  0.029 -0.038  0.057  0.013  0.073  0.000  0.006 -0.068     -0.059     -0.035  0.139  -0.226            -0.067            -0.082             0.073            -0.131           -0.022            -0.095 -0.103 -0.143             0.015           -0.027            0.002      -0.619       0.185  0.203 -0.011        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.29626482 -0.30824477  0.03584154  0.38131979  1.91841260 
    ## 
    ## Number of Observations: 23880
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    186                  23880

``` r
summary(modTfullHornexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   88747.37 89084.75 -44331.69
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.9519617 (Intr)
    ## temptrend_abs.sc 0.1579174 -0.001
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.620502 170.8519
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -4.810054 
    ## Fixed effects: log(Hornexp) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.7338419 0.4451620 22583   1.648483  0.0993
    ## temptrend_abs.sc                                  0.4434973 0.2420115 22583   1.832546  0.0669
    ## REALMMarine                                      -0.3829760 0.4723490   175  -0.810790  0.4186
    ## REALMTerrestrial                                  1.6423131 0.4287598   175   3.830381  0.0002
    ## tsign1                                           -0.0465806 0.0357630 22583  -1.302480  0.1928
    ## tempave_metab.sc                                 -0.7074789 0.0579184 22583 -12.215103  0.0000
    ## seas.sc                                          -0.2067921 0.0379832 22583  -5.444303  0.0000
    ## microclim.sc                                     -0.1533234 0.0198651 22583  -7.718223  0.0000
    ## mass.sc                                          -0.0403574 0.0290866 22583  -1.387491  0.1653
    ## speed.sc                                          0.2591566 0.0382719 22583   6.771455  0.0000
    ## consumerfrac.sc                                  -0.0147138 0.0228754 22583  -0.643213  0.5201
    ## nspp.sc                                          -0.3567665 0.0295943 22583 -12.055257  0.0000
    ## npp.sc                                           -0.0709740 0.0275581 22583  -2.575427  0.0100
    ## veg.sc                                            0.0495502 0.0713164 22583   0.694795  0.4872
    ## duration.sc                                       0.4616017 0.0227367 22583  20.302026  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.3926933 0.2452954 22583  -1.600899  0.1094
    ## temptrend_abs.sc:REALMTerrestrial                -0.2285477 0.2664713 22583  -0.857682  0.3911
    ## temptrend_abs.sc:tsign1                          -0.0622542 0.0397755 22583  -1.565138  0.1176
    ## temptrend_abs.sc:tempave_metab.sc                 0.1612868 0.0509504 22583   3.165568  0.0015
    ## temptrend_abs.sc:seas.sc                         -0.0206426 0.0321187 22583  -0.642697  0.5204
    ## temptrend_abs.sc:microclim.sc                     0.0609690 0.0225585 22583   2.702701  0.0069
    ## temptrend_abs.sc:mass.sc                         -0.0031066 0.0289562 22583  -0.107286  0.9146
    ## temptrend_abs.sc:speed.sc                        -0.0442968 0.0377725 22583  -1.172728  0.2409
    ## temptrend_abs.sc:consumerfrac.sc                  0.0663566 0.0266048 22583   2.494163  0.0126
    ## temptrend_abs.sc:nspp.sc                          0.1849798 0.0289103 22583   6.398410  0.0000
    ## tsign-1:thermal_bias.sc                          -0.0469505 0.0340943 22583  -1.377077  0.1685
    ## tsign1:thermal_bias.sc                           -0.0584104 0.0280626 22583  -2.081433  0.0374
    ## temptrend_abs.sc:npp.sc                           0.0926844 0.0318092 22583   2.913759  0.0036
    ## temptrend_abs.sc:veg.sc                          -0.1733157 0.0421955 22583  -4.107441  0.0000
    ## temptrend_abs.sc:duration.sc                      0.0707925 0.0289942 22583   2.441608  0.0146
    ## human_bowler.sc:REALM2TerrFresh                   0.1902327 0.0601915 22583   3.160457  0.0016
    ## human_bowler.sc:REALM2Marine                      0.0080473 0.0239417 22583   0.336123  0.7368
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0116061 0.0277014 22583   0.418971  0.6752
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0476509 0.0258639 22583   1.842373  0.0654
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.1091558 0.0274908 22583  -3.970627  0.0001
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0721648 0.0310577 22583  -2.323572  0.0202
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.323                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.958  0.310                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.757  0.199  0.706                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.057  0.060 -0.001  0.027                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.089 -0.010 -0.069 -0.181 -0.048                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.082  0.026  0.132 -0.107 -0.078  0.058                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.066  0.064  0.075 -0.003 -0.006 -0.110  0.087                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.043  0.026 -0.018 -0.016  0.008  0.143  0.057  0.005                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.077 -0.054 -0.028 -0.051 -0.051 -0.136  0.005  0.070 -0.434                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.019  0.029  0.011  0.090  0.058 -0.127 -0.073  0.010 -0.053 -0.069                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.012 -0.061 -0.046 -0.070  0.048 -0.121  0.014 -0.034 -0.068  0.172  0.075                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.049 -0.002 -0.067  0.073 -0.018 -0.002 -0.299 -0.299 -0.038  0.080 -0.019 -0.168                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.523  0.260  0.544  0.023 -0.012  0.004  0.112  0.039  0.008 -0.006 -0.004  0.014 -0.182                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.001  0.027  0.026  0.050 -0.160  0.028 -0.064 -0.021 -0.039  0.013  0.002 -0.286  0.045 -0.001                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.325 -0.956 -0.329 -0.191 -0.024  0.030 -0.063 -0.072 -0.032  0.035 -0.020  0.073  0.024 -0.289 -0.036                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.153 -0.629 -0.137 -0.301 -0.012  0.105  0.179  0.015 -0.065  0.088 -0.043  0.137 -0.097  0.050 -0.063  0.596                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.034 -0.164 -0.004 -0.012 -0.542  0.015 -0.007 -0.020  0.003  0.037 -0.026 -0.012 -0.003  0.016  0.095  0.077      0.031                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.014  0.038  0.008  0.084 -0.058 -0.325 -0.080  0.051 -0.136  0.192  0.078  0.019 -0.052  0.005  0.016 -0.029     -0.304      0.056                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.049 -0.112 -0.071  0.109  0.033 -0.032 -0.503 -0.094 -0.015 -0.008  0.037  0.004  0.134 -0.120 -0.019  0.154     -0.379      0.016  0.111                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.038 -0.081 -0.044 -0.002 -0.012  0.082 -0.065 -0.552  0.000 -0.037 -0.016  0.019  0.130 -0.011 -0.020  0.085     -0.070      0.048 -0.035   0.180                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                          0.011 -0.042 -0.009 -0.037  0.002 -0.106 -0.026  0.003 -0.547  0.312  0.092  0.015  0.006 -0.009  0.071  0.064      0.117     -0.019  0.269  -0.008            -0.040                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.020  0.102  0.017  0.057  0.038  0.144  0.016 -0.039  0.328 -0.508  0.018 -0.121 -0.017 -0.012  0.001 -0.082     -0.185     -0.062 -0.436   0.024             0.063            -0.580                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.001 -0.032  0.001 -0.017 -0.016  0.063  0.049 -0.007  0.078 -0.001 -0.590 -0.037  0.000  0.008 -0.009  0.032      0.073      0.028 -0.134  -0.056            -0.004            -0.176            -0.039                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.023  0.102  0.022  0.076 -0.006  0.010 -0.014  0.024 -0.009 -0.118 -0.026 -0.484  0.076 -0.029  0.090 -0.107     -0.244     -0.032  0.026  -0.007             0.000             0.030             0.200             0.061                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.028 -0.003 -0.026 -0.051 -0.050  0.273 -0.135 -0.033  0.003 -0.073 -0.015 -0.055 -0.025  0.013  0.057  0.009      0.010     -0.028 -0.021   0.111             0.011            -0.003             0.042            -0.003            0.043                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.072 -0.049 -0.072 -0.086  0.025  0.425 -0.207 -0.170  0.040 -0.118 -0.027 -0.045 -0.040 -0.009  0.057  0.035      0.023      0.083 -0.056   0.138             0.079            -0.013             0.063             0.005            0.014             0.366                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.003  0.007  0.007 -0.056  0.011 -0.019  0.102  0.125  0.004 -0.017  0.005  0.066 -0.527  0.107 -0.023 -0.043      0.152     -0.021  0.152  -0.187            -0.238             0.049            -0.042             0.006           -0.063             0.003  0.006                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.249 -0.455 -0.260  0.027 -0.013 -0.011 -0.143 -0.002 -0.009  0.000  0.015 -0.032  0.202 -0.515 -0.013  0.512     -0.154      0.058 -0.014   0.339            -0.021             0.013             0.035            -0.033            0.049             0.009  0.010 -0.374                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.024  0.172 -0.010 -0.026 -0.004  0.103 -0.052  0.011  0.053 -0.035  0.016  0.121  0.001  0.006 -0.350  0.013      0.056     -0.035 -0.024   0.020            -0.013            -0.097            -0.025             0.037           -0.132            -0.012 -0.065  0.003             0.007                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.079  0.027  0.084 -0.034  0.004  0.025  0.037  0.079  0.016 -0.016 -0.025  0.006 -0.100  0.156 -0.009 -0.035      0.054      0.000 -0.007  -0.067            -0.019            -0.007             0.018             0.015            0.015             0.035  0.050  0.050            -0.130            0.002                                                               
    ## human_bowler.sc:REALM2Marine                      0.029 -0.039 -0.033 -0.012 -0.003  0.071 -0.195 -0.078 -0.032  0.010 -0.007 -0.036 -0.136  0.010 -0.016  0.049      0.029      0.029 -0.070   0.158             0.031             0.049            -0.047             0.035            0.026             0.131  0.148  0.105            -0.002            0.033            0.006                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.006 -0.004 -0.001  0.010 -0.063 -0.057  0.123 -0.003 -0.015  0.075 -0.007  0.036 -0.003 -0.006 -0.047 -0.012     -0.016      0.242  0.115  -0.234            -0.012             0.056            -0.143             0.015           -0.075            -0.533 -0.247  0.023            -0.039            0.055           -0.013      -0.126                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.037  0.066  0.035  0.034  0.065 -0.116  0.137  0.061 -0.015  0.079  0.001  0.018  0.007  0.025 -0.072 -0.044     -0.074     -0.196  0.180  -0.197            -0.091             0.049            -0.148             0.006           -0.032            -0.231 -0.571  0.034            -0.024            0.063           -0.034      -0.150       0.471                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.076 -0.220 -0.083  0.050 -0.007 -0.030 -0.120 -0.048 -0.006  0.019  0.016  0.005  0.116 -0.188 -0.001  0.250     -0.171      0.044  0.051   0.334             0.040            -0.006            -0.028            -0.021           -0.032            -0.036 -0.039 -0.192             0.553           -0.003           -0.539       0.007       0.064  0.067               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.027  0.071  0.029  0.021  0.022 -0.055  0.110  0.046  0.024 -0.035  0.059  0.007  0.090  0.003  0.008 -0.079     -0.058     -0.031  0.126  -0.230            -0.063            -0.072             0.087            -0.122           -0.016            -0.088 -0.098 -0.188             0.015           -0.022            0.002      -0.610       0.164  0.190 -0.018        
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.356184e+00 -1.355293e-02  2.250806e-06  2.591982e-02  1.551790e+00 
    ## 
    ## Number of Observations: 22794
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    178                  22794

``` r
rsquared(modTfullJtuexp)
```

    ##   Response   family     link method   Marginal Conditional
    ## 1   Jtuexp gaussian identity   none 0.07921369   0.6052232

``` r
rsquared(modTfullJbetaexp)
```

    ##   Response   family     link method   Marginal Conditional
    ## 1 Jbetaexp gaussian identity   none 0.05292824   0.3635837

``` r
rsquared(modTfullHornexp)
```

    ##   Response   family     link method     Marginal  Conditional
    ## 1  Hornexp gaussian identity   none 1.452585e-05 0.0001364047

#### Summary mm

``` r
summary(modTfullJtumm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   103449.5 103805.2 -51682.77
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.6683633 (Intr)
    ## temptrend_abs.sc 0.1715579 -0.967
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:    1.042733 0.3150569
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -6.801563 
    ## Fixed effects: log(Jtumm + 1) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                                       1.0479985 0.28393583 34973   3.690970  0.0002
    ## temptrend_abs.sc                                  0.0266529 0.12667449 34973   0.210404  0.8334
    ## REALMMarine                                      -0.1423867 0.30037066   197  -0.474037  0.6360
    ## REALMTerrestrial                                  1.1357198 0.28161155   197   4.032930  0.0001
    ## tsign1                                           -0.0014113 0.01784461 34973  -0.079091  0.9370
    ## tempave_metab.sc                                 -0.3762870 0.03152904 34973 -11.934617  0.0000
    ## seas.sc                                          -0.1351314 0.01960131 34973  -6.894000  0.0000
    ## microclim.sc                                     -0.0343521 0.00997048 34973  -3.445384  0.0006
    ## mass.sc                                           0.0218029 0.01530315 34973   1.424731  0.1542
    ## speed.sc                                         -0.0058629 0.02001602 34973  -0.292911  0.7696
    ## consumerfrac.sc                                   0.0119982 0.01227613 34973   0.977357  0.3284
    ## nspp.sc                                          -0.3242936 0.01528028 34973 -21.223024  0.0000
    ## npp.sc                                            0.0140982 0.01390437 34973   1.013939  0.3106
    ## veg.sc                                           -0.0144594 0.03883551 34973  -0.372325  0.7097
    ## duration.sc                                       0.2354274 0.01146070 34973  20.542157  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.0106653 0.12927451 34973  -0.082501  0.9342
    ## temptrend_abs.sc:REALMTerrestrial                -0.2726116 0.13522219 34973  -2.016027  0.0438
    ## temptrend_abs.sc:tsign1                           0.0049447 0.01944852 34973   0.254244  0.7993
    ## temptrend_abs.sc:tempave_metab.sc                 0.0572086 0.02229436 34973   2.566056  0.0103
    ## temptrend_abs.sc:seas.sc                          0.0325672 0.01561363 34973   2.085820  0.0370
    ## temptrend_abs.sc:microclim.sc                     0.0119697 0.01062005 34973   1.127084  0.2597
    ## temptrend_abs.sc:mass.sc                         -0.0146844 0.01411803 34973  -1.040119  0.2983
    ## temptrend_abs.sc:speed.sc                         0.0139177 0.01821733 34973   0.763980  0.4449
    ## temptrend_abs.sc:consumerfrac.sc                  0.0014536 0.01181450 34973   0.123039  0.9021
    ## temptrend_abs.sc:nspp.sc                          0.0731462 0.01387321 34973   5.272475  0.0000
    ## tsign-1:thermal_bias.sc                          -0.0566822 0.01724546 34973  -3.286789  0.0010
    ## tsign1:thermal_bias.sc                            0.0030399 0.01432648 34973   0.212191  0.8320
    ## temptrend_abs.sc:npp.sc                           0.0010783 0.01508436 34973   0.071482  0.9430
    ## temptrend_abs.sc:veg.sc                          -0.0093046 0.02095555 34973  -0.444015  0.6570
    ## temptrend_abs.sc:duration.sc                      0.0101012 0.01433202 34973   0.704800  0.4809
    ## human_bowler.sc:REALM2TerrFresh                   0.0478780 0.03308661 34973   1.447050  0.1479
    ## human_bowler.sc:REALM2Marine                     -0.0005355 0.01204478 34973  -0.044462  0.9645
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0186701 0.01337529 34973   1.395863  0.1628
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0071569 0.01253755 34973   0.570839  0.5681
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0037303 0.01360540 34973  -0.274181  0.7839
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0075594 0.01504766 34973  -0.502365  0.6154
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.652                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.955  0.618                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.801  0.526  0.752                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.043  0.053 -0.003  0.020                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.093 -0.045 -0.068 -0.149 -0.050                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.091  0.067  0.133 -0.058 -0.079  0.076                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.054  0.060  0.060  0.000 -0.011 -0.119  0.086                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.055 -0.011 -0.026 -0.033  0.016  0.135  0.042  0.003                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.072 -0.064 -0.024 -0.040 -0.058 -0.036  0.019  0.064 -0.418                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.016  0.024  0.012  0.080  0.049 -0.107 -0.063  0.016 -0.031 -0.062                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                          -0.008 -0.020 -0.024 -0.042  0.044 -0.144  0.023 -0.025 -0.097  0.156  0.078                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.053 -0.021 -0.067  0.048 -0.006 -0.020 -0.310 -0.298 -0.040  0.076 -0.014 -0.171                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.448  0.278  0.468  0.014 -0.015  0.004  0.123  0.041  0.010 -0.005  0.000  0.022 -0.178                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.006  0.007  0.020  0.039 -0.176  0.011 -0.073 -0.008 -0.021  0.018 -0.007 -0.293  0.052 -0.011                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.642 -0.954 -0.667 -0.511 -0.015  0.053 -0.112 -0.065 -0.007  0.029 -0.015  0.045  0.044 -0.306 -0.028                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.482 -0.673 -0.452 -0.649 -0.022  0.123  0.117  0.011 -0.019  0.066 -0.054  0.100 -0.085  0.028 -0.045  0.643                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.019 -0.134  0.002 -0.005 -0.546  0.017 -0.002 -0.016  0.005  0.035 -0.025 -0.012 -0.012  0.025  0.099  0.051      0.020                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.039  0.045  0.027  0.089 -0.026 -0.524 -0.078  0.069 -0.153  0.144  0.099  0.065 -0.019 -0.001  0.020 -0.035     -0.239      0.035                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.052 -0.087 -0.066  0.076  0.041 -0.025 -0.549 -0.089 -0.028 -0.002  0.029 -0.010  0.156 -0.111 -0.008  0.125     -0.343      0.009  0.104                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.035 -0.078 -0.034 -0.003 -0.010  0.084 -0.066 -0.561 -0.004 -0.034 -0.020  0.013  0.133 -0.021 -0.031  0.077     -0.054      0.053 -0.058   0.165                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.005 -0.029 -0.003 -0.013  0.000 -0.103 -0.033  0.005 -0.589  0.312  0.068  0.044  0.009 -0.016  0.061  0.050      0.096     -0.043  0.343   0.018            -0.031                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.026  0.081  0.013  0.042  0.049  0.090  0.009 -0.041  0.325 -0.559 -0.001 -0.118 -0.018 -0.008 -0.023 -0.068     -0.161     -0.061 -0.469   0.005             0.060            -0.631                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.003 -0.021 -0.007 -0.022 -0.024  0.073  0.036 -0.014  0.053 -0.003 -0.613 -0.055 -0.002  0.005  0.007  0.024      0.051      0.052 -0.187  -0.047             0.010            -0.152             0.038                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.009  0.107  0.011  0.057  0.002  0.037 -0.019  0.018  0.019 -0.113 -0.046 -0.536  0.089 -0.035  0.106 -0.104     -0.264     -0.036  0.013   0.001             0.016             0.018             0.238             0.109                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.038 -0.022 -0.033 -0.047 -0.064  0.329 -0.119 -0.051  0.008 -0.043 -0.017 -0.058 -0.025  0.006  0.039  0.028      0.028     -0.025 -0.117   0.115             0.026             0.000             0.026             0.005            0.036                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.069 -0.067 -0.066 -0.078  0.035  0.478 -0.179 -0.177  0.053 -0.077 -0.026 -0.056 -0.053 -0.007  0.035  0.056      0.049      0.073 -0.185   0.145             0.083            -0.016             0.045             0.017            0.022             0.399                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.005 -0.017  0.008 -0.048  0.004 -0.003  0.114  0.125  0.005 -0.018 -0.005  0.071 -0.545  0.103 -0.017 -0.014      0.183     -0.019  0.090  -0.199            -0.216             0.047            -0.033             0.029           -0.069             0.006  0.016                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.233 -0.449 -0.244  0.017 -0.002 -0.009 -0.143 -0.011 -0.014 -0.003  0.012 -0.035  0.202 -0.545 -0.006  0.497     -0.125      0.040 -0.001   0.287            -0.007             0.011             0.028            -0.030            0.047             0.011  0.011 -0.357                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.018  0.168 -0.011 -0.022 -0.008  0.110 -0.047 -0.004  0.046 -0.034  0.031  0.122  0.005  0.010 -0.353  0.012      0.046     -0.032 -0.048   0.031             0.011            -0.098            -0.012             0.014           -0.116             0.005 -0.044 -0.024             0.011                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.071  0.058  0.074 -0.027  0.009  0.022  0.040  0.050  0.010 -0.011 -0.020 -0.004 -0.077  0.142 -0.006 -0.064      0.024     -0.005 -0.013  -0.060             0.016             0.002             0.014             0.012            0.023             0.021  0.038  0.028            -0.114           -0.004                                                               
    ## human_bowler.sc:REALM2Marine                      0.029 -0.039 -0.032 -0.010 -0.005  0.068 -0.199 -0.073 -0.021  0.013 -0.012 -0.034 -0.134  0.006 -0.015  0.050      0.025      0.022 -0.071   0.176             0.031             0.040            -0.040             0.056            0.016             0.133  0.145  0.115            -0.004            0.039            0.000                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.012  0.007  0.013  0.019 -0.054 -0.123  0.121  0.012 -0.010  0.052  0.004  0.035 -0.001  0.001 -0.036 -0.027     -0.024      0.245  0.115  -0.231            -0.032             0.032            -0.127             0.002           -0.059            -0.576 -0.285  0.014            -0.040            0.035            0.001      -0.137                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.039  0.067  0.036  0.037  0.064 -0.187  0.132  0.072 -0.019  0.055  0.015  0.027  0.016  0.020 -0.061 -0.050     -0.073     -0.207  0.189  -0.198            -0.103             0.037            -0.121            -0.027           -0.012            -0.269 -0.612  0.018            -0.022            0.043           -0.022      -0.153       0.460                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.069 -0.201 -0.073  0.035 -0.004 -0.025 -0.112 -0.037  0.002  0.013  0.014  0.004  0.106 -0.172  0.002  0.223     -0.131      0.039  0.040   0.273             0.023            -0.024            -0.020            -0.019           -0.032            -0.030 -0.026 -0.175             0.491            0.004           -0.585       0.008       0.050  0.042               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.025  0.065  0.026  0.013  0.019 -0.053  0.115  0.042  0.018 -0.030  0.059  0.006  0.098  0.004  0.007 -0.073     -0.049     -0.019  0.130  -0.268            -0.055            -0.062             0.070            -0.174           -0.001            -0.098 -0.102 -0.234             0.027           -0.031            0.010      -0.614       0.198  0.213 -0.015        
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -2.070110e-01 -5.206507e-06 -4.626918e-08  1.262174e-06  1.837784e-01 
    ## 
    ## Number of Observations: 35206
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    200                  35206

``` r
summary(modTfullJbetamm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC     BIC   logLik
    ##   72603.81 72959.5 -36259.9
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.7327520 (Intr)
    ## temptrend_abs.sc 0.1775951 -0.993
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:   0.4841602 0.7468424
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.2818268 
    ## Fixed effects: log(Jbetamm + 1) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.8796929 0.23502864 35007   3.742918  0.0002
    ## temptrend_abs.sc                                 -0.0544482 0.09119739 35007  -0.597036  0.5505
    ## REALMMarine                                      -0.4080962 0.25112911   197  -1.625045  0.1058
    ## REALMTerrestrial                                  0.7013623 0.24381465   197   2.876621  0.0045
    ## tsign1                                           -0.0361529 0.01164193 35007  -3.105402  0.0019
    ## tempave_metab.sc                                 -0.3499171 0.02159595 35007 -16.202901  0.0000
    ## seas.sc                                          -0.1010895 0.01282390 35007  -7.882896  0.0000
    ## microclim.sc                                     -0.0276120 0.00636313 35007  -4.339380  0.0000
    ## mass.sc                                          -0.0003189 0.01020841 35007  -0.031239  0.9751
    ## speed.sc                                         -0.0303364 0.01339044 35007  -2.265531  0.0235
    ## consumerfrac.sc                                  -0.0130083 0.00778214 35007  -1.671553  0.0946
    ## nspp.sc                                          -0.2352533 0.01006494 35007 -23.373554  0.0000
    ## npp.sc                                            0.0029952 0.00896352 35007   0.334156  0.7383
    ## veg.sc                                            0.0171835 0.02434648 35007   0.705791  0.4803
    ## duration.sc                                       0.1549009 0.00761045 35007  20.353717  0.0000
    ## temptrend_abs.sc:REALMMarine                      0.0723725 0.09372611 35007   0.772170  0.4400
    ## temptrend_abs.sc:REALMTerrestrial                -0.1507651 0.09680841 35007  -1.557355  0.1194
    ## temptrend_abs.sc:tsign1                           0.0173617 0.01315812 35007   1.319467  0.1870
    ## temptrend_abs.sc:tempave_metab.sc                 0.0669633 0.01384398 35007   4.836996  0.0000
    ## temptrend_abs.sc:seas.sc                          0.0165804 0.01038688 35007   1.596284  0.1104
    ## temptrend_abs.sc:microclim.sc                     0.0014075 0.00699911 35007   0.201097  0.8406
    ## temptrend_abs.sc:mass.sc                         -0.0012058 0.00938953 35007  -0.128421  0.8978
    ## temptrend_abs.sc:speed.sc                         0.0180913 0.01217257 35007   1.486231  0.1372
    ## temptrend_abs.sc:consumerfrac.sc                  0.0091304 0.00759663 35007   1.201897  0.2294
    ## temptrend_abs.sc:nspp.sc                          0.0599296 0.00912440 35007   6.568065  0.0000
    ## tsign-1:thermal_bias.sc                          -0.0598087 0.01144814 35007  -5.224313  0.0000
    ## tsign1:thermal_bias.sc                           -0.0133489 0.00947768 35007  -1.408459  0.1590
    ## temptrend_abs.sc:npp.sc                           0.0071515 0.00992763 35007   0.720368  0.4713
    ## temptrend_abs.sc:veg.sc                          -0.0134573 0.01408441 35007  -0.955474  0.3393
    ## temptrend_abs.sc:duration.sc                     -0.0010980 0.00950378 35007  -0.115531  0.9080
    ## human_bowler.sc:REALM2TerrFresh                   0.0520693 0.02057148 35007   2.531141  0.0114
    ## human_bowler.sc:REALM2Marine                     -0.0065215 0.00771954 35007  -0.844808  0.3982
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0196074 0.00889081 35007   2.205361  0.0274
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0049637 0.00824548 35007   0.601990  0.5472
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0142871 0.00902761 35007  -1.582601  0.1135
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0035806 0.01002671 35007  -0.357106  0.7210
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.705                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.938  0.659                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.850  0.593  0.791                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.035  0.047 -0.003  0.014                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.092 -0.050 -0.063 -0.113 -0.052                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.084  0.064  0.117 -0.026 -0.085  0.067                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.042  0.057  0.045  0.000 -0.009 -0.130  0.093                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.065 -0.024 -0.035 -0.042  0.016  0.139  0.031  0.004                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.066 -0.062 -0.020 -0.031 -0.066  0.037  0.038  0.057 -0.408                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.006  0.018  0.007  0.055  0.048 -0.091 -0.052  0.016 -0.034 -0.027                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                          -0.024 -0.008 -0.008 -0.019  0.049 -0.159  0.029 -0.027 -0.121  0.149  0.090                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.048 -0.019 -0.057  0.030 -0.007 -0.021 -0.307 -0.285 -0.035  0.074 -0.018 -0.182                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.335  0.248  0.347  0.014 -0.017 -0.001  0.116  0.043  0.012 -0.007 -0.005  0.020 -0.178                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.008  0.007  0.014  0.024 -0.176  0.022 -0.073 -0.007 -0.004  0.013 -0.018 -0.310  0.052 -0.007                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.686 -0.951 -0.723 -0.574 -0.011  0.054 -0.105 -0.061  0.005  0.026 -0.011  0.033  0.040 -0.271 -0.025                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.577 -0.715 -0.538 -0.709 -0.021  0.098  0.098  0.013 -0.008  0.060 -0.042  0.091 -0.086  0.024 -0.039  0.684                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.017 -0.125 -0.001 -0.004 -0.536  0.018 -0.004 -0.020  0.010  0.034 -0.031 -0.017 -0.010  0.021  0.096  0.046      0.016                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.039  0.042  0.025  0.063 -0.027 -0.526 -0.079  0.074 -0.171  0.132  0.107  0.074 -0.002 -0.001  0.011 -0.033     -0.182      0.033                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.037 -0.069 -0.046  0.057  0.048 -0.009 -0.539 -0.095 -0.023 -0.005  0.024 -0.017  0.157 -0.102 -0.007  0.099     -0.324      0.010  0.091                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.027 -0.080 -0.025 -0.003 -0.016  0.097 -0.069 -0.554 -0.005 -0.029 -0.017  0.012  0.113 -0.020 -0.040  0.076     -0.053      0.060 -0.067   0.177                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.009 -0.020  0.001 -0.010  0.003 -0.101 -0.031  0.006 -0.575  0.308  0.071  0.052  0.001 -0.014  0.056  0.037      0.092     -0.058  0.394   0.014            -0.034                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.022  0.071  0.012  0.036  0.052  0.071  0.008 -0.039  0.323 -0.544 -0.027 -0.122 -0.015 -0.005 -0.018 -0.063     -0.155     -0.054 -0.515   0.003             0.061            -0.655                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.001 -0.015 -0.005 -0.009 -0.029  0.071  0.032 -0.011  0.052 -0.027 -0.610 -0.071 -0.008  0.008  0.021  0.019      0.035      0.074 -0.214  -0.049             0.003            -0.153             0.082                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.006  0.098  0.008  0.045 -0.004  0.040 -0.021  0.015  0.019 -0.112 -0.062 -0.526  0.101 -0.031  0.107 -0.094     -0.257     -0.026  0.021   0.001             0.027             0.029             0.240             0.151                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.040 -0.024 -0.032 -0.039 -0.071  0.347 -0.126 -0.058  0.009 -0.012 -0.014 -0.069 -0.022  0.003  0.048  0.028      0.027     -0.022 -0.125   0.115             0.036             0.007             0.014             0.005            0.041                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.071 -0.073 -0.064 -0.065  0.032  0.502 -0.204 -0.186  0.057 -0.032 -0.021 -0.069 -0.051 -0.014  0.049  0.061      0.046      0.079 -0.193   0.160             0.087            -0.007             0.023             0.017            0.016             0.415                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.000 -0.022  0.000 -0.042  0.010  0.003  0.114  0.105  0.001 -0.016 -0.014  0.078 -0.533  0.105 -0.017 -0.006      0.204     -0.032  0.043  -0.216            -0.183             0.056            -0.034             0.057           -0.088             0.003  0.022                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.168 -0.418 -0.174  0.012  0.001 -0.006 -0.130 -0.014 -0.009 -0.004  0.014 -0.036  0.198 -0.530 -0.005  0.457     -0.114      0.042  0.004   0.264             0.005            -0.002             0.029            -0.034            0.048             0.009  0.016 -0.366                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.015  0.156 -0.010 -0.016 -0.009  0.112 -0.052 -0.007  0.044 -0.032  0.039  0.127  0.011  0.006 -0.344  0.009      0.041     -0.040 -0.049   0.034             0.022            -0.104            -0.014             0.006           -0.121             0.004 -0.043 -0.029             0.010                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.067  0.064  0.067 -0.003  0.011  0.020  0.037  0.061  0.007 -0.008 -0.015  0.001 -0.082  0.135 -0.014 -0.068      0.015     -0.009 -0.015  -0.067             0.011             0.001             0.011             0.010            0.020             0.020  0.043  0.031            -0.116           -0.001                                                               
    ## human_bowler.sc:REALM2Marine                      0.026 -0.036 -0.028 -0.010 -0.005  0.075 -0.196 -0.073 -0.016  0.017 -0.016 -0.030 -0.145  0.009 -0.011  0.045      0.023      0.020 -0.067   0.181             0.031             0.039            -0.044             0.074            0.012             0.133  0.148  0.132            -0.013            0.041            0.001                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.013  0.008  0.011  0.016 -0.050 -0.124  0.118  0.012 -0.003  0.031  0.003  0.039 -0.006  0.002 -0.044 -0.026     -0.026      0.246  0.095  -0.225            -0.040             0.003            -0.103             0.004           -0.056            -0.570 -0.281  0.012            -0.037            0.044           -0.001      -0.139                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.037  0.076  0.032  0.031  0.073 -0.187  0.141  0.072 -0.012  0.026  0.016  0.024  0.015  0.023 -0.074 -0.058     -0.073     -0.229  0.169  -0.212            -0.110             0.007            -0.086            -0.034            0.023            -0.273 -0.605  0.011            -0.029            0.050           -0.035      -0.159       0.443                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.056 -0.193 -0.057  0.020 -0.001 -0.020 -0.108 -0.052  0.003  0.010  0.013  0.000  0.112 -0.171  0.010  0.210     -0.118      0.041  0.036   0.271             0.054            -0.029            -0.015            -0.025           -0.028            -0.032 -0.026 -0.195             0.499           -0.004           -0.560       0.004       0.055  0.047               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.019  0.056  0.021  0.009  0.014 -0.052  0.118  0.041  0.013 -0.029  0.074  0.001  0.113 -0.001  0.007 -0.063     -0.039     -0.012  0.116  -0.282            -0.057            -0.052             0.071            -0.217            0.003            -0.098 -0.107 -0.271             0.044           -0.033            0.010      -0.610       0.207  0.229 -0.006        
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -3.6501993 -0.4314522 -0.2058641  0.2597143  4.5374871 
    ## 
    ## Number of Observations: 35240
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    200                  35240

``` r
summary(modTfullHornmm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   108865.1 109220.8 -54390.57
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.9597873 (Intr)
    ## temptrend_abs.sc 0.1570534 -0.944
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.072115 7.872225
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##    power 
    ## -2.28603 
    ## Fixed effects: log(Hornmm + 1) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF    t-value p-value
    ## (Intercept)                                       1.4576662 0.3441282 34968   4.235823  0.0000
    ## temptrend_abs.sc                                 -0.0571576 0.1364275 34968  -0.418960  0.6752
    ## REALMMarine                                      -0.5266352 0.3661929   197  -1.438136  0.1520
    ## REALMTerrestrial                                  1.3439205 0.3493411   197   3.847015  0.0002
    ## tsign1                                           -0.0234427 0.0194376 34968  -1.206047  0.2278
    ## tempave_metab.sc                                 -0.6770119 0.0347094 34968 -19.505154  0.0000
    ## seas.sc                                          -0.2176248 0.0214253 34968 -10.157361  0.0000
    ## microclim.sc                                     -0.0716429 0.0107327 34968  -6.675172  0.0000
    ## mass.sc                                          -0.0096868 0.0167216 34968  -0.579302  0.5624
    ## speed.sc                                          0.1475132 0.0219699 34968   6.714346  0.0000
    ## consumerfrac.sc                                  -0.0500810 0.0130479 34968  -3.838257  0.0001
    ## nspp.sc                                          -0.2111534 0.0167064 34968 -12.639077  0.0000
    ## npp.sc                                           -0.0460649 0.0150305 34968  -3.064767  0.0022
    ## veg.sc                                           -0.0118539 0.0418442 34968  -0.283286  0.7770
    ## duration.sc                                       0.1681249 0.0126872 34968  13.251578  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.0014287 0.1390716 34968  -0.010273  0.9918
    ## temptrend_abs.sc:REALMTerrestrial                -0.0607428 0.1465676 34968  -0.414436  0.6786
    ## temptrend_abs.sc:tsign1                          -0.0043388 0.0222418 34968  -0.195074  0.8453
    ## temptrend_abs.sc:tempave_metab.sc                 0.1128786 0.0248222 34968   4.547477  0.0000
    ## temptrend_abs.sc:seas.sc                          0.0013520 0.0179414 34968   0.075357  0.9399
    ## temptrend_abs.sc:microclim.sc                     0.0216506 0.0121711 34968   1.778849  0.0753
    ## temptrend_abs.sc:mass.sc                         -0.0050591 0.0160260 34968  -0.315683  0.7522
    ## temptrend_abs.sc:speed.sc                        -0.0142206 0.0209217 34968  -0.679707  0.4967
    ## temptrend_abs.sc:consumerfrac.sc                  0.0115986 0.0132202 34968   0.877338  0.3803
    ## temptrend_abs.sc:nspp.sc                          0.0604122 0.0159258 34968   3.793363  0.0001
    ## tsign-1:thermal_bias.sc                          -0.0513982 0.0189484 34968  -2.712537  0.0067
    ## tsign1:thermal_bias.sc                           -0.0304688 0.0157390 34968  -1.935878  0.0529
    ## temptrend_abs.sc:npp.sc                           0.0355402 0.0172430 34968   2.061137  0.0393
    ## temptrend_abs.sc:veg.sc                          -0.0636304 0.0241454 34968  -2.635300  0.0084
    ## temptrend_abs.sc:duration.sc                     -0.0253954 0.0159880 34968  -1.588408  0.1122
    ## human_bowler.sc:REALM2TerrFresh                   0.2070115 0.0356043 34968   5.814224  0.0000
    ## human_bowler.sc:REALM2Marine                     -0.0095233 0.0129679 34968  -0.734376  0.4627
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.0022829 0.0152864 34968  -0.149339  0.8813
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0075029 0.0141993 34968   0.528398  0.5972
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0736339 0.0154575 34968  -4.763646  0.0000
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0098880 0.0171110 34968  -0.577871  0.5634
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.541                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.946  0.511                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.827  0.422  0.772                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.039  0.054 -0.004  0.016                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.093 -0.036 -0.065 -0.130 -0.049                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.085  0.046  0.124 -0.045 -0.082  0.074                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.048  0.063  0.054  0.000 -0.010 -0.126  0.087                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.060 -0.002 -0.031 -0.034  0.016  0.140  0.035  0.003                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.069 -0.059 -0.023 -0.038 -0.061 -0.005  0.030  0.063 -0.419                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.010  0.024  0.007  0.064  0.050 -0.100 -0.060  0.016 -0.038 -0.037                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                          -0.015 -0.035 -0.017 -0.035  0.047 -0.149  0.028 -0.027 -0.108  0.158  0.088                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.048 -0.007 -0.060  0.041 -0.008 -0.025 -0.308 -0.287 -0.035  0.073 -0.016 -0.179                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.394  0.284  0.409  0.013 -0.015  0.001  0.115  0.040  0.010 -0.005 -0.003  0.019 -0.174                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.007  0.017  0.016  0.032 -0.176  0.020 -0.074 -0.011 -0.015  0.014 -0.015 -0.293  0.052 -0.007                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.535 -0.957 -0.557 -0.411 -0.016  0.048 -0.088 -0.068 -0.011  0.031 -0.015  0.055  0.030 -0.313 -0.030                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.389 -0.644 -0.361 -0.525 -0.018  0.106  0.148  0.011 -0.037  0.078 -0.041  0.126 -0.102  0.034 -0.053  0.614                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.019 -0.143  0.002 -0.005 -0.537  0.016 -0.004 -0.020  0.009  0.034 -0.029 -0.011 -0.009  0.020  0.096  0.054      0.017                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.027  0.042  0.017  0.066 -0.041 -0.422 -0.083  0.063 -0.160  0.169  0.099  0.047 -0.017  0.000  0.015 -0.033     -0.251      0.041                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.037 -0.071 -0.052  0.072  0.040 -0.021 -0.521 -0.083 -0.019 -0.007  0.027 -0.013  0.144 -0.101 -0.010  0.111     -0.387      0.012  0.113                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.029 -0.080 -0.032 -0.003 -0.012  0.095 -0.060 -0.552 -0.004 -0.030 -0.019  0.014  0.117 -0.015 -0.030  0.080     -0.059      0.056 -0.052   0.157                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.001 -0.034 -0.002 -0.019  0.004 -0.101 -0.028  0.005 -0.558  0.313  0.075  0.042  0.004 -0.013  0.062  0.052      0.113     -0.047  0.349   0.010            -0.033                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.018  0.082  0.013  0.042  0.045  0.099  0.010 -0.040  0.329 -0.523 -0.015 -0.129 -0.015 -0.008 -0.011 -0.073     -0.180     -0.054 -0.499   0.014             0.058            -0.631                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.002 -0.022  0.001 -0.010 -0.026  0.069  0.038 -0.012  0.060 -0.021 -0.591 -0.064 -0.007  0.008  0.013  0.024      0.047      0.055 -0.187  -0.048             0.009            -0.158             0.040                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.013  0.123  0.014  0.057 -0.001  0.026 -0.020  0.018  0.011 -0.119 -0.053 -0.505  0.089 -0.029  0.093 -0.119     -0.293     -0.037  0.029   0.004             0.018             0.020             0.242             0.119                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.039 -0.017 -0.032 -0.044 -0.071  0.334 -0.124 -0.049  0.007 -0.030 -0.015 -0.067 -0.026  0.003  0.048  0.023      0.022     -0.023 -0.074   0.114             0.026             0.003             0.027             0.004            0.041                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.070 -0.062 -0.065 -0.071  0.035  0.486 -0.192 -0.179  0.054 -0.058 -0.023 -0.062 -0.055 -0.011  0.047  0.051      0.038      0.075 -0.129   0.145             0.086            -0.010             0.044             0.014            0.016             0.407                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.001 -0.026  0.000 -0.046  0.008  0.001  0.106  0.109  0.001 -0.014 -0.008  0.072 -0.527  0.100 -0.018 -0.006      0.205     -0.025  0.086  -0.197            -0.207             0.050            -0.036             0.034           -0.074             0.007  0.020                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.198 -0.478 -0.206  0.013 -0.005 -0.009 -0.126 -0.006 -0.007 -0.004  0.012 -0.033  0.191 -0.532 -0.010  0.529     -0.131      0.046 -0.003   0.270            -0.015             0.004             0.032            -0.030            0.042             0.012  0.010 -0.348                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.016  0.166 -0.009 -0.020 -0.009  0.110 -0.049 -0.002  0.050 -0.030  0.036  0.120  0.007  0.006 -0.341  0.015      0.054     -0.035 -0.041   0.032             0.010            -0.103            -0.019             0.012           -0.125             0.000 -0.048 -0.017             0.012                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.069  0.044  0.071 -0.012  0.007  0.022  0.040  0.060  0.009 -0.007 -0.016  0.000 -0.083  0.136 -0.015 -0.049      0.031     -0.004 -0.011  -0.057            -0.002             0.002             0.010             0.009            0.017             0.025  0.041  0.037            -0.105            0.005                                                               
    ## human_bowler.sc:REALM2Marine                      0.027 -0.039 -0.030 -0.012 -0.007  0.073 -0.198 -0.070 -0.021  0.015 -0.020 -0.031 -0.146  0.008 -0.012  0.049      0.027      0.023 -0.068   0.176             0.025             0.043            -0.047             0.070            0.014             0.133  0.148  0.129            -0.010            0.038            0.004                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.009  0.006  0.009  0.016 -0.051 -0.097  0.117  0.005 -0.009  0.052  0.003  0.037 -0.001  0.002 -0.043 -0.025     -0.025      0.245  0.098  -0.225            -0.031             0.027            -0.127             0.004           -0.061            -0.554 -0.268  0.013            -0.041            0.047           -0.005      -0.136                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.034  0.069  0.030  0.031  0.064 -0.156  0.135  0.065 -0.017  0.052  0.013  0.021  0.017  0.022 -0.070 -0.048     -0.074     -0.204  0.169  -0.199            -0.105             0.028            -0.123            -0.023           -0.004            -0.259 -0.589  0.015            -0.020            0.055           -0.028      -0.155       0.458                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.061 -0.215 -0.065  0.027 -0.006 -0.025 -0.109 -0.043  0.001  0.010  0.011  0.002  0.107 -0.174  0.005  0.238     -0.146      0.044  0.044   0.281             0.035            -0.022            -0.018            -0.019           -0.027            -0.031 -0.032 -0.180             0.505           -0.003           -0.553       0.005       0.057  0.057               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.021  0.067  0.022  0.013  0.019 -0.051  0.117  0.037  0.018 -0.032  0.076  0.001  0.110  0.000  0.007 -0.075     -0.050     -0.018  0.122  -0.266            -0.051            -0.060             0.079            -0.194           -0.001            -0.097 -0.103 -0.248             0.034           -0.026            0.004      -0.610       0.192  0.210 -0.012        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -0.98055390 -0.11715029 -0.02237458  0.04219116  2.66049529 
    ## 
    ## Number of Observations: 35201
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    200                  35201

``` r
rsquared(modTfullJtumm)
```

    ##   Response   family     link method  Marginal Conditional
    ## 1    Jtumm gaussian identity   none 0.0768167   0.9448081

``` r
rsquared(modTfullJbetamm)
```

    ##   Response   family     link method   Marginal Conditional
    ## 1  Jbetamm gaussian identity   none 0.06571932   0.6163737

``` r
rsquared(modTfullHornmm)
```

    ##   Response   family     link method    Marginal Conditional
    ## 1   Hornmm gaussian identity   none 0.003396396  0.03596406

### Plots from the full models

#### Plot the coefficients

``` r
coefs1 <- summary(modTfullJturem0)$tTable
coefs2 <- summary(modTfullJbetarem0)$tTable
coefs3 <- summary(modTfullHornrem0)$tTable

varstoplot <- unique(c(rownames(coefs1), rownames(coefs2), rownames(coefs3)))
varstoplot <- varstoplot[which(!grepl('Intercept', varstoplot) | grepl(':', varstoplot))] # vars to plot

rows1_1 <- which(rownames(coefs1) %in% varstoplot) # rows in coefs
rows1_2 <- which(rownames(coefs2) %in% varstoplot)
rows1_3 <- which(rownames(coefs3) %in% varstoplot)
xlims <- range(c(coefs1[rows1_1,1] - coefs1[rows1_1,2], coefs1[rows1_1,1] + coefs1[rows1_1,2], 
                  coefs2[rows1_2,1] - coefs2[rows1_2,2], coefs2[rows1_2,1] + coefs2[rows1_2,2], 
                  coefs3[rows1_3,1] - coefs3[rows1_3,2], coefs3[rows1_3,1] + coefs3[rows1_3,2]))


cols <- brewer.pal(3, 'Dark2') # for Jtu, Jbeta and Horn models
pchs <- c(16, 16, 16)
offs <- c(0.1, 0, -0.1) # offset vertically for each model


par(las = 1, mai = c(0.5, 4, 0.1, 0.1))

plot(0,0, col = 'white', xlim = xlims, ylim = c(1,length(varstoplot)), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(varstoplot):1, labels = varstoplot, cex.axis = 0.7)
abline(v = 0, col = 'grey', lty = 2)
abline(h = 1:length(varstoplot), col = 'grey', lty = 3)
for(i in 1:length(varstoplot)){
  if(varstoplot[i] %in% rownames(coefs1)){
    x = coefs1[rownames(coefs1) == varstoplot[i], 1]
    se = coefs1[rownames(coefs1) == varstoplot[i], 2]
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
}
legend('topleft', col = cols, pch = 16, lwd = 1, legend = c('Jtu', 'Jbeta', 'Horn'), cex = 0.5)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20fullTmods-1.png)<!-- -->

#### Plot the main effects

##### Manually

``` r
# set up the variables to plot
# if variable is logged before scaling (see 'center and scale' above), then need to mark it here and express the limits on a log10 scale (even though log transforming is log)
vars <- data.frame(vars = c('temptrend_abs', 'tempave_metab', 'seas', 'microclim', 'mass', 'speed', 
                            'consumerfrac', 'nspp', 'thermal_bias', 'npp', 'veg', 'duration', 
                            'human_bowler', 'human_bowler'),
           min =      c(0,  0,   0.1, -2,  0,   0,   0,   0.3, -10, 1.9, 0,   0.5, 0,   0), 
           max =      c(2,  30,  16,  0.8, 8,   2,   1,   2.6, 10,  3.7, 0.3, 2,   1,   1),
           log =      c(F,  F,   F,   T,   T,   T,   F,   T,   F,   T,   T,   T,   T,   T),
           len =      c(100,  100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
           discrete = c(F,  F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F),
           plus =     c(0,  0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   1,   1), # what to add before log-scaling
           REALM = c(rep('Terrestrial', 13), 'Marine'),
           REALM2 = c(rep('TerrFresh', 13), 'Marine'),
           stringsAsFactors = FALSE)
baseall <- trends[, .(type = 'all', 
                      temptrend = -0.0001,
                      temptrend_abs.sc = 0.0001,
                      tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                      seas.sc = mean(seas.sc, na.rm=TRUE), 
                      microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                      speed.sc = mean(speed.sc, na.rm=TRUE), 
                      mass.sc = mean(mass.sc, na.rm=TRUE), 
                      nspp.sc = 0, 
                      thermal_bias.sc = mean(thermal_bias.sc, na.rm=TRUE), 
                      npp.sc = mean(npp.sc, na.rm=TRUE), 
                      human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                      veg.sc = mean(veg.sc, na.rm=TRUE), 
                      consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
baseterr <- trends[REALM == 'Terrestrial', 
                   .(type = 'Terrestrial', 
                     temptrend = -0.0001,
                     temptrend_abs.sc = 0.0001,
                     tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                     seas.sc = mean(seas.sc, na.rm=TRUE), 
                     microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                     speed.sc = mean(speed.sc, na.rm=TRUE), 
                     mass.sc = mean(mass.sc, na.rm=TRUE), 
                     nspp.sc = 0, 
                     thermal_bias.sc = 0, 
                     npp.sc = mean(npp.sc, na.rm=TRUE), 
                     human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                     veg.sc = mean(veg.sc, na.rm=TRUE), 
                     consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
basemar <- trends[REALM == 'Marine', 
                  .(type = 'Marine',
                    temptrend = -0.0001,
                    temptrend_abs.sc = 0.0001,
                    tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                    seas.sc = mean(seas.sc, na.rm=TRUE), 
                    microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                    speed.sc = mean(speed.sc, na.rm=TRUE), 
                    mass.sc = mean(mass.sc, na.rm=TRUE), 
                    nspp.sc = 0, 
                    thermal_bias.sc = 0, 
                    npp.sc = mean(npp.sc, na.rm=TRUE), 
                    human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                    veg.sc = mean(veg.sc, na.rm=TRUE), 
                    consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
basetab <- rbind(baseall, baseterr, basemar)
basetab[, ':='(duration.sc = 0, nyrBT = 20, STUDY_ID = 127L, rarefyID = '127_514668')]

# make the data frames for each interaction to plot                
for(j in 1:nrow(vars)){
  # set up the main effects
  if(vars$log[j]){
    thisdat <- data.frame(new = 10^seq(vars$min[j], vars$max[j], length.out = vars$len[j]),
                          var = vars$vars[j], stringsAsFactors = FALSE)
  } 
  if(!vars$log[j]){
    thisdat <- data.frame(new = seq(vars$min[j], vars$max[j], length.out = vars$len[j]),
                          var = vars$vars[j], stringsAsFactors = FALSE)
  }
  names(thisdat) <- c(vars$vars[j], 'var')

  # scale the variable
  cent <- attr(trends[[paste0(vars$vars[j], '.sc')]], 'scaled:center')
  scl <- attr(trends[[paste0(vars$vars[j], '.sc')]], 'scaled:scale')
  if(is.null(cent)) cent <- 0
  if(!is.null(cent) & !is.null(scl)){
    if(vars$log[j]) thisdat[[paste0(vars$var[j], '.sc')]] <- (log(thisdat[[vars$vars[j]]] + vars$plus[j]) - cent)/scl
    if(!vars$log[j]) thisdat[[paste0(vars$var[j], '.sc')]] <- (thisdat[[vars$var[j]]] - cent)/scl
  }

  # merge with the rest of the columns
  # use realm-specific averages for human impacts
  if(vars$vars[j] != 'tsign') colnamestouse <- setdiff(colnames(basetab), paste0(vars$vars[j], '.sc'))
  if(vars$vars[j] == 'tsign') colnamestouse <- setdiff(colnames(basetab), vars$var[j])
  if(vars$vars[j] != 'human_bowler'){
    thisdat <- cbind(thisdat, basetab[type == 'all', ..colnamestouse])
  }
  if(vars$vars[j] == 'human_bowler' & vars$REALM[j] == 'Terrestrial'){
    thisdat <- cbind(thisdat, basetab[type == 'Terrestrial', ..colnamestouse])
  }
  if(vars$vars[j] == 'human_bowler' & vars$REALM[j] == 'Marine'){
    thisdat <- cbind(thisdat, basetab[type == 'Marine', ..colnamestouse])
  }
  
  # add realm
  thisdat$REALM <- vars$REALM[j]
  thisdat$REALM2 <- vars$REALM2[j]
  
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
newdat <- rbind(newdat[1:6, ], newdat)
newdat$REALM[1:6] <- c('Marine', 'Marine', 'Freshwater', 'Freshwater', 'Terrestrial', 'Terrestrial')
newdat$REALM2[1:6] <- c('Marine', 'Marine', 'TerrFresh', 'TerrFresh', 'TerrFresh', 'TerrFresh')
newdat$temptrend_abs.sc[1:6] <- rep(0.0001, 6)
newdat$temptrend[1:6] <- rep(c(0.0001, -0.0001), 3)
newdat$var[1:6] <- 'test'

# trim to at least some temperature change (so that tsign is -1 or 1)
newdat <- newdat[newdat$temptrend_abs.sc != 0,]

# set up tsign
newdat$tsign <- factor(sign(newdat$temptrend))


# make predictions
newdat$predsJtu <- predict(object = modTfullJturem0, newdata = newdat, level = 0)
newdat$predsJbeta <- predict(object = modTfullJbetarem0, newdata = newdat, level = 0)
newdat$predsHorn <- predict(object = modTfullHornrem0, newdata = newdat, level = 0)

#compute standard error for predictions
# from https://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit
DesignmatJtu <- model.matrix(eval(eval(modTfullJturem0$call$fixed)[-2]), newdat)
DesignmatJbeta <- model.matrix(eval(eval(modTfullJbetarem0$call$fixed)[-2]), newdat)
DesignmatHorn <- model.matrix(eval(eval(modTfullHornrem0$call$fixed)[-2]), newdat)

predvarJtu <- diag(DesignmatJtu %*% modTfullJturem0$varFix %*% t(DesignmatJtu))
predvarJbeta <- diag(DesignmatJbeta %*% modTfullJbetarem0$varFix %*% t(DesignmatJbeta))
predvarHorn <- diag(DesignmatHorn %*% modTfullHornrem0$varFix %*% t(DesignmatHorn))

newdat$SE_Jtu <- sqrt(predvarJtu) 
newdat$SE_Jbeta <- sqrt(predvarJbeta) 
newdat$SE_Horn <- sqrt(predvarHorn) 

# prep the plots
varplots <- vector('list', nrow(vars))
for(j in 1:length(varplots)){
  subs <- newdat$var == vars$vars[j] # which rows of newdat
  xvar <- vars$vars[j]
  title <- vars$vars[j]
  if(vars$vars[j] %in% c('human_bowler')){
    subs <- newdat$var == vars$vars[j] & newdat$REALM2 == vars$REALM2[j]
    title <- paste0('human:', vars$REALM2[j])
  } 

  se <- 1
  thisplot <- ggplot(newdat[subs, ], 
                     aes_string(x = xvar, y = 'predsJtu')) +
    geom_line() +
    geom_ribbon(aes(ymin = predsJtu - se*SE_Jtu, ymax = predsJtu + se*SE_Jtu), alpha = 0.5, fill = "grey") +
    geom_line(aes(y = predsJbeta), color = 'red') +
    geom_ribbon(aes(ymin = predsJbeta - se*SE_Jbeta, ymax = predsJbeta + se*SE_Jbeta), alpha = 0.5, fill = "red") +
    geom_line(aes(y = predsHorn), color = 'blue') +
    geom_ribbon(aes(ymin = predsHorn - se*SE_Horn, ymax = predsHorn + se*SE_Horn), alpha = 0.5, fill = "blue") +
    #coord_cartesian(ylim = c(0, 0.4)) +
    theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm')) +
    labs(title = title) 
  varplots[[j]] <- thisplot
  if(vars$log[j] & !vars$discrete[j]){
    varplots[[j]] <- thisplot + scale_x_log10()
  }
}

grid.arrange(grobs = varplots, ncol = 3)
```

<img src="turnover_vs_temperature_MEmodels_files/figure-gfm/main effect plots-1.png" width="24" height="24" />

``` r
# write out the interactions
write.csv(newdat, file = 'output/maineffects.csv')
```

##### Using sjPlot

Doesn’t work now

``` r
require(sjPlot)
# p1 <- sjPlot::plot_model(modTfullJturem0, type = 'est', terms = c('temptrend_abs.sc'))
```

#### Plot interactions (Jaccard turnover)

``` r
# set up the interactions to plot
# if variable is logged before scaling (see 'center and scale' above), then need to mark it here and express the limits on a log10 scale (even though log transforming is log)
ints <- data.frame(vars = c('tsign', 'tempave_metab', 'seas', 'microclim', 'mass', 'speed', 
                            'consumerfrac', 'nspp', 'thermal_bias', 'npp', 'veg', 'duration', 
                            'human_bowler', 'human_bowler'),
           min =      c(1,  0,   0.1, -2,  0,   0,   0,   0.3, -10, 1.9, 0,   0.5, 0,   0), 
           max =      c(2,  30,  16,  0.8, 8,   2,   1,   2.6, 10,  3.7, 0.3, 2,   1,   1),
           log =      c(F,  F,   F,   T,   T,   T,   F,   T,   F,   T,   T,   T,   T,   T),
           len =      c(2,  100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
           discrete = c(T,  F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F),
           plus =     c(0,  0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   1,   1), # what to add before log-scaling
           REALM = c(rep('Terrestrial', 12), 'Terrestrial', 'Marine'),
           REALM2 = c(rep('TerrFresh', 13), 'Marine'),
           stringsAsFactors = FALSE)
baseall <- trends[, .(type = 'all', tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                      seas.sc = mean(seas.sc, na.rm=TRUE), 
                      microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                      speed.sc = mean(speed.sc, na.rm=TRUE), 
                      mass.sc = mean(mass.sc, na.rm=TRUE), 
                      nspp.sc = 0, 
                      thermal_bias.sc = mean(thermal_bias.sc, na.rm=TRUE), 
                      npp.sc = mean(npp.sc, na.rm=TRUE), 
                      human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                      veg.sc = mean(veg.sc, na.rm=TRUE), 
                      consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
baseterr <- trends[REALM == 'Terrestrial', 
                   .(type = 'Terrestrial', 
                     tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                     seas.sc = mean(seas.sc, na.rm=TRUE), 
                     microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                     speed.sc = mean(speed.sc, na.rm=TRUE), 
                     mass.sc = mean(mass.sc, na.rm=TRUE), 
                     nspp.sc = 0, 
                     thermal_bias.sc = 0, 
                     npp.sc = mean(npp.sc, na.rm=TRUE), 
                     human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                     veg.sc = mean(veg.sc, na.rm=TRUE), 
                     consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
basemar <- trends[REALM == 'Marine', 
                  .(type = 'Marine',
                    tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                    seas.sc = mean(seas.sc, na.rm=TRUE), 
                    microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                    speed.sc = mean(speed.sc, na.rm=TRUE), 
                    mass.sc = mean(mass.sc, na.rm=TRUE), 
                    nspp.sc = 0, 
                    thermal_bias.sc = 0, 
                    npp.sc = mean(npp.sc, na.rm=TRUE), 
                    human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                    veg.sc = mean(veg.sc, na.rm=TRUE), 
                    consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
basetab <- rbind(baseall, baseterr, basemar)
basetab[, ':='(duration.sc = 0, nyrBT = 20, STUDY_ID = 127L, rarefyID = '127_514668')]

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
    if(ints$log[j]) thisdat[[paste0(ints$var[j], '.sc')]] <- (log(thisdat[[ints$vars[j]]] + ints$plus[j]) - cent)/scl
    if(!ints$log[j]) thisdat[[paste0(ints$var[j], '.sc')]] <- (thisdat[[ints$var[j]]] - cent)/scl
  }

  # merge with the rest of the columns
  # use realm-specific averages for human impacts
  if(ints$vars[j] != 'tsign') colnamestouse <- setdiff(colnames(basetab), paste0(ints$var[j], '.sc'))
  if(ints$vars[j] == 'tsign') colnamestouse <- setdiff(colnames(basetab), ints$var[j])
  if(ints$vars[j] != 'human_bowler'){
    thisdat <- cbind(thisdat, basetab[type == 'all', ..colnamestouse])
  }
  if(ints$vars[j] == 'human_bowler' & ints$REALM[j] == 'Terrestrial'){
    thisdat <- cbind(thisdat, basetab[type == 'Terrestrial', ..colnamestouse])
  }
  if(ints$vars[j] == 'human_bowler' & ints$REALM[j] == 'Marine'){
    thisdat <- cbind(thisdat, basetab[type == 'Marine', ..colnamestouse])
  }
  
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
newdat <- rbind(newdat[1:6, ], newdat)
newdat$REALM[1:6] <- c('Marine', 'Marine', 'Freshwater', 'Freshwater', 'Terrestrial', 'Terrestrial')
newdat$REALM2[1:6] <- c('Marine', 'Marine', 'TerrFresh', 'TerrFresh', 'TerrFresh', 'TerrFresh')
newdat$temptrend[1:6] <- c(-1, 1, -1, 1, -1, 1)

# trim to at least some temperature change (so that tsign is -1 or 1)
newdat <- newdat[newdat$temptrend != 0,]

# scale the temperature vars
newdat$temptrend.sc <- newdat$temptrend/attr(trends$temptrend.sc, 'scaled:scale') 
newdat$temptrend_abs <- abs(newdat$temptrend)
newdat$temptrend_abs.sc <- (newdat$temptrend_abs)/attr(trends$temptrend_abs.sc, 'scaled:scale')
newdat$tsign <- factor(sign(newdat$temptrend))

# make predictions
newdat$preds <- predict(object = modTfullJturem0, newdata = newdat, level = 0)

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
    coord_cartesian(ylim = c(-0.2, 0.4)) +
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

![](turnover_vs_temperature_MEmodels_files/figure-gfm/interaction%20plots%20modTfullJturem0-1.png)<!-- -->

``` r
# write out the interactions
write.csv(newdat, file = 'temp/interactions.csv')
```

#### Plot residuals against each predictor (Jaccard turnover)

``` r
resids <- resid(modTfullJturem0)
preds <- getData(modTfullJturem0)
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

``` r
AICnas <- function(x){
  if(class(x) == 'NULL'){
    return(NA)
  } else {
    return(AIC(x))
  }
}

if(file.exists('output/aics_from_full.csv')){
  aicsfromfull <- read.csv('output/aics_from_full.csv')
  
  if('dAIC_Jtu' %in% colnames(aicsfromfull)){
    runJtu <- FALSE
  } else {
    runJtu <- TRUE
  }
  
  if('dAIC_Jbeta' %in% colnames(aicsfromfull)){
    runJbeta <- FALSE
  } else {
    runJbeta <- TRUE
  }
  
  if('dAIC_Horn' %in% colnames(aicsfromfull)){
    runHorn <- FALSE
  } else {
    runHorn <- TRUE
  }
  
} else {
  runJtu <- TRUE
  runJbeta <- TRUE
  runHorn <- TRUE
}

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

terms <- c('temptrend_abs.sc*REALM', 
           'temptrend_abs.sc*tsign',
           'temptrend_abs.sc*tempave_metab.sc',
           'temptrend_abs.sc*seas.sc',
           'temptrend_abs.sc*microclim.sc',
           'temptrend_abs.sc*mass.sc',
           'temptrend_abs.sc*speed.sc', 
           'temptrend_abs.sc*consumerfrac.sc',
           'temptrend_abs.sc*nspp.sc',
           'temptrend_abs.sc*thermal_bias.sc:tsign',
           'temptrend_abs.sc*npp.sc',
           'temptrend_abs.sc*veg.sc',
           'temptrend_abs.sc*duration.sc',
           'temptrend_abs.sc*human_bowler.sc:REALM2')


if(runJtu){
  i <- trends[, complete.cases(Jtutrendrem0, REALM, tempave_metab.sc, seas.sc, microclim.sc, 
                               temptrend_abs.sc, mass.sc, speed.sc, 
                               consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                               veg.sc, duration.sc, human_bowler.sc)]
  
  modTdrops <- vector('list', length(terms)+2)
  names(modTdrops) <- c('full', '-temptrend_abs.sc', paste0('-', terms))
  
  # fit full model with ML for model comparison
  modTdrops[[1]] <- lme(formula(paste0('Jtutrendrem0 ~ ', paste(terms, collapse = ' + '))),
                        random = randef, weights = varef, data = trends[i,], method = 'ML')
  
  # w/out temptrend
  modTdrops[[2]] <- lme(formula(paste0('Jtutrendrem0 ~ ', paste(gsub('temptrend_abs.sc\\*', '', terms), collapse = ' + '))),
                        random = list(STUDY_ID = ~ 1, rarefyID = ~1), weights = varef, data = trends[i,], method = 'ML')
  
  for(j in 1:length(terms)){
    print(j)
    tryCatch({
      modTdrops[[j+2]] <- lme(formula(paste0('Jtutrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                              random = randef, weights = varef, data = trends[i,], method = 'ML')
      
    }, error = function(e){
      print('going to optim (Jtu)')
      tryCatch({
        modTdrops[[j+2]] <- lme(formula(paste0('Jtutrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                                random = randef, weights = varef, data = trends[i,], method = 'ML',
                                control = lmeControl(opt = 'optim'))
        
      }, error = function(e){
        print('going to more iters (Jtu)') 
        tryCatch({
          modTdrops[[j+2]] <- lme(formula(paste0('Jtutrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                                  random = randef, weights = varef, data = trends[i,], method = 'ML',
                                  control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
          
        }, error= function(e){
          print('giving up on this one')
          modTdrops[[j+2]] <- NA
        })
      })
    })
  }
  
  aicsJtu <- sapply(modTdrops, AICnas)
}


if(runJbeta){
  i <- trends[, complete.cases(Jbetatrendrem0, REALM, tempave_metab.sc, seas.sc, microclim.sc, 
                               temptrend_abs.sc, mass.sc, speed.sc, 
                               consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                               veg.sc, duration.sc, human_bowler.sc)]
  
  modTJbetadrops <- vector('list', length(terms)+2)
  names(modTJbetadrops) <- c('full', '-temptrend_abs.sc', paste0('-', terms))
  
  # fit full model with ML for model comparison
  modTJbetadrops[[1]] <- lme(formula(paste0('Jbetatrendrem0 ~ ', paste(terms, collapse = ' + '))),
                             random = randef, weights = varef, data = trends[i,], method = 'ML')
  
  # w/out temptrend
  modTJbetadrops[[2]] <- lme(formula(paste0('Jbetatrendrem0 ~ ', paste(gsub('temptrend_abs.sc\\*', '', terms), collapse = ' + '))),
                             random = list(STUDY_ID = ~ 1, rarefyID = ~1), weights = varef, data = trends[i,], method = 'ML')
  
  for(j in 1:length(terms)){
    print(j)
    tryCatch({
      modTJbetadrops[[j+2]] <- lme(formula(paste0('Jbetatrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                                   random = randef, weights = varef, data = trends[i,], method = 'ML')
    }, error = function(e){
      print('going to optim (Jbeta)')
      tryCatch({
        modTJbetadrops[[j+2]] <- lme(formula(paste0('Jbetatrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                                     random = randef, weights = varef, data = trends[i,], method = 'ML',
                                     control = lmeControl(opt = 'optim'))
        
      }, error = function(e){
        print('going to more iters (Jbeta)') 
        tryCatch({
          modTJbetadrops[[j+2]] <- lme(formula(paste0('Jbetatrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                                       random = randef, weights = varef, data = trends[i,], method = 'ML',
                                       control = lmeControl(maxIter = 100, msMaxIter = 100, 
                                                            niterEM = 50, msMaxEval = 500))
        }, error= function(e){
          print('giving up on this one (Jbeta)')
          modTJbetadrops[[j+2]] <- NA
        })
      }
      )
    }
    )
  }
  aicsJbeta <- sapply(modTJbetadrops, AICnas)
}

if(runHorn){
  i2 <- trends[, complete.cases(Horntrendrem0, REALM, tempave_metab.sc, seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]
  
  modTHorndrops <- vector('list', length(terms)+2)
  names(modTHorndrops) <- c('full', '-temptrend_abs.sc', paste0('-', terms))
  modTHorndrops[[1]] <- lme(formula(paste0('Horntrendrem0 ~ ', paste(terms, collapse = ' + '))),
                            random = randef, weights = varef, data = trends[i2,], method = 'ML')
  modTHorndrops[[2]] <- lme(formula(paste0('Horntrendrem0 ~ ', paste(gsub('temptrend_abs.sc\\*', '', terms), collapse = ' + '))),
                            random = list(STUDY_ID = ~ 1, rarefyID = ~1), weights = varef, data = trends[i2,], method = 'ML')
  
  for(j in 1:length(terms)){
    print(j)
    tryCatch({
      modTHorndrops[[j+2]] <- lme(formula(paste0('Horntrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                                  random = randef, weights = varef, data = trends[i2,], method = 'ML')
    }, error = function(e){
      print('going to optim (Horn)')
      tryCatch({
        modTHorndrops[[j+2]] <- lme(formula(paste0('Horntrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                                    random = randef, weights = varef, data = trends[i2,], method = 'ML',
                                    control = lmeControl(opt = 'optim'))
        
      }, error = function(e){
        print('going to more iters (Horn)') 
        tryCatch({
          modTHorndrops[[j+2]] <- lme(formula(paste0('Horntrendrem0 ~ ', paste(terms[-j], collapse = ' + '))),
                                      random = randef, weights = varef, data = trends[i2,], method = 'ML',
                                      control = lmeControl(maxIter = 100, msMaxIter = 100, 
                                                           niterEM = 50, msMaxEval = 500))
          
        }, error= function(e){
          print('giving up on this one (Horn)')
          modTHorndrops[[j+2]] <- NA
        })
      })
    })
  }
  aicsHorn <- sapply(modTHorndrops, AICnas)
}

# if there was anything new
if(runJtu | runJbeta | runHorn){
  if(!exists('aicsfromfull')){
    aicsfromfull <- data.frame(mod = names(aicsJtu))
  }
  
  # subtract full from each model AIC. Negative means term removal is supported. Positive means full is the better model.
  if(runJtu){
    aicsfromfull$dAIC_Jtu <- aicsJtu - aicsJtu[1]
  }
  if(runJbeta){
    aicsfromfull$dAIC_Jbeta <- aicsJbeta - aicsJbeta[1]
  }
  if(runHorn){
    aicsfromfull$dAIC_Horn <- aicsHorn - aicsHorn[1]
  }
  
  # write out
  write.csv(aicsfromfull, file = 'output/aics_from_full.csv', row.names = FALSE)
}

aicsfromfull
```

    ##                                         mod   dAIC_Jtu  dAIC_Jbeta   dAIC_Horn
    ## 1                                      full   0.000000   0.0000000   0.0000000
    ## 2                         -temptrend_abs.sc 153.675151 148.3806423 192.7915870
    ## 3                   -temptrend_abs.sc*REALM  -2.910674   5.6796861   3.0212006
    ## 4                   -temptrend_abs.sc*tsign  17.489508   9.8748407  24.7129173
    ## 5        -temptrend_abs.sc*tempave_metab.sc  12.841486  25.9051193  45.6279310
    ## 6                 -temptrend_abs.sc*seas.sc  -2.619386  -3.6269879  10.1424377
    ## 7            -temptrend_abs.sc*microclim.sc   9.208844   6.2335885  -3.6834624
    ## 8                 -temptrend_abs.sc*mass.sc   7.402343  12.4013310  -0.5844054
    ## 9                -temptrend_abs.sc*speed.sc   1.437238  -1.2770636  -0.5801870
    ## 10        -temptrend_abs.sc*consumerfrac.sc  41.459152   3.4897487 112.3183306
    ## 11                -temptrend_abs.sc*nspp.sc  -3.859265  -1.5612518  -0.6699811
    ## 12  -temptrend_abs.sc*thermal_bias.sc:tsign  -1.667301   4.5298601   6.0876346
    ## 13                 -temptrend_abs.sc*npp.sc  15.013137  23.0866513 110.2403343
    ## 14                 -temptrend_abs.sc*veg.sc   1.541348   8.1718838   2.9372433
    ## 15            -temptrend_abs.sc*duration.sc   1.830280  24.9672469  11.1598650
    ## 16 -temptrend_abs.sc*human_bowler.sc:REALM2  16.200259   0.9073633   8.1089261

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
xlims <- range(aicsfromfulllong$dAIC_tr, na.rm = TRUE)
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
black is for Morisita-Horn. Clear that removing temperature trend makes
the model quite a bit worse and has the biggest effect.

## Simplify the full models

This takes a couple days on a laptop to run if temp/ files not
available.

``` r
i1 <- trends[, complete.cases(Jtutrendrem0, REALM, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc)]
i2 <- trends[, complete.cases(Jbetatrendrem0, REALM, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc)]
i3 <- trends[, complete.cases(Horntrendrem0, REALM, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# simplify the full models
if(file.exists('temp/modTsimpJturem0.rds')){
  modTsimpJturem0 <- readRDS('temp/modTsimpJturem0.rds')
} else {
  modTfullJturem0ML <- lme(Jtutrendrem0 ~ temptrend_abs.sc*REALM + 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*consumerfrac.sc +
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*veg.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc:REALM2,
                       random = randef, weights = varef, data = trends[i1,], method = 'ML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  modTsimpJturem0 <- stepAIC(modTfullJturem0ML, direction = 'backward')
  saveRDS(modTsimpJturem0, file = 'temp/modTsimpJturem0.rds')
}

if(file.exists('temp/modTsimpJbetarem0.rds')){
  modTsimpJbetarem0 <- readRDS('temp/modTsimpJbetarem0.rds')
} else {
  modTfullJbetarem0ML <- lme(Jbetatrendrem0 ~ temptrend_abs.sc*REALM + 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*consumerfrac.sc +
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*veg.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc:REALM2,
                       random = randef, weights = varef, data = trends[i2,], method = 'ML')
  modTsimpJbetarem0 <- stepAIC(modTfullJbetarem0ML, direction = 'backward')
  saveRDS(modTsimpJbetarem0, file = 'temp/modTsimpJbetarem0.rds')
}

if(file.exists('temp/modTsimpHornrem0.rds')){
  modTsimpHornrem0 <- readRDS('temp/modTsimpHornrem0.rds')
} else {
  modTfullHornrem0ML <- lme(Horntrendrem0 ~ temptrend_abs.sc*REALM + 
                        temptrend_abs.sc*tsign +
                        temptrend_abs.sc*tempave_metab.sc + 
                        temptrend_abs.sc*seas.sc + 
                        temptrend_abs.sc*microclim.sc + 
                        temptrend_abs.sc*mass.sc + 
                        temptrend_abs.sc*speed.sc + 
                        temptrend_abs.sc*consumerfrac.sc +
                        temptrend_abs.sc*nspp.sc +
                        temptrend_abs.sc*thermal_bias.sc:tsign +
                        temptrend_abs.sc*npp.sc +
                        temptrend_abs.sc*veg.sc +
                        temptrend_abs.sc*duration.sc +
                        temptrend_abs.sc*human_bowler.sc:REALM2,
                      random = randef, weights = varef, data = trends[i3,], method = 'ML')
  modTsimpHornrem0 <- stepAIC(modTfullHornrem0ML, direction = 'backward')
  saveRDS(modTsimpHornrem0, file = 'temp/modTsimpHornrem0.rds')
}

summary(modTsimpJturem0)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: trends[i1, ] 
    ##         AIC       BIC   logLik
    ##   -103051.6 -102788.3 51556.79
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.009359204 (Intr)
    ## temptrend_abs.sc 0.024191157 -0.98 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01090463 2.006731
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.129313 
    ## Fixed effects: Jtutrendrem0 ~ temptrend_abs.sc + tsign + tempave_metab.sc +      seas.sc + microclim.sc + mass.sc + speed.sc + consumerfrac.sc +      npp.sc + veg.sc + duration.sc + temptrend_abs.sc:tsign +      temptrend_abs.sc:tempave_metab.sc + temptrend_abs.sc:microclim.sc +      temptrend_abs.sc:consumerfrac.sc + tsign:thermal_bias.sc +      temptrend_abs.sc:npp.sc + human_bowler.sc:REALM2 + temptrend_abs.sc:tsign:thermal_bias.sc +      temptrend_abs.sc:human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.007626715 0.001313336 35762  5.807133  0.0000
    ## temptrend_abs.sc                                  0.002242230 0.003323388 35762  0.674682  0.4999
    ## tsign1                                           -0.002213744 0.000635908 35762 -3.481233  0.0005
    ## tempave_metab.sc                                 -0.003425682 0.000837651 35762 -4.089631  0.0000
    ## seas.sc                                           0.000600012 0.000371436 35762  1.615385  0.1062
    ## microclim.sc                                      0.000580912 0.000329202 35762  1.764610  0.0776
    ## mass.sc                                          -0.001354388 0.000369086 35762 -3.669576  0.0002
    ## speed.sc                                          0.000923336 0.000415938 35762  2.219886  0.0264
    ## consumerfrac.sc                                  -0.000515151 0.000311974 35762 -1.651262  0.0987
    ## npp.sc                                           -0.000251530 0.000441304 35762 -0.569969  0.5687
    ## veg.sc                                           -0.001119135 0.000356679 35762 -3.137658  0.0017
    ## duration.sc                                      -0.000872264 0.000405010 35762 -2.153686  0.0313
    ## temptrend_abs.sc:tsign1                          -0.001478051 0.001813657 35762 -0.814956  0.4151
    ## temptrend_abs.sc:tempave_metab.sc                 0.006151630 0.002242899 35762  2.742714  0.0061
    ## temptrend_abs.sc:microclim.sc                    -0.003113480 0.000932424 35762 -3.339123  0.0008
    ## temptrend_abs.sc:consumerfrac.sc                  0.003884322 0.000730491 35762  5.317409  0.0000
    ## tsign-1:thermal_bias.sc                          -0.000933793 0.000645813 35762 -1.445918  0.1482
    ## tsign1:thermal_bias.sc                           -0.000382418 0.000441970 35762 -0.865258  0.3869
    ## temptrend_abs.sc:npp.sc                           0.004160100 0.001272634 35762  3.268891  0.0011
    ## human_bowler.sc:REALM2TerrFresh                   0.000083579 0.000686413 35762  0.121762  0.9031
    ## human_bowler.sc:REALM2Marine                      0.000377053 0.000393159 35762  0.959035  0.3375
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.000686557 0.001335537 35762  0.514069  0.6072
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.002279928 0.001097409 35762  2.077555  0.0378
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.001631622 0.001271849 35762 -1.282874  0.1995
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.003722025 0.001074965 35762 -3.462462  0.0005
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. npp.sc veg.sc drtn.s tm_.:1 tm_.:_. tmptrnd_bs.sc:m. tmptrnd_bs.sc:c. t-1:_. ts1:_. tmptrnd_bs.sc:n. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.860                                                                                                                                                                                                                  
    ## tsign1                                           -0.318  0.167                                                                                                                                                                                                           
    ## tempave_metab.sc                                 -0.095  0.079  0.081                                                                                                                                                                                                    
    ## seas.sc                                           0.099 -0.049 -0.051  0.071                                                                                                                                                                                             
    ## microclim.sc                                     -0.054  0.062  0.028 -0.203  0.032                                                                                                                                                                                      
    ## mass.sc                                          -0.012  0.059 -0.021  0.072  0.147  0.013                                                                                                                                                                               
    ## speed.sc                                          0.096 -0.056 -0.036 -0.113 -0.045  0.075 -0.560                                                                                                                                                                        
    ## consumerfrac.sc                                   0.085 -0.084  0.000 -0.096 -0.041  0.047  0.040 -0.123                                                                                                                                                                 
    ## npp.sc                                           -0.134  0.131  0.059  0.050 -0.053 -0.222 -0.045  0.130 -0.023                                                                                                                                                          
    ## veg.sc                                           -0.134  0.076 -0.031 -0.074 -0.424 -0.074  0.120 -0.043  0.052 -0.064                                                                                                                                                   
    ## duration.sc                                      -0.174  0.145 -0.163  0.103 -0.095 -0.049 -0.010 -0.025  0.019  0.023 -0.034                                                                                                                                            
    ## temptrend_abs.sc:tsign1                           0.157 -0.361 -0.516 -0.057 -0.024 -0.009 -0.009  0.017 -0.001 -0.050  0.004 -0.049                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                 0.074 -0.076 -0.047 -0.820 -0.072  0.148  0.070 -0.182  0.112 -0.049 -0.110 -0.013  0.046                                                                                                                              
    ## temptrend_abs.sc:microclim.sc                     0.057 -0.056 -0.024  0.131 -0.003 -0.763  0.006 -0.048 -0.023  0.210 -0.025  0.056  0.024 -0.112                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.089  0.092  0.011  0.072 -0.006 -0.023 -0.073  0.061 -0.783  0.001 -0.014 -0.041  0.034 -0.121   0.003                                                                                                               
    ## tsign-1:thermal_bias.sc                           0.010  0.041 -0.089  0.133 -0.155 -0.005 -0.025  0.006  0.006 -0.025  0.035  0.076 -0.044 -0.116   0.002           -0.002                                                                                              
    ## tsign1:thermal_bias.sc                           -0.052  0.030  0.055  0.228 -0.272 -0.089 -0.041  0.006 -0.013 -0.091  0.079  0.047 -0.003 -0.221   0.045            0.021            0.323                                                                             
    ## temptrend_abs.sc:npp.sc                           0.125 -0.173 -0.037 -0.017 -0.008  0.191  0.032 -0.048 -0.014 -0.743 -0.004 -0.038  0.018  0.070  -0.316           -0.002            0.024  0.081                                                                      
    ## human_bowler.sc:REALM2TerrFresh                  -0.029  0.031  0.021 -0.007 -0.092  0.090  0.046 -0.042 -0.031 -0.059  0.009 -0.024 -0.015 -0.009   0.005            0.019            0.074  0.118 -0.018                                                               
    ## human_bowler.sc:REALM2Marine                      0.009 -0.004 -0.014  0.092 -0.124  0.030 -0.016  0.032 -0.069 -0.231  0.062  0.019  0.006 -0.080  -0.027            0.069            0.085  0.112  0.176            0.033                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.058 -0.129 -0.059 -0.154  0.016  0.043  0.025 -0.044  0.009  0.010 -0.057 -0.048  0.241  0.178  -0.083            0.005           -0.488 -0.281  0.006           -0.085      -0.053                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.010 -0.006  0.017 -0.234  0.003  0.065  0.034 -0.116  0.010  0.083 -0.064 -0.032 -0.026  0.304  -0.179            0.003           -0.242 -0.714 -0.061           -0.091      -0.085       0.369                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.018 -0.031 -0.010 -0.025 -0.006 -0.042  0.005 -0.016  0.023  0.017  0.003  0.014  0.024  0.031   0.036           -0.020           -0.071 -0.086 -0.013           -0.630      -0.014       0.137  0.106               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.003  0.009  0.002 -0.090 -0.003 -0.016  0.031 -0.040  0.128  0.180  0.022 -0.020 -0.025  0.126  -0.020           -0.139           -0.046 -0.080 -0.231           -0.014      -0.771       0.044  0.085  0.011        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.63842790 -0.23658104 -0.02245775  0.26774537  5.41648169 
    ## 
    ## Number of Observations: 36017
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    231                  36017

``` r
summary(modTsimpJbetarem0)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -133460.2 -133154.5 66766.11
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.006875615 (Intr)
    ## temptrend_abs.sc 0.016512162 -0.015
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:  0.00367393 0.9559692
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.889706 
    ## Fixed effects: Jbetatrendrem0 ~ temptrend_abs.sc + REALM + tsign + tempave_metab.sc +      microclim.sc + mass.sc + consumerfrac.sc + nspp.sc + npp.sc +      veg.sc + duration.sc + temptrend_abs.sc:REALM + temptrend_abs.sc:tsign +      temptrend_abs.sc:tempave_metab.sc + temptrend_abs.sc:microclim.sc +      temptrend_abs.sc:consumerfrac.sc + temptrend_abs.sc:nspp.sc +      tsign:thermal_bias.sc + temptrend_abs.sc:npp.sc + temptrend_abs.sc:veg.sc +      human_bowler.sc:REALM2 + temptrend_abs.sc:tsign:thermal_bias.sc +      temptrend_abs.sc:human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.001316067 0.003379887 35759  0.389382  0.6970
    ## temptrend_abs.sc                                  0.024519339 0.009152657 35759  2.678931  0.0074
    ## REALMMarine                                       0.007140466 0.003560174   228  2.005651  0.0461
    ## REALMTerrestrial                                  0.007939402 0.003480317   228  2.281229  0.0235
    ## tsign1                                           -0.001059073 0.000390989 35759 -2.708705  0.0068
    ## tempave_metab.sc                                 -0.003003806 0.000575397 35759 -5.220409  0.0000
    ## microclim.sc                                      0.000392647 0.000200379 35759  1.959523  0.0501
    ## mass.sc                                          -0.000980677 0.000251343 35759 -3.901743  0.0001
    ## consumerfrac.sc                                  -0.000014445 0.000197364 35759 -0.073190  0.9417
    ## nspp.sc                                          -0.000483653 0.000324584 35759 -1.490072  0.1362
    ## npp.sc                                           -0.000498440 0.000272470 35759 -1.829340  0.0674
    ## veg.sc                                            0.000553301 0.000476348 35759  1.161548  0.2454
    ## duration.sc                                      -0.001515978 0.000287331 35759 -5.276070  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.025612882 0.009618098 35759 -2.662988  0.0077
    ## temptrend_abs.sc:REALMTerrestrial                -0.011515664 0.009404729 35759 -1.224455  0.2208
    ## temptrend_abs.sc:tsign1                          -0.001243098 0.001233498 35759 -1.007782  0.3136
    ## temptrend_abs.sc:tempave_metab.sc                 0.007153231 0.001723894 35759  4.149461  0.0000
    ## temptrend_abs.sc:microclim.sc                    -0.001842043 0.000601700 35759 -3.061396  0.0022
    ## temptrend_abs.sc:consumerfrac.sc                  0.000993789 0.000489573 35759  2.029908  0.0424
    ## temptrend_abs.sc:nspp.sc                          0.001312372 0.000882727 35759  1.486724  0.1371
    ## tsign-1:thermal_bias.sc                          -0.000435198 0.000410740 35759 -1.059545  0.2894
    ## tsign1:thermal_bias.sc                           -0.000143220 0.000250881 35759 -0.570868  0.5681
    ## temptrend_abs.sc:npp.sc                           0.004095385 0.000859939 35759  4.762412  0.0000
    ## temptrend_abs.sc:veg.sc                          -0.003800729 0.001337859 35759 -2.840904  0.0045
    ## human_bowler.sc:REALM2TerrFresh                   0.000607178 0.000362912 35759  1.673071  0.0943
    ## human_bowler.sc:REALM2Marine                     -0.000184572 0.000237244 35759 -0.777983  0.4366
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.000327224 0.000902334 35759 -0.362641  0.7169
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.001799921 0.000714435 35759  2.519363  0.0118
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.001846528 0.000891197 35759 -2.071964  0.0383
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.000521369 0.000697667 35759 -0.747305  0.4549
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. mcrcl. mss.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:m. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.565                                                                                                                                                                                                                                                                                   
    ## REALMMarine                                      -0.949  0.549                                                                                                                                                                                                                                                                            
    ## REALMTerrestrial                                 -0.765  0.391  0.722                                                                                                                                                                                                                                                                     
    ## tsign1                                           -0.088  0.041  0.013 -0.012                                                                                                                                                                                                                                                              
    ## tempave_metab.sc                                  0.123 -0.049 -0.106 -0.257  0.095                                                                                                                                                                                                                                                       
    ## microclim.sc                                     -0.066  0.050  0.071  0.039  0.013 -0.195                                                                                                                                                                                                                                                
    ## mass.sc                                           0.168  0.015 -0.118 -0.096 -0.025  0.032  0.041                                                                                                                                                                                                                                         
    ## consumerfrac.sc                                  -0.006  0.009  0.009  0.088 -0.036 -0.118  0.028 -0.001                                                                                                                                                                                                                                  
    ## nspp.sc                                          -0.013 -0.021 -0.045 -0.016  0.062 -0.188 -0.145  0.026  0.129                                                                                                                                                                                                                           
    ## npp.sc                                           -0.001  0.014 -0.022  0.038  0.080  0.104 -0.138  0.047 -0.039 -0.245                                                                                                                                                                                                                    
    ## veg.sc                                           -0.449  0.370  0.471  0.019 -0.011 -0.014 -0.013 -0.004 -0.017  0.029 -0.197                                                                                                                                                                                                             
    ## duration.sc                                      -0.048  0.019  0.020 -0.026 -0.150  0.160 -0.034 -0.025  0.008 -0.141  0.074 -0.011                                                                                                                                                                                                      
    ## temptrend_abs.sc:REALMMarine                      0.551 -0.958 -0.566 -0.372 -0.012  0.041 -0.052  0.008 -0.009  0.043  0.020 -0.389  0.019                                                                                                                                                                                               
    ## temptrend_abs.sc:REALMTerrestrial                 0.385 -0.774 -0.366 -0.531  0.018  0.151  0.002  0.001 -0.062  0.034 -0.053  0.008  0.025  0.730                                                                                                                                                                                        
    ## temptrend_abs.sc:tsign1                           0.053 -0.126 -0.024  0.002 -0.440 -0.049 -0.010  0.017  0.007 -0.003 -0.068 -0.025 -0.037  0.047      0.004                                                                                                                                                                             
    ## temptrend_abs.sc:tempave_metab.sc                -0.035  0.101  0.029  0.158 -0.101 -0.521  0.030 -0.007  0.109  0.081 -0.046  0.000 -0.051 -0.072     -0.281      0.044                                                                                                                                                                  
    ## temptrend_abs.sc:microclim.sc                     0.048 -0.064 -0.052 -0.018 -0.017  0.124 -0.781 -0.027 -0.014  0.084  0.162  0.003  0.028  0.070     -0.029      0.015  0.028                                                                                                                                                           
    ## temptrend_abs.sc:consumerfrac.sc                 -0.002 -0.020 -0.001 -0.057  0.031  0.090  0.005 -0.035 -0.715 -0.100  0.017  0.003 -0.008  0.021      0.086      0.015 -0.152  -0.021                                                                                                                                                   
    ## temptrend_abs.sc:nspp.sc                         -0.007  0.009  0.026  0.020 -0.031  0.112  0.073 -0.020 -0.070 -0.592  0.133 -0.031 -0.010 -0.058     -0.050     -0.016 -0.083  -0.048            0.133                                                                                                                                  
    ## tsign-1:thermal_bias.sc                           0.033 -0.005 -0.020 -0.066 -0.088  0.157  0.011  0.000  0.000 -0.053 -0.019  0.005  0.077  0.011      0.048     -0.053 -0.072  -0.009           -0.002            0.013                                                                                                                 
    ## tsign1:thermal_bias.sc                            0.067 -0.023 -0.052 -0.144  0.041  0.318 -0.089  0.063 -0.049 -0.085 -0.088  0.030  0.075  0.014      0.107      0.005 -0.198   0.037            0.031            0.004             0.258                                                                                               
    ## temptrend_abs.sc:npp.sc                           0.014  0.023  0.005 -0.044 -0.050 -0.043  0.148 -0.006  0.009  0.144 -0.737  0.153 -0.069 -0.063      0.061      0.029  0.076  -0.283           -0.025           -0.143             0.018  0.094                                                                                        
    ## temptrend_abs.sc:veg.sc                           0.358 -0.449 -0.374  0.005 -0.008 -0.005  0.006  0.003  0.009 -0.035  0.174 -0.815  0.021  0.474     -0.006      0.051 -0.023   0.035           -0.004            0.018            -0.006 -0.055 -0.264                                                                                 
    ## human_bowler.sc:REALM2TerrFresh                  -0.107  0.124  0.114 -0.027  0.016  0.003  0.116  0.011 -0.039  0.005 -0.130  0.256 -0.029 -0.134      0.023     -0.024 -0.008  -0.031            0.022            0.025             0.054  0.129  0.059            -0.302                                                               
    ## human_bowler.sc:REALM2Marine                      0.023 -0.010 -0.013 -0.055 -0.019  0.134  0.051  0.035 -0.087 -0.045 -0.235  0.049  0.024  0.005      0.042      0.008 -0.087  -0.040            0.084            0.056             0.072  0.069  0.177            -0.044            0.045                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.013  0.012  0.012  0.071 -0.069 -0.103  0.002 -0.016  0.013  0.003  0.012 -0.005 -0.064 -0.026     -0.112      0.274  0.225  -0.029           -0.010            0.014            -0.414 -0.258 -0.009            -0.001           -0.092      -0.044                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.021  0.053  0.013  0.099 -0.006 -0.161  0.018 -0.045  0.025  0.037  0.090 -0.047 -0.058 -0.031     -0.166     -0.018  0.374  -0.096           -0.027           -0.033            -0.186 -0.712 -0.085             0.049           -0.104      -0.070       0.384                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.145 -0.206 -0.155  0.024 -0.015 -0.017 -0.068 -0.010  0.025  0.006  0.100 -0.354  0.025  0.220     -0.026      0.044  0.037   0.091           -0.026           -0.019            -0.054 -0.109 -0.143             0.463           -0.752      -0.036       0.124  0.123               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.015  0.016  0.009  0.048  0.006 -0.101 -0.024 -0.016  0.134  0.037  0.181 -0.040 -0.014 -0.008     -0.054     -0.015  0.126  -0.020           -0.137           -0.063            -0.048 -0.066 -0.214             0.059           -0.040      -0.759       0.041  0.074  0.037        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.39440359 -0.32029428 -0.03087913  0.31695654  8.33262939 
    ## 
    ## Number of Observations: 36017
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    231                  36017

``` r
summary(modTsimpHornrem0)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -98480.02 -98191.96 49274.01
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.01429668 (Intr)
    ## temptrend_abs.sc 0.02508293 0.046 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01877495  2.44238
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.326836 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc + REALM + tsign + tempave_metab.sc +      seas.sc + mass.sc + speed.sc + consumerfrac.sc + nspp.sc +      npp.sc + veg.sc + duration.sc + temptrend_abs.sc:REALM +      temptrend_abs.sc:tsign + temptrend_abs.sc:consumerfrac.sc +      tsign:thermal_bias.sc + temptrend_abs.sc:npp.sc + temptrend_abs.sc:veg.sc +      human_bowler.sc:REALM2 + temptrend_abs.sc:tsign:thermal_bias.sc +      temptrend_abs.sc:human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.006250838 0.007649521 35102  0.817154  0.4138
    ## temptrend_abs.sc                                  0.028516450 0.014915822 35102  1.911826  0.0559
    ## REALMMarine                                       0.004883179 0.008116819   197  0.601612  0.5481
    ## REALMTerrestrial                                  0.011158345 0.007317293   197  1.524928  0.1289
    ## tsign1                                           -0.001980450 0.000764540 35102 -2.590382  0.0096
    ## tempave_metab.sc                                 -0.007298768 0.001006061 35102 -7.254799  0.0000
    ## seas.sc                                          -0.002380184 0.000651974 35102 -3.650735  0.0003
    ## mass.sc                                          -0.000988552 0.000542207 35102 -1.823203  0.0683
    ## speed.sc                                          0.001328395 0.000645638 35102  2.057491  0.0396
    ## consumerfrac.sc                                   0.000616883 0.000403809 35102  1.527661  0.1266
    ## nspp.sc                                           0.000822459 0.000531277 35102  1.548080  0.1216
    ## npp.sc                                           -0.001374283 0.000553700 35102 -2.481997  0.0131
    ## veg.sc                                            0.000952492 0.001250726 35102  0.761552  0.4463
    ## duration.sc                                      -0.002132503 0.000555488 35102 -3.838970  0.0001
    ## temptrend_abs.sc:REALMMarine                     -0.018527495 0.015710874 35102 -1.179278  0.2383
    ## temptrend_abs.sc:REALMTerrestrial                 0.004776132 0.014827907 35102  0.322104  0.7474
    ## temptrend_abs.sc:tsign1                          -0.005311137 0.001911520 35102 -2.778488  0.0055
    ## temptrend_abs.sc:consumerfrac.sc                  0.005433770 0.000828454 35102  6.558929  0.0000
    ## tsign-1:thermal_bias.sc                          -0.001596230 0.000766622 35102 -2.082159  0.0373
    ## tsign1:thermal_bias.sc                           -0.000464288 0.000572131 35102 -0.811506  0.4171
    ## temptrend_abs.sc:npp.sc                           0.006912130 0.001402662 35102  4.927865  0.0000
    ## temptrend_abs.sc:veg.sc                          -0.005418716 0.002196751 35102 -2.466695  0.0136
    ## human_bowler.sc:REALM2TerrFresh                   0.002095166 0.001021467 35102  2.051134  0.0403
    ## human_bowler.sc:REALM2Marine                      0.000318239 0.000482834 35102  0.659107  0.5098
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.000950297 0.001350700 35102 -0.703559  0.4817
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.002012589 0.001098795 35102  1.831632  0.0670
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.003525913 0.001483837 35102 -2.376213  0.0175
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.003184388 0.001202124 35102 -2.648968  0.0081
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tmptrnd_bs.sc:c. t-1:_. ts1:_. tmptrnd_bs.sc:n. tmptrnd_bs.sc:v. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.472                                                                                                                                                                                                                                              
    ## REALMMarine                                      -0.957  0.456                                                                                                                                                                                                                                       
    ## REALMTerrestrial                                 -0.747  0.335  0.698                                                                                                                                                                                                                                
    ## tsign1                                           -0.069  0.049  0.003  0.012                                                                                                                                                                                                                         
    ## tempave_metab.sc                                  0.100 -0.003 -0.081 -0.205  0.009                                                                                                                                                                                                                  
    ## seas.sc                                          -0.085 -0.012  0.130 -0.038 -0.076  0.083                                                                                                                                                                                                           
    ## mass.sc                                           0.078  0.006 -0.054 -0.043  0.002  0.081  0.078                                                                                                                                                                                                    
    ## speed.sc                                          0.073  0.004 -0.027 -0.054 -0.057 -0.058  0.024 -0.428                                                                                                                                                                                             
    ## consumerfrac.sc                                  -0.013  0.009  0.009  0.072  0.031 -0.087 -0.058  0.004 -0.063                                                                                                                                                                                      
    ## nspp.sc                                          -0.011 -0.025 -0.031 -0.030  0.054 -0.155  0.037 -0.065  0.132  0.089                                                                                                                                                                               
    ## npp.sc                                            0.000  0.037 -0.013  0.067  0.014 -0.041 -0.137 -0.046  0.136 -0.018 -0.203                                                                                                                                                                        
    ## veg.sc                                           -0.539  0.287  0.557  0.044 -0.007 -0.003  0.035  0.006 -0.016  0.002  0.011 -0.159                                                                                                                                                                 
    ## duration.sc                                      -0.027  0.038  0.016  0.010 -0.157  0.127 -0.080  0.002 -0.005  0.016 -0.201  0.039 -0.005                                                                                                                                                          
    ## temptrend_abs.sc:REALMMarine                      0.460 -0.958 -0.466 -0.319 -0.023  0.014  0.019  0.008  0.004 -0.017  0.016 -0.006 -0.301 -0.001                                                                                                                                                   
    ## temptrend_abs.sc:REALMTerrestrial                 0.319 -0.781 -0.301 -0.438 -0.009  0.020 -0.017 -0.009  0.010 -0.025  0.025 -0.054  0.008 -0.004  0.738                                                                                                                                            
    ## temptrend_abs.sc:tsign1                           0.043 -0.116 -0.014  0.000 -0.488 -0.021 -0.015  0.005  0.022 -0.013 -0.016 -0.038 -0.011 -0.007  0.043      0.012                                                                                                                                 
    ## temptrend_abs.sc:consumerfrac.sc                  0.008 -0.018 -0.012 -0.034 -0.011  0.028  0.003 -0.027  0.019 -0.712 -0.017 -0.010 -0.009 -0.008  0.029      0.040      0.038                                                                                                                      
    ## tsign-1:thermal_bias.sc                           0.049 -0.002 -0.044 -0.067 -0.099  0.201 -0.173 -0.037 -0.004 -0.003 -0.047 -0.017  0.007  0.071  0.006      0.040     -0.025  0.003                                                                                                               
    ## tsign1:thermal_bias.sc                            0.083 -0.006 -0.086 -0.114  0.060  0.324 -0.303  0.000 -0.031 -0.008 -0.073 -0.113  0.008  0.056  0.001      0.057      0.010  0.013            0.370                                                                                              
    ## temptrend_abs.sc:npp.sc                           0.033 -0.026 -0.028 -0.057 -0.004  0.060 -0.044  0.024 -0.065 -0.017  0.064 -0.650  0.098 -0.017 -0.015      0.076      0.001  0.014            0.018  0.095                                                                                       
    ## temptrend_abs.sc:veg.sc                           0.317 -0.445 -0.326  0.003 -0.022 -0.013  0.009 -0.003  0.010 -0.003 -0.016  0.148 -0.616  0.006  0.467     -0.024      0.061  0.007           -0.007 -0.039 -0.231                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.056  0.056  0.057 -0.026  0.018  0.034 -0.048  0.031 -0.033 -0.029  0.018 -0.082  0.138 -0.009 -0.062      0.035     -0.013  0.011            0.057  0.115  0.043           -0.194                                                               
    ## human_bowler.sc:REALM2Marine                      0.032 -0.007 -0.036 -0.037 -0.015  0.092 -0.127 -0.010  0.017 -0.054 -0.025 -0.230  0.036  0.034  0.013      0.025      0.013  0.068            0.113  0.140  0.174           -0.041            0.035                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.007 -0.013  0.008  0.036 -0.038  0.016  0.042  0.007 -0.007 -0.014  0.000  0.009 -0.009 -0.033 -0.003     -0.060      0.262  0.026           -0.452 -0.220 -0.033            0.018           -0.057      -0.046                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.017  0.009  0.015  0.042  0.023  0.040  0.075  0.003 -0.022 -0.028 -0.012  0.100 -0.027 -0.028  0.002     -0.078     -0.029  0.034           -0.210 -0.614 -0.159            0.067           -0.062      -0.085       0.361                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.116 -0.209 -0.121  0.017 -0.016 -0.005  0.027 -0.009  0.003  0.004 -0.009  0.071 -0.262  0.017  0.220     -0.031      0.049  0.002           -0.063 -0.091 -0.112            0.492           -0.537      -0.031       0.142  0.132               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.019  0.009  0.024  0.029  0.007 -0.039  0.047  0.011 -0.016  0.127  0.002  0.174 -0.024 -0.014 -0.015     -0.028     -0.009 -0.147           -0.062 -0.094 -0.265            0.060           -0.021      -0.728       0.021  0.052  0.033        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.66375573 -0.23464298 -0.02584232  0.23186565  5.76940752 
    ## 
    ## Number of Observations: 35327
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    200                  35327

## Make realm-specific models

``` r
i1 <- trends[, REALM == 'Terrestrial' & complete.cases(Horntrendrem0, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc)]
i2 <- trends[, REALM == 'Freshwater' & complete.cases(Horntrendrem0, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc)] # no consumerfrac
i3 <- trends[, REALM == 'Marine' & complete.cases(Horntrendrem0, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             duration.sc, human_bowler.sc)] # no veg

print(paste('Terrestrial', sum(i1)))
```

    ## [1] "Terrestrial 2299"

``` r
print(paste('Freshwater', sum(i2)))
```

    ## [1] "Freshwater 608"

``` r
print(paste('Marine', sum(i3)))
```

    ## [1] "Marine 32420"

``` r
randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# land
if(file.exists('temp/modTfullHornTerr.rds')){
  modTfullHornTerr <- readRDS('temp/modTfullHornTerr.rds')
} else {
  modTfullHornTerr <- lme(Horntrendrem0 ~ 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*consumerfrac.sc +
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*veg.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc,
                       random = randef, weights = varef, data = trends[i1,], method = 'REML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullHornTerr, file = 'temp/modTfullHornTerr.rds')
}

# freshwater
if(file.exists('temp/modTfullHornFresh.rds')){
  modTfullHornFresh <- readRDS('temp/modTfullHornFresh.rds')
} else {
  modTfullHornFresh <- lme(Horntrendrem0 ~ 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*veg.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc,
                       random = randef, weights = varef, data = trends[i2,], method = 'REML')
  saveRDS(modTfullHornFresh, file = 'temp/modTfullHornFresh.rds')
}

# marine
if(file.exists('temp/modTfullHornMar.rds')){
  modTfullHornMar <- readRDS('temp/modTfullHornMar.rds')
} else {
  modTfullHornMar <- lme(Horntrendrem0 ~ 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*consumerfrac.sc +
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc,
                       random = randef, weights = varef, data = trends[i3,], method = 'REML')
  saveRDS(modTfullHornMar, file = 'temp/modTfullHornMar.rds')
}

summary(modTfullHornTerr)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i1, ] 
    ##         AIC       BIC logLik
    ##   -7791.801 -7585.625 3931.9
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.009555472 (Intr)
    ## temptrend_abs.sc 0.022365121 0.313 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev: 0.004029061 1.695873
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.036154 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * tsign + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * duration.sc +      temptrend_abs.sc * human_bowler.sc 
    ##                                                 Value   Std.Error   DF    t-value p-value
    ## (Intercept)                               0.005370186 0.005866953 2180  0.9153280  0.3601
    ## temptrend_abs.sc                          0.011801970 0.016761061 2180  0.7041302  0.4814
    ## tsign1                                    0.008025381 0.003159965 2180  2.5397053  0.0112
    ## tempave_metab.sc                         -0.007369070 0.003513868 2180 -2.0971388  0.0361
    ## seas.sc                                   0.000744444 0.000907070 2180  0.8207127  0.4119
    ## microclim.sc                             -0.000226097 0.000618360 2180 -0.3656389  0.7147
    ## mass.sc                                  -0.001156008 0.001905833 2180 -0.6065630  0.5442
    ## speed.sc                                  0.003241949 0.002815481 2180  1.1514722  0.2497
    ## consumerfrac.sc                          -0.000081225 0.001661787 2180 -0.0488782  0.9610
    ## nspp.sc                                   0.000640954 0.001193980 2180  0.5368215  0.5914
    ## npp.sc                                    0.000513481 0.001205463 2180  0.4259615  0.6702
    ## veg.sc                                    0.000624416 0.000755592 2180  0.8263929  0.4087
    ## duration.sc                              -0.001406395 0.001190151 2180 -1.1816945  0.2375
    ## human_bowler.sc                           0.001125545 0.000516963 2180  2.1772259  0.0296
    ## temptrend_abs.sc:tsign1                  -0.001874351 0.007269033 2180 -0.2578542  0.7965
    ## temptrend_abs.sc:tempave_metab.sc         0.009049264 0.010197675 2180  0.8873850  0.3750
    ## temptrend_abs.sc:seas.sc                  0.001216727 0.002516110 2180  0.4835746  0.6287
    ## temptrend_abs.sc:microclim.sc            -0.000331152 0.001800633 2180 -0.1839083  0.8541
    ## temptrend_abs.sc:mass.sc                 -0.000069130 0.004769385 2180 -0.0144946  0.9884
    ## temptrend_abs.sc:speed.sc                 0.005064810 0.008597226 2180  0.5891215  0.5558
    ## temptrend_abs.sc:consumerfrac.sc         -0.005701039 0.005241337 2180 -1.0877069  0.2768
    ## temptrend_abs.sc:nspp.sc                 -0.002537190 0.002808060 2180 -0.9035382  0.3663
    ## tsign-1:thermal_bias.sc                   0.003304500 0.002433135 2180  1.3581245  0.1746
    ## tsign1:thermal_bias.sc                    0.000319074 0.000955843 2180  0.3338139  0.7386
    ## temptrend_abs.sc:npp.sc                  -0.001260031 0.003719838 2180 -0.3387328  0.7348
    ## temptrend_abs.sc:veg.sc                  -0.001837064 0.002319106 2180 -0.7921431  0.4284
    ## temptrend_abs.sc:duration.sc             -0.000243755 0.002176334 2180 -0.1120025  0.9108
    ## temptrend_abs.sc:human_bowler.sc         -0.002178483 0.001342933 2180 -1.6221829  0.1049
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.003666404 0.003634288 2180 -1.0088369  0.3132
    ## temptrend_abs.sc:tsign1:thermal_bias.sc  -0.006111275 0.002880406 2180 -2.1216715  0.0340
    ##  Correlation: 
    ##                                          (Intr) tmpt_. tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s hmn_b. tm_.:1 tmptrnd_bs.sc:t_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. tmptrnd_bs.sc:h_. t_.:-1
    ## temptrend_abs.sc                         -0.629                                                                                                                                                                                                                                                                                                                          
    ## tsign1                                   -0.215  0.019                                                                                                                                                                                                                                                                                                                   
    ## tempave_metab.sc                         -0.486  0.307 -0.209                                                                                                                                                                                                                                                                                                            
    ## seas.sc                                  -0.258  0.267 -0.078  0.049                                                                                                                                                                                                                                                                                                     
    ## microclim.sc                             -0.182  0.232  0.094 -0.018  0.473                                                                                                                                                                                                                                                                                              
    ## mass.sc                                   0.515 -0.246 -0.010 -0.141  0.008  0.039                                                                                                                                                                                                                                                                                       
    ## speed.sc                                  0.102  0.048 -0.033 -0.595 -0.025 -0.031 -0.141                                                                                                                                                                                                                                                                                
    ## consumerfrac.sc                           0.356 -0.340  0.181  0.004 -0.048  0.072  0.424 -0.654                                                                                                                                                                                                                                                                         
    ## nspp.sc                                  -0.064  0.025 -0.016 -0.061 -0.118 -0.161 -0.012  0.032  0.000                                                                                                                                                                                                                                                                  
    ## npp.sc                                   -0.044  0.007  0.118 -0.045  0.204  0.149  0.011  0.013 -0.007 -0.064                                                                                                                                                                                                                                                           
    ## veg.sc                                   -0.213  0.262 -0.078 -0.082 -0.069 -0.106  0.001  0.198 -0.160 -0.038 -0.553                                                                                                                                                                                                                                                    
    ## duration.sc                              -0.098  0.042 -0.186  0.067  0.069 -0.016 -0.013 -0.051  0.138 -0.012 -0.020 -0.022                                                                                                                                                                                                                                             
    ## human_bowler.sc                          -0.229  0.250  0.042  0.057 -0.059  0.150 -0.120  0.007 -0.147  0.019 -0.190  0.287 -0.127                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:tsign1                   0.169 -0.306 -0.490  0.093 -0.109 -0.195  0.015  0.050 -0.111 -0.013 -0.123  0.012  0.122 -0.120                                                                                                                                                                                                                               
    ## temptrend_abs.sc:tempave_metab.sc         0.296 -0.503  0.138 -0.662 -0.039 -0.003  0.002  0.427 -0.032  0.065  0.034  0.059 -0.090 -0.025 -0.077                                                                                                                                                                                                                        
    ## temptrend_abs.sc:seas.sc                  0.219 -0.380  0.028 -0.025 -0.748 -0.378 -0.012  0.008  0.045  0.073 -0.106  0.012 -0.028 -0.077  0.146  0.033                                                                                                                                                                                                                 
    ## temptrend_abs.sc:microclim.sc             0.174 -0.290 -0.097  0.007 -0.367 -0.845 -0.034  0.038 -0.071  0.090 -0.089  0.054  0.026 -0.139  0.232  0.001             0.453                                                                                                                                                                                               
    ## temptrend_abs.sc:mass.sc                 -0.313  0.361  0.005 -0.015 -0.041 -0.038 -0.642  0.280 -0.433  0.018 -0.026  0.020  0.031  0.088  0.053  0.063             0.028             0.039                                                                                                                                                                             
    ## temptrend_abs.sc:speed.sc                 0.039  0.012  0.007  0.427  0.048  0.058  0.213 -0.732  0.560 -0.016 -0.007 -0.133  0.034  0.004 -0.073 -0.598            -0.029            -0.063            -0.439                                                                                                                                                           
    ## temptrend_abs.sc:consumerfrac.sc         -0.283  0.407 -0.096 -0.047 -0.009 -0.066 -0.336  0.552 -0.784  0.017  0.007  0.087 -0.118  0.068  0.101  0.017            -0.017             0.081             0.549            -0.699                                                                                                                                         
    ## temptrend_abs.sc:nspp.sc                  0.025 -0.086  0.006  0.050  0.053  0.060 -0.005  0.003 -0.014 -0.603  0.011  0.028 -0.007  0.055  0.019 -0.015            -0.055            -0.058             0.002            -0.032             0.012                                                                                                                       
    ## tsign-1:thermal_bias.sc                  -0.091 -0.062  0.667 -0.179 -0.286 -0.090  0.013 -0.024  0.140  0.004  0.115 -0.079 -0.062  0.055 -0.214  0.119             0.198             0.050             0.011            -0.014            -0.061            0.022                                                                                                      
    ## tsign1:thermal_bias.sc                    0.212 -0.208 -0.065 -0.033 -0.666 -0.403  0.049  0.039 -0.026 -0.028  0.143 -0.091 -0.017  0.084  0.301  0.012             0.496             0.298             0.021            -0.027             0.039            0.049             0.281                                                                                    
    ## temptrend_abs.sc:npp.sc                   0.027  0.023 -0.116  0.050 -0.131 -0.072 -0.023 -0.018  0.004  0.003 -0.880  0.494  0.021  0.147  0.149 -0.046             0.078             0.074             0.011             0.033            -0.034           -0.024            -0.123 -0.152                                                                             
    ## temptrend_abs.sc:veg.sc                   0.225 -0.308  0.055  0.044  0.027  0.048  0.031 -0.133  0.116  0.013  0.504 -0.847  0.026 -0.301 -0.005 -0.071             0.016            -0.034            -0.028             0.147            -0.097           -0.050             0.076  0.111 -0.609                                                                      
    ## temptrend_abs.sc:duration.sc             -0.064  0.147  0.020 -0.086 -0.033  0.033  0.001  0.020 -0.059  0.152  0.001  0.010 -0.250  0.105 -0.223  0.094             0.012            -0.038            -0.011            -0.015             0.043           -0.259             0.026 -0.064 -0.005             0.004                                                    
    ## temptrend_abs.sc:human_bowler.sc          0.242 -0.331 -0.037 -0.055 -0.037 -0.159  0.096  0.004  0.106  0.032  0.135 -0.334  0.116 -0.798  0.151  0.051             0.066             0.201            -0.080            -0.030            -0.066           -0.053            -0.040 -0.044 -0.180             0.418           -0.105                                   
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc  0.049 -0.049 -0.379  0.107  0.218  0.025 -0.004  0.036 -0.095 -0.015 -0.268  0.129  0.041 -0.116  0.635 -0.107            -0.305            -0.022             0.054            -0.045             0.083           -0.009            -0.457 -0.267  0.311            -0.156           -0.092            0.169                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc  -0.140  0.284  0.118 -0.001  0.450  0.297 -0.010 -0.009  0.032  0.026 -0.186  0.122 -0.028 -0.006 -0.369 -0.014            -0.615            -0.339            -0.018             0.016            -0.024           -0.077            -0.196 -0.831  0.198            -0.155            0.085            0.038             0.319
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.32506581 -0.38770014 -0.04491269  0.36618931  5.19217938 
    ## 
    ## Number of Observations: 2299
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                     90                   2299

``` r
summary(modTfullHornFresh)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -972.9369 -824.5939 520.4684
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev       Corr  
    ## (Intercept)      3.088553e-08 (Intr)
    ## temptrend_abs.sc 2.055026e-02 0.025 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.02085277  2.44783
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.345198 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * tsign + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc 
    ##                                                Value  Std.Error  DF   t-value p-value
    ## (Intercept)                               0.01342964 0.02774716 563  0.484000  0.6286
    ## temptrend_abs.sc                         -0.08728721 0.07498855 563 -1.164007  0.2449
    ## tsign1                                   -0.00129613 0.00936826 563 -0.138353  0.8900
    ## tempave_metab.sc                          0.00229079 0.02341751 563  0.097824  0.9221
    ## seas.sc                                   0.00887804 0.00676282 563  1.312771  0.1898
    ## microclim.sc                             -0.00054835 0.00518003 563 -0.105859  0.9157
    ## mass.sc                                  -0.00321190 0.00502149 563 -0.639630  0.5227
    ## speed.sc                                  0.00441196 0.00662144 563  0.666314  0.5055
    ## nspp.sc                                  -0.01612930 0.00473283 563 -3.407963  0.0007
    ## npp.sc                                    0.01007588 0.00640617 563  1.572840  0.1163
    ## veg.sc                                   -0.00617104 0.00587586 563 -1.050236  0.2941
    ## duration.sc                               0.00417228 0.00497137 563  0.839263  0.4017
    ## human_bowler.sc                           0.00593081 0.00586290 563  1.011583  0.3122
    ## temptrend_abs.sc:tsign1                   0.01502038 0.02346993 563  0.639984  0.5224
    ## temptrend_abs.sc:tempave_metab.sc         0.01023082 0.06251469 563  0.163655  0.8701
    ## temptrend_abs.sc:seas.sc                 -0.02056624 0.02029530 563 -1.013350  0.3113
    ## temptrend_abs.sc:microclim.sc            -0.00172573 0.01617519 563 -0.106690  0.9151
    ## temptrend_abs.sc:mass.sc                  0.00811168 0.01574522 563  0.515184  0.6066
    ## temptrend_abs.sc:speed.sc                -0.01972287 0.02063490 563 -0.955802  0.3396
    ## temptrend_abs.sc:nspp.sc                  0.04782508 0.01402193 563  3.410735  0.0007
    ## tsign-1:thermal_bias.sc                   0.01601591 0.01541784 563  1.038790  0.2993
    ## tsign1:thermal_bias.sc                    0.01482133 0.00977478 563  1.516282  0.1300
    ## temptrend_abs.sc:npp.sc                  -0.02032586 0.01914174 563 -1.061860  0.2888
    ## temptrend_abs.sc:veg.sc                   0.03337712 0.02431163 563  1.372887  0.1703
    ## temptrend_abs.sc:duration.sc             -0.01968358 0.00979032 563 -2.010514  0.0449
    ## temptrend_abs.sc:human_bowler.sc         -0.00875220 0.01936998 563 -0.451844  0.6516
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.04534739 0.04437608 563 -1.021888  0.3073
    ## temptrend_abs.sc:tsign1:thermal_bias.sc  -0.02252824 0.03273324 563 -0.688238  0.4916
    ##  Correlation: 
    ##                                          (Intr) tmpt_. tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc nspp.s npp.sc veg.sc drtn.s hmn_b. tm_.:1 tmptrnd_bs.sc:t_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. tmptrnd_bs.sc:h_. t_.:-1
    ## temptrend_abs.sc                         -0.602                                                                                                                                                                                                                                                                                                  
    ## tsign1                                   -0.231  0.096                                                                                                                                                                                                                                                                                           
    ## tempave_metab.sc                          0.484 -0.185  0.121                                                                                                                                                                                                                                                                                    
    ## seas.sc                                   0.008 -0.055  0.022  0.354                                                                                                                                                                                                                                                                             
    ## microclim.sc                             -0.132  0.056  0.107  0.381  0.099                                                                                                                                                                                                                                                                      
    ## mass.sc                                  -0.034  0.042 -0.073 -0.301  0.196 -0.393                                                                                                                                                                                                                                                               
    ## speed.sc                                  0.241 -0.155  0.004  0.232  0.111  0.078 -0.670                                                                                                                                                                                                                                                        
    ## nspp.sc                                   0.125 -0.007 -0.163 -0.343 -0.442 -0.166  0.205 -0.078                                                                                                                                                                                                                                                 
    ## npp.sc                                   -0.272 -0.025  0.019 -0.132  0.452  0.147  0.096 -0.007 -0.032                                                                                                                                                                                                                                          
    ## veg.sc                                   -0.543  0.574 -0.027 -0.010 -0.068 -0.128 -0.108  0.065 -0.146 -0.358                                                                                                                                                                                                                                   
    ## duration.sc                              -0.134  0.121 -0.187  0.102 -0.026  0.072 -0.112  0.023 -0.208  0.041  0.119                                                                                                                                                                                                                            
    ## human_bowler.sc                           0.065  0.062 -0.078 -0.106 -0.267 -0.011  0.165 -0.105 -0.048 -0.123  0.057  0.041                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tsign1                   0.171 -0.362 -0.534 -0.016 -0.013 -0.007 -0.046  0.047  0.081 -0.018 -0.047  0.054 -0.032                                                                                                                                                                                                              
    ## temptrend_abs.sc:tempave_metab.sc        -0.159  0.098 -0.029 -0.614 -0.427 -0.311  0.204 -0.189  0.196 -0.151 -0.024 -0.045  0.298 -0.040                                                                                                                                                                                                       
    ## temptrend_abs.sc:seas.sc                 -0.021  0.063  0.027 -0.352 -0.697 -0.152 -0.026 -0.180  0.187 -0.353  0.010  0.000  0.322 -0.050  0.779                                                                                                                                                                                                
    ## temptrend_abs.sc:microclim.sc             0.039  0.000  0.013 -0.258 -0.174 -0.622  0.213 -0.074  0.034 -0.108  0.098 -0.008  0.118 -0.070  0.524             0.397                                                                                                                                                                              
    ## temptrend_abs.sc:mass.sc                 -0.016 -0.069 -0.010  0.150 -0.015  0.239 -0.674  0.446 -0.105  0.021  0.040  0.076 -0.298  0.097 -0.376            -0.201            -0.358                                                                                                                                                            
    ## temptrend_abs.sc:speed.sc                -0.130  0.244  0.058 -0.123 -0.179 -0.049  0.420 -0.644  0.022 -0.042 -0.028  0.001  0.200 -0.093  0.226             0.319             0.166            -0.674                                                                                                                                          
    ## temptrend_abs.sc:nspp.sc                 -0.072 -0.009  0.157  0.177  0.134  0.083 -0.160 -0.005 -0.624 -0.031  0.046  0.116 -0.059 -0.111 -0.251            -0.150            -0.005             0.180             0.016                                                                                                                        
    ## tsign-1:thermal_bias.sc                   0.429 -0.242 -0.411  0.253  0.043  0.065 -0.035  0.048  0.091  0.126 -0.256  0.003  0.055  0.244 -0.183            -0.117            -0.061             0.052            -0.035            -0.028                                                                                                      
    ## tsign1:thermal_bias.sc                    0.218 -0.069  0.138  0.581  0.367  0.124  0.040  0.015 -0.204  0.124  0.060  0.027  0.188 -0.140 -0.399            -0.320            -0.052            -0.104             0.034             0.041             0.196                                                                                    
    ## temptrend_abs.sc:npp.sc                   0.014  0.271  0.009 -0.067 -0.351 -0.074 -0.028 -0.023  0.022 -0.690  0.340 -0.066  0.117  0.033  0.341             0.459             0.196            -0.075             0.013             0.050            -0.173 -0.199                                                                             
    ## temptrend_abs.sc:veg.sc                   0.354 -0.817  0.020 -0.011  0.016  0.060  0.023 -0.036  0.006  0.244 -0.692 -0.075 -0.037  0.106  0.071             0.032            -0.138            -0.023             0.016            -0.022             0.144 -0.042 -0.513                                                                      
    ## temptrend_abs.sc:duration.sc              0.043 -0.011 -0.037  0.044  0.153  0.060  0.045  0.061  0.117  0.044 -0.090 -0.342 -0.061  0.021 -0.044            -0.210            -0.189            -0.003            -0.108            -0.279             0.071  0.059 -0.047             0.116                                                    
    ## temptrend_abs.sc:human_bowler.sc         -0.018 -0.038  0.025  0.193  0.275  0.108 -0.232  0.190  0.007  0.116  0.005 -0.030 -0.709  0.056 -0.499            -0.501            -0.182             0.480            -0.329             0.004            -0.006 -0.029 -0.145            -0.075            0.053                                   
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.223  0.440  0.187 -0.198 -0.170 -0.066  0.098 -0.065  0.040 -0.246  0.236 -0.021  0.087 -0.477  0.423             0.374             0.182            -0.228             0.072            -0.085            -0.534 -0.177  0.494            -0.339           -0.057           -0.153                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc  -0.099  0.022 -0.106 -0.344 -0.342 -0.044 -0.075 -0.026  0.075 -0.216 -0.030  0.076 -0.004  0.173  0.585             0.520             0.087             0.003            -0.045            -0.045            -0.159 -0.711  0.411             0.068           -0.056           -0.145             0.337
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -4.25987257 -0.29069598 -0.02883191  0.25023757  5.23849247 
    ## 
    ## Number of Observations: 608
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                     18                    608

``` r
summary(modTfullHornMar)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -89811.54 -89526.43 44939.77
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.02167793 (Intr)
    ## temptrend_abs.sc 0.02656924 0.035 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.02055808 2.782391
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.432839 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * tsign + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc 
    ##                                                 Value   Std.Error    DF   t-value p-value
    ## (Intercept)                               0.009417700 0.003221245 32301  2.923622  0.0035
    ## temptrend_abs.sc                          0.013723984 0.005278704 32301  2.599878  0.0093
    ## tsign1                                   -0.002771867 0.000847672 32301 -3.269977  0.0011
    ## tempave_metab.sc                         -0.011144629 0.001466704 32301 -7.598416  0.0000
    ## seas.sc                                  -0.004877220 0.001167078 32301 -4.179002  0.0000
    ## microclim.sc                              0.000168668 0.000458639 32301  0.367758  0.7131
    ## mass.sc                                  -0.001088684 0.000744248 32301 -1.462798  0.1435
    ## speed.sc                                  0.000080710 0.000934878 32301  0.086333  0.9312
    ## consumerfrac.sc                           0.000591355 0.000443181 32301  1.334341  0.1821
    ## nspp.sc                                   0.000499063 0.000787859 32301  0.633442  0.5264
    ## npp.sc                                   -0.000922892 0.000674078 32301 -1.369118  0.1710
    ## duration.sc                              -0.002215334 0.000694507 32301 -3.189791  0.0014
    ## human_bowler.sc                           0.000295339 0.000511846 32301  0.577007  0.5639
    ## temptrend_abs.sc:tsign1                  -0.003528588 0.002184310 32301 -1.615425  0.1062
    ## temptrend_abs.sc:tempave_metab.sc         0.005423715 0.003717238 32301  1.459071  0.1446
    ## temptrend_abs.sc:seas.sc                  0.002887307 0.002581617 32301  1.118410  0.2634
    ## temptrend_abs.sc:microclim.sc            -0.000714943 0.001251111 32301 -0.571446  0.5677
    ## temptrend_abs.sc:mass.sc                  0.000152209 0.001668004 32301  0.091252  0.9273
    ## temptrend_abs.sc:speed.sc                 0.002394258 0.002266011 32301  1.056596  0.2907
    ## temptrend_abs.sc:consumerfrac.sc          0.005764073 0.000900718 32301  6.399419  0.0000
    ## temptrend_abs.sc:nspp.sc                  0.002102625 0.001930497 32301  1.089162  0.2761
    ## tsign-1:thermal_bias.sc                  -0.001834972 0.000867779 32301 -2.114563  0.0345
    ## tsign1:thermal_bias.sc                   -0.001864918 0.000692132 32301 -2.694455  0.0071
    ## temptrend_abs.sc:npp.sc                   0.006761866 0.001833005 32301  3.688952  0.0002
    ## temptrend_abs.sc:duration.sc              0.000693145 0.001279466 32301  0.541746  0.5880
    ## temptrend_abs.sc:human_bowler.sc         -0.003061024 0.001265012 32301 -2.419760  0.0155
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.001725531 0.001785383 32301 -0.966476  0.3338
    ## temptrend_abs.sc:tsign1:thermal_bias.sc   0.004528860 0.001490020 32301  3.039462  0.0024
    ##  Correlation: 
    ##                                          (Intr) tmpt_. tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc drtn.s hmn_b. tm_.:1 tmptrnd_bs.sc:t_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:d. tmptrnd_bs.sc:h_. t_.:-1
    ## temptrend_abs.sc                         -0.263                                                                                                                                                                                                                                                                                                  
    ## tsign1                                   -0.188  0.115                                                                                                                                                                                                                                                                                           
    ## tempave_metab.sc                          0.102 -0.003  0.022                                                                                                                                                                                                                                                                                    
    ## seas.sc                                   0.185 -0.084 -0.091  0.221                                                                                                                                                                                                                                                                             
    ## microclim.sc                              0.037 -0.045  0.008 -0.204  0.067                                                                                                                                                                                                                                                                      
    ## mass.sc                                   0.058 -0.033  0.000  0.062  0.037  0.012                                                                                                                                                                                                                                                               
    ## speed.sc                                  0.153 -0.080 -0.079  0.072  0.073  0.039 -0.434                                                                                                                                                                                                                                                        
    ## consumerfrac.sc                           0.000  0.012  0.029 -0.068 -0.012  0.012 -0.055  0.000                                                                                                                                                                                                                                                 
    ## nspp.sc                                  -0.131  0.096  0.058 -0.114  0.092 -0.084 -0.105  0.204  0.112                                                                                                                                                                                                                                          
    ## npp.sc                                   -0.043  0.057 -0.002 -0.053 -0.345 -0.238 -0.024  0.081 -0.040 -0.237                                                                                                                                                                                                                                   
    ## duration.sc                               0.041 -0.022 -0.137  0.053 -0.064 -0.016 -0.001 -0.016 -0.016 -0.279  0.065                                                                                                                                                                                                                            
    ## human_bowler.sc                          -0.030  0.030  0.005  0.098 -0.193  0.005 -0.026  0.040 -0.055 -0.051 -0.158  0.018                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tsign1                   0.089 -0.249 -0.520 -0.048 -0.009 -0.003  0.008  0.032 -0.017 -0.036 -0.002  0.018 -0.006                                                                                                                                                                                                              
    ## temptrend_abs.sc:tempave_metab.sc        -0.031  0.078 -0.081 -0.437 -0.161  0.032 -0.019  0.058  0.040  0.038  0.018  0.022 -0.044  0.092                                                                                                                                                                                                       
    ## temptrend_abs.sc:seas.sc                 -0.078  0.140  0.048 -0.144 -0.689 -0.051 -0.018 -0.038 -0.016 -0.057  0.261  0.016  0.155  0.013  0.196                                                                                                                                                                                                
    ## temptrend_abs.sc:microclim.sc            -0.022  0.066 -0.016  0.131 -0.010 -0.723  0.014 -0.028 -0.003  0.052  0.180  0.004 -0.028 -0.019  0.029             0.053                                                                                                                                                                              
    ## temptrend_abs.sc:mass.sc                  0.002  0.093  0.012 -0.038 -0.016  0.026 -0.573  0.252  0.060  0.028 -0.019  0.033  0.031 -0.013  0.041             0.005            -0.063                                                                                                                                                            
    ## temptrend_abs.sc:speed.sc                -0.031  0.176  0.047  0.026 -0.018 -0.039  0.237 -0.575  0.016 -0.167 -0.017  0.027 -0.043 -0.019 -0.156             0.067             0.051            -0.412                                                                                                                                          
    ## temptrend_abs.sc:consumerfrac.sc         -0.013 -0.027 -0.018  0.051  0.006  0.017  0.052  0.007 -0.719 -0.087  0.016 -0.004  0.075  0.046 -0.049            -0.029            -0.039            -0.112            -0.042                                                                                                                        
    ## temptrend_abs.sc:nspp.sc                  0.038 -0.174 -0.024  0.031 -0.065  0.046  0.026 -0.170 -0.064 -0.643  0.144  0.179  0.044  0.032 -0.064             0.061            -0.014            -0.009             0.212             0.115                                                                                                      
    ## tsign-1:thermal_bias.sc                   0.037  0.000 -0.166  0.241 -0.140 -0.031 -0.046  0.047 -0.004 -0.055 -0.011  0.077  0.141  0.042 -0.054             0.091             0.011             0.023            -0.015            -0.001            0.009                                                                                     
    ## tsign1:thermal_bias.sc                    0.006 -0.004  0.129  0.432 -0.171 -0.128 -0.002  0.027 -0.005 -0.059 -0.118  0.093  0.170 -0.076 -0.125             0.065             0.049             0.017            -0.004             0.002            0.000             0.358                                                                   
    ## temptrend_abs.sc:npp.sc                   0.007 -0.082  0.014  0.031  0.208  0.182 -0.009 -0.036  0.014  0.133 -0.656 -0.036  0.097 -0.050 -0.027            -0.379            -0.338             0.065            -0.012            -0.003           -0.152            -0.005  0.086                                                            
    ## temptrend_abs.sc:duration.sc             -0.104  0.249  0.028  0.111  0.011 -0.038  0.036  0.016  0.030  0.253 -0.038 -0.484  0.034 -0.079 -0.042            -0.046             0.044            -0.068            -0.044             0.024           -0.356            -0.028 -0.055  0.054                                                     
    ## temptrend_abs.sc:human_bowler.sc          0.029 -0.029 -0.017 -0.067  0.145 -0.007  0.023 -0.033  0.133  0.039  0.105 -0.005 -0.735  0.018  0.089            -0.171            -0.022            -0.025             0.038            -0.153           -0.057            -0.092 -0.124 -0.140            -0.020                                   
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.012  0.010  0.028 -0.091  0.105  0.023  0.013  0.009 -0.003  0.008 -0.003 -0.030 -0.083  0.059  0.277            -0.127            -0.034            -0.014            -0.047             0.029           -0.009            -0.482 -0.174  0.011             0.066            0.087                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc   0.006 -0.004 -0.104 -0.174  0.072  0.038  0.012  0.042 -0.019  0.015  0.082 -0.046 -0.122  0.164  0.469            -0.069            -0.123            -0.016            -0.124             0.050           -0.033            -0.177 -0.571 -0.084             0.037            0.137             0.341
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.59847381 -0.21790065 -0.02487884  0.21515504  5.76491820 
    ## 
    ## Number of Observations: 32420
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                     92                  32420

### Plot the realm-specific coefficients

Also uses the full models across all realms

``` r
coefs1 <- summary(modTfullHornrem0)$tTable
coefs2 <- summary(modTfullHornTerr)$tTable
coefs3 <- summary(modTfullHornFresh)$tTable
coefs4 <- summary(modTfullHornMar)$tTable

varstoplot <- unique(c(rownames(coefs1), rownames(coefs2), rownames(coefs3), rownames(coefs4)))
varstoplot <- varstoplot[which(!grepl('Intercept', varstoplot) | grepl(':', varstoplot))] # vars to plot

rows1_1 <- which(rownames(coefs1) %in% varstoplot) # rows in coefs
rows1_2 <- which(rownames(coefs2) %in% varstoplot)
rows1_3 <- which(rownames(coefs3) %in% varstoplot)
rows1_4 <- which(rownames(coefs4) %in% varstoplot)
xlims <- range(c(coefs1[rows1_1,1] - coefs1[rows1_1,2], coefs1[rows1_1,1] + coefs1[rows1_1,2], 
                 coefs2[rows1_2,1] - coefs2[rows1_2,2], coefs2[rows1_2,1] + coefs2[rows1_2,2], 
                 coefs3[rows1_3,1] - coefs3[rows1_3,2], coefs3[rows1_3,1] + coefs3[rows1_3,2],
                 coefs4[rows1_4,1] - coefs4[rows1_4,2], coefs4[rows1_4,1] + coefs4[rows1_4,2]))


cols <- brewer.pal(4, 'Dark2') # for full, terr, fresh, mar
pchs <- c(1, 16, 16, 16)
offs <- c(0.1, 0, -0.1, -0.2) # offset vertically for each model


par(las = 1, mai = c(0.5, 4, 0.1, 0.1))

plot(0,0, col = 'white', xlim = xlims, ylim = c(1,length(varstoplot)), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(varstoplot):1, labels = varstoplot, cex.axis = 0.7)
abline(v = 0, col = 'grey', lty = 2)
abline(h = 1:length(varstoplot), col = 'grey', lty = 3)
for(i in 1:length(varstoplot)){
  if(varstoplot[i] %in% rownames(coefs1)){
    x = coefs1[rownames(coefs1) == varstoplot[i], 1]
    se = coefs1[rownames(coefs1) == varstoplot[i], 2]
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
}
legend('bottomleft', col = cols, pch = pchs, lwd = 1, legend = c('All', 'Terestrial', 'Freshwater', 'Marine'))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20realm%20mods-1.png)<!-- -->

\[End text in hopes this helps the last figure show up when knitted\]

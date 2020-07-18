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
library(beanplot) # for beanplots
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

# realm that combined Terrestrial and Freshwater, for interacting with human impact
trends[, REALM2 := REALM]
levels(trends$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")

# group Marine invertebrates/plants in with All
trends[, taxa_mod2 := taxa_mod]
trends[taxa_mod == 'Marine invertebrates/plants', taxa_mod2 := 'All']
```

Log-transform some variables, then center and scale.

``` r
trends[, tempave.sc := scale(tempave)]
trends[, tempave_metab.sc := scale(tempave_metab)]
trends[, seas.sc := scale(seas)]
trends[, microclim.sc := scale(log(microclim))]
trends[, temptrend.sc := scale(temptrend)]
trends[, temptrend_abs.sc := scale(abs(temptrend))]
trends[, npp.sc := scale(log(npp))]
trends[, mass.sc := scale(log(mass_mean_weight))]
trends[, speed.sc := scale(log(speed_mean_weight+1))]
trends[, lifespan.sc := scale(log(lifespan_mean_weight))]
trends[, thermal_bias.sc := scale(thermal_bias)]
trends[, consumerfrac.sc := scale(consfrac)]
trends[, endothermfrac.sc := scale(endofrac)]
trends[, nspp.sc := scale(log(Nspp))]
trends[, human.sc := scale(human)]
```

### Do the variables look ok?

``` r
# histograms to examine
cexmain = 0.6
par(mfrow = c(3,5))
invisible(trends[, hist(tempave.sc, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(trends[, hist(tempave_metab.sc, main = 'Metabolic temperature (°C)', cex.main = cexmain)])
invisible(trends[, hist(seas.sc, main = 'Seasonality (°C)', cex.main = cexmain)])
invisible(trends[, hist(microclim.sc, main = 'log Microclimates (°C)', cex.main = cexmain)])
invisible(trends[, hist(temptrend.sc, main = 'Temperature trend (°C/yr)', cex.main = cexmain)])
invisible(trends[, hist(temptrend_abs.sc, main = 'log abs(Temperature trend) (°C/yr)', cex.main = cexmain)])
invisible(trends[, hist(mass.sc, main = 'log Mass (g)', cex.main = cexmain)])
invisible(trends[, hist(speed.sc, main = 'log Speed (km/hr)', cex.main = cexmain)])
invisible(trends[, hist(lifespan.sc, main = 'log Lifespan (yr)', cex.main = cexmain)])
invisible(trends[, hist(consumerfrac.sc, main = 'Consumers (fraction)', cex.main = cexmain)])
invisible(trends[, hist(endothermfrac.sc, main = 'Endotherms (fraction)', cex.main = cexmain)])
invisible(trends[, hist(nspp.sc, main = 'log Species richness', cex.main = cexmain)])
invisible(trends[, hist(thermal_bias.sc, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(trends[, hist(npp.sc, main = 'log Net primary productivity', cex.main = cexmain)])
invisible(trends[, hist(human.sc, main = 'Human impact score', cex.main = cexmain)])
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/histograms-1.png)<!-- -->

### Check correlations among variables. Pearson’s r is on the lower diagonal.

``` r
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = 'pairwise.complete.obs')
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt) #, cex = cex.cor * r)
}
pairs(formula = ~ REALM + tempave.sc + tempave_metab.sc + seas.sc + microclim.sc + temptrend.sc + temptrend_abs.sc +  mass.sc + speed.sc + lifespan.sc + consumerfrac.sc + endothermfrac.sc + nspp.sc + thermal_bias.sc + npp.sc + human.sc, data = trends, gap = 1/10, cex = 0.2, col = '#00000022', lower.panel = panel.cor)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/pairs-1.png)<!-- -->

Mass and lifespan look tightly correlated, but r only 0.56…?
Tempave\_metab and lifespan don’t look tightly correlated, but r= -0.81
Tempave\_metab and speed don’t look tightly correlated, but r= -0.83
Lifespan and speed don’t look tightly correlated, but r = 0.73

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

With all covariates

``` r
# the cases we can compare
apply(trends[, .(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, human.sc)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##         Jtutrend            REALM       tempave.sc tempave_metab.sc          seas.sc     microclim.sc     temptrend.sc          mass.sc         speed.sc      lifespan.sc  consumerfrac.sc endothermfrac.sc          nspp.sc  thermal_bias.sc           npp.sc         human.sc 
    ##            53013            53013            49916            49916            49916            51834            49916            52820            52689            51540            47534            53013            53013            49371            52863            53013

``` r
i <- trends[, complete.cases(Jtutrend, temptrend.sc, tempave_metab.sc, REALM, seas.sc, microclim.sc, npp.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, thermal_bias.sc)]
cat('Overall # time-series: ', sum(i), '\n')
```

    ## Overall # time-series:  43585

``` r
cat('# studies: ', trends[i, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  250

``` r
cat('Data points: ', trends[i, sum(nyrBT)], '\n')
```

    ## Data points:  222824

``` r
trends[i, table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##        1008       39735        2842

``` r
trends[i, table(taxa_mod)]
```

    ## taxa_mod
    ##           All    Amphibians       Benthos         Birds          Fish Invertebrates       Mammals         Plant      Reptiles 
    ##           521            12           590         11803         27372          2567           518           200             2

``` r
trends[i, table(taxa_mod, REALM)]
```

    ##                REALM
    ## taxa_mod        Freshwater Marine Terrestrial
    ##   All                    0    520           1
    ##   Amphibians             2      0          10
    ##   Benthos                0    590           0
    ##   Birds                  0   9221        2582
    ##   Fish                 993  26379           0
    ##   Invertebrates         12   2484          71
    ##   Mammals                0    477          41
    ##   Plant                  1     64         135
    ##   Reptiles               0      0           2

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

## Plot turnover vs. explanatory variables

Lines are ggplot smoother fits by realm.
![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-1.png)<!-- -->

Strong trends with temperature change, but trends are pretty symmetric
around no trend in temperature, which implies warming or cooling drives
similar degree of community turnover. Some indication of less turnover
for larger organisms (mass) Higher turnover on land with higher
seasonality? More turnover for shorter-lived organisms? No really clear
differences among realms.

Average rates of
turnover

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

## Compare covariates across realms

``` r
i <- trends[, !duplicated(rarefyID)]; sum(i)
```

    ## [1] 53013

``` r
par(mfrow=c(5,3))
beanplot(rarefyID_y ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Latitude (degN)', ll = 0.05)
beanplot(tempave ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature (degC)', ll = 0.05)
beanplot(tempave_metab ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Metabolic Temperature (degC)', ll = 0.05, bw = 'nrd0') # nrd0 bandwidth to calculation gap
beanplot(seas ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Seasonality (degC)', ll = 0.05)
beanplot(microclim ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Microclimates (degC)', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(temptrend ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature trend (degC/yr)', ll = 0.05)
beanplot(mass_mean_weight ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Mass (g)', ll = 0.05, log = 'y')
beanplot(speed_mean_weight +1 ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Speed (km/hr)', ll = 0.05, log = 'y')
beanplot(lifespan_mean_weight ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Lifespan (yr)', ll = 0.05, log = 'y')
#beanplot(consfrac ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Consumers (fraction)', ll = 0.05, log = '') # too sparse
#beanplot(endofrac ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Endotherms (fraction)', ll = 0.05, log = '') # too sparse
beanplot(Nspp ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Number of species', ll = 0.05, log = 'y')
beanplot(thermal_bias ~ REALM, data = trends[i & !is.na(thermal_bias),], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Thermal bias (degC)', ll = 0.05)
beanplot(npp ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'NPP', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(human ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Human impact score', ll = 0.05)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->

Marine are in generally warmer locations (seawater doesn’t freeze)
Marine have much lower seasonality. Marine and freshwater have some very
small masses (plankton), but much of dataset is similar to terrestrial.
Marine has a lot of slow, crawling organisms, but land has plants. Land
also has birds (fast).

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
Scatterplot, violin plots, and coefficient plots all
together

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
  - environmental temperature
  - average metabolic temperature
  - seasonality
  - microclimates
  - NPP
  - speed
  - mass
  - lifespan
  - consumer vs. producer
  - thermal bias

Except for thermal bias: interact with temperature trend (not
abs)

### Full model for Jaccard total

``` r
i <- trends[, complete.cases(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, human.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modTfull1.rds')){
  modTfull1 <- readRDS('temp/modTfull1.rds')
} else {
  modTfull1 <- lme(Jtutrend ~ temptrend_abs.sc*REALM + 
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
                     temptrend.sc*thermal_bias.sc +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*human.sc:REALM2,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modTfull1, file = 'temp/modTfull1.rds')
}

summary(modTfull1)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC  logLik
    ##   -134191.8 -133835.8 67136.9
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.04834247 (Intr)
    ## temptrend_abs.sc 0.03064449 0.335 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 0.0008908894 0.2861661
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.235108 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tempave.sc +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      lifespan.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      endothermfrac.sc + temptrend_abs.sc * nspp.sc + temptrend.sc *      thermal_bias.sc + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      human.sc:REALM2 
    ##                                                 Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                0.03485684 0.012488014 43303   2.791224  0.0053
    ## temptrend_abs.sc                           0.06059353 0.013950373 43303   4.343506  0.0000
    ## REALMMarine                                0.04353310 0.013342238   247   3.262803  0.0013
    ## REALMTerrestrial                           0.01283571 0.014736639   247   0.871007  0.3846
    ## tempave.sc                                 0.00132970 0.000771340 43303   1.723878  0.0847
    ## tempave_metab.sc                           0.00623790 0.002070164 43303   3.013241  0.0026
    ## seas.sc                                    0.00106939 0.000438407 43303   2.439258  0.0147
    ## microclim.sc                              -0.00108978 0.000256948 43303  -4.241269  0.0000
    ## mass.sc                                   -0.00266149 0.000859922 43303  -3.095035  0.0020
    ## speed.sc                                  -0.00010731 0.000521053 43303  -0.205952  0.8368
    ## lifespan.sc                               -0.00006994 0.001723658 43303  -0.040578  0.9676
    ## consumerfrac.sc                            0.00081968 0.000850291 43303   0.964002  0.3351
    ## endothermfrac.sc                          -0.01009966 0.004056904 43303  -2.489499  0.0128
    ## nspp.sc                                   -0.00898207 0.000466114 43303 -19.270110  0.0000
    ## temptrend.sc                              -0.00038927 0.000462965 43303  -0.840817  0.4005
    ## thermal_bias.sc                           -0.00032265 0.000442144 43303  -0.729735  0.4656
    ## npp.sc                                    -0.00102418 0.000394709 43303  -2.594769  0.0095
    ## temptrend_abs.sc:REALMMarine              -0.01603037 0.014395697 43303  -1.113553  0.2655
    ## temptrend_abs.sc:REALMTerrestrial         -0.03673714 0.017428303 43303  -2.107901  0.0350
    ## temptrend_abs.sc:tempave.sc                0.00874920 0.001170711 43303   7.473408  0.0000
    ## temptrend_abs.sc:tempave_metab.sc          0.00182583 0.003932507 43303   0.464291  0.6424
    ## temptrend_abs.sc:seas.sc                   0.00369032 0.000831083 43303   4.440373  0.0000
    ## temptrend_abs.sc:microclim.sc             -0.00248228 0.000541074 43303  -4.587685  0.0000
    ## temptrend_abs.sc:mass.sc                  -0.00501285 0.001517754 43303  -3.302811  0.0010
    ## temptrend_abs.sc:speed.sc                  0.00168939 0.000762144 43303   2.216621  0.0267
    ## temptrend_abs.sc:lifespan.sc               0.00815228 0.003096863 43303   2.632433  0.0085
    ## temptrend_abs.sc:consumerfrac.sc          -0.00150838 0.001128378 43303  -1.336769  0.1813
    ## temptrend_abs.sc:endothermfrac.sc          0.01087517 0.005149856 43303   2.111742  0.0347
    ## temptrend_abs.sc:nspp.sc                  -0.00472719 0.000790837 43303  -5.977454  0.0000
    ## temptrend.sc:thermal_bias.sc              -0.00056051 0.000335586 43303  -1.670231  0.0949
    ## temptrend_abs.sc:npp.sc                   -0.00299048 0.000800678 43303  -3.734933  0.0002
    ## human.sc:REALM2TerrFresh                  -0.00062625 0.000274700 43303  -2.279776  0.0226
    ## human.sc:REALM2Marine                     -0.00024969 0.000353901 43303  -0.705526  0.4805
    ## temptrend_abs.sc:human.sc:REALM2TerrFresh -0.00143439 0.000536378 43303  -2.674225  0.0075
    ## temptrend_abs.sc:human.sc:REALM2Marine    -0.00146468 0.000715136 43303  -2.048112  0.0406
    ##  Correlation: 
    ##                                           (Intr) tmpt_. REALMM REALMT tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s tmptr. thrm_. npp.sc t_.:REALMM t_.:REALMT tmptrnd_bs.sc:t. t_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. tm.:_. tmptrnd_bs.sc:np. h.:REALM2T h.:REALM2M t_.:.:REALM2T
    ## temptrend_abs.sc                           0.344                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                               -0.917 -0.314                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                          -0.878 -0.314  0.802                                                                                                                                                                                                                                                                                                                                                           
    ## tempave.sc                                 0.002 -0.008 -0.015  0.002                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                           0.047  0.030 -0.035 -0.050 -0.338                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                   -0.036 -0.019  0.051 -0.020  0.181  0.007                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                              -0.012 -0.019  0.015 -0.016 -0.044  0.133  0.135                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                    0.012 -0.012  0.000  0.003 -0.042 -0.515 -0.052 -0.027                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                   0.019  0.002 -0.032 -0.032 -0.011  0.121 -0.059 -0.015 -0.090                                                                                                                                                                                                                                                                                                                 
    ## lifespan.sc                                0.046  0.047 -0.028 -0.048 -0.023  0.713  0.075  0.031 -0.813  0.190                                                                                                                                                                                                                                                                                                          
    ## consumerfrac.sc                           -0.072 -0.057  0.081  0.342 -0.001 -0.109  0.002  0.004  0.060 -0.164 -0.131                                                                                                                                                                                                                                                                                                   
    ## endothermfrac.sc                           0.172  0.080 -0.091 -0.347  0.134 -0.234  0.004 -0.044  0.011  0.085 -0.024 -0.415                                                                                                                                                                                                                                                                                            
    ## nspp.sc                                   -0.013  0.021 -0.010 -0.023 -0.013  0.032 -0.033 -0.024 -0.169 -0.062  0.138  0.012  0.053                                                                                                                                                                                                                                                                                     
    ## temptrend.sc                              -0.008 -0.015  0.008  0.008 -0.087  0.004 -0.034  0.018  0.017 -0.004 -0.014  0.005 -0.007 -0.019                                                                                                                                                                                                                                                                              
    ## thermal_bias.sc                            0.029 -0.011 -0.045 -0.008  0.671 -0.024 -0.113 -0.112  0.033  0.028 -0.059 -0.002 -0.012 -0.069  0.000                                                                                                                                                                                                                                                                       
    ## npp.sc                                    -0.007 -0.033 -0.001  0.008 -0.221  0.176 -0.229 -0.291  0.004 -0.012  0.044 -0.043 -0.070 -0.119  0.022 -0.085                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:REALMMarine              -0.327 -0.952  0.337  0.296  0.006 -0.021  0.031  0.020  0.014 -0.006 -0.032  0.060 -0.055 -0.026  0.013  0.009  0.019                                                                                                                                                                                                                                                         
    ## temptrend_abs.sc:REALMTerrestrial         -0.297 -0.852  0.269  0.369  0.014 -0.048 -0.051 -0.025  0.015 -0.010 -0.040  0.151 -0.156 -0.041  0.010  0.003  0.023  0.801                                                                                                                                                                                                                                                  
    ## temptrend_abs.sc:tempave.sc                0.001 -0.008 -0.004  0.013  0.489 -0.305  0.006 -0.082 -0.083 -0.040  0.014  0.026  0.151  0.021 -0.085  0.049 -0.066  0.006      0.027                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:tempave_metab.sc          0.014  0.050 -0.008 -0.024 -0.183  0.668  0.078  0.126 -0.299  0.033  0.434 -0.077 -0.191  0.030  0.013  0.056  0.072 -0.038     -0.067     -0.502                                                                                                                                                                                                                            
    ## temptrend_abs.sc:seas.sc                  -0.013 -0.038  0.018 -0.027  0.109  0.079  0.648  0.080 -0.045 -0.018  0.036  0.006 -0.034 -0.015 -0.005  0.079 -0.136  0.054     -0.064      0.137            0.073                                                                                                                                                                                                           
    ## temptrend_abs.sc:microclim.sc             -0.012 -0.030  0.009 -0.012 -0.042  0.131  0.086  0.768 -0.033  0.010  0.013 -0.001 -0.056  0.013  0.029 -0.003 -0.256  0.034     -0.026     -0.060            0.106  0.113                                                                                                                                                                                                    
    ## temptrend_abs.sc:mass.sc                  -0.006 -0.014  0.008  0.009 -0.085 -0.330 -0.051 -0.028  0.605  0.007 -0.528  0.031  0.001 -0.105  0.020 -0.011  0.019  0.020      0.023     -0.117           -0.479 -0.047            -0.018                                                                                                                                                                                  
    ## temptrend_abs.sc:speed.sc                  0.001  0.010 -0.003 -0.010 -0.025  0.015 -0.007  0.007  0.018  0.286  0.005 -0.057  0.039 -0.040 -0.010  0.002  0.004 -0.026     -0.031     -0.043            0.063 -0.016            -0.003            -0.026                                                                                                                                                                
    ## temptrend_abs.sc:lifespan.sc               0.026  0.071 -0.013 -0.023  0.026  0.465  0.045  0.014 -0.515  0.021  0.645 -0.076 -0.034  0.112 -0.010  0.004  0.020 -0.052     -0.064      0.028            0.667  0.042             0.010            -0.831             0.093                                                                                                                                              
    ## temptrend_abs.sc:consumerfrac.sc          -0.054 -0.094  0.054  0.135  0.032 -0.122 -0.011 -0.002  0.039 -0.061 -0.111  0.360 -0.118 -0.002  0.003 -0.003 -0.034  0.103      0.309      0.046           -0.164 -0.026             0.000             0.061            -0.188            -0.175                                                                                                                            
    ## temptrend_abs.sc:endothermfrac.sc          0.072  0.164 -0.046 -0.145  0.140 -0.270 -0.050 -0.077 -0.002  0.023 -0.035 -0.135  0.443  0.054 -0.018 -0.035 -0.033 -0.111     -0.365      0.369           -0.413 -0.048            -0.065            -0.018             0.126            -0.032           -0.292                                                                                                           
    ## temptrend_abs.sc:nspp.sc                   0.012  0.036 -0.019 -0.027 -0.005  0.066 -0.030  0.010 -0.110 -0.029  0.125 -0.001  0.021  0.583 -0.029 -0.032 -0.037 -0.051     -0.071      0.034            0.040 -0.050            -0.022            -0.208            -0.083             0.201            0.001            0.103                                                                                          
    ## temptrend.sc:thermal_bias.sc               0.006  0.006 -0.007 -0.004  0.076  0.038  0.031 -0.054 -0.010  0.013  0.006 -0.004 -0.013 -0.023 -0.354 -0.023 -0.008 -0.007     -0.004      0.044            0.046  0.047            -0.040            -0.016             0.028             0.002           -0.003           -0.031            -0.029                                                                        
    ## temptrend_abs.sc:npp.sc                   -0.011 -0.052  0.001  0.009 -0.077  0.105 -0.138 -0.265  0.023  0.004  0.017 -0.026 -0.046 -0.030  0.050  0.013  0.750  0.033      0.029     -0.133            0.142 -0.132            -0.302             0.009            -0.006             0.031           -0.042           -0.078            -0.076            -0.005                                                      
    ## human.sc:REALM2TerrFresh                  -0.014 -0.009  0.006 -0.006 -0.230  0.092 -0.041  0.137  0.057  0.057 -0.062 -0.040 -0.063  0.022  0.044 -0.052 -0.008  0.008     -0.014     -0.172            0.084  0.097             0.056             0.038             0.031            -0.043           -0.012           -0.074            -0.015            -0.011  0.034                                               
    ## human.sc:REALM2Marine                      0.016  0.020 -0.017  0.003 -0.083  0.034 -0.187 -0.054  0.004 -0.001 -0.007  0.003 -0.015 -0.054 -0.003  0.001 -0.178 -0.022      0.005     -0.098            0.037 -0.176            -0.018             0.002            -0.024            -0.005           -0.005           -0.027            -0.054            -0.003 -0.135             0.007                             
    ## temptrend_abs.sc:human.sc:REALM2TerrFresh -0.009 -0.011  0.009 -0.008 -0.125  0.094  0.172  0.038  0.029  0.020 -0.032 -0.001 -0.051 -0.024  0.044  0.011  0.023  0.010     -0.024     -0.325            0.155  0.180             0.061             0.049             0.061            -0.058           -0.024           -0.133            -0.052             0.011  0.047             0.525     -0.020                  
    ## temptrend_abs.sc:human.sc:REALM2Marine     0.012  0.026 -0.013  0.001 -0.102  0.041 -0.161 -0.016  0.000 -0.003 -0.004 -0.004 -0.016 -0.056 -0.016 -0.035 -0.133 -0.029      0.001     -0.120            0.037 -0.200            -0.005            -0.002            -0.031            -0.004           -0.008           -0.022            -0.060             0.007 -0.181            -0.001      0.809     -0.016       
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.34115139 -0.34664125  0.05370777  0.53577626  6.61354550 
    ## 
    ## Number of Observations: 43585
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  43585

#### Plot the coefficients from the full model

``` r
coefs <- summary(modTfull1)$tTable
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
}
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20fullTmod1-1.png)<!-- -->

#### Plot interactions

``` r
# nspp effect
# trends[, range(Nspp)] # 2:1427
# trends[, range(temptrend, na.rm = TRUE)] # -1.3 to 1.9
# which.min(abs(ranef(modTfull1simpreml)$STUDY_ID$`(Intercept)`)) # study 127 has ranefs near 0
# trends[STUDY_ID == 127, sort(unique(rarefyID))] # ranefs near zero for 127_514668
# 

newdat <- cbind(expand.grid(temptrend = seq(-1.3, 2, length.out = 100), nspp = 2^seq(1, 9, length.out=500)), 
                data.table(REALM = 'Freshwater', REALM2 = 'TerrFresh', tempave.sc = 0, tempave_metab.sc = 0, 
                           seas.sc = 0, microclim.sc = 0, mass.sc = 0, 
                           speed.sc = 0, lifespan.sc = 0, endothermfrac.sc = 0, 
                           thermal_bias.sc = 0, npp.sc = 0, human.sc = 0,
                           consumerfrac.sc = 0,
                           nyrBT = 20, STUDY_ID = 127L, rarefyID = '127_514668'))
newdat <- rbind(newdat[1:2, ], newdat) # add two extra rows so that all factor levels are represented (for predict.lme to work)
newdat$REALM[1:2] <- c('Marine', 'Terrestrial')
newdat$REALM2[1:2] <- c('Marine', 'Marine')

# scale the vars
newdat$nspp.sc <- (log(newdat$nspp) - attr(trends$nspp.sc, 'scaled:center'))/attr(trends$nspp.sc, 'scaled:scale')
newdat$temptrend.sc <- (newdat$temptrend - attr(trends$temptrend.sc, 'scaled:center'))/attr(trends$temptrend.sc, 'scaled:scale') 
newdat$temptrend_abs <- abs(newdat$temptrend)
newdat$temptrend_abs.sc <- ((newdat$temptrend_abs) - attr(trends$temptrend_abs.sc, 'scaled:center'))/attr(trends$temptrend_abs.sc, 'scaled:scale')

newdat <- newdat[is.finite(newdat$temptrend_abs.sc),]

# make predictions
newdat$preds <- predict(object = modTfull1, newdata = newdat, level = 0)

#remove the extra rows
newdat <- newdat[newdat$REALM == 'Freshwater', ]

# plot
ggplot(newdat, aes(temptrend, preds, group = nspp, color = nspp)) +
  geom_line() +
  scale_color_distiller(palette = "YlGnBu", trans = 'log', breaks = c(2, 10, 50, 500))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/interaction%20plots%20modTfull1-1.png)<!-- -->

#### Plot residuals against each predictor

``` r
resids <- resid(modTfull1)
preds <- getData(modTfull1)
col = '#00000033'
cex = 0.5
par(mfrow = c(4,4))
boxplot(resids ~ preds$REALM, cex = cex, col = col)
plot(preds$temptrend_abs.sc, resids, cex = cex, col = col)
plot(preds$temptrend.sc, resids, cex = cex, col = col)
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
plot(preds$human.sc, resids, cex = cex, col = col)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/resids%20modTfull1-1.png)<!-- -->

### Sensitivity analysis: total turnover and Horn-Morisita models

#### Fit full models for total and HM

``` r
i2 <- trends[, complete.cases(Jbetatrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, human.sc)]
i3 <- trends[, complete.cases(Horntrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend.sc, temptrend_abs.sc, mass.sc, speed.sc, lifespan.sc, 
                             consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, human.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# full models
if(file.exists('temp/modTfullJbeta.rds')){
  modTfullJbeta <- readRDS('temp/modTfullJbeta.rds')
} else {
  modTfullJbeta <- lme(Jbetatrend ~ temptrend_abs.sc*REALM + 
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
                     temptrend.sc*thermal_bias.sc +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*human.sc:REALM2,
                   random = randef, weights = varef, data = trends[i2,], method = 'REML')
  saveRDS(modTfullJbeta, file = 'temp/modTfullJbeta.rds')
}

if(file.exists('temp/modTfullHorn.rds')){
  modTfullHorn <- readRDS('temp/modTfullHorn.rds')
} else {
  modTfullHorn <- lme(Horntrend ~ temptrend_abs.sc*REALM + 
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
                     temptrend.sc*thermal_bias.sc +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*human.sc:REALM2,
                   random = randef, weights = varef, data = trends[i3,], method = 'REML')
  saveRDS(modTfullHorn, file = 'temp/modTfullHorn.rds')
}

summary(modTfullJbeta)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##       AIC     BIC   logLik
    ##   -138984 -138628 69532.99
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06303875 (Intr)
    ## temptrend_abs.sc 0.03799046 0.262 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 0.0009214091 0.3174282
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.350847 
    ## Fixed effects: Jbetatrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tempave.sc +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      lifespan.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      endothermfrac.sc + temptrend_abs.sc * nspp.sc + temptrend.sc *      thermal_bias.sc + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      human.sc:REALM2 
    ##                                                 Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                0.06091987 0.015752189 43303   3.86739  0.0001
    ## temptrend_abs.sc                           0.07721276 0.016684119 43303   4.62792  0.0000
    ## REALMMarine                                0.04870999 0.016909043   247   2.88071  0.0043
    ## REALMTerrestrial                           0.02682193 0.018325686   247   1.46362  0.1446
    ## tempave.sc                                -0.00362418 0.000677129 43303  -5.35227  0.0000
    ## tempave_metab.sc                           0.01386063 0.001904526 43303   7.27774  0.0000
    ## seas.sc                                   -0.00286018 0.000391790 43303  -7.30029  0.0000
    ## microclim.sc                              -0.00209132 0.000231229 43303  -9.04437  0.0000
    ## mass.sc                                   -0.00229640 0.000791640 43303  -2.90082  0.0037
    ## speed.sc                                   0.00020316 0.000474085 43303   0.42854  0.6683
    ## lifespan.sc                                0.00108710 0.001592260 43303   0.68274  0.4948
    ## consumerfrac.sc                            0.00223456 0.000997386 43303   2.24041  0.0251
    ## endothermfrac.sc                          -0.01464435 0.004443529 43303  -3.29566  0.0010
    ## nspp.sc                                   -0.01638199 0.000425981 43303 -38.45710  0.0000
    ## temptrend.sc                              -0.00091281 0.000444927 43303  -2.05159  0.0402
    ## thermal_bias.sc                           -0.00135763 0.000381993 43303  -3.55408  0.0004
    ## npp.sc                                    -0.00016470 0.000356541 43303  -0.46194  0.6441
    ## temptrend_abs.sc:REALMMarine              -0.02150086 0.017228779 43303  -1.24796  0.2121
    ## temptrend_abs.sc:REALMTerrestrial         -0.02069394 0.020621387 43303  -1.00352  0.3156
    ## temptrend_abs.sc:tempave.sc                0.00539845 0.001107960 43303   4.87243  0.0000
    ## temptrend_abs.sc:tempave_metab.sc          0.00169231 0.003747805 43303   0.45155  0.6516
    ## temptrend_abs.sc:seas.sc                  -0.00163423 0.000781315 43303  -2.09164  0.0365
    ## temptrend_abs.sc:microclim.sc             -0.00467534 0.000504111 43303  -9.27442  0.0000
    ## temptrend_abs.sc:mass.sc                  -0.00457315 0.001433053 43303  -3.19119  0.0014
    ## temptrend_abs.sc:speed.sc                  0.00220093 0.000737068 43303   2.98606  0.0028
    ## temptrend_abs.sc:lifespan.sc               0.00710206 0.002948531 43303   2.40868  0.0160
    ## temptrend_abs.sc:consumerfrac.sc          -0.00085940 0.001290687 43303  -0.66585  0.5055
    ## temptrend_abs.sc:endothermfrac.sc          0.00658979 0.005663210 43303   1.16361  0.2446
    ## temptrend_abs.sc:nspp.sc                  -0.01172193 0.000755049 43303 -15.52474  0.0000
    ## temptrend.sc:thermal_bias.sc              -0.00094126 0.000323792 43303  -2.90699  0.0037
    ## temptrend_abs.sc:npp.sc                   -0.00190402 0.000753418 43303  -2.52718  0.0115
    ## human.sc:REALM2TerrFresh                  -0.00033521 0.000231701 43303  -1.44673  0.1480
    ## human.sc:REALM2Marine                      0.00039608 0.000321827 43303   1.23072  0.2184
    ## temptrend_abs.sc:human.sc:REALM2TerrFresh -0.00084918 0.000504605 43303  -1.68287  0.0924
    ## temptrend_abs.sc:human.sc:REALM2Marine    -0.00039613 0.000663029 43303  -0.59745  0.5502
    ##  Correlation: 
    ##                                           (Intr) tmpt_. REALMM REALMT tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s tmptr. thrm_. npp.sc t_.:REALMM t_.:REALMT tmptrnd_bs.sc:t. t_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. tm.:_. tmptrnd_bs.sc:np. h.:REALM2T h.:REALM2M t_.:.:REALM2T
    ## temptrend_abs.sc                           0.281                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                               -0.917 -0.256                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                          -0.883 -0.262  0.810                                                                                                                                                                                                                                                                                                                                                           
    ## tempave.sc                                -0.003 -0.007 -0.009 -0.003                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                           0.036  0.029 -0.026 -0.038 -0.316                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                   -0.024 -0.018  0.035 -0.015  0.205  0.011                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                              -0.008 -0.017  0.010 -0.011 -0.026  0.142  0.149                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                    0.011 -0.010 -0.001  0.000 -0.046 -0.516 -0.054 -0.033                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                   0.014  0.003 -0.023 -0.022 -0.011  0.123 -0.060 -0.015 -0.087                                                                                                                                                                                                                                                                                                                 
    ## lifespan.sc                                0.034  0.042 -0.021 -0.037 -0.021  0.720  0.075  0.034 -0.807  0.194                                                                                                                                                                                                                                                                                                          
    ## consumerfrac.sc                           -0.054 -0.054  0.068  0.324 -0.027 -0.080  0.007  0.009  0.051 -0.123 -0.109                                                                                                                                                                                                                                                                                                   
    ## endothermfrac.sc                           0.157  0.081 -0.083 -0.318  0.099 -0.190  0.008 -0.036  0.022  0.065 -0.027 -0.410                                                                                                                                                                                                                                                                                            
    ## nspp.sc                                   -0.010  0.017 -0.007 -0.018 -0.016  0.027 -0.029 -0.022 -0.162 -0.062  0.130  0.003  0.040                                                                                                                                                                                                                                                                                     
    ## temptrend.sc                              -0.006 -0.012  0.006  0.006 -0.099  0.004 -0.031  0.014  0.018 -0.004 -0.013  0.006 -0.005 -0.020                                                                                                                                                                                                                                                                              
    ## thermal_bias.sc                            0.019 -0.008 -0.030 -0.009  0.654 -0.003 -0.093 -0.093  0.030  0.028 -0.055 -0.012 -0.015 -0.075  0.001                                                                                                                                                                                                                                                                       
    ## npp.sc                                    -0.005 -0.027 -0.002  0.001 -0.215  0.172 -0.237 -0.290  0.009 -0.017  0.044 -0.051 -0.062 -0.120  0.024 -0.077                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:REALMMarine              -0.267 -0.954  0.274  0.248  0.006 -0.022  0.030  0.018  0.012 -0.008 -0.030  0.056 -0.056 -0.023  0.011  0.008  0.015                                                                                                                                                                                                                                                         
    ## temptrend_abs.sc:REALMTerrestrial         -0.248 -0.858  0.224  0.316  0.010 -0.042 -0.039 -0.020  0.013 -0.009 -0.035  0.144 -0.166 -0.035  0.008  0.001  0.019  0.809                                                                                                                                                                                                                                                  
    ## temptrend_abs.sc:tempave.sc               -0.002 -0.009 -0.001  0.009  0.549 -0.328  0.055 -0.066 -0.084 -0.041  0.012  0.017  0.132  0.020 -0.097  0.036 -0.083  0.007      0.021                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:tempave_metab.sc          0.013  0.044 -0.009 -0.019 -0.204  0.724  0.064  0.127 -0.325  0.043  0.477 -0.061 -0.174  0.024  0.012  0.070  0.086 -0.034     -0.055     -0.498                                                                                                                                                                                                                            
    ## temptrend_abs.sc:seas.sc                  -0.010 -0.033  0.016 -0.016  0.154  0.065  0.715  0.092 -0.049 -0.024  0.039  0.021 -0.017 -0.016 -0.004  0.086 -0.158  0.047     -0.045      0.169            0.056                                                                                                                                                                                                           
    ## temptrend_abs.sc:microclim.sc             -0.009 -0.025  0.007 -0.008 -0.028  0.135  0.097  0.810 -0.039  0.007  0.017  0.003 -0.047  0.012  0.025  0.002 -0.264  0.028     -0.020     -0.043            0.103  0.119                                                                                                                                                                                                    
    ## temptrend_abs.sc:mass.sc                  -0.003 -0.010  0.005  0.007 -0.087 -0.361 -0.056 -0.031  0.661 -0.003 -0.569  0.029  0.003 -0.113  0.022 -0.009  0.018  0.015      0.017     -0.113           -0.474 -0.049            -0.021                                                                                                                                                                                  
    ## temptrend_abs.sc:speed.sc                  0.001  0.007 -0.004 -0.009 -0.025  0.017 -0.014  0.003  0.013  0.363  0.014 -0.053  0.039 -0.047 -0.008  0.004  0.003 -0.022     -0.023     -0.041            0.056 -0.021            -0.006            -0.023                                                                                                                                                                
    ## temptrend_abs.sc:lifespan.sc               0.021  0.059 -0.011 -0.018  0.023  0.513  0.048  0.018 -0.553  0.038  0.703 -0.063 -0.040  0.114 -0.010  0.002  0.024 -0.043     -0.053      0.022            0.673  0.041             0.013            -0.820             0.088                                                                                                                                              
    ## temptrend_abs.sc:consumerfrac.sc          -0.049 -0.085  0.048  0.125  0.032 -0.109 -0.010 -0.001  0.033 -0.065 -0.100  0.324 -0.130  0.002  0.002  0.001 -0.026  0.097      0.302      0.035           -0.137 -0.024             0.002             0.048            -0.156            -0.148                                                                                                                            
    ## temptrend_abs.sc:endothermfrac.sc          0.068  0.167 -0.043 -0.141  0.133 -0.252 -0.036 -0.066  0.000  0.028 -0.038 -0.149  0.452  0.054 -0.016 -0.038 -0.035 -0.112     -0.373      0.312           -0.355 -0.034            -0.053            -0.010             0.109            -0.035           -0.313                                                                                                           
    ## temptrend_abs.sc:nspp.sc                   0.007  0.026 -0.013 -0.019 -0.009  0.059 -0.031  0.004 -0.117 -0.032  0.127  0.001  0.019  0.638 -0.028 -0.039 -0.050 -0.039     -0.055      0.030            0.029 -0.051            -0.026            -0.203            -0.088             0.189            0.003            0.091                                                                                          
    ## temptrend.sc:thermal_bias.sc               0.005  0.006 -0.006 -0.002  0.093  0.039  0.021 -0.055 -0.010  0.013  0.005 -0.002 -0.011 -0.028 -0.341 -0.027 -0.010 -0.007     -0.003      0.061            0.046  0.031            -0.043            -0.018             0.027             0.002           -0.003           -0.026            -0.033                                                                        
    ## temptrend_abs.sc:npp.sc                   -0.008 -0.040 -0.001  0.005 -0.096  0.117 -0.157 -0.270  0.026  0.002  0.018 -0.025 -0.046 -0.041  0.048  0.010  0.797  0.025      0.021     -0.142            0.145 -0.140            -0.296             0.010            -0.006             0.031           -0.036           -0.070            -0.083            -0.005                                                      
    ## human.sc:REALM2TerrFresh                  -0.012 -0.007  0.005 -0.006 -0.230  0.086 -0.038  0.134  0.059  0.060 -0.068 -0.044 -0.058  0.024  0.048 -0.044 -0.002  0.005     -0.011     -0.201            0.090  0.075             0.070             0.044             0.037            -0.052           -0.010           -0.071            -0.014            -0.018  0.037                                               
    ## human.sc:REALM2Marine                      0.012  0.017 -0.012  0.003 -0.092  0.036 -0.184 -0.060  0.003 -0.001 -0.005  0.007 -0.012 -0.058 -0.002 -0.008 -0.172 -0.018      0.004     -0.096            0.035 -0.172            -0.025             0.003            -0.023            -0.003           -0.004           -0.022            -0.059            -0.003 -0.138             0.006                             
    ## temptrend_abs.sc:human.sc:REALM2TerrFresh -0.006 -0.009  0.006 -0.004 -0.133  0.095  0.154  0.049  0.031  0.022 -0.036  0.005 -0.039 -0.022  0.046  0.021  0.020  0.008     -0.017     -0.319            0.147  0.149             0.069             0.048             0.060            -0.059           -0.022           -0.110            -0.048             0.001  0.042             0.615     -0.016                  
    ## temptrend_abs.sc:human.sc:REALM2Marine     0.009  0.021 -0.009  0.001 -0.108  0.040 -0.161 -0.024  0.000 -0.004 -0.002 -0.002 -0.013 -0.062 -0.013 -0.040 -0.139 -0.022      0.001     -0.111            0.034 -0.189            -0.012            -0.002            -0.030            -0.002           -0.004           -0.017            -0.063             0.007 -0.179             0.001      0.841     -0.012       
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -7.6056540 -0.2375536  0.1357522  0.6187386  6.5688252 
    ## 
    ## Number of Observations: 43585
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  43585

``` r
summary(modTfullHorn)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -125510.9 -125155.9 62796.43
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06204110 (Intr)
    ## temptrend_abs.sc 0.03226359 0.305 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.004394388 0.2966333
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.223201 
    ## Fixed effects: Horntrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tempave.sc +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      lifespan.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      endothermfrac.sc + temptrend_abs.sc * nspp.sc + temptrend.sc *      thermal_bias.sc + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      human.sc:REALM2 
    ##                                                 Value   Std.Error    DF    t-value p-value
    ## (Intercept)                                0.05107769 0.015619568 42336   3.270109  0.0011
    ## temptrend_abs.sc                           0.07705446 0.015193284 42336   5.071613  0.0000
    ## REALMMarine                                0.04486173 0.016961571   215   2.644904  0.0088
    ## REALMTerrestrial                           0.00872105 0.018497638   215   0.471468  0.6378
    ## tempave.sc                                 0.00132923 0.000898701 42336   1.479058  0.1391
    ## tempave_metab.sc                           0.00625370 0.002303117 42336   2.715322  0.0066
    ## seas.sc                                   -0.00006127 0.000507162 42336  -0.120802  0.9038
    ## microclim.sc                               0.00010304 0.000291187 42336   0.353859  0.7234
    ## mass.sc                                   -0.00067629 0.000930442 42336  -0.726847  0.4673
    ## speed.sc                                   0.00066257 0.000592004 42336   1.119194  0.2631
    ## lifespan.sc                               -0.00284719 0.001877082 42336  -1.516818  0.1293
    ## consumerfrac.sc                            0.00172298 0.001127264 42336   1.528465  0.1264
    ## endothermfrac.sc                          -0.01150114 0.004893353 42336  -2.350360  0.0188
    ## nspp.sc                                   -0.01365832 0.000516249 42336 -26.456832  0.0000
    ## temptrend.sc                              -0.00095985 0.000493696 42336  -1.944211  0.0519
    ## thermal_bias.sc                           -0.00032141 0.000502138 42336  -0.640085  0.5221
    ## npp.sc                                     0.00128657 0.000451059 42336   2.852338  0.0043
    ## temptrend_abs.sc:REALMMarine              -0.02349670 0.015758558 42336  -1.491044  0.1360
    ## temptrend_abs.sc:REALMTerrestrial         -0.05912389 0.019258933 42336  -3.069946  0.0021
    ## temptrend_abs.sc:tempave.sc                0.00745542 0.001273805 42336   5.852877  0.0000
    ## temptrend_abs.sc:tempave_metab.sc          0.00502134 0.004259192 42336   1.178943  0.2384
    ## temptrend_abs.sc:seas.sc                   0.00523937 0.000914639 42336   5.728342  0.0000
    ## temptrend_abs.sc:microclim.sc             -0.00103517 0.000601104 42336  -1.722118  0.0851
    ## temptrend_abs.sc:mass.sc                  -0.00132842 0.001621600 42336  -0.819202  0.4127
    ## temptrend_abs.sc:speed.sc                  0.00041821 0.000821272 42336   0.509224  0.6106
    ## temptrend_abs.sc:lifespan.sc               0.00329795 0.003322219 42336   0.992695  0.3209
    ## temptrend_abs.sc:consumerfrac.sc          -0.00064235 0.001408870 42336  -0.455933  0.6484
    ## temptrend_abs.sc:endothermfrac.sc          0.00635519 0.005711003 42336   1.112798  0.2658
    ## temptrend_abs.sc:nspp.sc                  -0.00755289 0.000853988 42336  -8.844257  0.0000
    ## temptrend.sc:thermal_bias.sc              -0.00030957 0.000356571 42336  -0.868189  0.3853
    ## temptrend_abs.sc:npp.sc                    0.00167656 0.000881725 42336   1.901455  0.0572
    ## human.sc:REALM2TerrFresh                   0.00065743 0.000351184 42336   1.872049  0.0612
    ## human.sc:REALM2Marine                      0.00094452 0.000391701 42336   2.411319  0.0159
    ## temptrend_abs.sc:human.sc:REALM2TerrFresh  0.00093774 0.000580983 42336   1.614049  0.1065
    ## temptrend_abs.sc:human.sc:REALM2Marine    -0.00008815 0.000809242 42336  -0.108928  0.9133
    ##  Correlation: 
    ##                                           (Intr) tmpt_. REALMM REALMT tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s tmptr. thrm_. npp.sc t_.:REALMM t_.:REALMT tmptrnd_bs.sc:t. t_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. tm.:_. tmptrnd_bs.sc:np. h.:REALM2T h.:REALM2M t_.:.:REALM2T
    ## temptrend_abs.sc                           0.295                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                               -0.904 -0.264                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                          -0.872 -0.271  0.782                                                                                                                                                                                                                                                                                                                                                           
    ## tempave.sc                                 0.001 -0.007 -0.012  0.002                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                           0.045  0.031 -0.032 -0.052 -0.374                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                   -0.038 -0.015  0.050 -0.012  0.161 -0.001                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                              -0.010 -0.016  0.012 -0.015 -0.058  0.126  0.118                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                    0.012 -0.013 -0.001  0.005 -0.043 -0.498 -0.048 -0.019                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                   0.016  0.000 -0.030 -0.027 -0.019  0.120 -0.069 -0.022 -0.093                                                                                                                                                                                                                                                                                                                 
    ## lifespan.sc                                0.042  0.046 -0.025 -0.045 -0.025  0.698  0.073  0.024 -0.808  0.190                                                                                                                                                                                                                                                                                                          
    ## consumerfrac.sc                           -0.057 -0.052  0.053  0.336 -0.001 -0.086 -0.007  0.003  0.051 -0.132 -0.112                                                                                                                                                                                                                                                                                                   
    ## endothermfrac.sc                           0.167  0.077 -0.085 -0.343  0.140 -0.213  0.002 -0.038  0.009  0.070 -0.020 -0.401                                                                                                                                                                                                                                                                                            
    ## nspp.sc                                   -0.014  0.021 -0.006 -0.017 -0.008  0.027 -0.029 -0.021 -0.173 -0.059  0.137  0.008  0.048                                                                                                                                                                                                                                                                                     
    ## temptrend.sc                              -0.006 -0.015  0.006  0.008 -0.082  0.003 -0.038  0.018  0.019 -0.002 -0.015  0.007 -0.007 -0.016                                                                                                                                                                                                                                                                              
    ## thermal_bias.sc                            0.029 -0.010 -0.041 -0.013  0.661 -0.032 -0.131 -0.131  0.032  0.028 -0.061  0.006  0.001 -0.066  0.000                                                                                                                                                                                                                                                                       
    ## npp.sc                                    -0.004 -0.031 -0.002  0.000 -0.224  0.183 -0.236 -0.305  0.000 -0.003  0.045 -0.029 -0.059 -0.121  0.015 -0.084                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:REALMMarine              -0.279 -0.948  0.288  0.254  0.003 -0.020  0.023  0.017  0.014 -0.003 -0.031  0.055 -0.055 -0.023  0.013  0.008  0.021                                                                                                                                                                                                                                                         
    ## temptrend_abs.sc:REALMTerrestrial         -0.255 -0.847  0.226  0.321  0.015 -0.056 -0.045 -0.025  0.018 -0.006 -0.039  0.127 -0.153 -0.037  0.011  0.000  0.016  0.791                                                                                                                                                                                                                                                  
    ## temptrend_abs.sc:tempave.sc                0.003 -0.008 -0.006  0.013  0.417 -0.273 -0.034 -0.092 -0.081 -0.034  0.013  0.025  0.125  0.018 -0.079  0.059 -0.049  0.006      0.030                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:tempave_metab.sc          0.012  0.053 -0.006 -0.026 -0.152  0.623  0.081  0.127 -0.281  0.025  0.410 -0.064 -0.163  0.031  0.010  0.055  0.062 -0.038     -0.077     -0.507                                                                                                                                                                                                                            
    ## temptrend_abs.sc:seas.sc                  -0.009 -0.033  0.011 -0.022  0.060  0.093  0.560  0.067 -0.038 -0.010  0.030  0.002 -0.042 -0.019 -0.003  0.070 -0.112  0.046     -0.068      0.120            0.080                                                                                                                                                                                                           
    ## temptrend_abs.sc:microclim.sc             -0.009 -0.030  0.006 -0.012 -0.054  0.128  0.071  0.700 -0.025  0.015  0.005 -0.005 -0.052  0.013  0.031  0.002 -0.236  0.033     -0.027     -0.065            0.103  0.116                                                                                                                                                                                                    
    ## temptrend_abs.sc:mass.sc                  -0.006 -0.017  0.007  0.010 -0.085 -0.305 -0.045 -0.023  0.577  0.013 -0.499  0.024 -0.003 -0.095  0.022 -0.017  0.019  0.022      0.031     -0.118           -0.475 -0.043            -0.016                                                                                                                                                                                  
    ## temptrend_abs.sc:speed.sc                  0.000  0.010  0.000 -0.008 -0.023  0.007  0.001  0.010  0.025  0.202 -0.007 -0.039  0.029 -0.038 -0.012 -0.001  0.002 -0.027     -0.028     -0.043            0.060 -0.016            -0.006            -0.027                                                                                                                                                                
    ## temptrend_abs.sc:lifespan.sc               0.023  0.072 -0.010 -0.021  0.028  0.434  0.040  0.008 -0.487  0.008  0.614 -0.064 -0.032  0.105 -0.010  0.008  0.017 -0.052     -0.066      0.027            0.663  0.039             0.007            -0.828             0.094                                                                                                                                              
    ## temptrend_abs.sc:consumerfrac.sc          -0.049 -0.094  0.054  0.105  0.031 -0.113 -0.015 -0.006  0.036 -0.039 -0.106  0.302 -0.094 -0.006  0.005  0.003 -0.012  0.102      0.286      0.041           -0.159 -0.032            -0.003             0.059            -0.165            -0.174                                                                                                                            
    ## temptrend_abs.sc:endothermfrac.sc          0.067  0.171 -0.043 -0.134  0.114 -0.239 -0.056 -0.075 -0.006  0.011 -0.032 -0.116  0.406  0.046 -0.017 -0.029 -0.020 -0.120     -0.387      0.360           -0.397 -0.053            -0.061            -0.025             0.122            -0.028           -0.273                                                                                                           
    ## temptrend_abs.sc:nspp.sc                   0.011  0.039 -0.015 -0.021 -0.008  0.064 -0.033  0.014 -0.102 -0.025  0.119 -0.001  0.013  0.545 -0.027 -0.031 -0.025 -0.050     -0.071      0.038            0.038 -0.043            -0.019            -0.208            -0.079             0.202           -0.003            0.099                                                                                          
    ## temptrend.sc:thermal_bias.sc               0.005  0.006 -0.006 -0.006  0.064  0.043  0.037 -0.054 -0.010  0.012  0.006 -0.005 -0.013 -0.022 -0.361 -0.018 -0.003 -0.007     -0.008      0.033            0.052  0.058            -0.037            -0.017             0.030             0.002           -0.006           -0.033            -0.028                                                                        
    ## temptrend_abs.sc:npp.sc                   -0.010 -0.053  0.001  0.005 -0.047  0.095 -0.116 -0.251  0.018  0.005  0.018 -0.008 -0.032 -0.019  0.048  0.022  0.686  0.036      0.025     -0.132            0.146 -0.140            -0.320             0.007            -0.010             0.031           -0.014           -0.073            -0.078            -0.005                                                      
    ## human.sc:REALM2TerrFresh                  -0.013 -0.007  0.007 -0.008 -0.223  0.100 -0.042  0.135  0.049  0.052 -0.052 -0.043 -0.057  0.018  0.039 -0.055 -0.014  0.006     -0.016     -0.134            0.072  0.105             0.036             0.030             0.022            -0.035           -0.006           -0.058            -0.010            -0.008  0.030                                               
    ## human.sc:REALM2Marine                      0.016  0.020 -0.018  0.000 -0.072  0.029 -0.187 -0.054  0.006  0.005 -0.011  0.000 -0.008 -0.058  0.001  0.012 -0.182 -0.025      0.001     -0.087            0.034 -0.162            -0.009             0.002            -0.025            -0.004           -0.011           -0.016            -0.053            -0.007 -0.132             0.007                             
    ## temptrend_abs.sc:human.sc:REALM2TerrFresh -0.008 -0.009  0.007 -0.008 -0.117  0.091  0.178  0.029  0.027  0.017 -0.031  0.002 -0.045 -0.024  0.043  0.002  0.024  0.008     -0.029     -0.331            0.163  0.204             0.055             0.049             0.060            -0.058           -0.019           -0.132            -0.056             0.018  0.049             0.399     -0.023                  
    ## temptrend_abs.sc:human.sc:REALM2Marine     0.009  0.028 -0.011 -0.001 -0.085  0.034 -0.132  0.001 -0.002 -0.008 -0.002 -0.008 -0.006 -0.046 -0.018 -0.038 -0.126 -0.036     -0.005     -0.135            0.041 -0.216            -0.001             0.000            -0.028            -0.006           -0.017           -0.011            -0.072             0.009 -0.168            -0.007      0.721     -0.016       
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.83865023 -0.37581189  0.04109756  0.54195029  6.16746015 
    ## 
    ## Number of Observations: 42586
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    218                  42586

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

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20simplified%20total%20and%20HM%20mod%20coefs-1.png)<!-- -->

Black is for Jaccard total turnover (pres/abs), grey is for
Morisita-Horn turnover (considers
abundance)

### Plot coefficients from all full models

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

# add base temperature:human effect (freshwater) to realm-specific temperature:human effects
allcoefsfull$Value[allcoefsfull$var == 'temptrend_abs.sc:REALMMarine:human.sc'] <- with(allcoefsfull, Value[var == 'temptrend_abs.sc:REALMMarine:human.sc'] + Value[var == 'temptrend_abs.sc:human.sc'])
allcoefsfull$Value[allcoefsfull$var == 'temptrend_abs.sc:REALMTerrestrial:human.sc'] <- with(allcoefsfull, Value[var == 'temptrend_abs.sc:REALMTerrestrial:human.sc'] + Value[var == 'temptrend_abs.sc:human.sc'])


# remove non-temperature effects
allcoefsfull <- allcoefsfull[grepl(':', allcoefsfull$var) & grepl('temptrend', allcoefsfull$var), ]

# add info for plotting
allcoefsfull$lCI <- allcoefsfull$Value - allcoefsfull$Std.Error # lower confidence interval
allcoefsfull$uCI <- allcoefsfull$Value + allcoefsfull$Std.Error
nvar <- nrow(allcoefsfull)/3
allcoefsfull$y <- 1:nvar + rep(c(0, 0.1, 0.2), c(nvar, nvar, nvar)) # y-values

allcoefsfull$varname <- gsub('temptrend_abs.sc:|temptrend.sc:', '', allcoefsfull$var)
allcoefsfull$varname <- gsub('REALM', '', allcoefsfull$varname)
allcoefsfull$varname <- gsub('.sc', '', allcoefsfull$varname)
allcoefsfull$varname <- gsub('^human$', 'Freshwater:human', allcoefsfull$varname)

xlims1 <- c(0.0, 0.05) # for realms
xlims2 <- c(-0.01, 0.015) # for traits
xlims3 <- c(-0.004, 0.0025) # for environment
xlims4 <- c(-0.016, 0.005) # for community
xlims5 <- c(-0.01, 0.015) # for human

ddg <- 0.5 # dodge for each model

set1 <- c('Terrestrial', 'Marine', 'Freshwater')
set2 <- c('mass', 'speed', 'lifespan', 'consumerfrac', 'endothermfrac', 'tempave_metab')
set3 <- c('seas', 'microclim', 'tempave')
set4 <- c('npp', 'nspp', 'thermal_bias')
set5 <- c('Terrestrial:human', 'Marine:human', 'Freshwater:human')

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

# To do

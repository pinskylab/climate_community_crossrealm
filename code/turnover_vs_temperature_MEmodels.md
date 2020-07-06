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
library(gridExtra)

options(width=500) # turn off most text wrapping

# tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) 
```

``` r
# Turnover and covariates assembled by turnover_vs_temperature_prep.Rmd
trends <- fread('output/turnover_w_covariates.csv.gz')

# set realm order
trends[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = FALSE)]

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
trends[, temptrend_abs.sc := scale(log(abs(temptrend)))]
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

Do the variables look ok?

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

Check correlations among variables. Pearson’s r is on the lower
diagonal.

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

Examine how many data points are available

``` r
# the cases we can compare
apply(trends[, .(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, human.sc)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##         Jtutrend            REALM       tempave.sc tempave_metab.sc          seas.sc     microclim.sc     temptrend.sc          mass.sc         speed.sc      lifespan.sc  consumerfrac.sc endothermfrac.sc          nspp.sc  thermal_bias.sc           npp.sc         human.sc 
    ##            53013            53013            49916            49916            49916            51834            49916            52820            52689            51540            47534            53013            53013            49371            52863            53013

``` r
i <- trends[, complete.cases(Jtutrend, temptrend.sc, tempave_metab.sc, REALM, seas.sc, microclim.sc, npp.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, thermal_bias.sc)]
cat('Overall:\n')
```

    ## Overall:

``` r
sum(i)
```

    ## [1] 43585

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
    scale_color_brewer(palette="Set1")
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

## Jaccard turnover temporal trend

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

### Full model

``` r
i <- trends[, complete.cases(Jtutrend, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, temptrend_abs.sc,  
                             mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, human.sc)]

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
                     temptrend_abs.sc*human.sc*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modTfull1, file = 'temp/modTfull1.rds')
}

summary(modTfull1)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##       AIC       BIC   logLik
    ##   -130336 -129962.7 65211.02
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.05028541 (Intr)
    ## temptrend_abs.sc 0.01934781 0.555 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.001189033 0.3037859
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.246996 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tempave.sc +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      lifespan.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      endothermfrac.sc + temptrend_abs.sc * nspp.sc + temptrend.sc *      thermal_bias.sc + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      human.sc * REALM 
    ##                                                  Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                 0.02193634 0.012919362 43301   1.69794  0.0895
    ## temptrend_abs.sc                            0.02533530 0.008335922 43301   3.03929  0.0024
    ## REALMMarine                                 0.05429560 0.013804566   247   3.93316  0.0001
    ## REALMTerrestrial                            0.02584205 0.015139481   247   1.70693  0.0891
    ## tempave.sc                                 -0.00144234 0.000735841 43301  -1.96013  0.0500
    ## tempave_metab.sc                            0.00962407 0.001761682 43301   5.46300  0.0000
    ## seas.sc                                     0.00048937 0.000356628 43301   1.37221  0.1700
    ## microclim.sc                               -0.00032055 0.000183570 43301  -1.74622  0.0808
    ## mass.sc                                    -0.00410866 0.000800371 43301  -5.13344  0.0000
    ## speed.sc                                   -0.00080881 0.000529985 43301  -1.52610  0.1270
    ## lifespan.sc                                 0.00289214 0.001551058 43301   1.86463  0.0622
    ## consumerfrac.sc                             0.00067647 0.000858159 43301   0.78828  0.4305
    ## endothermfrac.sc                           -0.01381349 0.004043987 43301  -3.41581  0.0006
    ## nspp.sc                                    -0.01320044 0.000416290 43301 -31.70976  0.0000
    ## temptrend.sc                               -0.00006458 0.000473430 43301  -0.13641  0.8915
    ## thermal_bias.sc                             0.00005235 0.000479384 43301   0.10920  0.9130
    ## npp.sc                                      0.00123480 0.000295708 43301   4.17575  0.0000
    ## human.sc                                    0.00835896 0.004222447 43301   1.97965  0.0477
    ## temptrend_abs.sc:REALMMarine               -0.00463118 0.008642421 43301  -0.53587  0.5921
    ## temptrend_abs.sc:REALMTerrestrial          -0.00594687 0.009874128 43301  -0.60227  0.5470
    ## temptrend_abs.sc:tempave.sc                -0.00036731 0.000817673 43301  -0.44922  0.6533
    ## temptrend_abs.sc:tempave_metab.sc           0.00736415 0.002169688 43301   3.39410  0.0007
    ## temptrend_abs.sc:seas.sc                    0.00124640 0.000530485 43301   2.34955  0.0188
    ## temptrend_abs.sc:microclim.sc              -0.00115860 0.000263202 43301  -4.40193  0.0000
    ## temptrend_abs.sc:mass.sc                   -0.00482207 0.000826466 43301  -5.83457  0.0000
    ## temptrend_abs.sc:speed.sc                   0.00070334 0.000580200 43301   1.21223  0.2254
    ## temptrend_abs.sc:lifespan.sc                0.00657667 0.001593947 43301   4.12603  0.0000
    ## temptrend_abs.sc:consumerfrac.sc           -0.00072374 0.000640651 43301  -1.12969  0.2586
    ## temptrend_abs.sc:endothermfrac.sc          -0.00187699 0.002817693 43301  -0.66614  0.5053
    ## temptrend_abs.sc:nspp.sc                   -0.00970260 0.000474539 43301 -20.44637  0.0000
    ## temptrend.sc:thermal_bias.sc               -0.00073599 0.000341854 43301  -2.15294  0.0313
    ## temptrend_abs.sc:npp.sc                     0.00091619 0.000394971 43301   2.31964  0.0204
    ## temptrend_abs.sc:human.sc                   0.00672778 0.004646579 43301   1.44790  0.1477
    ## REALMMarine:human.sc                       -0.00863342 0.004228477 43301  -2.04173  0.0412
    ## REALMTerrestrial:human.sc                  -0.00904734 0.004229618 43301  -2.13904  0.0324
    ## temptrend_abs.sc:REALMMarine:human.sc      -0.00815689 0.004652338 43301  -1.75329  0.0796
    ## temptrend_abs.sc:REALMTerrestrial:human.sc -0.00887481 0.004677152 43301  -1.89748  0.0578
    ##  Correlation: 
    ##                                            (Intr) tmpt_. REALMMr REALMTr tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s tmptr. thrm_. npp.sc hmn.sc tm_.:REALMM tm_.:REALMT tmptrnd_bs.sc:t. t_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. tm.:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:h. REALMM: REALMT: t_.:REALMM:
    ## temptrend_abs.sc                            0.368                                                                                                                                                                                                                                                                                                                                                                                             
    ## REALMMarine                                -0.918 -0.337                                                                                                                                                                                                                                                                                                                                                                                      
    ## REALMTerrestrial                           -0.882 -0.333  0.807                                                                                                                                                                                                                                                                                                                                                                               
    ## tempave.sc                                  0.000 -0.007 -0.013  -0.003                                                                                                                                                                                                                                                                                                                                                                       
    ## tempave_metab.sc                            0.053  0.015 -0.042  -0.045  -0.278                                                                                                                                                                                                                                                                                                                                                               
    ## seas.sc                                    -0.037 -0.004  0.053  -0.003   0.271 -0.083                                                                                                                                                                                                                                                                                                                                                        
    ## microclim.sc                               -0.008  0.002  0.013  -0.010   0.049  0.067  0.170                                                                                                                                                                                                                                                                                                                                                 
    ## mass.sc                                     0.019  0.000 -0.004   0.001  -0.002 -0.550 -0.046 -0.015                                                                                                                                                                                                                                                                                                                                          
    ## speed.sc                                    0.021  0.007 -0.034  -0.032  -0.003  0.170 -0.075 -0.032 -0.135                                                                                                                                                                                                                                                                                                                                   
    ## lifespan.sc                                 0.037  0.021 -0.024  -0.043  -0.037  0.752  0.094  0.047 -0.806  0.265                                                                                                                                                                                                                                                                                                                            
    ## consumerfrac.sc                            -0.067 -0.057  0.074   0.342  -0.014 -0.090  0.004  0.007  0.063 -0.163 -0.123                                                                                                                                                                                                                                                                                                                     
    ## endothermfrac.sc                            0.169  0.077 -0.087  -0.349   0.080 -0.167  0.048 -0.002  0.014  0.083 -0.014 -0.440                                                                                                                                                                                                                                                                                                              
    ## nspp.sc                                    -0.023  0.000 -0.002  -0.014  -0.047 -0.021 -0.033 -0.056 -0.156 -0.050  0.098  0.014  0.057                                                                                                                                                                                                                                                                                                       
    ## temptrend.sc                               -0.005 -0.009  0.006   0.008  -0.060 -0.018 -0.044 -0.013  0.015 -0.011 -0.012  0.011  0.002  0.015                                                                                                                                                                                                                                                                                                
    ## thermal_bias.sc                             0.030 -0.006 -0.044  -0.008   0.759 -0.044 -0.152 -0.131  0.039  0.019 -0.067  0.001 -0.006 -0.106 -0.001                                                                                                                                                                                                                                                                                         
    ## npp.sc                                      0.004  0.000 -0.006  -0.001  -0.320  0.241 -0.256 -0.220 -0.011 -0.021  0.059 -0.047 -0.090 -0.162 -0.012 -0.137                                                                                                                                                                                                                                                                                  
    ## human.sc                                   -0.166 -0.023  0.157   0.146   0.017 -0.005  0.008  0.017  0.006  0.004  0.007  0.013 -0.002 -0.015  0.008  0.016 -0.029                                                                                                                                                                                                                                                                           
    ## temptrend_abs.sc:REALMMarine               -0.349 -0.953  0.369   0.315   0.007 -0.009  0.014  0.001  0.006 -0.011 -0.015  0.060 -0.051 -0.007  0.006  0.005 -0.008  0.023                                                                                                                                                                                                                                                                    
    ## temptrend_abs.sc:REALMTerrestrial          -0.327 -0.874  0.301   0.367  -0.007 -0.016 -0.020 -0.012  0.012 -0.010 -0.022  0.139 -0.130 -0.015  0.009  0.001  0.013  0.020  0.829                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave.sc                 0.001 -0.004 -0.002  -0.003   0.099 -0.059  0.036 -0.028 -0.035 -0.031  0.005 -0.005  0.024 -0.002 -0.071  0.016 -0.031 -0.004  0.003       0.003                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:tempave_metab.sc           0.009  0.032 -0.004  -0.005   0.033  0.330  0.018  0.085 -0.207  0.074  0.300 -0.034 -0.060 -0.025  0.001  0.090  0.027  0.003 -0.023      -0.030      -0.639                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:seas.sc                    0.012 -0.029 -0.008  -0.023   0.153  0.024  0.227  0.027 -0.039 -0.029  0.042  0.016  0.002 -0.006 -0.029  0.098 -0.178 -0.003  0.054      -0.050       0.206            0.017                                                                                                                                                                                                                    
    ## temptrend_abs.sc:microclim.sc              -0.004 -0.032  0.001  -0.002   0.012  0.099  0.059  0.312 -0.033  0.008  0.022  0.002 -0.032  0.005  0.003 -0.029 -0.141  0.002  0.042      -0.002      -0.039            0.010  0.048                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                   -0.001 -0.002  0.007   0.009  -0.016 -0.272 -0.038 -0.024  0.459 -0.043 -0.385  0.037  0.006 -0.056  0.019  0.014  0.007  0.009  0.015       0.027      -0.059           -0.424 -0.029             0.014                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                  -0.001 -0.004  0.000  -0.005  -0.009  0.028 -0.028  0.002 -0.014  0.180  0.051 -0.041  0.027 -0.024 -0.006  0.010  0.012  0.012 -0.024      -0.008      -0.042            0.093 -0.037            -0.038            -0.059                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                0.009  0.040 -0.003  -0.010   0.000  0.378  0.049  0.026 -0.385  0.095  0.475 -0.057 -0.035  0.040 -0.015 -0.014  0.016 -0.001 -0.033      -0.040      -0.003            0.633  0.055             0.013            -0.805             0.150                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc           -0.050 -0.058  0.050   0.124  -0.008 -0.067 -0.005 -0.003  0.037 -0.055 -0.077  0.278 -0.127  0.008  0.009 -0.009 -0.018 -0.011  0.071       0.267       0.023           -0.114 -0.029             0.008             0.062            -0.234            -0.147                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc           0.065  0.136 -0.038  -0.120  -0.027 -0.082 -0.008 -0.044 -0.003  0.008 -0.026 -0.144  0.346  0.048 -0.011 -0.064 -0.014 -0.001 -0.097      -0.290       0.450           -0.452  0.001            -0.011            -0.038             0.150            -0.037           -0.282                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                    0.001  0.040 -0.006  -0.012  -0.035 -0.005  0.002  0.006 -0.064 -0.001  0.057 -0.003  0.022  0.308  0.018 -0.063 -0.020 -0.013 -0.060      -0.071       0.041           -0.012 -0.056            -0.164            -0.211            -0.139             0.186            0.024            0.112                                                                                                   
    ## temptrend.sc:thermal_bias.sc                0.006  0.012 -0.008  -0.005   0.081  0.032  0.005 -0.033 -0.008  0.010  0.009 -0.007 -0.004 -0.014 -0.407 -0.008 -0.005 -0.006 -0.013      -0.006       0.100            0.001 -0.003            -0.041            -0.012             0.012             0.005           -0.004            0.007            -0.003                                                                                 
    ## temptrend_abs.sc:npp.sc                     0.003 -0.051 -0.008   0.005  -0.127  0.106 -0.164 -0.156  0.016  0.006  0.013 -0.021 -0.051 -0.034 -0.005 -0.044  0.367 -0.012  0.037       0.044      -0.158            0.138 -0.192            -0.203            -0.023            -0.006             0.043           -0.019           -0.079            -0.189            -0.014                                                               
    ## temptrend_abs.sc:human.sc                  -0.068  0.149  0.064   0.048  -0.005 -0.007  0.004  0.006  0.007  0.004 -0.004 -0.035  0.017 -0.003  0.007 -0.009 -0.010  0.408 -0.142      -0.119      -0.007           -0.018 -0.021             0.038             0.048             0.022            -0.025            0.016            0.001            -0.034            -0.005 -0.033                                                        
    ## REALMMarine:human.sc                        0.166  0.024 -0.157  -0.146  -0.020  0.004 -0.017 -0.024 -0.005 -0.003 -0.008 -0.013  0.002  0.013 -0.008 -0.015  0.019 -0.998 -0.023      -0.019       0.002           -0.002 -0.002            -0.005            -0.008            -0.011             0.001            0.011            0.000             0.012             0.005  0.008            -0.407                                      
    ## REALMTerrestrial:human.sc                   0.166  0.023 -0.157  -0.147  -0.031  0.008 -0.022 -0.006 -0.002  0.000 -0.012 -0.016 -0.001  0.018 -0.008 -0.021  0.028 -0.998 -0.023      -0.019       0.005           -0.005 -0.002             0.001            -0.008            -0.011            -0.001            0.011            0.001             0.014             0.004  0.012            -0.407            0.997                     
    ## temptrend_abs.sc:REALMMarine:human.sc       0.068 -0.147 -0.064  -0.048   0.001  0.007 -0.008 -0.008 -0.007 -0.003  0.004  0.035 -0.017  0.001 -0.006  0.006  0.007 -0.407  0.140       0.119       0.001            0.022  0.011            -0.037            -0.047            -0.022             0.025           -0.016           -0.003             0.029             0.005  0.024            -0.998            0.408   0.406             
    ## temptrend_abs.sc:REALMTerrestrial:human.sc  0.068 -0.148 -0.064  -0.047   0.010  0.006 -0.002 -0.004 -0.006 -0.004  0.004  0.037 -0.016  0.004 -0.005  0.011  0.013 -0.405  0.141       0.116      -0.019            0.030  0.013            -0.025            -0.043            -0.019             0.020           -0.022           -0.010             0.033             0.004  0.033            -0.993            0.404   0.403   0.991     
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -7.58461466 -0.34521721  0.05131203  0.55413461  6.93144917 
    ## 
    ## Number of Observations: 43585
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  43585

#### Try simplifying the model

``` r
if(file.exists('temp/modTfull1simp.rds')){
  modTfull1simp <- readRDS('temp/modTfull1simp.rds')
} else {
  require(MASS) # for stepAIC
  modTfull1ml <- update(modTfull1, method = 'ML')
  modTfull1simp <- stepAIC(modTfull1ml, direction = 'backward')
  saveRDS(modTfull1simp, file = 'temp/modTfull1simp.rds')
}
summary(modTfull1simp)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -130797.6 -130467.7 65436.79
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.04935417 (Intr)
    ## temptrend_abs.sc 0.01813999 0.582 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.001168172 0.3039957
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.247572 
    ## Fixed effects: Jtutrend ~ temptrend_abs.sc + REALM + tempave.sc + tempave_metab.sc +      seas.sc + microclim.sc + mass.sc + speed.sc + lifespan.sc +      endothermfrac.sc + nspp.sc + temptrend.sc + thermal_bias.sc +      npp.sc + human.sc + temptrend_abs.sc:REALM + temptrend_abs.sc:tempave_metab.sc +      temptrend_abs.sc:seas.sc + temptrend_abs.sc:microclim.sc +      temptrend_abs.sc:mass.sc + temptrend_abs.sc:lifespan.sc +      temptrend_abs.sc:nspp.sc + temptrend.sc:thermal_bias.sc +      temptrend_abs.sc:npp.sc + temptrend_abs.sc:human.sc + REALM:human.sc +      temptrend_abs.sc:REALM:human.sc 
    ##                                                  Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                 0.02291872 0.012617466 43306   1.81643  0.0693
    ## temptrend_abs.sc                            0.02582900 0.007804995 43306   3.30929  0.0009
    ## REALMMarine                                 0.05343691 0.013495830   247   3.95951  0.0001
    ## REALMTerrestrial                            0.02083818 0.013905374   247   1.49857  0.1353
    ## tempave.sc                                 -0.00146784 0.000727891 43306  -2.01657  0.0437
    ## tempave_metab.sc                            0.00936746 0.001738319 43306   5.38880  0.0000
    ## seas.sc                                     0.00049069 0.000355592 43306   1.37992  0.1676
    ## microclim.sc                               -0.00033312 0.000183109 43306  -1.81923  0.0689
    ## mass.sc                                    -0.00409834 0.000797222 43306  -5.14077  0.0000
    ## speed.sc                                   -0.00086515 0.000513200 43306  -1.68579  0.0918
    ## lifespan.sc                                 0.00273645 0.001530552 43306   1.78788  0.0738
    ## endothermfrac.sc                           -0.01106519 0.003328370 43306  -3.32451  0.0009
    ## nspp.sc                                    -0.01315130 0.000414867 43306 -31.70005  0.0000
    ## temptrend.sc                               -0.00006325 0.000471918 43306  -0.13402  0.8934
    ## thermal_bias.sc                             0.00000442 0.000476728 43306   0.00926  0.9926
    ## npp.sc                                      0.00123123 0.000294793 43306   4.17659  0.0000
    ## human.sc                                    0.00817197 0.004206767 43306   1.94258  0.0521
    ## temptrend_abs.sc:REALMMarine               -0.00445514 0.008134219 43306  -0.54770  0.5839
    ## temptrend_abs.sc:REALMTerrestrial          -0.00634443 0.008666955 43306  -0.73203  0.4642
    ## temptrend_abs.sc:tempave_metab.sc           0.00582453 0.001556759 43306   3.74145  0.0002
    ## temptrend_abs.sc:seas.sc                    0.00122237 0.000512005 43306   2.38742  0.0170
    ## temptrend_abs.sc:microclim.sc              -0.00114277 0.000262252 43306  -4.35753  0.0000
    ## temptrend_abs.sc:mass.sc                   -0.00476283 0.000821190 43306  -5.79991  0.0000
    ## temptrend_abs.sc:lifespan.sc                0.00594367 0.001550134 43306   3.83429  0.0001
    ## temptrend_abs.sc:nspp.sc                   -0.00952124 0.000463453 43306 -20.54412  0.0000
    ## temptrend.sc:thermal_bias.sc               -0.00073827 0.000339622 43306  -2.17381  0.0297
    ## temptrend_abs.sc:npp.sc                     0.00088836 0.000388931 43306   2.28410  0.0224
    ## temptrend_abs.sc:human.sc                   0.00696016 0.004605571 43306   1.51125  0.1307
    ## REALMMarine:human.sc                       -0.00846061 0.004212891 43306  -2.00827  0.0446
    ## REALMTerrestrial:human.sc                  -0.00884525 0.004213690 43306  -2.09917  0.0358
    ## temptrend_abs.sc:REALMMarine:human.sc      -0.00839988 0.004611409 43306  -1.82154  0.0685
    ## temptrend_abs.sc:REALMTerrestrial:human.sc -0.00920552 0.004634476 43306  -1.98631  0.0470
    ##  Correlation: 
    ##                                            (Intr) tmpt_. REALMMr REALMTr tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. endth. nspp.s tmptr. thrm_. npp.sc hmn.sc tm_.:REALMM tm_.:REALMT t_.:_. tmptrnd_bs.sc:s. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:l. tmptrnd_bs.sc:ns. tm.:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:h. REALMM: REALMT: t_.:REALMM:
    ## temptrend_abs.sc                            0.372                                                                                                                                                                                                                                                                                                               
    ## REALMMarine                                -0.918 -0.344                                                                                                                                                                                                                                                                                                        
    ## REALMTerrestrial                           -0.916 -0.336  0.834                                                                                                                                                                                                                                                                                                 
    ## tempave.sc                                  0.004  0.006 -0.015  -0.006                                                                                                                                                                                                                                                                                         
    ## tempave_metab.sc                            0.053  0.023 -0.038  -0.023  -0.289                                                                                                                                                                                                                                                                                 
    ## seas.sc                                    -0.037 -0.001  0.053  -0.006   0.268 -0.084                                                                                                                                                                                                                                                                          
    ## microclim.sc                               -0.005  0.007  0.012  -0.017   0.049  0.063  0.171                                                                                                                                                                                                                                                                   
    ## mass.sc                                     0.022  0.000 -0.009  -0.020   0.006 -0.551 -0.045 -0.015                                                                                                                                                                                                                                                            
    ## speed.sc                                    0.014  0.005 -0.025   0.022  -0.003  0.153 -0.071 -0.034 -0.128                                                                                                                                                                                                                                                     
    ## lifespan.sc                                 0.033  0.024 -0.016  -0.006  -0.048  0.750  0.095  0.046 -0.807  0.246                                                                                                                                                                                                                                              
    ## endothermfrac.sc                            0.145  0.009 -0.055  -0.222   0.126 -0.213  0.071  0.018  0.037  0.024 -0.054                                                                                                                                                                                                                                       
    ## nspp.sc                                    -0.027 -0.009 -0.002  -0.014  -0.041 -0.012 -0.032 -0.053 -0.161 -0.043  0.108  0.051                                                                                                                                                                                                                                
    ## temptrend.sc                               -0.007 -0.013  0.006   0.007  -0.050 -0.017 -0.041 -0.014  0.010 -0.009 -0.007 -0.003  0.013                                                                                                                                                                                                                         
    ## thermal_bias.sc                             0.036  0.008 -0.047  -0.017   0.760 -0.054 -0.156 -0.135  0.043  0.016 -0.076  0.029 -0.100  0.003                                                                                                                                                                                                                  
    ## npp.sc                                      0.001 -0.002 -0.003   0.016  -0.322  0.237 -0.255 -0.221 -0.009 -0.032  0.053 -0.130 -0.161 -0.013 -0.138                                                                                                                                                                                                           
    ## human.sc                                   -0.168 -0.019  0.159   0.154   0.018 -0.004  0.008  0.017  0.006  0.004  0.008  0.006 -0.015  0.008  0.015 -0.029                                                                                                                                                                                                    
    ## temptrend_abs.sc:REALMMarine               -0.353 -0.955  0.377   0.318   0.000 -0.010  0.012 -0.003  0.004 -0.002 -0.012  0.005 -0.003  0.008 -0.003 -0.006  0.020                                                                                                                                                                                             
    ## temptrend_abs.sc:REALMTerrestrial          -0.331 -0.896  0.309   0.349  -0.030 -0.022 -0.028 -0.025  0.009 -0.002 -0.018  0.018  0.000  0.015 -0.024  0.020  0.019  0.855                                                                                                                                                                                      
    ## temptrend_abs.sc:tempave_metab.sc           0.026  0.083 -0.009  -0.018   0.102  0.364  0.049  0.079 -0.303  0.047  0.385  0.034 -0.009 -0.049  0.107  0.003 -0.001 -0.050      -0.101                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                    0.019 -0.013 -0.011  -0.041   0.126  0.025  0.223  0.029 -0.028 -0.017  0.035  0.055  0.001 -0.011  0.086 -0.178 -0.004  0.046      -0.086       0.168                                                                                                                                                              
    ## temptrend_abs.sc:microclim.sc              -0.005 -0.038  0.002  -0.002   0.017  0.100  0.060  0.314 -0.036  0.015  0.025 -0.043  0.003  0.000 -0.026 -0.142  0.002  0.045       0.002      -0.012  0.058                                                                                                                                                       
    ## temptrend_abs.sc:mass.sc                    0.002 -0.002  0.003  -0.004  -0.008 -0.274 -0.037 -0.025  0.457 -0.031 -0.382  0.024 -0.060  0.013  0.017  0.008  0.010  0.012       0.018      -0.630 -0.016            0.008                                                                                                                                      
    ## temptrend_abs.sc:lifespan.sc                0.010  0.051  0.001   0.003  -0.011  0.367  0.049  0.022 -0.385  0.063  0.466 -0.029  0.054 -0.010 -0.027  0.011 -0.004 -0.033      -0.039       0.829  0.044            0.020            -0.815                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                   -0.008  0.023 -0.002   0.002  -0.027  0.015  0.001  0.013 -0.071  0.029  0.076 -0.031  0.303  0.015 -0.051 -0.017 -0.012 -0.056      -0.041       0.083 -0.056           -0.173            -0.224             0.233                                                                                                  
    ## temptrend.sc:thermal_bias.sc                0.009  0.019 -0.009  -0.007   0.067  0.032  0.000 -0.032 -0.002  0.008  0.003  0.009 -0.009 -0.402 -0.015 -0.003 -0.006 -0.017      -0.020       0.071 -0.031           -0.036            -0.004            -0.004            0.003                                                                                 
    ## temptrend_abs.sc:npp.sc                     0.002 -0.057 -0.007   0.013  -0.116  0.095 -0.162 -0.164  0.013  0.000  0.011 -0.063 -0.034 -0.016 -0.044  0.367 -0.012  0.042       0.056       0.046 -0.168           -0.212            -0.032             0.042           -0.188             0.001                                                               
    ## temptrend_abs.sc:human.sc                  -0.069  0.174  0.065   0.062  -0.005 -0.009  0.005  0.006  0.009 -0.007 -0.009  0.000 -0.003  0.006 -0.009 -0.012  0.410 -0.165      -0.151      -0.030 -0.019            0.039             0.049            -0.026           -0.033            -0.004 -0.033                                                        
    ## REALMMarine:human.sc                        0.169  0.020 -0.160  -0.154  -0.020  0.004 -0.016 -0.024 -0.005 -0.004 -0.009 -0.006  0.013 -0.008 -0.014  0.018 -0.998 -0.020      -0.019       0.001 -0.001           -0.005            -0.009             0.003            0.010             0.005  0.008            -0.409                                      
    ## REALMTerrestrial:human.sc                   0.168  0.019 -0.159  -0.153  -0.032  0.007 -0.022 -0.006 -0.002 -0.001 -0.013 -0.011  0.018 -0.007 -0.020  0.027 -0.998 -0.020      -0.018      -0.001 -0.002            0.001            -0.008             0.002            0.013             0.004  0.012            -0.409            0.997                     
    ## temptrend_abs.sc:REALMMarine:human.sc       0.069 -0.172 -0.065  -0.062   0.001  0.009 -0.009 -0.009 -0.009  0.007  0.008  0.000  0.000 -0.006  0.006  0.008 -0.409  0.163       0.151       0.030  0.010           -0.038            -0.049             0.026            0.029             0.005  0.023            -0.998            0.410   0.408             
    ## temptrend_abs.sc:REALMTerrestrial:human.sc  0.069 -0.174 -0.065  -0.062   0.012  0.007 -0.003 -0.005 -0.008  0.007  0.008  0.002  0.004 -0.006  0.011  0.013 -0.407  0.164       0.150       0.023  0.016           -0.027            -0.045             0.021            0.034             0.006  0.029            -0.993            0.406   0.406   0.991     
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -7.59895515 -0.34513651  0.05160027  0.55556600  6.91372100 
    ## 
    ## Number of Observations: 43585
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  43585

#### Plot the coefficients

``` r
modTfull1simpreml <- update(modTfull1simp, method = 'REML')

# fails at insight::model_info() in plot_model()... ugh
# require(sjPlot)
# p1 <- plot_model(modTfull1simpreml, type = 'est', terms = 'temptrends_abs.sc') + ylim(-0.011, 0.02)
# p2<- plot_model(modTfull1simpreml, type = 'est', terms = 'REALM [Terrestrial,Marine]') + ylim(-0.0, 0.1)
# grid.arrange(p1, p2, nrow = 1)

coefs <- summary(modTfull1simpreml)$tTable
par(mfrow=c(1,2), las = 1, mai = c(0.5, 2, 0.1, 0.1))
rows1 <- which(!grepl('Intercept|REALM', rownames(coefs)))
plot(0,0, col = 'white', xlim=c(-0.02, 0.035), ylim = c(1,length(rows1)), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(rows1):1, labels = rownames(coefs)[rows1], cex.axis = 0.7)
abline(v = 0, col = 'grey')
for(i in 1:length(rows1)){
  x = coefs[rows1[i], 1]
  se = coefs[rows1[i], 2]
  points(x, length(rows1) + 1 - i, pch = 16)
  lines(x = c(x-se, x+se), y = c(length(rows1) + 1 - i, length(rows1) + 1 - i))
}

rows2 <- which(grepl('REALM', rownames(coefs)) & !grepl(':', rownames(coefs)))
plot(0,0, col = 'white', xlim=c(-0.0, 0.1), ylim = c(1,length(rows2)), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(rows2):1, labels = rownames(coefs)[rows2], cex.axis = 0.7)
abline(v = 0, col = 'grey')
for(i in 1:length(rows2)){
  x = coefs[rows2[i], 1]
  se = coefs[rows2[i], 2]
  points(x, length(rows2) + 1 - i, pch = 16)
  lines(x = c(x-se, x+se), y = c(length(rows2) + 1 - i, length(rows2) + 1 - i))
}
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20fullTmod1simp-1.png)<!-- -->

#### Plot residuals against each predictor

``` r
resids <- resid(modTfull1simpreml)
preds <- getData(modTfull1simpreml)
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

![](turnover_vs_temperature_MEmodels_files/figure-gfm/resids%20modTfull1simp-1.png)<!-- -->

## Sensitivity analysis: total turnover and Morisita-Horn models

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
                     temptrend_abs.sc*human.sc*REALM,
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
                     temptrend_abs.sc*human.sc*REALM,
                   random = randef, weights = varef, data = trends[i3,], method = 'REML')
  saveRDS(modTfullHorn, file = 'temp/modTfullHorn.rds')
}

summary(modTfullJbeta)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -131614.5 -131241.2 65850.23
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06574154 (Intr)
    ## temptrend_abs.sc 0.02529340 0.519 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.002266703 0.3555714
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.375398 
    ## Fixed effects: Jbetatrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tempave.sc +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      lifespan.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      endothermfrac.sc + temptrend_abs.sc * nspp.sc + temptrend.sc *      thermal_bias.sc + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      human.sc * REALM 
    ##                                                  Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                 0.04555299 0.016406397 43301   2.77654  0.0055
    ## temptrend_abs.sc                            0.03060307 0.010472914 43301   2.92212  0.0035
    ## REALMMarine                                 0.05972756 0.017600211   247   3.39357  0.0008
    ## REALMTerrestrial                            0.03638054 0.019004592   247   1.91430  0.0567
    ## tempave.sc                                 -0.00611321 0.000670313 43301  -9.11994  0.0000
    ## tempave_metab.sc                            0.01745540 0.001615230 43301  10.80676  0.0000
    ## seas.sc                                    -0.00266674 0.000323931 43301  -8.23245  0.0000
    ## microclim.sc                               -0.00060300 0.000165642 43301  -3.64040  0.0003
    ## mass.sc                                    -0.00436866 0.000742090 43301  -5.88697  0.0000
    ## speed.sc                                   -0.00073747 0.000496290 43301  -1.48596  0.1373
    ## lifespan.sc                                 0.00468815 0.001428641 43301   3.28154  0.0010
    ## consumerfrac.sc                             0.00164041 0.001025948 43301   1.59892  0.1098
    ## endothermfrac.sc                           -0.01699366 0.004557241 43301  -3.72893  0.0002
    ## nspp.sc                                    -0.02091884 0.000383981 43301 -54.47885  0.0000
    ## temptrend.sc                               -0.00011408 0.000473054 43301  -0.24115  0.8094
    ## thermal_bias.sc                            -0.00059755 0.000437607 43301  -1.36549  0.1721
    ## npp.sc                                      0.00243359 0.000267518 43301   9.09690  0.0000
    ## human.sc                                    0.00871132 0.004162944 43301   2.09259  0.0364
    ## temptrend_abs.sc:REALMMarine               -0.00442901 0.010865455 43301  -0.40762  0.6836
    ## temptrend_abs.sc:REALMTerrestrial           0.00254177 0.012350191 43301   0.20581  0.8369
    ## temptrend_abs.sc:tempave.sc                -0.00300186 0.000778780 43301  -3.85456  0.0001
    ## temptrend_abs.sc:tempave_metab.sc           0.00825935 0.002032677 43301   4.06329  0.0000
    ## temptrend_abs.sc:seas.sc                   -0.00138996 0.000497841 43301  -2.79197  0.0052
    ## temptrend_abs.sc:microclim.sc              -0.00103294 0.000246707 43301  -4.18691  0.0000
    ## temptrend_abs.sc:mass.sc                   -0.00459360 0.000765149 43301  -6.00354  0.0000
    ## temptrend_abs.sc:speed.sc                   0.00127432 0.000564231 43301   2.25851  0.0239
    ## temptrend_abs.sc:lifespan.sc                0.00531446 0.001469905 43301   3.61551  0.0003
    ## temptrend_abs.sc:consumerfrac.sc           -0.00037198 0.000768092 43301  -0.48429  0.6282
    ## temptrend_abs.sc:endothermfrac.sc          -0.00446006 0.003190560 43301  -1.39789  0.1622
    ## temptrend_abs.sc:nspp.sc                   -0.01443895 0.000448313 43301 -32.20728  0.0000
    ## temptrend.sc:thermal_bias.sc               -0.00064687 0.000343583 43301  -1.88272  0.0597
    ## temptrend_abs.sc:npp.sc                     0.00256923 0.000366002 43301   7.01971  0.0000
    ## temptrend_abs.sc:human.sc                   0.00845448 0.004689852 43301   1.80272  0.0714
    ## REALMMarine:human.sc                       -0.00826935 0.004168192 43301  -1.98392  0.0473
    ## REALMTerrestrial:human.sc                  -0.00918762 0.004168909 43301  -2.20384  0.0275
    ## temptrend_abs.sc:REALMMarine:human.sc      -0.00918357 0.004695218 43301  -1.95594  0.0505
    ## temptrend_abs.sc:REALMTerrestrial:human.sc -0.00999338 0.004714460 43301  -2.11973  0.0340
    ##  Correlation: 
    ##                                            (Intr) tmpt_. REALMMr REALMTr tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s tmptr. thrm_. npp.sc hmn.sc tm_.:REALMM tm_.:REALMT tmptrnd_bs.sc:t. t_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. tm.:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:h. REALMM: REALMT: t_.:REALMM:
    ## temptrend_abs.sc                            0.344                                                                                                                                                                                                                                                                                                                                                                                             
    ## REALMMarine                                -0.918 -0.315                                                                                                                                                                                                                                                                                                                                                                                      
    ## REALMTerrestrial                           -0.886 -0.314  0.812                                                                                                                                                                                                                                                                                                                                                                               
    ## tempave.sc                                 -0.002 -0.006 -0.009  -0.006                                                                                                                                                                                                                                                                                                                                                                       
    ## tempave_metab.sc                            0.042  0.013 -0.032  -0.036  -0.264                                                                                                                                                                                                                                                                                                                                                               
    ## seas.sc                                    -0.026 -0.002  0.038  -0.002   0.282 -0.083                                                                                                                                                                                                                                                                                                                                                        
    ## microclim.sc                               -0.004  0.002  0.010  -0.007   0.062  0.066  0.188                                                                                                                                                                                                                                                                                                                                                 
    ## mass.sc                                     0.018  0.002 -0.005  -0.003   0.000 -0.550 -0.041 -0.013                                                                                                                                                                                                                                                                                                                                          
    ## speed.sc                                    0.016  0.006 -0.025  -0.022  -0.005  0.180 -0.080 -0.037 -0.138                                                                                                                                                                                                                                                                                                                                   
    ## lifespan.sc                                 0.027  0.016 -0.017  -0.033  -0.036  0.757  0.095  0.049 -0.801  0.280                                                                                                                                                                                                                                                                                                                            
    ## consumerfrac.sc                            -0.050 -0.051  0.062   0.323  -0.034 -0.067  0.003  0.013  0.053 -0.123 -0.103                                                                                                                                                                                                                                                                                                                     
    ## endothermfrac.sc                            0.157  0.075 -0.079  -0.320   0.057 -0.125  0.044  0.004  0.023  0.061 -0.012 -0.427                                                                                                                                                                                                                                                                                                              
    ## nspp.sc                                    -0.020  0.001  0.000  -0.009  -0.048 -0.034 -0.026 -0.059 -0.150 -0.048  0.085  0.006  0.043                                                                                                                                                                                                                                                                                                       
    ## temptrend.sc                               -0.004 -0.008  0.005   0.007  -0.065 -0.019 -0.048 -0.017  0.015 -0.011 -0.012  0.012  0.003  0.020                                                                                                                                                                                                                                                                                                
    ## thermal_bias.sc                             0.022 -0.006 -0.033  -0.008   0.757 -0.033 -0.152 -0.127  0.038  0.018 -0.065 -0.006 -0.003 -0.112 -0.001                                                                                                                                                                                                                                                                                         
    ## npp.sc                                      0.003  0.001 -0.007  -0.007  -0.322  0.246 -0.265 -0.223 -0.013 -0.023  0.063 -0.059 -0.080 -0.165 -0.015 -0.135                                                                                                                                                                                                                                                                                  
    ## human.sc                                   -0.131 -0.012  0.123   0.116   0.020 -0.009  0.013  0.017  0.003  0.001  0.005  0.010  0.001 -0.009  0.009  0.015 -0.034                                                                                                                                                                                                                                                                           
    ## temptrend_abs.sc:REALMMarine               -0.327 -0.953  0.347   0.299   0.006 -0.008  0.010  0.000  0.004 -0.010 -0.012  0.054 -0.049 -0.007  0.006  0.005 -0.007  0.012                                                                                                                                                                                                                                                                    
    ## temptrend_abs.sc:REALMTerrestrial          -0.307 -0.878  0.282   0.346  -0.004 -0.014 -0.015 -0.009  0.007 -0.007 -0.017  0.134 -0.132 -0.013  0.008  0.002  0.008  0.010  0.833                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave.sc                -0.001 -0.006 -0.001  -0.003   0.108 -0.063  0.050 -0.026 -0.033 -0.030  0.003 -0.011  0.017 -0.003 -0.077  0.008 -0.033 -0.003  0.005       0.002                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:tempave_metab.sc           0.008  0.027 -0.005  -0.004   0.031  0.353  0.014  0.086 -0.218  0.088  0.319 -0.022 -0.053 -0.028  0.002  0.098  0.030 -0.001 -0.020      -0.024      -0.649                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:seas.sc                    0.010 -0.024 -0.005  -0.014   0.158  0.026  0.232  0.025 -0.036 -0.029  0.041  0.028  0.007 -0.010 -0.033  0.101 -0.180  0.000  0.044      -0.035       0.245           -0.011                                                                                                                                                                                                                    
    ## temptrend_abs.sc:microclim.sc              -0.002 -0.023  0.000  -0.001   0.010  0.101  0.054  0.277 -0.029  0.009  0.018  0.003 -0.028 -0.001  0.001 -0.022 -0.126  0.001  0.032      -0.003      -0.024           -0.001  0.059                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                    0.001  0.004  0.004   0.006  -0.013 -0.290 -0.034 -0.020  0.496 -0.054 -0.409  0.034  0.007 -0.060  0.019  0.015  0.005  0.011  0.009       0.014      -0.053           -0.407 -0.027             0.015                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                   0.000 -0.002 -0.001  -0.005  -0.006  0.038 -0.029  0.004 -0.023  0.222  0.070 -0.036  0.024 -0.025 -0.004  0.012  0.012  0.011 -0.020      -0.005      -0.044            0.090 -0.046            -0.047            -0.055                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                0.007  0.031 -0.004  -0.008  -0.002  0.409  0.050  0.025 -0.409  0.122  0.510 -0.047 -0.038  0.038 -0.016 -0.014  0.021 -0.005 -0.025      -0.031      -0.007            0.628  0.053             0.016            -0.790             0.150                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc           -0.049 -0.055  0.047   0.121  -0.007 -0.057 -0.003 -0.003  0.030 -0.053 -0.067  0.251 -0.134  0.005  0.008 -0.005 -0.019 -0.007  0.072       0.258       0.013           -0.090 -0.017             0.010             0.041            -0.188            -0.120                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc           0.064  0.145 -0.036  -0.120  -0.027 -0.072 -0.002 -0.035  0.002  0.009 -0.026 -0.158  0.375  0.042 -0.010 -0.060 -0.015  0.000 -0.099      -0.299       0.371           -0.373  0.017             0.001            -0.020             0.124            -0.033           -0.293                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                   -0.001  0.028 -0.004  -0.008  -0.037 -0.013  0.003  0.004 -0.066 -0.001  0.054 -0.002  0.016  0.333  0.022 -0.069 -0.028 -0.013 -0.044      -0.053       0.045           -0.023 -0.045            -0.177            -0.207            -0.147             0.172            0.022            0.095                                                                                                   
    ## temptrend.sc:thermal_bias.sc                0.005  0.011 -0.006  -0.004   0.088  0.036  0.002 -0.033 -0.007  0.009  0.007 -0.009 -0.003 -0.016 -0.395 -0.011 -0.004 -0.006 -0.012      -0.006       0.111           -0.004 -0.013            -0.042            -0.011             0.011             0.003           -0.004            0.010            -0.003                                                                                 
    ## temptrend_abs.sc:npp.sc                     0.003 -0.036 -0.007   0.001  -0.124  0.107 -0.162 -0.143  0.013  0.005  0.018 -0.026 -0.045 -0.042 -0.009 -0.044  0.362 -0.015  0.025       0.032      -0.164            0.149 -0.201            -0.206            -0.023            -0.005             0.046           -0.020           -0.074            -0.197            -0.013                                                               
    ## temptrend_abs.sc:human.sc                  -0.062  0.104  0.058   0.046  -0.003 -0.010  0.006  0.005  0.009  0.004 -0.007 -0.025  0.012 -0.002  0.008 -0.007 -0.014  0.439 -0.099      -0.084      -0.001           -0.028 -0.008             0.039             0.046             0.018            -0.033            0.015            0.005            -0.028            -0.005 -0.045                                                        
    ## REALMMarine:human.sc                        0.131  0.012 -0.123  -0.115  -0.022  0.009 -0.021 -0.023 -0.002  0.000 -0.006 -0.009 -0.001  0.007 -0.008 -0.014  0.024 -0.998 -0.012      -0.010       0.002            0.001 -0.003            -0.004            -0.010            -0.011             0.005            0.007            0.000             0.012             0.005  0.011            -0.438                                      
    ## REALMTerrestrial:human.sc                   0.130  0.012 -0.123  -0.116  -0.032  0.011 -0.026 -0.007  0.000  0.002 -0.009 -0.013 -0.004  0.012 -0.008 -0.020  0.033 -0.999 -0.012      -0.010       0.003           -0.001 -0.005             0.001            -0.009            -0.010             0.003            0.006            0.000             0.014             0.005  0.015            -0.438            0.997                     
    ## temptrend_abs.sc:REALMMarine:human.sc       0.062 -0.103 -0.057  -0.046   0.000  0.010 -0.009 -0.007 -0.008 -0.003  0.007  0.025 -0.012  0.000 -0.007  0.005  0.010 -0.438  0.098       0.084      -0.005            0.032 -0.002            -0.037            -0.045            -0.017             0.033           -0.015           -0.007             0.022             0.006  0.038            -0.998            0.439   0.437             
    ## temptrend_abs.sc:REALMTerrestrial:human.sc  0.062 -0.104 -0.057  -0.046   0.007  0.009 -0.005 -0.003 -0.007 -0.003  0.006  0.026 -0.012  0.003 -0.006  0.010  0.016 -0.437  0.098       0.083      -0.023            0.039 -0.002            -0.027            -0.041            -0.014             0.028           -0.019           -0.012             0.027             0.005  0.045            -0.994            0.436   0.435   0.993     
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -6.8341663 -0.2345669  0.1396972  0.6505144  7.5323730 
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
    ##       AIC       BIC  logLik
    ##   -121401 -121028.7 60743.5
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06615737 (Intr)
    ## temptrend_abs.sc 0.02037063 0.557 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:  0.00442701 0.3182813
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.239152 
    ## Fixed effects: Horntrend ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tempave.sc +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      lifespan.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      endothermfrac.sc + temptrend_abs.sc * nspp.sc + temptrend.sc *      thermal_bias.sc + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      human.sc * REALM 
    ##                                                  Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                 0.03630111 0.016537186 42334   2.19512  0.0282
    ## temptrend_abs.sc                            0.02566124 0.009133866 42334   2.80946  0.0050
    ## REALMMarine                                 0.05779454 0.017961469   215   3.21770  0.0015
    ## REALMTerrestrial                            0.02820967 0.019458153   215   1.44976  0.1486
    ## tempave.sc                                 -0.00094889 0.000878233 42334  -1.08045  0.2799
    ## tempave_metab.sc                            0.00886991 0.002018731 42334   4.39381  0.0000
    ## seas.sc                                    -0.00074525 0.000438496 42334  -1.69956  0.0892
    ## microclim.sc                                0.00048662 0.000224756 42334   2.16510  0.0304
    ## mass.sc                                    -0.00341515 0.000879926 42334  -3.88118  0.0001
    ## speed.sc                                    0.00015063 0.000610684 42334   0.24666  0.8052
    ## lifespan.sc                                 0.00265407 0.001722861 42334   1.54050  0.1234
    ## consumerfrac.sc                             0.00155926 0.001162095 42334   1.34176  0.1797
    ## endothermfrac.sc                           -0.01391882 0.004986453 42334  -2.79133  0.0053
    ## nspp.sc                                    -0.01859915 0.000471489 42334 -39.44772  0.0000
    ## temptrend.sc                               -0.00035296 0.000508826 42334  -0.69367  0.4879
    ## thermal_bias.sc                            -0.00006384 0.000543004 42334  -0.11757  0.9064
    ## npp.sc                                      0.00218807 0.000362539 42334   6.03542  0.0000
    ## human.sc                                    0.00767713 0.004658957 42334   1.64782  0.0994
    ## temptrend_abs.sc:REALMMarine               -0.00118061 0.009512313 42334  -0.12411  0.9012
    ## temptrend_abs.sc:REALMTerrestrial          -0.00170601 0.010934635 42334  -0.15602  0.8760
    ## temptrend_abs.sc:tempave.sc                -0.00203056 0.000943077 42334  -2.15312  0.0313
    ## temptrend_abs.sc:tempave_metab.sc           0.01085011 0.002461522 42334   4.40789  0.0000
    ## temptrend_abs.sc:seas.sc                    0.00127544 0.000630698 42334   2.02227  0.0432
    ## temptrend_abs.sc:microclim.sc              -0.00005515 0.000315399 42334  -0.17485  0.8612
    ## temptrend_abs.sc:mass.sc                   -0.00395703 0.000905259 42334  -4.37116  0.0000
    ## temptrend_abs.sc:speed.sc                   0.00036505 0.000639893 42334   0.57049  0.5683
    ## temptrend_abs.sc:lifespan.sc                0.00672732 0.001760210 42334   3.82189  0.0001
    ## temptrend_abs.sc:consumerfrac.sc           -0.00006437 0.000790691 42334  -0.08140  0.9351
    ## temptrend_abs.sc:endothermfrac.sc          -0.00704634 0.003145965 42334  -2.23980  0.0251
    ## temptrend_abs.sc:nspp.sc                   -0.01351167 0.000531750 42334 -25.40981  0.0000
    ## temptrend.sc:thermal_bias.sc               -0.00032558 0.000364962 42334  -0.89211  0.3723
    ## temptrend_abs.sc:npp.sc                     0.00196009 0.000467162 42334   4.19574  0.0000
    ## temptrend_abs.sc:human.sc                   0.00191049 0.005111880 42334   0.37374  0.7086
    ## REALMMarine:human.sc                       -0.00689867 0.004666934 42334  -1.47820  0.1394
    ## REALMTerrestrial:human.sc                  -0.00784904 0.004670654 42334  -1.68050  0.0929
    ## temptrend_abs.sc:REALMMarine:human.sc      -0.00258119 0.005118777 42334  -0.50426  0.6141
    ## temptrend_abs.sc:REALMTerrestrial:human.sc -0.00217049 0.005154709 42334  -0.42107  0.6737
    ##  Correlation: 
    ##                                            (Intr) tmpt_. REALMMr REALMTr tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s tmptr. thrm_. npp.sc hmn.sc tm_.:REALMM tm_.:REALMT tmptrnd_bs.sc:t. t_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:c. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. tm.:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:h. REALMM: REALMT: t_.:REALMM:
    ## temptrend_abs.sc                            0.336                                                                                                                                                                                                                                                                                                                                                                                             
    ## REALMMarine                                -0.905 -0.302                                                                                                                                                                                                                                                                                                                                                                                      
    ## REALMTerrestrial                           -0.876 -0.303  0.786                                                                                                                                                                                                                                                                                                                                                                               
    ## tempave.sc                                 -0.002 -0.005 -0.008  -0.003                                                                                                                                                                                                                                                                                                                                                                       
    ## tempave_metab.sc                            0.050  0.014 -0.037  -0.047  -0.337                                                                                                                                                                                                                                                                                                                                                               
    ## seas.sc                                    -0.039 -0.003  0.051   0.003   0.258 -0.087                                                                                                                                                                                                                                                                                                                                                        
    ## microclim.sc                               -0.007  0.005  0.012  -0.008   0.027  0.057  0.138                                                                                                                                                                                                                                                                                                                                                 
    ## mass.sc                                     0.017 -0.001 -0.004   0.002  -0.010 -0.525 -0.046 -0.009                                                                                                                                                                                                                                                                                                                                          
    ## speed.sc                                    0.017  0.005 -0.030  -0.026  -0.017  0.165 -0.089 -0.041 -0.133                                                                                                                                                                                                                                                                                                                                   
    ## lifespan.sc                                 0.034  0.021 -0.021  -0.040  -0.035  0.730  0.092  0.037 -0.802  0.256                                                                                                                                                                                                                                                                                                                            
    ## consumerfrac.sc                            -0.051 -0.051  0.046   0.335  -0.014 -0.066 -0.003  0.010  0.052 -0.129 -0.101                                                                                                                                                                                                                                                                                                                     
    ## endothermfrac.sc                            0.163  0.075 -0.079  -0.340   0.100 -0.159  0.042  0.000  0.013  0.065 -0.012 -0.421                                                                                                                                                                                                                                                                                                              
    ## nspp.sc                                    -0.023 -0.001  0.000  -0.008  -0.044 -0.021 -0.028 -0.044 -0.163 -0.043  0.100  0.010  0.049                                                                                                                                                                                                                                                                                                       
    ## temptrend.sc                               -0.004 -0.007  0.005   0.008  -0.060 -0.016 -0.043 -0.013  0.016 -0.008 -0.012  0.010  0.001  0.018                                                                                                                                                                                                                                                                                                
    ## thermal_bias.sc                             0.029 -0.004 -0.039  -0.013   0.726 -0.049 -0.151 -0.144  0.037  0.016 -0.067  0.007  0.006 -0.107 -0.005                                                                                                                                                                                                                                                                                         
    ## npp.sc                                      0.007 -0.002 -0.008  -0.007  -0.313  0.241 -0.273 -0.259 -0.010 -0.006  0.056 -0.040 -0.077 -0.154 -0.016 -0.130                                                                                                                                                                                                                                                                                  
    ## human.sc                                   -0.144  0.005  0.134   0.127   0.019 -0.009  0.014  0.022  0.005  0.000  0.006  0.013  0.000 -0.012  0.009  0.014 -0.038                                                                                                                                                                                                                                                                           
    ## temptrend_abs.sc:REALMMarine               -0.316 -0.951  0.339   0.285   0.004 -0.007  0.010 -0.003  0.006 -0.008 -0.014  0.052 -0.050 -0.005  0.005  0.003 -0.005 -0.005                                                                                                                                                                                                                                                                    
    ## temptrend_abs.sc:REALMTerrestrial          -0.297 -0.869  0.268   0.338  -0.006 -0.021 -0.018 -0.014  0.014 -0.007 -0.023  0.123 -0.125 -0.012  0.009 -0.003  0.013 -0.004  0.820                                                                                                                                                                                                                                                             
    ## temptrend_abs.sc:tempave.sc                 0.001  0.002 -0.003  -0.001   0.087 -0.056  0.006 -0.035 -0.038 -0.029  0.004 -0.005  0.024 -0.007 -0.068  0.028 -0.014 -0.003 -0.005       0.004                                                                                                                                                                                                                                                 
    ## temptrend_abs.sc:tempave_metab.sc           0.007  0.026 -0.001  -0.008   0.036  0.301  0.030  0.087 -0.190  0.063  0.281 -0.023 -0.054 -0.023  0.000  0.084  0.012  0.002 -0.015      -0.032      -0.660                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:seas.sc                    0.013 -0.019 -0.011  -0.020   0.117  0.031  0.181  0.007 -0.032 -0.024  0.037  0.014 -0.007 -0.015 -0.025  0.091 -0.158 -0.003  0.044      -0.075       0.135            0.055                                                                                                                                                                                                                    
    ## temptrend_abs.sc:microclim.sc              -0.003 -0.033  0.001  -0.003  -0.007  0.098  0.043  0.248 -0.026  0.017  0.012 -0.002 -0.032  0.016  0.009 -0.034 -0.134  0.000  0.042      -0.002      -0.061            0.014  0.037                                                                                                                                                                                                             
    ## temptrend_abs.sc:mass.sc                   -0.001 -0.005  0.006   0.009  -0.025 -0.249 -0.036 -0.018  0.443 -0.034 -0.368  0.028  0.002 -0.051  0.021  0.007  0.009  0.008  0.016       0.032      -0.064           -0.408 -0.024             0.016                                                                                                                                                                                           
    ## temptrend_abs.sc:speed.sc                  -0.002 -0.006  0.003  -0.003  -0.008  0.019 -0.019  0.010 -0.006  0.120  0.036 -0.026  0.018 -0.022 -0.005  0.009  0.016  0.011 -0.026      -0.005      -0.046            0.092 -0.039            -0.043            -0.059                                                                                                                                                                         
    ## temptrend_abs.sc:lifespan.sc                0.008  0.037  0.000  -0.009   0.006  0.350  0.046  0.016 -0.368  0.078  0.456 -0.043 -0.034  0.036 -0.016 -0.008  0.012 -0.001 -0.031      -0.038      -0.006            0.620  0.051             0.011            -0.804             0.148                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc           -0.046 -0.059  0.048   0.104  -0.006 -0.057 -0.013 -0.006  0.035 -0.034 -0.069  0.225 -0.103  0.005  0.010 -0.003 -0.007 -0.009  0.064       0.275       0.024           -0.107 -0.060             0.011             0.068            -0.205            -0.136                                                                                                                                     
    ## temptrend_abs.sc:endothermfrac.sc           0.060  0.137 -0.034  -0.109  -0.026 -0.075 -0.017 -0.047 -0.006 -0.004 -0.027 -0.124  0.329  0.043 -0.009 -0.054 -0.002  0.000 -0.106      -0.292       0.472           -0.471 -0.024            -0.016            -0.047             0.148            -0.040           -0.271                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                    0.001  0.047 -0.006  -0.008  -0.038 -0.004 -0.008  0.016 -0.059  0.005  0.054 -0.003  0.014  0.281  0.019 -0.061 -0.003 -0.011 -0.063      -0.077       0.057           -0.018 -0.029            -0.171            -0.211            -0.132             0.185            0.008            0.116                                                                                                   
    ## temptrend.sc:thermal_bias.sc                0.005  0.011 -0.007  -0.006   0.077  0.032  0.011 -0.034 -0.009  0.007  0.010 -0.007 -0.003 -0.015 -0.414 -0.001 -0.002 -0.006 -0.011      -0.008       0.087            0.009  0.006            -0.044            -0.013             0.012             0.006           -0.006            0.002            -0.003                                                                                 
    ## temptrend_abs.sc:npp.sc                     0.002 -0.060 -0.006   0.003  -0.099  0.091 -0.139 -0.147  0.013  0.012  0.015 -0.012 -0.037 -0.015 -0.009 -0.037  0.335 -0.012  0.045       0.057      -0.150            0.133 -0.230            -0.238            -0.023            -0.007             0.042            0.000           -0.078            -0.196            -0.013                                                               
    ## temptrend_abs.sc:human.sc                  -0.059  0.200  0.056   0.040  -0.005 -0.007  0.006  0.005  0.007  0.003 -0.005 -0.035  0.016 -0.002  0.007 -0.010 -0.012  0.391 -0.191      -0.160      -0.005           -0.022 -0.018             0.042             0.046             0.020            -0.027            0.017            0.002            -0.031            -0.005 -0.040                                                        
    ## REALMMarine:human.sc                        0.144 -0.005 -0.134  -0.126  -0.021  0.008 -0.024 -0.029 -0.003  0.001 -0.008 -0.013  0.000  0.009 -0.008 -0.012  0.027 -0.998  0.004       0.004       0.003           -0.002  0.000            -0.003            -0.007            -0.011             0.001            0.009            0.000             0.010             0.005  0.008            -0.390                                      
    ## REALMTerrestrial:human.sc                   0.143 -0.005 -0.134  -0.127  -0.035  0.014 -0.028 -0.008 -0.001  0.004 -0.011 -0.017 -0.003  0.015 -0.008 -0.020  0.036 -0.997  0.004       0.004       0.004           -0.005 -0.002             0.003            -0.006            -0.010            -0.001            0.009            0.001             0.013             0.004  0.012            -0.389            0.996                     
    ## temptrend_abs.sc:REALMMarine:human.sc       0.059 -0.198 -0.055  -0.040   0.002  0.007 -0.008 -0.006 -0.007 -0.002  0.005  0.035 -0.015  0.000 -0.006  0.007  0.008 -0.390  0.189       0.160      -0.003            0.026  0.005            -0.039            -0.045            -0.019             0.027           -0.017           -0.005             0.024             0.005  0.032            -0.998            0.391   0.389             
    ## temptrend_abs.sc:REALMTerrestrial:human.sc  0.059 -0.199 -0.055  -0.039   0.010  0.006 -0.003 -0.004 -0.006 -0.002  0.005  0.038 -0.015  0.003 -0.005  0.011  0.014 -0.387  0.190       0.156      -0.021            0.036  0.016            -0.031            -0.042            -0.016             0.023           -0.025           -0.013             0.029             0.005  0.039            -0.991            0.386   0.385   0.989     
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.78006527 -0.36573591  0.04098243  0.55898457  6.73978237 
    ## 
    ## Number of Observations: 42586
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    218                  42586

``` r
# simplify
if(file.exists('temp/modTfullJbetasimp.rds')){
  modTfullJbetasimp <- readRDS('temp/modTfullJbetasimp.rds')
} else {
  require(MASS) # for stepAIC
  modTfullJbetaml <- update(modTfullJbeta, method = 'ML')
  modTfullJbetasimp <- stepAIC(modTfullJbetaml, direction = 'backward')
  saveRDS(modTfullJbetasimp, file = 'temp/modTfullJbetasimp.rds')
}

if(file.exists('temp/modTfullHornsimp.rds')){
  modTfullHornsimp <- readRDS('temp/modTfullHornsimp.rds')
} else {
  require(MASS) # for stepAIC
  modTfullHornml <- update(modTfullHorn, method = 'ML')
  modTfullHornsimp <- stepAIC(modTfullHornml, direction = 'backward')
  saveRDS(modTfullHornsimp, file = 'temp/modTfullHornsimp.rds')
}

summary(modTfullJbetasimp)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -132071.5 -131706.9 66077.77
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06485034 (Intr)
    ## temptrend_abs.sc 0.02424001 0.542 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.002184697 0.3556416
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.375443 
    ## Fixed effects: Jbetatrend ~ temptrend_abs.sc + REALM + tempave.sc + tempave_metab.sc +      seas.sc + microclim.sc + mass.sc + speed.sc + lifespan.sc +      consumerfrac.sc + endothermfrac.sc + nspp.sc + temptrend.sc +      thermal_bias.sc + npp.sc + human.sc + temptrend_abs.sc:REALM +      temptrend_abs.sc:tempave.sc + temptrend_abs.sc:tempave_metab.sc +      temptrend_abs.sc:seas.sc + temptrend_abs.sc:microclim.sc +      temptrend_abs.sc:mass.sc + temptrend_abs.sc:speed.sc + temptrend_abs.sc:lifespan.sc +      temptrend_abs.sc:endothermfrac.sc + temptrend_abs.sc:nspp.sc +      temptrend.sc:thermal_bias.sc + temptrend_abs.sc:npp.sc +      temptrend_abs.sc:human.sc + REALM:human.sc + temptrend_abs.sc:REALM:human.sc 
    ##                                                  Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                 0.04499864 0.016157352 43302   2.78503  0.0054
    ## temptrend_abs.sc                            0.02982205 0.010058655 43302   2.96482  0.0030
    ## REALMMarine                                 0.06005560 0.017334232   247   3.46457  0.0006
    ## REALMTerrestrial                            0.03725759 0.018605543   247   2.00250  0.0463
    ## tempave.sc                                 -0.00607685 0.000667488 43302  -9.10406  0.0000
    ## tempave_metab.sc                            0.01735752 0.001608762 43302  10.78936  0.0000
    ## seas.sc                                    -0.00263548 0.000322099 43302  -8.18220  0.0000
    ## microclim.sc                               -0.00059669 0.000164821 43302  -3.62026  0.0003
    ## mass.sc                                    -0.00436676 0.000741161 43302  -5.89178  0.0000
    ## speed.sc                                   -0.00075888 0.000495051 43302  -1.53294  0.1253
    ## lifespan.sc                                 0.00462503 0.001423553 43302   3.24893  0.0012
    ## consumerfrac.sc                             0.00173575 0.000978662 43302   1.77359  0.0761
    ## endothermfrac.sc                           -0.01724298 0.004480859 43302  -3.84814  0.0001
    ## nspp.sc                                    -0.02086332 0.000383313 43302 -54.42893  0.0000
    ## temptrend.sc                               -0.00011671 0.000473108 43302  -0.24668  0.8052
    ## thermal_bias.sc                            -0.00059469 0.000436666 43302  -1.36189  0.1732
    ## npp.sc                                      0.00241963 0.000266048 43302   9.09473  0.0000
    ## human.sc                                    0.00872861 0.004156462 43302   2.10001  0.0357
    ## temptrend_abs.sc:REALMMarine               -0.00356687 0.010424068 43302  -0.34218  0.7322
    ## temptrend_abs.sc:REALMTerrestrial           0.00423260 0.011487792 43302   0.36844  0.7125
    ## temptrend_abs.sc:tempave.sc                -0.00299429 0.000775237 43302  -3.86242  0.0001
    ## temptrend_abs.sc:tempave_metab.sc           0.00815789 0.002016940 43302   4.04469  0.0001
    ## temptrend_abs.sc:seas.sc                   -0.00139542 0.000495156 43302  -2.81814  0.0048
    ## temptrend_abs.sc:microclim.sc              -0.00102866 0.000245890 43302  -4.18343  0.0000
    ## temptrend_abs.sc:mass.sc                   -0.00458888 0.000763766 43302  -6.00823  0.0000
    ## temptrend_abs.sc:speed.sc                   0.00122244 0.000551503 43302   2.21656  0.0267
    ## temptrend_abs.sc:lifespan.sc                0.00521989 0.001456312 43302   3.58432  0.0003
    ## temptrend_abs.sc:endothermfrac.sc          -0.00487181 0.002968801 43302  -1.64100  0.1008
    ## temptrend_abs.sc:nspp.sc                   -0.01439868 0.000447256 43302 -32.19336  0.0000
    ## temptrend.sc:thermal_bias.sc               -0.00064170 0.000343647 43302  -1.86733  0.0619
    ## temptrend_abs.sc:npp.sc                     0.00255762 0.000364834 43302   7.01036  0.0000
    ## temptrend_abs.sc:human.sc                   0.00856396 0.004673233 43302   1.83256  0.0669
    ## REALMMarine:human.sc                       -0.00829722 0.004161671 43302  -1.99372  0.0462
    ## REALMTerrestrial:human.sc                  -0.00920412 0.004162340 43302  -2.21129  0.0270
    ## temptrend_abs.sc:REALMMarine:human.sc      -0.00928940 0.004678524 43302  -1.98554  0.0471
    ## temptrend_abs.sc:REALMTerrestrial:human.sc -0.01010289 0.004697204 43302  -2.15083  0.0315
    ##  Correlation: 
    ##                                            (Intr) tmpt_. REALMMr REALMTr tmpv.s tmpv_. ses.sc mcrcl. mss.sc spd.sc lfspn. cnsmr. endth. nspp.s tmptr. thrm_. npp.sc hmn.sc tm_.:REALMM tm_.:REALMT tmptrnd_bs.sc:t. t_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:l. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. tm.:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:h. REALMM: REALMT: t_.:REALMM:
    ## temptrend_abs.sc                            0.356                                                                                                                                                                                                                                                                                                                                                                            
    ## REALMMarine                                -0.917 -0.326                                                                                                                                                                                                                                                                                                                                                                     
    ## REALMTerrestrial                           -0.887 -0.321  0.813                                                                                                                                                                                                                                                                                                                                                              
    ## tempave.sc                                 -0.002 -0.007 -0.009  -0.005                                                                                                                                                                                                                                                                                                                                                      
    ## tempave_metab.sc                            0.039  0.011 -0.030  -0.029  -0.262                                                                                                                                                                                                                                                                                                                                              
    ## seas.sc                                    -0.027 -0.003  0.039  -0.002   0.282 -0.083                                                                                                                                                                                                                                                                                                                                       
    ## microclim.sc                               -0.004  0.002  0.010  -0.007   0.063  0.067  0.189                                                                                                                                                                                                                                                                                                                                
    ## mass.sc                                     0.019  0.004 -0.006  -0.007   0.001 -0.550 -0.041 -0.013                                                                                                                                                                                                                                                                                                                         
    ## speed.sc                                    0.013  0.003 -0.023  -0.016  -0.004  0.178 -0.080 -0.037 -0.137                                                                                                                                                                                                                                                                                                                  
    ## lifespan.sc                                 0.024  0.013 -0.014  -0.026  -0.037  0.757  0.095  0.049 -0.801  0.278                                                                                                                                                                                                                                                                                                           
    ## consumerfrac.sc                            -0.039 -0.038  0.052   0.304  -0.032 -0.055  0.005  0.014  0.047 -0.115 -0.090                                                                                                                                                                                                                                                                                                    
    ## endothermfrac.sc                            0.153  0.070 -0.074  -0.310   0.056 -0.135  0.044  0.004  0.027  0.055 -0.022 -0.411                                                                                                                                                                                                                                                                                             
    ## nspp.sc                                    -0.020  0.001 -0.001  -0.010  -0.048 -0.034 -0.026 -0.060 -0.150 -0.048  0.086  0.005  0.045                                                                                                                                                                                                                                                                                      
    ## temptrend.sc                               -0.004 -0.007  0.004   0.006  -0.065 -0.019 -0.048 -0.017  0.015 -0.010 -0.012  0.010  0.004  0.020                                                                                                                                                                                                                                                                               
    ## thermal_bias.sc                             0.022 -0.006 -0.033  -0.008   0.759 -0.033 -0.152 -0.126  0.038  0.018 -0.066 -0.005 -0.004 -0.112 -0.001                                                                                                                                                                                                                                                                        
    ## npp.sc                                      0.002 -0.001 -0.006  -0.004  -0.322  0.245 -0.264 -0.222 -0.012 -0.025  0.062 -0.056 -0.084 -0.166 -0.014 -0.135                                                                                                                                                                                                                                                                 
    ## human.sc                                   -0.133 -0.011  0.125   0.119   0.019 -0.009  0.013  0.017  0.003  0.001  0.005  0.012  0.000 -0.009  0.009  0.015 -0.034                                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:REALMMarine               -0.337 -0.953  0.358   0.303   0.006 -0.004  0.011  0.000  0.002 -0.006 -0.008  0.037 -0.040 -0.008  0.005  0.005 -0.005  0.011                                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:REALMTerrestrial          -0.317 -0.895  0.289   0.341  -0.003  0.001 -0.015 -0.009 -0.001  0.007  0.000  0.074 -0.104 -0.016  0.006  0.003  0.014  0.011  0.845                                                                                                                                                                                                                                            
    ## temptrend_abs.sc:tempave.sc                 0.000 -0.005 -0.002  -0.005   0.108 -0.063  0.051 -0.025 -0.033 -0.030  0.004 -0.015  0.019 -0.003 -0.077  0.008 -0.034 -0.003  0.004      -0.002                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tempave_metab.sc           0.004  0.022  0.000   0.008   0.031  0.350  0.013  0.086 -0.217  0.084  0.316  0.001 -0.067 -0.028  0.002  0.098  0.029 -0.001 -0.013       0.000      -0.649                                                                                                                                                                                                                    
    ## temptrend_abs.sc:seas.sc                    0.009 -0.025 -0.005  -0.012   0.159  0.025  0.234  0.026 -0.035 -0.031  0.040  0.034  0.005 -0.009 -0.033  0.101 -0.181 -0.001  0.045      -0.034       0.248           -0.014                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc              -0.002 -0.023  0.000  -0.002   0.011  0.101  0.055  0.278 -0.029  0.009  0.019  0.001 -0.027 -0.001  0.000 -0.022 -0.126  0.001  0.033      -0.005      -0.023            0.000  0.060                                                                                                                                                                                            
    ## temptrend_abs.sc:mass.sc                    0.003  0.005  0.002   0.001  -0.012 -0.289 -0.034 -0.021  0.496 -0.052 -0.408  0.024  0.013 -0.060  0.019  0.016  0.006  0.011  0.006       0.004      -0.054           -0.407 -0.027             0.015                                                                                                                                                                          
    ## temptrend_abs.sc:speed.sc                  -0.009 -0.014  0.008   0.019  -0.008  0.027 -0.030  0.004 -0.017  0.217  0.059  0.015 -0.002 -0.025 -0.003  0.011  0.008  0.011 -0.007       0.048      -0.042            0.074 -0.050            -0.046            -0.048                                                                                                                                                        
    ## temptrend_abs.sc:lifespan.sc                0.002  0.024  0.002   0.008  -0.003  0.406  0.050  0.025 -0.409  0.117  0.507 -0.016 -0.056  0.039 -0.015 -0.015  0.019 -0.006 -0.017       0.000      -0.005            0.624  0.052             0.017            -0.792             0.130                                                                                                                                      
    ## temptrend_abs.sc:endothermfrac.sc           0.053  0.134 -0.024  -0.091  -0.031 -0.095 -0.003 -0.039  0.012 -0.006 -0.048 -0.090  0.359  0.047 -0.008 -0.066 -0.022 -0.002 -0.082      -0.242       0.401           -0.430  0.014             0.003            -0.010             0.076            -0.073                                                                                                                    
    ## temptrend_abs.sc:nspp.sc                    0.000  0.031 -0.005  -0.011  -0.037 -0.012  0.004  0.003 -0.067  0.000  0.056 -0.008  0.020  0.334  0.022 -0.069 -0.028 -0.013 -0.048      -0.063       0.044           -0.020 -0.046            -0.177            -0.208            -0.145             0.177            0.108                                                                                                   
    ## temptrend.sc:thermal_bias.sc                0.005  0.011 -0.006  -0.004   0.088  0.036  0.002 -0.033 -0.006  0.009  0.007 -0.008 -0.004 -0.016 -0.395 -0.011 -0.004 -0.006 -0.012      -0.005       0.112           -0.004 -0.013            -0.042            -0.011             0.010             0.003            0.009            -0.003                                                                                 
    ## temptrend_abs.sc:npp.sc                     0.002 -0.040 -0.006   0.004  -0.125  0.106 -0.163 -0.143  0.014  0.004  0.017 -0.021 -0.049 -0.042 -0.008 -0.044  0.362 -0.015  0.029       0.040      -0.164            0.148 -0.200            -0.205            -0.022            -0.009             0.044           -0.085            -0.196            -0.013                                                               
    ## temptrend_abs.sc:human.sc                  -0.061  0.117  0.057   0.044  -0.003 -0.009  0.006  0.005  0.008  0.004 -0.006 -0.031  0.015 -0.002  0.008 -0.007 -0.013  0.440 -0.112      -0.102      -0.002           -0.026 -0.009             0.039             0.046             0.021            -0.031            0.009            -0.029            -0.005 -0.044                                                        
    ## REALMMarine:human.sc                        0.133  0.011 -0.125  -0.119  -0.022  0.009 -0.020 -0.023 -0.002  0.000 -0.006 -0.011  0.000  0.007 -0.008 -0.014  0.024 -0.998 -0.011      -0.011       0.002            0.002 -0.003            -0.004            -0.010            -0.010             0.006            0.002             0.012             0.005  0.011            -0.439                                      
    ## REALMTerrestrial:human.sc                   0.132  0.011 -0.125  -0.119  -0.032  0.011 -0.025 -0.007  0.000  0.003 -0.009 -0.015 -0.003  0.012 -0.008 -0.020  0.032 -0.999 -0.011      -0.011       0.003           -0.001 -0.004             0.001            -0.009            -0.010             0.004            0.002             0.014             0.005  0.015            -0.439            0.997                     
    ## temptrend_abs.sc:REALMMarine:human.sc       0.061 -0.116 -0.056  -0.044   0.000  0.009 -0.009 -0.007 -0.008 -0.004  0.006  0.031 -0.015  0.000 -0.006  0.005  0.010 -0.439  0.111       0.102      -0.004            0.030 -0.001            -0.037            -0.045            -0.021             0.031           -0.011             0.023             0.005  0.036            -0.998            0.440   0.438             
    ## temptrend_abs.sc:REALMTerrestrial:human.sc  0.060 -0.117 -0.056  -0.043   0.007  0.008 -0.005 -0.003 -0.007 -0.004  0.005  0.033 -0.015  0.003 -0.006  0.010  0.015 -0.437  0.111       0.101      -0.022            0.037 -0.002            -0.027            -0.041            -0.019             0.025           -0.019             0.028             0.005  0.043            -0.994            0.436   0.436   0.993     
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -6.8573436 -0.2341121  0.1401190  0.6510143  7.5263164 
    ## 
    ## Number of Observations: 43585
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  43585

``` r
summary(modTfullHornsimp)
```

    ## Linear mixed-effects model fit by maximum likelihood
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -121863.5 -121595.1 60962.77
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.06476501 (Intr)
    ## temptrend_abs.sc 0.01910837 0.582 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev: 0.004350937 0.3182529
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.238924 
    ## Fixed effects: Horntrend ~ temptrend_abs.sc + REALM + tempave.sc + tempave_metab.sc +      seas.sc + microclim.sc + mass.sc + lifespan.sc + endothermfrac.sc +      nspp.sc + npp.sc + human.sc + temptrend_abs.sc:tempave.sc +      temptrend_abs.sc:tempave_metab.sc + temptrend_abs.sc:seas.sc +      temptrend_abs.sc:mass.sc + temptrend_abs.sc:lifespan.sc +      temptrend_abs.sc:endothermfrac.sc + temptrend_abs.sc:nspp.sc +      temptrend_abs.sc:npp.sc + temptrend_abs.sc:human.sc + REALM:human.sc 
    ##                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                        0.03771688 0.015106577 42346   2.49672  0.0125
    ## temptrend_abs.sc                   0.02414876 0.002234205 42346  10.80866  0.0000
    ## REALMMarine                        0.05577564 0.016260433   215   3.43015  0.0007
    ## REALMTerrestrial                   0.01940087 0.016793098   215   1.15529  0.2493
    ## tempave.sc                        -0.00080677 0.000595208 42346  -1.35545  0.1753
    ## tempave_metab.sc                   0.00897650 0.001969741 42346   4.55720  0.0000
    ## seas.sc                           -0.00075315 0.000428779 42346  -1.75649  0.0790
    ## microclim.sc                       0.00047323 0.000214488 42346   2.20632  0.0274
    ## mass.sc                           -0.00344038 0.000869747 42346  -3.95561  0.0001
    ## lifespan.sc                        0.00267540 0.001652630 42346   1.61888  0.1055
    ## endothermfrac.sc                  -0.01136049 0.004431156 42346  -2.56378  0.0104
    ## nspp.sc                           -0.01858053 0.000467436 42346 -39.74990  0.0000
    ## npp.sc                             0.00216514 0.000353795 42346   6.11974  0.0000
    ## human.sc                           0.00653169 0.004258003 42346   1.53398  0.1250
    ## temptrend_abs.sc:tempave.sc       -0.00186809 0.000924142 42346  -2.02143  0.0432
    ## temptrend_abs.sc:tempave_metab.sc  0.01060019 0.002412658 42346   4.39357  0.0000
    ## temptrend_abs.sc:seas.sc           0.00122693 0.000608701 42346   2.01565  0.0438
    ## temptrend_abs.sc:mass.sc          -0.00398191 0.000898981 42346  -4.42936  0.0000
    ## temptrend_abs.sc:lifespan.sc       0.00660469 0.001719320 42346   3.84145  0.0001
    ## temptrend_abs.sc:endothermfrac.sc -0.00706338 0.002739933 42346  -2.57794  0.0099
    ## temptrend_abs.sc:nspp.sc          -0.01346559 0.000512652 42346 -26.26651  0.0000
    ## temptrend_abs.sc:npp.sc            0.00191164 0.000449475 42346   4.25305  0.0000
    ## temptrend_abs.sc:human.sc         -0.00056433 0.000316598 42346  -1.78248  0.0747
    ## REALMMarine:human.sc              -0.00573273 0.004264581 42346  -1.34427  0.1789
    ## REALMTerrestrial:human.sc         -0.00666778 0.004270004 42346  -1.56154  0.1184
    ##  Correlation: 
    ##                                   (Intr) tmpt_. REALMMr REALMTr tmpv.s tmpv_. ses.sc mcrcl. mss.sc lfspn. endth. nspp.s npp.sc hmn.sc tmptrnd_bs.sc:t. t_.:_. tmptrnd_bs.sc:s. tmptrnd_bs.sc:m. tmptrnd_bs.sc:l. tmptrnd_bs.sc:nd. tmptrnd_bs.sc:ns. tmptrnd_bs.sc:np. tmptrnd_bs.sc:h. REALMM:
    ## temptrend_abs.sc                   0.115                                                                                                                                                                                                                                                       
    ## REALMMarine                       -0.895  0.059                                                                                                                                                                                                                                                
    ## REALMTerrestrial                  -0.906 -0.043  0.815                                                                                                                                                                                                                                         
    ## tempave.sc                        -0.038 -0.010  0.032   0.022                                                                                                                                                                                                                                 
    ## tempave_metab.sc                   0.047  0.033 -0.030  -0.028  -0.453                                                                                                                                                                                                                         
    ## seas.sc                           -0.036 -0.002  0.043   0.008   0.542 -0.089                                                                                                                                                                                                                  
    ## microclim.sc                      -0.008  0.000  0.012  -0.008   0.197  0.037  0.110                                                                                                                                                                                                           
    ## mass.sc                            0.026  0.021 -0.012  -0.018  -0.057 -0.513 -0.051 -0.004                                                                                                                                                                                                    
    ## lifespan.sc                        0.027  0.037 -0.008  -0.011   0.027  0.722  0.108  0.039 -0.802                                                                                                                                                                                             
    ## endothermfrac.sc                   0.150  0.049 -0.064  -0.225   0.147 -0.213  0.046  0.014  0.039 -0.065                                                                                                                                                                                      
    ## nspp.sc                           -0.024 -0.035 -0.004  -0.009   0.049 -0.020 -0.051 -0.068 -0.168  0.109  0.061                                                                                                                                                                               
    ## npp.sc                             0.012 -0.002 -0.013   0.000  -0.328  0.258 -0.301 -0.261 -0.008  0.049 -0.108 -0.172                                                                                                                                                                        
    ## human.sc                          -0.112 -0.023  0.102   0.100   0.014 -0.003  0.015  0.030  0.000  0.013  0.009 -0.010 -0.037                                                                                                                                                                 
    ## temptrend_abs.sc:tempave.sc        0.001 -0.013 -0.003   0.000   0.093 -0.049  0.013 -0.016 -0.044  0.015  0.026  0.000 -0.013 -0.006                                                                                                                                                          
    ## temptrend_abs.sc:tempave_metab.sc -0.007  0.056  0.016   0.009  -0.038  0.299  0.045  0.100 -0.188  0.280 -0.076 -0.012  0.022  0.013 -0.667                                                                                                                                                   
    ## temptrend_abs.sc:seas.sc           0.012 -0.035 -0.030  -0.007   0.071  0.035  0.188  0.009 -0.038  0.051 -0.023 -0.011 -0.145 -0.005  0.150            0.043                                                                                                                                  
    ## temptrend_abs.sc:mass.sc           0.014  0.043 -0.008  -0.013  -0.045 -0.247 -0.039 -0.021  0.444 -0.372  0.019 -0.054  0.013 -0.012 -0.065           -0.408 -0.023                                                                                                                           
    ## temptrend_abs.sc:lifespan.sc      -0.008  0.060  0.021   0.015   0.022  0.343  0.053  0.011 -0.366  0.454 -0.060  0.043  0.011  0.011  0.001            0.616  0.060           -0.808                                                                                                          
    ## temptrend_abs.sc:endothermfrac.sc  0.004 -0.041 -0.005   0.005   0.023 -0.107 -0.043 -0.069  0.005 -0.055  0.304  0.043 -0.006  0.006  0.551           -0.577 -0.083           -0.020           -0.099                                                                                         
    ## temptrend_abs.sc:nspp.sc          -0.030 -0.098  0.020   0.025   0.013  0.009 -0.014  0.057 -0.061  0.053  0.000  0.283 -0.036  0.003  0.053           -0.004 -0.036           -0.221            0.217            0.128                                                                        
    ## temptrend_abs.sc:npp.sc            0.023 -0.049 -0.028  -0.015  -0.105  0.117 -0.145 -0.101  0.011  0.010 -0.054 -0.018  0.308 -0.007 -0.163            0.145 -0.233           -0.023            0.051           -0.078            -0.258                                                      
    ## temptrend_abs.sc:human.sc         -0.013  0.008  0.020   0.012   0.001 -0.019 -0.022 -0.042  0.009 -0.006  0.019 -0.029 -0.039  0.029 -0.176            0.100 -0.175            0.034           -0.033           -0.076            -0.093            -0.093                                    
    ## REALMMarine:human.sc               0.113  0.024 -0.103  -0.100  -0.020  0.003 -0.025 -0.037  0.001 -0.015 -0.009  0.008  0.024 -0.998  0.004           -0.012  0.002            0.013           -0.011           -0.006            -0.004             0.002            -0.007                  
    ## REALMTerrestrial:human.sc          0.112  0.023 -0.102  -0.100  -0.033  0.007 -0.031 -0.016  0.005 -0.019 -0.014  0.013  0.034 -0.997  0.007           -0.015  0.001            0.014           -0.013           -0.005            -0.001             0.008            -0.033            0.995 
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.77663334 -0.36634606  0.04052091  0.55819364  6.71061228 
    ## 
    ## Number of Observations: 42586
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    218                  42586

``` r
if(file.exists('temp/modTfullJbetasimpreml.rds')){
  modTfullJbetasimpreml <- readRDS('temp/modTfullJbetasimpreml.rds')
} else {
  modTfullJbetasimpreml <- update(modTfullJbetasimp, method = 'REML')
  saveRDS(modTfullJbetasimpreml, file = 'temp/modTfullJbetasimpreml.rds')
}
if(file.exists('temp/modTfullHornsimpreml.rds')){
  modTfullHornsimpreml <- readRDS('temp/modTfullHornsimpreml.rds')
} else {
  modTfullHornsimpreml <- update(modTfullHornsimp, method = 'REML')
  saveRDS(modTfullHornsimpreml, file = 'temp/modTfullHornsimpreml.rds')
}

# plot coefs
coefs2 <- summary(modTfullJbetasimpreml)$tTable
coefs3 <- summary(modTfullHornsimpreml)$tTable
varstoplot <- unique(c(rownames(coefs2), rownames(coefs3)))

rows1 <- which(!grepl('Intercept|REALM', varstoplot) | grepl(':', varstoplot)) # vars to plot in first graph
rows1_2 <- which(rownames(coefs2) %in% varstoplot[rows1]) # rows in coefs2
rows1_3 <- which(rownames(coefs3) %in% varstoplot[rows1]) # rows in coefs3
xlims1 <- range(c(coefs2[rows1_2,1] - coefs2[rows1_2,2], 
                  coefs2[rows1_2,1] + coefs2[rows1_2,2], 
                  coefs3[rows1_3,1] - coefs3[rows1_3,2], 
                  coefs3[rows1_3,1] + coefs3[rows1_3,2]))

rows2 <- which(grepl('REALM', varstoplot) & !grepl(':', varstoplot)) # vars to plot in 2nd graph
rows2_2 <- which(rownames(coefs2) %in% varstoplot[rows2]) # rows in coefs2
rows2_3 <- which(rownames(coefs3) %in% varstoplot[rows2]) # rows in coefs3
xlims2 <- range(c(coefs2[rows2_2,1] - coefs2[rows2_2,2], 
                  coefs2[rows2_2,1] + coefs2[rows2_2,2], 
                  coefs3[rows2_3,1] - coefs3[rows2_3,2], 
                  coefs3[rows2_3,1] + coefs3[rows2_3,2]))

cols <- c('black', 'grey') # for Jbeta and Horn models, respectively
offs1 <- 0.1 # offset vertically for the two models
offs2 <- 0.01 # offset vertically for the two models (plot 2)

par(mfrow=c(1,2), las = 1, mai = c(0.5, 3, 0.1, 0.1))

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

plot(0,0, col = 'white', xlim=xlims2, ylim = c(0.9, length(rows2) + 0.1), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(rows2):1, labels = varstoplot[rows2], cex.axis = 0.7)
abline(v = 0, col = 'grey')
for(i in 1:length(rows2)){
  if(varstoplot[rows2[i]] %in% rownames(coefs2)){
    x = coefs2[rownames(coefs2) == varstoplot[rows2[i]], 1]
    se = coefs2[rownames(coefs2) == varstoplot[rows2[i]], 2]
    points(x, length(rows2) + 1 - i + offs2, pch = 16, col = cols[1])
    lines(x = c(x-se, x+se), y = c(length(rows2) + 1 - i + offs2, length(rows2) + 1 - i + offs2), col = cols[1])
  }
  if(varstoplot[rows2[i]] %in% rownames(coefs3)){
    x = coefs3[rownames(coefs3) == varstoplot[rows2[i]], 1]
    se = coefs3[rownames(coefs3) == varstoplot[rows2[i]], 2]
    points(x, length(rows2) + 1 - i - offs2, pch = 16, col = cols[2])
    lines(x = c(x-se, x+se), y = c(length(rows2) + 1 - i - offs2, length(rows2) + 1 - i - offs2), col = cols[2])
  }
}
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/LME%20Jacard%20total%20and%20MH%20models-1.png)<!-- -->

Black is for Jaccard total turnover (pres/abs), grey is for
Morisita-Horn turnover (considers abundance)

# To do

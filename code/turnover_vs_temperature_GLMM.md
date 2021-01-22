Temperature change and the community response to temperature change
across realms
================

# Prep

First run turnover\_vs\_temperature\_GLMM.R to fit the models.

<details>

<summary>Click to expand code</summary>

``` r
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(glmmTMB) # for ME models
library(bbmle) # for AICtab
```

    ## Loading required package: stats4

``` r
library(gridExtra) # to combine ggplots together
library(grid) # to combine ggplots together
library(RColorBrewer)
library(scales) # for defining custom scales in ggplot

options(width=500) # turn off most text wrapping

signedsqrt = function(x) sign(x)*sqrt(abs(x))
signedsq = function(x) sign(x) * x^2
signedsqrttrans <- trans_new(name = 'signedsqrt', transform = signedsqrt, inverse = signedsq)


# tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) 
```

``` r
# Turnover and covariates assembled by turnover_vs_temperature_prep.Rmd
trends <- fread('output/turnover_w_covariates.csv.gz')

# trim to only data with some temperature change
# important since sign of temperature change is a variable
trends[tempchange == 0, .N] # number to remove
```

    ## [1] 95

``` r
trends <- trends[tempchange != 0, ] # also removes any NA values

# set realm order
trends[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = FALSE)]

# set up sign of temperature change
trends[, tsign := factor(sign(tempchange))]

# realm that combines Terrestrial and Freshwater, for interacting with human impact
trends[, REALM2 := REALM]
levels(trends$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")

# group Marine invertebrates/plants in with All
trends[, taxa_mod2 := taxa_mod]
trends[taxa_mod == 'Marine invertebrates/plants', taxa_mod2 := 'All']

# calculate duration
trends[, duration := year2 - year1]

#add a comparison id
trends[, compID := paste0(rarefyID, '_', year1, '_', year2)]

# transformation for 2 categories. Eq. 1 in Douma & Weedon 2019
transform01 <- function(x) (x * (length(x) - 1) + 0.5) / (length(x))

trends[, Jtu.sc := transform01(Jtu)]
trends[, Jbeta.sc := transform01(Jbeta)]
trends[, Horn.sc := transform01(Horn)]

# Log-transform some variables, then center and scale. 
trends[, tempave.sc := scale(tempave)]
trends[, tempave_metab.sc := scale(tempave_metab)]
trends[, seas.sc := scale(seas)]
trends[, microclim.sc := scale(log(microclim))]
trends[, tempchange.sc := scale(tempchange, center = FALSE)] # do not center
trends[, tempchange_abs.sc := scale(abs(tempchange), center = FALSE)] # do not center, so that 0 is still 0 temperature change
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

</details>

# Summary stats

## Examine how many data points are available

### Just turnover and non-zero tempchange

``` r
cat('Overall # dissimilarities: ', nrow(trends), '\n')
```

    ## Overall # dissimilarities:  1288067

``` r
cat('# studies: ', trends[, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  330

``` r
cat('# timeseries: ', trends[, length(unique(rarefyID))], '\n')
```

    ## # timeseries:  52623

``` r
trends[!duplicated(rarefyID), table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##        1054       48195        3374

``` r
trends[!duplicated(rarefyID), table(taxa_mod)]
```

    ## taxa_mod
    ##                         All                  Amphibians                     Benthos                       Birds                        Fish               Invertebrates                     Mammals Marine invertebrates/plants                       Plant                    Reptiles 
    ##                        1715                         381                        4672                       13617                       28276                        2906                         543                         206                         303                           4

``` r
trends[!duplicated(rarefyID), table(taxa_mod, REALM)]
```

    ##                              REALM
    ## taxa_mod                      Freshwater Marine Terrestrial
    ##   All                                  0   1712           3
    ##   Amphibians                           2      0         379
    ##   Benthos                              0   4672           0
    ##   Birds                                0  10945        2672
    ##   Fish                              1037  27239           0
    ##   Invertebrates                       13   2813          80
    ##   Mammals                              0    496          47
    ##   Marine invertebrates/plants          0    206           0
    ##   Plant                                1    112         190
    ##   Reptiles                             1      0           3

### With all covariates

Use Bowler for human impact

``` r
# the cases we can compare
apply(trends[, .(Jtu, REALM, tempave_metab.sc, seas.sc, microclim.sc, tempchange.sc, mass.sc, speed.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##              Jtu            REALM tempave_metab.sc          seas.sc     microclim.sc    tempchange.sc          mass.sc         speed.sc  consumerfrac.sc endothermfrac.sc          nspp.sc  thermal_bias.sc           npp.sc           veg.sc  human_bowler.sc 
    ##          1288067          1288067          1288067          1288067          1274691          1288067          1285632          1285874          1288067          1288067          1277351          1178645          1284409          1276521          1288067

``` r
i <- trends[, complete.cases(Jtu, tempave_metab.sc, seas.sc, microclim.sc, tempchange.sc, mass.sc, speed.sc, consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
cat('Overall # dissimilarities: ', sum(i), '\n')
```

    ## Overall # dissimilarities:  1160127

``` r
cat('# studies: ', trends[i, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  231

``` r
cat('# timeseries: ', trends[i, length(unique(rarefyID))], '\n')
```

    ## # timeseries:  48863

``` r
trends[i & !duplicated(rarefyID), table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##         996       45083        2784

``` r
trends[i & !duplicated(rarefyID), table(taxa_mod)]
```

    ## taxa_mod
    ##                         All                  Amphibians                     Benthos                       Birds                        Fish               Invertebrates                     Mammals Marine invertebrates/plants                       Plant                    Reptiles 
    ##                        1683                          12                        4643                       11888                       27152                        2569                         533                         206                         175                           2

``` r
trends[i & !duplicated(rarefyID), table(taxa_mod, REALM)]
```

    ##                              REALM
    ## taxa_mod                      Freshwater Marine Terrestrial
    ##   All                                  0   1681           2
    ##   Amphibians                           2      0          10
    ##   Benthos                              0   4643           0
    ##   Birds                                0   9338        2550
    ##   Fish                               984  26168           0
    ##   Invertebrates                        9   2494          66
    ##   Mammals                              0    495          38
    ##   Marine invertebrates/plants          0    206           0
    ##   Plant                                1     58         116
    ##   Reptiles                             0      0           2

# Models

## Choose the variance structure

Try combinations of

  - beta errors
  - variance scaled to Nspp
  - random intercept for taxa\_mod2, STUDY\_ID, rarefyID
  - random slopes (abs temperature trend)
  - random intercept for compID (for overdispersion)?

And choose the one with lowest AIC

``` r
if(file.exists('temp/modRFgauss.rds')) modRFgauss <- readRDS('temp/modRFgauss.rds')
if(file.exists('temp/modRFbeta.rds')) modRFbeta <- readRDS('temp/modRFbeta.rds')
if(file.exists('temp/modRFrID.rds')) modRFrID <- readRDS('temp/modRFrID.rds')
if(file.exists('temp/modRFnestedRE.rds')) modRFnestedRE <- readRDS('temp/modRFnestedRE.rds')
if(file.exists('temp/modRFslopeRE.rds')) modRFslopeRE <- readRDS('temp/modRFslopeRE.rds')
if(file.exists('temp/modRFdisp.rds')) modRFdisp <- readRDS('temp/modRFdisp.rds')


AICtab(modRFgauss, modRFbeta, modRFrID, modRFnestedRE, modRFslopeRE, modRFdisp)
```

    ##               dAIC      df
    ## modRFdisp           0.0 21
    ## modRFslopeRE   193325.1 20
    ## modRFnestedRE  201811.7 14
    ## modRFrID       244070.8 12
    ## modRFbeta      795419.3 11
    ## modRFgauss    9037050.0 11

``` r
summary(modRFdisp)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jtu.sc ~ tempchange_abs.sc * REALM + tempchange_abs.sc * tempave_metab.sc +      tempchange_abs.sc * duration.sc + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:              ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -8326779 -8326526  4163410 -8326821  1277330 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance  Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.2589213 0.50884        
    ##                                tempchange_abs.sc 0.1295919 0.35999  -0.42 
    ##  STUDY_ID:taxa_mod2            (Intercept)       0.7962878 0.89235        
    ##                                tempchange_abs.sc 0.0207984 0.14422  -0.16 
    ##  taxa_mod2                     (Intercept)       0.0007081 0.02661        
    ##                                tempchange_abs.sc 0.0045132 0.06718  1.00  
    ## Number of obs: 1277351, groups:  rarefyID:(STUDY_ID:taxa_mod2), 52322; STUDY_ID:taxa_mod2, 305; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                        -0.417113   0.225285  -1.851   0.0641 .  
    ## tempchange_abs.sc                   0.200808   0.078579   2.555   0.0106 *  
    ## REALMMarine                         0.355581   0.238757   1.489   0.1364    
    ## REALMTerrestrial                   -0.184977   0.242119  -0.764   0.4449    
    ## tempave_metab.sc                    0.161907   0.015866  10.205  < 2e-16 ***
    ## duration.sc                         0.030800   0.001962  15.695  < 2e-16 ***
    ## tempchange_abs.sc:REALMMarine      -0.082965   0.072916  -1.138   0.2552    
    ## tempchange_abs.sc:REALMTerrestrial -0.204213   0.083687  -2.440   0.0147 *  
    ## tempchange_abs.sc:tempave_metab.sc  0.065653   0.015706   4.180 2.91e-05 ***
    ## tempchange_abs.sc:duration.sc      -0.004202   0.002018  -2.082   0.0374 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.404660   0.001124  -360.0   <2e-16 ***
    ## nspp.sc      0.573810   0.001143   502.2   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Chooses beta errors, random slopes (tempchange\_abs.sc) & intercepts for
taxa\_mod2, STUDY\_ID, rarefyID, and variance scaled to nspp. We haven’t
dealt with potential testing on the boundary issues here yet.

## Temperature-only models

### Load the models

``` r
# all years
modonlyTchangeJtu <- readRDS('temp/modonlyTchangeJtu.rds')
modonlyTchangeJbeta <- readRDS('temp/modonlyTchangeJbeta.rds')
modonlyTchangeHorn <- readRDS('temp/modonlyTchangeHorn.rds')

# 1 year
modonlyTchange1yrJtu <- readRDS('temp/modonlyTchange1yrJtu.rds')
modonlyTchange1yrJbeta <- readRDS('temp/modonlyTchange1yrJbeta.rds')
modonlyTchange1yrHorn <- readRDS('temp/modonlyTchange1yrHorn.rds')

# 10 year
modonlyTchange10yrJtu <- readRDS('temp/modonlyTchange10yrJtu.rds')
modonlyTchange10yrJbeta <- readRDS('temp/modonlyTchange10yrJbeta.rds')
modonlyTchange10yrHorn <- readRDS('temp/modonlyTchange10yrHorn.rds')
```

### Summary (all years)

``` r
summary(modonlyTchangeJtu)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jtu.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:              ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -8326073 -8325916  4163049 -8326099  1277338 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.261364 0.51124        
    ##                                tempchange_abs.sc 0.129587 0.35998  -0.41 
    ##  STUDY_ID:taxa_mod2            (Intercept)       0.801747 0.89540        
    ##                                tempchange_abs.sc 0.025212 0.15878  -0.07 
    ##  taxa_mod2                     (Intercept)       0.002956 0.05437        
    ##                                tempchange_abs.sc 0.006248 0.07905  0.33  
    ## Number of obs: 1277351, groups:  rarefyID:(STUDY_ID:taxa_mod2), 52322; STUDY_ID:taxa_mod2, 305; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -0.35819    0.06754  -5.304 1.13e-07 ***
    ## abs(tempchange)  0.09818    0.04939   1.988   0.0468 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.405277   0.001123  -360.7   <2e-16 ***
    ## nspp.sc      0.574250   0.001142   502.9   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modonlyTchangeJbeta)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jbeta.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:                ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -5446438 -5446281  2723232 -5446464  1277338 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.329225 0.57378        
    ##                                tempchange_abs.sc 0.076171 0.27599  -0.27 
    ##  STUDY_ID:taxa_mod2            (Intercept)       2.389237 1.54572        
    ##                                tempchange_abs.sc 0.045406 0.21309  -0.08 
    ##  taxa_mod2                     (Intercept)       0.084495 0.29068        
    ##                                tempchange_abs.sc 0.003781 0.06149  0.67  
    ## Number of obs: 1277351, groups:  rarefyID:(STUDY_ID:taxa_mod2), 52322; STUDY_ID:taxa_mod2, 305; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                 Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)      0.35838    0.15142   2.367   0.0179 *
    ## abs(tempchange)  0.07697    0.04600   1.673   0.0943 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 1.447742   0.001300  1113.5   <2e-16 ***
    ## nspp.sc     1.023839   0.001145   894.1   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modonlyTchangeHorn)
```

    ##  Family: beta  ( logit )
    ## Formula:          Horn.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:               ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -5466135 -5465979  2733081 -5466161  1248248 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance  Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.3432690 0.58589        
    ##                                tempchange_abs.sc 0.0474461 0.21782  -0.27 
    ##  STUDY_ID:taxa_mod2            (Intercept)       1.8234191 1.35034        
    ##                                tempchange_abs.sc 0.0502827 0.22424  -0.07 
    ##  taxa_mod2                     (Intercept)       0.0873317 0.29552        
    ##                                tempchange_abs.sc 0.0009186 0.03031  1.00  
    ## Number of obs: 1248261, groups:  rarefyID:(STUDY_ID:taxa_mod2), 51273; STUDY_ID:taxa_mod2, 285; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                 Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)      0.09895    0.14704   0.673   0.5010  
    ## abs(tempchange)  0.07443    0.03323   2.240   0.0251 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.515472   0.001276   403.9   <2e-16 ***
    ## nspp.sc     0.547730   0.001193   459.0   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Summary (1 year)

``` r
summary(modonlyTchange1yrJtu)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jtu.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:              ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##       AIC       BIC    logLik  deviance  df.resid 
    ## -941822.5 -941695.9  470924.3 -941848.5    125573 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.133034 0.36474        
    ##                                tempchange_abs.sc 0.001666 0.04082  1.00  
    ##  STUDY_ID:taxa_mod2            (Intercept)       0.819956 0.90551        
    ##                                tempchange_abs.sc 0.003414 0.05843  -0.57 
    ##  taxa_mod2                     (Intercept)       0.030364 0.17425        
    ##                                tempchange_abs.sc 0.001411 0.03756  -1.00 
    ## Number of obs: 125586, groups:  rarefyID:(STUDY_ID:taxa_mod2), 30311; STUDY_ID:taxa_mod2, 249; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -0.32225    0.09484  -3.398 0.000679 ***
    ## abs(tempchange) -0.01983    0.02709  -0.732 0.464284    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.459115   0.004016  -114.3   <2e-16 ***
    ## nspp.sc      0.417117   0.003327   125.4   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modonlyTchange1yrJbeta)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jbeta.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:                ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##       NA       NA       NA       NA   125573 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.260796 0.51068        
    ##                                tempchange_abs.sc 0.002705 0.05201  -1.00 
    ##  STUDY_ID:taxa_mod2            (Intercept)       2.380236 1.54280        
    ##                                tempchange_abs.sc 0.039148 0.19786  -0.24 
    ##  taxa_mod2                     (Intercept)       0.064088 0.25316        
    ##                                tempchange_abs.sc 0.001661 0.04075  -1.00 
    ## Number of obs: 125586, groups:  rarefyID:(STUDY_ID:taxa_mod2), 30311; STUDY_ID:taxa_mod2, 249; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                  Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)      0.265131   0.144215   1.838    0.066 .
    ## abs(tempchange) -0.003164   0.037369  -0.085    0.933  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 1.220390   0.004477   272.6   <2e-16 ***
    ## nspp.sc     0.822600   0.003526   233.3   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modonlyTchange1yrHorn)
```

    ##  Family: beta  ( logit )
    ## Formula:          Horn.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:               ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##       NA       NA       NA       NA   122716 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance  Std.Dev.  Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       1.491e-01 3.861e-01       
    ##                                tempchange_abs.sc 1.529e-14 1.236e-07 1.00  
    ##  STUDY_ID:taxa_mod2            (Intercept)       1.605e+00 1.267e+00       
    ##                                tempchange_abs.sc 3.229e-02 1.797e-01 -0.21 
    ##  taxa_mod2                     (Intercept)       1.123e-02 1.060e-01       
    ##                                tempchange_abs.sc 2.155e-10 1.468e-05 -0.96 
    ## Number of obs: 122729, groups:  rarefyID:(STUDY_ID:taxa_mod2), 29717; STUDY_ID:taxa_mod2, 232; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)     -0.01320    0.10047  -0.131    0.895
    ## abs(tempchange) -0.02805    0.03287  -0.853    0.393
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.259762   0.004446   58.43   <2e-16 ***
    ## nspp.sc     0.327831   0.003445   95.16   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Summary (10 year)

``` r
summary(modonlyTchange10yrJtu)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jtu.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:              ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##       NA       NA       NA       NA    49461 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance   Std.Dev.   Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)        3.063e-01  5.534e-01       
    ##                                tempchange_abs.sc  9.565e-03  9.780e-02 1.00  
    ##  STUDY_ID:taxa_mod2            (Intercept)        9.199e-01  9.591e-01       
    ##                                tempchange_abs.sc  4.203e-02  2.050e-01 -0.15 
    ##  taxa_mod2                     (Intercept)       1.315e-268 1.147e-134       
    ##                                tempchange_abs.sc 1.583e-184  1.258e-92 -1.00 
    ## Number of obs: 49474, groups:  rarefyID:(STUDY_ID:taxa_mod2), 17246; STUDY_ID:taxa_mod2, 156; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                  Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)     -0.282290   0.093422  -3.022  0.00251 **
    ## abs(tempchange)  0.009551   0.045092   0.212  0.83225   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.209160   0.007181  -29.13   <2e-16 ***
    ## nspp.sc      0.708102   0.006408  110.50   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modonlyTchange10yrJbeta)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jbeta.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:                ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##       NA       NA       NA       NA    49461 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance  Std.Dev.  Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       2.465e-01 4.965e-01       
    ##                                tempchange_abs.sc 5.922e-33 7.696e-17 1.00  
    ##  STUDY_ID:taxa_mod2            (Intercept)       2.797e+00 1.673e+00       
    ##                                tempchange_abs.sc 2.433e-02 1.560e-01 -0.11 
    ##  taxa_mod2                     (Intercept)       1.383e-02 1.176e-01       
    ##                                tempchange_abs.sc 2.546e-03 5.046e-02 -1.00 
    ## Number of obs: 49474, groups:  rarefyID:(STUDY_ID:taxa_mod2), 17246; STUDY_ID:taxa_mod2, 156; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)      0.52725    0.14975   3.521  0.00043 ***
    ## abs(tempchange)  0.02565    0.04397   0.583  0.55968    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 1.670844   0.007819   213.7   <2e-16 ***
    ## nspp.sc     1.141039   0.006250   182.6   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modonlyTchange10yrHorn)
```

    ##  Family: beta  ( logit )
    ## Formula:          Horn.sc ~ abs(tempchange) + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:               ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##       AIC       BIC    logLik  deviance  df.resid 
    ## -222426.0 -222311.8  111226.0 -222452.0     48429 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance  Std.Dev.  Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       1.007e-01 3.173e-01       
    ##                                tempchange_abs.sc 1.752e-11 4.186e-06 1.00  
    ##  STUDY_ID:taxa_mod2            (Intercept)       1.929e+00 1.389e+00       
    ##                                tempchange_abs.sc 1.997e-02 1.413e-01 -0.21 
    ##  taxa_mod2                     (Intercept)       2.773e-03 5.266e-02       
    ##                                tempchange_abs.sc 1.143e-03 3.380e-02 -0.99 
    ## Number of obs: 48442, groups:  rarefyID:(STUDY_ID:taxa_mod2), 17015; STUDY_ID:taxa_mod2, 148; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                 Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)      0.17022    0.12666   1.344    0.179
    ## abs(tempchange)  0.04021    0.04125   0.975    0.330
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.506410   0.007810   64.84   <2e-16 ***
    ## nspp.sc     0.608834   0.006213   97.99   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Plot the temp-only coefficients

``` r
# make table of coefficients
coefs1 <- as.data.frame(summary(modonlyTchangeJtu)$coefficients$cond)
coefs2 <- as.data.frame(summary(modonlyTchangeJbeta)$coefficients$cond)
coefs3 <- as.data.frame(summary(modonlyTchangeHorn)$coefficients$cond)
coefs4 <- as.data.frame(summary(modonlyTchange1yrJtu)$coefficients$cond)
coefs5 <- as.data.frame(summary(modonlyTchange1yrJbeta)$coefficients$cond)
coefs6 <- as.data.frame(summary(modonlyTchange1yrHorn)$coefficients$cond)
coefs7 <- as.data.frame(summary(modonlyTchange10yrJtu)$coefficients$cond)
coefs8 <- as.data.frame(summary(modonlyTchange10yrJbeta)$coefficients$cond)
coefs9 <- as.data.frame(summary(modonlyTchange10yrHorn)$coefficients$cond)
coefs1$response <- coefs4$response <- coefs7$response <- 'Jtu'
coefs2$response <- coefs5$response <- coefs8$response <- 'Jbeta'
coefs3$response <- coefs6$response <- coefs9$response <- 'Horn'
coefs1$fit <- coefs2$fit <- coefs3$fit <- 'all'
coefs4$fit <- coefs5$fit <- coefs6$fit <- '1yr'
coefs7$fit <- coefs8$fit <- coefs9$fit <- '10yr'
rows1 <- which(grepl('tempchange', rownames(coefs1))) # extract temperature effect
cols <- c('Estimate', 'Std. Error', 'response', 'fit')
allcoefs <- rbind(coefs1[rows1, cols], coefs2[rows1, cols], coefs3[rows1, cols],
                  coefs4[rows1, cols], coefs5[rows1, cols], coefs6[rows1, cols],
                  coefs7[rows1, cols], coefs8[rows1, cols], coefs9[rows1, cols])

allcoefs$lCI <- allcoefs$Estimate - 1.96 * allcoefs$`Std. Error` # lower confidence interval
allcoefs$uCI <- allcoefs$Estimate + 1.96 * allcoefs$`Std. Error`
allcoefs$y <- c(3, 2, 1) + rep(c(0, -0.1, -0.2), c(3, 3, 3)) # y-values
allcoefs$group <- paste0(allcoefs$response, allcoefs$fit)

pd <- position_dodge(width = 0.5)
ggplot(allcoefs, aes(x = response, y = Estimate, color = fit, group = group)) +
  geom_point(position = pd) + 
  geom_errorbar(aes(ymin=lCI, ymax=uCI), width=.1, position = pd) +
  labs(y = 'Dissimilarity per |°C change|') #+ 
```

![](turnover_vs_temperature_GLMM_files/figure-gfm/modonlyT%20coefs-1.png)<!-- -->

``` r
  #coord_cartesian(ylim = c(-0.01, 0.05))
```

CIs are 1.96\*SE

## Temperature\&duration models

### Load the models

``` r
# all years
modTDJtu <- readRDS('temp/modTDJtu.rds')
modTDJbeta <- readRDS('temp/modTDJbeta.rds')
modTDHorn <- readRDS('temp/modTDHorn.rds')
```

### Summary

``` r
summary(modTDJtu)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jtu.sc ~ abs(tempchange) * duration + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:              ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -8326075 -8325894  4163052 -8326105  1277336 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.261364 0.51124        
    ##                                tempchange_abs.sc 0.129587 0.35998  -0.41 
    ##  STUDY_ID:taxa_mod2            (Intercept)       0.802169 0.89564        
    ##                                tempchange_abs.sc 0.025021 0.15818  -0.08 
    ##  taxa_mod2                     (Intercept)       0.003199 0.05656        
    ##                                tempchange_abs.sc 0.006189 0.07867  0.30  
    ## Number of obs: 1277351, groups:  rarefyID:(STUDY_ID:taxa_mod2), 52322; STUDY_ID:taxa_mod2, 305; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                            Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -0.3614083  0.0678327  -5.328 9.93e-08 ***
    ## abs(tempchange)           0.1007895  0.0492199   2.048   0.0406 *  
    ## duration                  0.0005195  0.0002138   2.430   0.0151 *  
    ## abs(tempchange):duration -0.0004244  0.0002271  -1.868   0.0617 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.405234   0.001124  -360.6   <2e-16 ***
    ## nspp.sc      0.574185   0.001142   502.7   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modTDJbeta)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jbeta.sc ~ abs(tempchange) * duration + (tempchange_abs.sc |      taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:                ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -5472783 -5472602  2736406 -5472813  1277336 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.327830 0.57256        
    ##                                tempchange_abs.sc 0.079556 0.28206  -0.27 
    ##  STUDY_ID:taxa_mod2            (Intercept)       2.447754 1.56453        
    ##                                tempchange_abs.sc 0.033751 0.18371  -0.09 
    ##  taxa_mod2                     (Intercept)       0.098591 0.31399        
    ##                                tempchange_abs.sc 0.002793 0.05285  0.60  
    ## Number of obs: 1277351, groups:  rarefyID:(STUDY_ID:taxa_mod2), 52322; STUDY_ID:taxa_mod2, 305; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                            Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               0.2769920  0.1587945    1.74   0.0811 .  
    ## abs(tempchange)           0.0702443  0.0420753    1.67   0.0950 .  
    ## duration                  0.0141588  0.0001303  108.66   <2e-16 ***
    ## abs(tempchange):duration -0.0014920  0.0001311  -11.38   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 1.477618   0.001304  1133.3   <2e-16 ***
    ## nspp.sc     1.039118   0.001151   903.2   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modTDHorn)
```

    ##  Family: beta  ( logit )
    ## Formula:          Horn.sc ~ abs(tempchange) * duration + (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:               ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -5480772 -5480592  2740401 -5480802  1248246 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance  Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.3411517 0.5841         
    ##                                tempchange_abs.sc 0.0506017 0.2249   -0.27 
    ##  STUDY_ID:taxa_mod2            (Intercept)       1.8689780 1.3671         
    ##                                tempchange_abs.sc 0.0363508 0.1907   -0.08 
    ##  taxa_mod2                     (Intercept)       0.1030181 0.3210         
    ##                                tempchange_abs.sc 0.0003921 0.0198   0.99  
    ## Number of obs: 1248261, groups:  rarefyID:(STUDY_ID:taxa_mod2), 51273; STUDY_ID:taxa_mod2, 285; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                            Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               0.0088979  0.1556696    0.06   0.9544    
    ## abs(tempchange)           0.0671778  0.0283034    2.37   0.0176 *  
    ## duration                  0.0152950  0.0001873   81.66  < 2e-16 ***
    ## abs(tempchange):duration -0.0015492  0.0001949   -7.95 1.88e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.533141   0.001282   416.0   <2e-16 ***
    ## nspp.sc     0.553580   0.001195   463.4   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Plot the coefficients

``` r
coefs1 <- summary(modTDJtu)$coefficients$cond
coefs2 <- summary(modTDJbeta)$coefficients$cond
coefs3 <- summary(modTDHorn)$coefficients$cond

varstoplot <- unique(c(rownames(coefs1), rownames(coefs2), rownames(coefs3)))
varstoplot <- varstoplot[which(!grepl('Intercept', varstoplot) | grepl(':', varstoplot))] # vars to plot

rows1_1 <- which(rownames(coefs1) %in% varstoplot) # rows in coefs
rows1_2 <- which(rownames(coefs2) %in% varstoplot)
rows1_3 <- which(rownames(coefs3) %in% varstoplot)
xlims <- range(c(coefs1[rows1_1,1] - 1.96*coefs1[rows1_1,2], coefs1[rows1_1,1] + 1.96*coefs1[rows1_1,2], 
                  coefs2[rows1_2,1] - 1.96*coefs2[rows1_2,2], coefs2[rows1_2,1] + 1.96*coefs2[rows1_2,2], 
                  coefs3[rows1_3,1] - 1.96*coefs3[rows1_3,2], coefs3[rows1_3,1] + 1.96*coefs3[rows1_3,2]))

cols <- brewer.pal(3, 'Dark2') # for Jtu, Jbeta and Horn models
pchs <- c(16, 16, 16)
offs <- c(0.1, 0, -0.1) # offset vertically for each model

par(las = 1, mai = c(0.5, 4, 0.1, 0.1))

plot(0,0, col = 'white', xlim = xlims, ylim = c(0.9,length(varstoplot)+0.1), yaxt='n', xlab = '', ylab ='')
axis(2, at = length(varstoplot):1, labels = varstoplot, cex.axis = 0.7)
abline(v = 0, col = 'grey', lty = 2)
abline(h = 1:length(varstoplot), col = 'grey', lty = 3)
for(i in 1:length(varstoplot)){
  if(varstoplot[i] %in% rownames(coefs1)){
    x = coefs1[rownames(coefs1) == varstoplot[i], 1] # the coef in col 1
    se = coefs1[rownames(coefs1) == varstoplot[i], 2] # the SE in col 2
    points(x, length(varstoplot) + 1 - i + offs[1], pch = pchs[1], col = cols[1])
    lines(x = c(x-1.96*se, x+1.96*se), y = c(length(varstoplot) + 1 - i + offs[1], length(varstoplot) + 1 - i + offs[1]), col = cols[1])
  }
  if(varstoplot[i] %in% rownames(coefs2)){
    x = coefs2[rownames(coefs2) == varstoplot[i], 1]
    se = coefs2[rownames(coefs2) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[2], pch = pchs[2], col = cols[2])
    lines(x = c(x-1.96*se, x+1.96*se), y = c(length(varstoplot) + 1 - i + offs[2], length(varstoplot) + 1 - i + offs[2]), col = cols[2])
  }
  if(varstoplot[i] %in% rownames(coefs3)){
    x = coefs3[rownames(coefs3) == varstoplot[i], 1]
    se = coefs3[rownames(coefs3) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[3], pch = pchs[3], col = cols[3])
    lines(x = c(x-1.96*se, x+1.96*se), y = c(length(varstoplot) + 1 - i + offs[3], length(varstoplot) + 1 - i + offs[3]), col = cols[3])
  }
}
legend('bottomright', col = cols, pch = 16, lwd = 1, legend = c('Jtu', 'Jbeta', 'Horn'), cex = 0.5)
```

![](turnover_vs_temperature_GLMM_files/figure-gfm/modTD%20coefs-1.png)<!-- -->
abs(tempchange) is in degC duration is in years

### Plot the response

``` r
if(!file.exists('temp/modTDJtu_preds.rds')) {
  newdat <- expand.grid(duration = c(1, 50), tempchange = seq(0.01, 5, length.out = 5), 
                      taxa_mod2 = NA, STUDY_ID = NA, rarefyID = NA,
                      Nspp = 60) # no REs, ave number of species
  scl <- attr(trends$tempchange_abs.sc, 'scaled:scale') # scaling factor for tempchange
  newdat$tempchange_abs.sc <- abs(newdat$tempchange)/scl
  cent <- attr(trends$nspp.sc, 'scaled:center') 
  scl <- attr(trends$nspp.sc, 'scaled:scale') 
  newdat$nspp.sc <- (log(newdat$Nspp) - cent)/scl
  Jtu <- predict(modTDJtu, newdata = newdat, se.fit = TRUE, re.form = NA, type = 'response')
  newdat$durationfac <- as.factor(newdat$duration)
  newdat$Jtu <- Jtu$fit
  newdat$Jtu.se <- Jtu$se.fit
  
  saveRDS(newdat, file = 'temp/modTDJtu_preds.rds')
} else {
  newdat <- readRDS('temp/modTDJtu_preds.rds')
}

ggplot(newdat, aes(tempchange, Jtu, color = durationfac, group = durationfac)) +
  geom_line() +
  geom_ribbon(aes(ymin=Jtu - 1.96*Jtu.se, ymax=Jtu + 1.96*Jtu.se), alpha = 0.1, linetype = 'blank')
```

![](turnover_vs_temperature_GLMM_files/figure-gfm/plot%20TD%20response-1.png)<!-- -->

## Temperature and duration by REALM

### Load the models

``` r
if(file.exists('temp/modTDrealmJtu.rds')) modTDrealmJtu <- readRDS('temp/modTDrealmJtu.rds')
if(file.exists('temp/modTDrealmJbeta.rds')) modTDrealmJbeta <- readRDS('temp/modTDrealmJbeta.rds')
if(file.exists('temp/modTDrealmHorn.rds')) modTDrealmHorn <- readRDS('temp/modTDrealmHorn.rds')
```

### Summary

``` r
summary(modTDrealmJtu)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jtu.sc ~ abs(tempchange) * duration + abs(tempchange) * REALM +      (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:              ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -8326093 -8325864  4163066 -8326131  1277332 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance  Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       2.614e-01 0.511229       
    ##                                tempchange_abs.sc 1.296e-01 0.359965 -0.41 
    ##  STUDY_ID:taxa_mod2            (Intercept)       7.557e-01 0.869312       
    ##                                tempchange_abs.sc 2.114e-02 0.145403 -0.18 
    ##  taxa_mod2                     (Intercept)       9.803e-05 0.009901       
    ##                                tempchange_abs.sc 1.163e-02 0.107855 0.78  
    ## Number of obs: 1277351, groups:  rarefyID:(STUDY_ID:taxa_mod2), 52322; STUDY_ID:taxa_mod2, 305; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                                    Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)                      -0.5644194  0.2432059  -2.321   0.0203 *
    ## abs(tempchange)                   0.2669777  0.1108265   2.409   0.0160 *
    ## duration                          0.0005213  0.0002138   2.439   0.0147 *
    ## REALMMarine                       0.4503284  0.2350978   1.916   0.0554 .
    ## REALMTerrestrial                 -0.0138683  0.2672532  -0.052   0.9586  
    ## abs(tempchange):duration         -0.0004244  0.0002271  -1.869   0.0617 .
    ## abs(tempchange):REALMMarine      -0.0908671  0.0964062  -0.942   0.3459  
    ## abs(tempchange):REALMTerrestrial -0.2682967  0.1139343  -2.355   0.0185 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.405232   0.001124  -360.6   <2e-16 ***
    ## nspp.sc      0.574188   0.001142   502.7   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modTDrealmJbeta)
```

    ##  Family: beta  ( logit )
    ## Formula:          Jbeta.sc ~ abs(tempchange) * duration + abs(tempchange) * REALM +      (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:                ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -5472805 -5472576  2736422 -5472843  1277332 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.327827 0.57256        
    ##                                tempchange_abs.sc 0.079604 0.28214  -0.27 
    ##  STUDY_ID:taxa_mod2            (Intercept)       2.292589 1.51413        
    ##                                tempchange_abs.sc 0.029300 0.17117  -0.16 
    ##  taxa_mod2                     (Intercept)       0.010721 0.10354        
    ##                                tempchange_abs.sc 0.006834 0.08267  1.00  
    ## Number of obs: 1277351, groups:  rarefyID:(STUDY_ID:taxa_mod2), 52322; STUDY_ID:taxa_mod2, 305; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                                    Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                       0.1513262  0.3594605    0.42   0.6738    
    ## abs(tempchange)                   0.2056933  0.1026803    2.00   0.0452 *  
    ## duration                          0.0141593  0.0001303  108.66   <2e-16 ***
    ## REALMMarine                       0.6669576  0.3767038    1.77   0.0766 .  
    ## REALMTerrestrial                 -0.3896767  0.3807822   -1.02   0.3061    
    ## abs(tempchange):duration         -0.0014921  0.0001311  -11.38   <2e-16 ***
    ## abs(tempchange):REALMMarine      -0.0718025  0.0933696   -0.77   0.4419    
    ## abs(tempchange):REALMTerrestrial -0.2053292  0.1069832   -1.92   0.0550 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 1.477622   0.001304  1133.3   <2e-16 ***
    ## nspp.sc     1.039121   0.001151   903.2   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(modTDrealmHorn)
```

    ##  Family: beta  ( logit )
    ## Formula:          Horn.sc ~ abs(tempchange) * duration + abs(tempchange) * REALM +      (tempchange_abs.sc | taxa_mod2/STUDY_ID/rarefyID)
    ## Dispersion:               ~nspp.sc
    ## Data: trends[i, ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -5480787 -5480558  2740412 -5480825  1248242 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                        Name              Variance Std.Dev. Corr  
    ##  rarefyID:(STUDY_ID:taxa_mod2) (Intercept)       0.341110 0.58405        
    ##                                tempchange_abs.sc 0.050586 0.22491  -0.27 
    ##  STUDY_ID:taxa_mod2            (Intercept)       1.764545 1.32836        
    ##                                tempchange_abs.sc 0.030492 0.17462  -0.13 
    ##  taxa_mod2                     (Intercept)       0.016944 0.13017        
    ##                                tempchange_abs.sc 0.004509 0.06715  1.00  
    ## Number of obs: 1248261, groups:  rarefyID:(STUDY_ID:taxa_mod2), 51273; STUDY_ID:taxa_mod2, 285; taxa_mod2, 9
    ## 
    ## Conditional model:
    ##                                    Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                       0.0552813  0.3224437    0.17   0.8639    
    ## abs(tempchange)                   0.1823237  0.0966884    1.89   0.0593 .  
    ## duration                          0.0152958  0.0001873   81.67  < 2e-16 ***
    ## REALMMarine                       0.4299906  0.3401971    1.26   0.2062    
    ## REALMTerrestrial                 -0.5214365  0.3415231   -1.53   0.1268    
    ## abs(tempchange):duration         -0.0015498  0.0001949   -7.95 1.84e-15 ***
    ## abs(tempchange):REALMMarine      -0.0589580  0.0909988   -0.65   0.5171    
    ## abs(tempchange):REALMTerrestrial -0.1524211  0.1033472   -1.47   0.1403    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.533129   0.001282   416.0   <2e-16 ***
    ## nspp.sc     0.553562   0.001194   463.4   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Plot the temp coefficients from TD realm models

``` r
coefs1 <- summary(modTDrealmJtu)$coefficients$cond
coefs2 <- summary(modTDrealmJbeta)$coefficients$cond
coefs3 <- summary(modTDrealmHorn)$coefficients$cond

varstoplot <- unique(c(rownames(coefs1), rownames(coefs2), rownames(coefs3)))
varstoplot <- varstoplot[which(!grepl('Intercept', varstoplot) | grepl(':', varstoplot))] # vars to plot

rows1_1 <- which(rownames(coefs1) %in% varstoplot) # rows in coefs
rows1_2 <- which(rownames(coefs2) %in% varstoplot)
rows1_3 <- which(rownames(coefs3) %in% varstoplot)
xlims <- range(c(coefs1[rows1_1,1] - 1.96*coefs1[rows1_1,2], coefs1[rows1_1,1] + 1.96*coefs1[rows1_1,2], 
                  coefs2[rows1_2,1] - 1.96*coefs2[rows1_2,2], coefs2[rows1_2,1] + 1.96*coefs2[rows1_2,2], 
                  coefs3[rows1_3,1] - 1.96*coefs3[rows1_3,2], coefs3[rows1_3,1] + 1.96*coefs3[rows1_3,2]))

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
    lines(x = c(x-1.96*se, x+1.96*se), y = c(length(varstoplot) + 1 - i + offs[1], length(varstoplot) + 1 - i + offs[1]), col = cols[1])
  }
  if(varstoplot[i] %in% rownames(coefs2)){
    x = coefs2[rownames(coefs2) == varstoplot[i], 1]
    se = coefs2[rownames(coefs2) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[2], pch = pchs[2], col = cols[2])
    lines(x = c(x-1.96*se, x+1.96*se), y = c(length(varstoplot) + 1 - i + offs[2], length(varstoplot) + 1 - i + offs[2]), col = cols[2])
  }
  if(varstoplot[i] %in% rownames(coefs3)){
    x = coefs3[rownames(coefs3) == varstoplot[i], 1]
    se = coefs3[rownames(coefs3) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[3], pch = pchs[3], col = cols[3])
    lines(x = c(x-1.96*se, x+1.96*se), y = c(length(varstoplot) + 1 - i + offs[3], length(varstoplot) + 1 - i + offs[3]), col = cols[3])
  }
}
legend('topleft', col = cols, pch = 16, lwd = 1, legend = c('Jtu', 'Jbeta', 'Horn'), cex = 0.5)
```

![](turnover_vs_temperature_GLMM_files/figure-gfm/modTD%20realm%20coefs-1.png)<!-- -->

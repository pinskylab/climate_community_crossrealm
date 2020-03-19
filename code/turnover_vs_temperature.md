Turnover vs. temperature exploration
================

``` r
library(data.table)
library(ggplot2)
library(lme4)
```

    ## Loading required package: Matrix

``` r
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.
```

### Load data

``` r
# BioTime
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin)

# Temperature trends from CRU TS (on land)
tmptrendbyyr <- fread('output/cruts_tmptrend_byyear.csv.gz', drop = 1)
tmptrendallyr <- fread('output/cruts_tmptrend_allyear.csv.gz', drop = 1)

# Temperature trends from ERSST (ocean)
ssttrendbyyr <- fread('output/ersst_ssttrend_byyear.csv.gz', drop = 1)
ssttrendallyr <- fread('output/ersst_ssttrend_allyear.csv.gz', drop = 1)

# Seasonality
seascru <- fread('output/cruts_seas.csv.gz', drop = 1)
seasersst <- fread('output/ersst_seas.csv.gz', drop = 1)

# Average temperature
tavecru <- fread('output/cruts_tmpave_allyear.csv.gz', drop = 1)
taveersst <- fread('output/ersst_sstave_allyear.csv.gz', drop = 1)

# NPP
npp <- fread('output/npplandocean.csv.gz')

# Body size
bs <- fread('output/mass_byrarefyid.csv.gz', drop = 1)
bs[, ':='(STUDY_ID = NULL, REALM = NULL, taxa_mod = NULL)] # remove unnecessary columns 
```

### Assemble dataset of trends (turnover and temperature) and covariates

``` r
# calculate temporal turnover
calctrend <- function(y, YEAR, nm = 'y'){ # function to calc trends
    mod <- lm(y ~ YEAR)
    out <- list(y = coef(mod)[2], # coef for the slope
         y_se = sqrt(diag(vcov(mod)))[2]) # SE
    names(out) <- c(nm, paste0(nm, '_se'))
    return(out)
}

setkey(bt, STUDY_ID, rarefyID, YEAR)
trends <- bt[, calctrend(Jtu_base, YEAR, 'Jtutrend'), 
    by = .(REALM, Biome, taxa_mod, STUDY_ID, rarefyID, rarefyID_x, rarefyID_y)] # calculate trend in Jaccard turnover from first year, plus SEs
trends2 <- bt[, calctrend(Jbeta_base, YEAR, 'Jbetatrend'), 
    by = .(rarefyID)] # calculate trend in total Jaccard' beta diversity's from first year, 
```

    ## Warning in summary.lm(object, ...): essentially perfect fit: summary may be
    ## unreliable
    
    ## Warning in summary.lm(object, ...): essentially perfect fit: summary may be
    ## unreliable

``` r
trends3 <- bt[, calctrend(1-Horn_base, YEAR, 'Horntrend'), 
    by = .(rarefyID)] # calculate trend in Horn-Morisita from first year. Convert to dissimilarity.
trends4 <- bt[, .(Strend = coef(lm(I(log(S)) ~ YEAR))[2]), by = .(rarefyID)] # trend in log(S)
nyrBT <-  bt[, .(nyrBT = length(YEAR)), by = .(rarefyID)] # number of years in time-series

trends <- merge(trends, trends2) # merge in total J and Horn-Morisita
trends <- merge(trends, trends3)
trends <- merge(trends, trends4)
trends <- merge(trends, nyrBT)
rm(trends2, trends3, trends4, nyrBT)
```

Find a minimum trend SE to use in cases where
missing

``` r
trends[Jtutrend_se>0, hist(log10(Jtutrend_se))] # plot of SE on the trends
```

![](turnover_vs_temperature_files/figure-gfm/find%20minse-1.png)<!-- -->

    ## $breaks
    ##  [1] -17 -16 -15 -14 -13 -12 -11 -10  -9  -8  -7  -6  -5  -4  -3  -2  -1
    ## [18]   0
    ## 
    ## $counts
    ##  [1]     1     5    20    21     2     0     0     0     0     0     0
    ## [12]     0     1   135  9017 23995  2747
    ## 
    ## $density
    ##  [1] 2.782105e-05 1.391053e-04 5.564211e-04 5.842422e-04 5.564211e-05
    ##  [6] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
    ## [11] 0.000000e+00 0.000000e+00 2.782105e-05 3.755842e-03 2.508625e-01
    ## [16] 6.675662e-01 7.642444e-02
    ## 
    ## $mids
    ##  [1] -16.5 -15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5
    ## [12]  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5
    ## 
    ## $xname
    ## [1] "log10(Jtutrend_se)"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
trends[Jtutrend_se>0, plot(nyrBT, Jtutrend_se, log = 'y')]
```

![](turnover_vs_temperature_files/figure-gfm/find%20minse-2.png)<!-- -->

    ## NULL

``` r
fillse <- trends[Jtutrend_se > 1e-10, mean(Jtutrend_se)] # find smallest non-zero se, avoiding the perfect fits
trends[, Jtutrend_se_nz := Jtutrend_se] # create an SE that is not zero or NA
trends[is.na(Jtutrend_se_nz), Jtutrend_se_nz := fillse] 
trends[Jtutrend_se_nz < 1e-10, Jtutrend_se_nz := fillse] 

trends[, hist(log10(1/(Jtutrend_se_nz^2)))] # inverse variance, used for weighting
```

![](turnover_vs_temperature_files/figure-gfm/find%20minse-3.png)<!-- -->

    ## $breaks
    ##  [1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0
    ## [18] 8.5
    ## 
    ## $counts
    ##  [1]    30    91   746  1884  3741 23351  7649  6822  5102  2635   808
    ## [12]   472    93    33     7     2     1
    ## 
    ## $density
    ##  [1] 1.122188e-03 3.403969e-03 2.790506e-02 7.047338e-02 1.399368e-01
    ##  [6] 8.734734e-01 2.861204e-01 2.551854e-01 1.908467e-01 9.856547e-02
    ## [11] 3.022425e-02 1.765575e-02 3.478781e-03 1.234406e-03 2.618438e-04
    ## [16] 7.481250e-05 3.740625e-05
    ## 
    ## $mids
    ##  [1] 0.25 0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75
    ## [15] 7.25 7.75 8.25
    ## 
    ## $xname
    ## [1] "log10(1/(Jtutrend_se_nz^2))"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
trends[, plot(nyrBT, 1/(Jtutrend_se_nz^2), log = 'y')]
```

![](turnover_vs_temperature_files/figure-gfm/find%20minse-4.png)<!-- -->

    ## NULL

``` r
# add temperature trends
trends <- merge(trends, tmptrendbyyr, all.x = TRUE, by = 'rarefyID') # add temperature trend (years that match biotime)
trends <- merge(trends, tmptrendallyr[, .(rarefyID, tmptrendallyr)], all.x = TRUE, by = 'rarefyID') # add temperature trend (all years from first to last of each time-series)
trends <- merge(trends, ssttrendbyyr[, .(rarefyID, ssttrendbyyr)], all.x = TRUE, by = 'rarefyID') # sst
trends <- merge(trends, ssttrendallyr[, .(rarefyID, ssttrendallyr)], all.x = TRUE, by = 'rarefyID')

# add covariates
trends <- merge(trends, seascru[, .(rarefyID, seas_cruts)], all.x = TRUE, by = 'rarefyID') # seasonality from CRU TS
trends <- merge(trends, seasersst[, .(rarefyID, seas_ersst)], all.x = TRUE, by = 'rarefyID') # seasonality from ERSST
trends <- merge(trends, tavecru[, .(rarefyID, tmpaveallyr)], all.x = TRUE, by = 'rarefyID')
trends <- merge(trends, taveersst[, .(rarefyID, sstaveallyr)], all.x = TRUE, by = 'rarefyID')
trends <- merge(trends, npp, all.x = TRUE, by = 'rarefyID') # npp, centered and scale
trends <- merge(trends, bs, all.x = TRUE) # npp, centered and scale
```

Do some basic checks

``` r
# basic checks
trends
```

    ##           rarefyID  REALM                      Biome      taxa_mod
    ##     1:  100_606491 Marine     Northern_European_Seas          Fish
    ##     2:  101_606491 Marine     Northern_European_Seas Invertebrates
    ##     3: 108_3933165 Marine Continental_High_Antarctic         Birds
    ##     4: 108_3941181 Marine Continental_High_Antarctic         Birds
    ##     5: 108_3941182 Marine Continental_High_Antarctic         Birds
    ##    ---                                                            
    ## 53463:  99_4377155 Marine Northwest_Australian_Shelf          Fish
    ## 53464:  99_4383724 Marine Northwest_Australian_Shelf          Fish
    ## 53465:  99_4386651 Marine Northwest_Australian_Shelf          Fish
    ## 53466:  99_4390299 Marine Northwest_Australian_Shelf          Fish
    ## 53467:  99_4394671 Marine Northwest_Australian_Shelf          Fish
    ##        STUDY_ID rarefyID_x rarefyID_y     Jtutrend Jtutrend_se  Jbetatrend
    ##     1:      100    -3.0800   51.14000  0.004763952 0.001601483 0.004886186
    ##     2:      101    -3.0800   51.14000  0.000000000 0.000000000 0.005573571
    ##     3:      108    57.9650  -65.28500  0.071428571         NaN 0.071428571
    ##     4:      108    59.9275  -66.29250  0.100000000         NaN 0.107142857
    ##     5:      108    59.9700  -66.19500  0.062500000         NaN 0.089285714
    ##    ---                                                                    
    ## 53463:       99   116.8010  -19.84227 -0.010268562 0.020521929 0.043500339
    ## 53464:       99   117.5515  -19.61040  0.000000000         NaN 0.000000000
    ## 53465:       99   117.8600  -18.80625  0.000000000         NaN 0.250000000
    ## 53466:       99   118.3213  -18.79645  0.088888889 0.071860741 0.092380952
    ## 53467:       99   118.8110  -19.35420  0.000000000         NaN 0.080000000
    ##        Jbetatrend_se    Horntrend Horntrend_se       Strend nyrBT
    ##     1:   0.001622134 0.0031009875 0.0015680649  0.001734836    31
    ##     2:   0.002121344 0.0004638415 0.0008664591  0.010650816    31
    ##     3:           NaN 0.0714285714          NaN -0.078472306     2
    ##     4:           NaN 0.1147144190          NaN  0.063853203     2
    ##     5:           NaN 0.0878275862          NaN -0.086643398     2
    ##    ---                                                           
    ## 53463:   0.027264849           NA           NA -0.123358425     3
    ## 53464:           NaN           NA           NA  0.000000000     2
    ## 53465:           NaN           NA           NA -0.346573590     2
    ## 53466:   0.072315133           NA           NA  0.042711681     4
    ## 53467:           NaN           NA           NA -0.102165125     2
    ##        Jtutrend_se_nz tmptrendbyyr YEAR_min YEAR_max tmptrendallyr
    ##     1:    0.001601483    0.0252789     1981     2011     0.0252789
    ##     2:    0.036812583    0.0252789     1981     2011     0.0252789
    ##     3:    0.036812583           NA       NA       NA            NA
    ##     4:    0.036812583           NA       NA       NA            NA
    ##     5:    0.036812583           NA       NA       NA            NA
    ##    ---                                                            
    ## 53463:    0.020521929           NA       NA       NA            NA
    ## 53464:    0.036812583           NA       NA       NA            NA
    ## 53465:    0.036812583           NA       NA       NA            NA
    ## 53466:    0.071860741           NA       NA       NA            NA
    ## 53467:    0.036812583           NA       NA       NA            NA
    ##        ssttrendbyyr ssttrendallyr seas_cruts seas_ersst tmpaveallyr
    ##     1:  0.041039992   0.041039992   4.308796  3.0779262     10.2379
    ##     2:  0.041039992   0.041039992   4.308796  3.0779262     10.2379
    ##     3: -0.013489248  -0.006945833         NA  0.4781110          NA
    ##     4: -0.009487191  -0.004529737         NA  0.5235794          NA
    ##     5: -0.009487191  -0.004529737         NA  0.5235794          NA
    ##    ---                                                             
    ## 53463: -0.006082084  -0.003552321         NA  2.0016941          NA
    ## 53464:           NA            NA         NA         NA          NA
    ## 53465:  0.620516857   0.620516857         NA  1.8848515          NA
    ## 53466: -0.018872044  -0.002693142         NA  1.8848515          NA
    ## 53467:           NA            NA         NA         NA          NA
    ##        sstaveallyr       npp   mass_mean     mass_sd mass_geomean
    ##     1:   12.050894 1685.9467  7875.25419 19167.24413   479.944935
    ##     2:   12.050894 1685.9467    16.10904    39.50364     1.520777
    ##     3:   -1.320041  126.3856  1375.26250  1226.77458  1038.456233
    ##     4:   -1.271692  152.3814   507.55429   422.07819   307.942688
    ##     5:   -1.271692  150.4591  1524.17714  1894.15945   617.583072
    ##    ---                                                           
    ## 53463:   26.938956  552.9129 13765.15579 29554.30442  3279.363673
    ## 53464:          NA  543.8439  4447.13533  9074.39246  1333.573456
    ## 53465:   27.370296  366.2237  4594.27542 12629.71915   349.741494
    ## 53466:   27.549504  384.6634 15454.40947 48499.40539  1155.380521
    ## 53467:          NA  694.7807  5858.67208  9466.09115  2221.402814
    ##        mass_geosd nspp nspp_wdata
    ##     1:  27.865437   83         76
    ##     2:   8.097404   15          9
    ##     3:   2.345711    4          4
    ##     4:   3.640893    7          7
    ##     5:   5.469016    7          7
    ##    ---                           
    ## 53463:   6.375602   37         19
    ## 53464:   4.394475   23         15
    ## 53465:  47.531678   46         21
    ## 53466:  19.294515   50         30
    ## 53467:   4.467042   36         24

``` r
trends[, .(minJtu = min(Jtutrend), maxJtu = max(Jtutrend), minJbe = min(Jbetatrend), maxJbe = max(Jbetatrend), 
           minHo = min(Horntrend, na.rm = TRUE), maxHo = max(Horntrend, na.rm = TRUE), minS = min(Strend), maxS = max(Strend)), by = REALM]
```

    ##          REALM      minJtu maxJtu      minJbe maxJbe       minHo maxHo
    ## 1:      Marine -0.10000000    0.5 -0.08000000    0.5 -0.09421855   0.5
    ## 2: Terrestrial -0.10000000    0.5 -0.07142857    0.5 -0.09373434   0.5
    ## 3:  Freshwater -0.06666667    0.5 -0.04166667    0.5 -0.06129032   0.5
    ##         minS     maxS
    ## 1: -2.124248 1.812170
    ## 2: -1.522261 1.510212
    ## 3: -1.039721 0.804719

``` r
trends[, .(nJtu = sum(Jtutrend < 0), nS = sum(Strend < 0)), by = REALM] # why are some turnover trends < 0?
```

    ##          REALM nJtu    nS
    ## 1:      Marine 2891 20037
    ## 2: Terrestrial  223  1513
    ## 3:  Freshwater   71   462

``` r
trends[, .(nJ = sum(Jtutrend > 0), nS = sum(Strend > 0)), by = REALM]
```

    ##          REALM    nJ    nS
    ## 1:      Marine 39219 25505
    ## 2: Terrestrial  2854  1768
    ## 3:  Freshwater   693   524

``` r
# are turnover calculations correlated?
ggplot(trends, aes(Jbetatrend, Jtutrend)) +
    geom_point() +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](turnover_vs_temperature_files/figure-gfm/basic%20graphs%20of%20turnover-1.png)<!-- -->

``` r
ggplot(trends, aes(Jbetatrend, Horntrend)) +
    geom_point() +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 1451 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1451 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/basic%20graphs%20of%20turnover-2.png)<!-- -->

``` r
# are the temperature trends correlated?
# compare trends calculated from all sampled years in BioTime vs. all years from min to max year in BioTime
ggplot(trends, aes(tmptrendallyr, tmptrendbyyr)) +
    geom_point() +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 40885 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 40885 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/basic%20graphs%20of%20temperatures-1.png)<!-- -->

``` r
ggplot(trends, aes(ssttrendallyr, ssttrendbyyr)) +
    geom_point() +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 5199 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5199 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/basic%20graphs%20of%20temperatures-2.png)<!-- -->

``` r
# compare temperature trends from CRU TS vs ERSST
ggplot(trends, aes(tmptrendallyr, ssttrendbyyr)) + # ERSST more muted at fastest rates of change in CRU TS
    geom_point() +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 45240 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 45240 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/basic%20graphs%20of%20temperatures-3.png)<!-- -->

``` r
# average temperature from CRU TS and ERSST
ggplot(trends, aes(tmpaveallyr, sstaveallyr, color = rarefyID_y)) + # ERSST more muted at fastest rates of change in CRU TS
    geom_point() +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 45240 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 45240 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/basic%20graphs%20of%20temperatures-4.png)<!-- -->

``` r
# compare seasonality from CRU TS and ERSST
ggplot(trends, aes(seas_cruts, seas_ersst, color = rarefyID_y)) + # ERSST more muted at fastest rates of change in CRU TS
    geom_point() +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 45240 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 45240 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/basic%20graphs%20of%20temperatures-5.png)<!-- -->

### Combine CRU TS and ERSST

``` r
# temperature trend
trends[ , temptrend_comb_allyr := tmptrendallyr] # CRU TS for land and freshwater
trends[REALM == 'Marine', temptrend_comb_allyr := ssttrendallyr] # ERSST for marine

# average temperature
trends[, tempave_comb := tmpaveallyr]
trends[REALM == 'Marine', tempave_comb := sstaveallyr]

# seasonality
trends[, seas_comb := seas_cruts]
trends[REALM == 'Marine', seas_comb := seas_ersst]
```

### Center and scale covariates

``` r
trends[, temptrend_comb_allyr.sc := scale(temptrend_comb_allyr)]
trends[, temptrend_comb_allyr_abs.sc := scale(abs(temptrend_comb_allyr))]
trends[, seas_comb.sc := scale(seas_comb)]
trends[, tempave_comb.sc := scale(tempave_comb)]
trends[, npp.sc := scale(npp)]
trends[, mass_geomean.sc := scale(mass_geomean)]
```

### Plot turnover vs. temperature trends

Trends are pretty symmetric around no trend in temperature. Instead,
abs(trend) looks better.

``` r
# Jaccard turnover trend vs. temperature trend (across all years)
ggplot(trends, aes(temptrend_comb_allyr, Jtutrend, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/plot%20turnover%20v%20temp%20trend-1.png)<!-- -->

``` r
# Jaccard turnover trend vs. abs temperature trend (across all years)
ggplot(trends, aes(abs(temptrend_comb_allyr), Jtutrend, group = REALM, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    geom_smooth(method = 'lm')
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/plot%20turnover%20v%20temp%20trend-2.png)<!-- -->

``` r
# Jaccard total trend vs. temperature trend (across all years)
ggplot(trends, aes(temptrend_comb_allyr, Jbetatrend, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/plot%20turnover%20v%20temp%20trend-3.png)<!-- -->

``` r
# Jaccard turnover trend vs. abs temperature trend (across all years)
ggplot(trends, aes(abs(temptrend_comb_allyr), Jbetatrend, group = REALM, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    geom_smooth(method = 'lm')
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/plot%20turnover%20v%20temp%20trend-4.png)<!-- -->

``` r
# Jaccard turnover trend vs. temperature trend (across all years)
ggplot(trends, aes(temptrend_comb_allyr, Jtutrend, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/plot%20turnover%20v%20temp%20trend-5.png)<!-- -->

``` r
# Horn-Morisita turnover trend vs. abs temperature trend (across all years)
ggplot(trends, aes(abs(temptrend_comb_allyr), Horntrend, group = REALM, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    geom_smooth(method = 'lm')
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 4262 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 4262 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 4262 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/plot%20turnover%20v%20temp%20trend-6.png)<!-- -->

### Plot Richness trend vs. temperature trends

``` r
# Richness trend vs. temperature trend (across all years)
ggplot(trends, aes(temptrend_comb_allyr, Strend, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/richness%20trend%20v%20temp%20trend-1.png)<!-- -->

``` r
# Richness trend vs. abs temperature trend (across all years)
ggplot(trends, aes(abs(temptrend_comb_allyr), Strend, group = REALM, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_files/figure-gfm/richness%20trend%20v%20temp%20trend-2.png)<!-- -->

### Model choice: Jaccard turnover trend vs. temperature trend (across all years)

Baseline turnover is slowest in terrestrial, fastest in marine In
relation to temperature, freshwater turnover is fastest, terrestrial
slowest Seasonality, baseline temperature, NPP, and organismal mass are
not better explanatory factors than REALM

``` r
# the cases we can compare
cat('Number of data points available:\n')
```

    ## Number of data points available:

``` r
apply(trends[, .(Jtutrend, temptrend_comb_allyr, REALM, seas_comb, npp, mass_geomean)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##             Jtutrend temptrend_comb_allyr                REALM 
    ##                53467                50335                53467 
    ##            seas_comb                  npp         mass_geomean 
    ##                50335                53314                53092

``` r
i <- trends[, complete.cases(Jtutrend, temptrend_comb_allyr, REALM, seas_comb, npp, tempave_comb, mass_geomean)]
cat('Overall:\n')
```

    ## Overall:

``` r
sum(i) 
```

    ## [1] 49865

``` r
# fit models
mods <- vector('list', 0)
mods[[1]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr.sc + REALM + (1|STUDY_ID))]#, weights = 1/(Jtutrend_se_nz^2))]
mods[[2]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc + REALM + (1|STUDY_ID))] 
mods[[3]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*REALM + (1|STUDY_ID))]
mods[[4]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*seas_comb.sc + (1|STUDY_ID))]
mods[[5]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*tempave_comb.sc + (1|STUDY_ID))]
mods[[6]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*npp.sc + (1|STUDY_ID))]
mods[[7]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*mass_geomean.sc + (1|STUDY_ID))]

# examine models
cat('AIC:\n')
```

    ## AIC:

``` r
aics <- sapply(mods, AIC) - min(sapply(mods, AIC)); aics
```

    ## [1] 7117.49284  744.03614    0.00000   91.65972  148.30268  734.64240
    ## [7]  784.58964

``` r
cat('\nModel terms:\n')
```

    ## 
    ## Model terms:

``` r
summary(mods[[which.min(aics)]])
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: Jtutrend ~ temptrend_comb_allyr_abs.sc * REALM + (1 | STUDY_ID)
    ## 
    ## REML criterion at convergence: -116829.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.7603 -0.4124 -0.1307  0.1926  6.5007 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  STUDY_ID (Intercept) 0.003077 0.05547 
    ##  Residual             0.005562 0.07458 
    ## Number of obs: 49865, groups:  STUDY_ID, 293
    ## 
    ## Fixed effects:
    ##                                               Estimate Std. Error t value
    ## (Intercept)                                   0.049894   0.017368   2.873
    ## temptrend_comb_allyr_abs.sc                   0.048069   0.003226  14.900
    ## REALMMarine                                   0.038079   0.018260   2.085
    ## REALMTerrestrial                             -0.009009   0.018640  -0.483
    ## temptrend_comb_allyr_abs.sc:REALMMarine      -0.007453   0.003264  -2.284
    ## temptrend_comb_allyr_abs.sc:REALMTerrestrial -0.030044   0.003292  -9.125
    ## 
    ## Correlation of Fixed Effects:
    ##              (Intr) tm___. REALMM REALMT t___.:REALMM
    ## tmptrnd___.   0.021                                  
    ## REALMMarine  -0.951 -0.020                           
    ## REALMTrrstr  -0.932 -0.019  0.886                    
    ## t___.:REALMM -0.020 -0.988  0.019  0.019             
    ## t___.:REALMT -0.020 -0.980  0.019  0.017  0.968

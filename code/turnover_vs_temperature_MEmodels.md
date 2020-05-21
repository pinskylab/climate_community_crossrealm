Turnover vs. temperature ME models
================

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
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.
```

### Load data

``` r
# Turnover and covariates assembled by turnover_vs_temperature_prep.Rmd
trends <- fread('output/turnover_w_covariates.csv.gz')

# set realm order
trends[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = TRUE)]
```

### Where do we have data?

Mostly northern hemisphere, but spread all over. Little in Africa.

``` r
world <- map_data('world')
ggplot(world, aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = 'lightgray', color = 'white') +
    geom_point(data = trends, aes(rarefyID_x, rarefyID_y, group = REALM, color = REALM), size = 0.5, alpha = 0.4)  +
    scale_color_brewer(palette="Set1")
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/map-1.png)<!-- -->

### Plot turnover vs. temperature trends

Trends are pretty symmetric around no trend in temperature. Instead,
abs(trend) looks better.

``` r
# Jaccard turnover trend vs. temperature trend (across all years)
ggplot(trends, aes(temptrend_comb_allyr, Jtutrend, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    scale_color_brewer(palette="Set1")
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-1.png)<!-- -->

``` r
# Jaccard turnover trend vs. abs temperature trend (across all years)
ggplot(trends, aes(abs(temptrend_comb_allyr), Jtutrend, group = REALM, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    geom_smooth(method = 'lm') +
    scale_color_brewer(palette="Set1")
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-2.png)<!-- -->

``` r
# Jaccard total trend vs. temperature trend (across all years)
ggplot(trends, aes(temptrend_comb_allyr, Jbetatrend, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    scale_color_brewer(palette="Set1")
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-3.png)<!-- -->

``` r
# Jaccard turnover trend vs. abs temperature trend (across all years)
ggplot(trends, aes(abs(temptrend_comb_allyr), Jbetatrend, group = REALM, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    geom_smooth(method = 'lm') +
    scale_color_brewer(palette="Set1")
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-4.png)<!-- -->

``` r
# Horn-Morisita turnover trend vs. temperature trend (across all years)
ggplot(trends, aes(temptrend_comb_allyr, Horntrend, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    scale_color_brewer(palette="Set1")
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 4262 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 4262 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-5.png)<!-- -->

``` r
# Horn-Morisita turnover trend vs. abs temperature trend (across all years)
ggplot(trends, aes(abs(temptrend_comb_allyr), Horntrend, group = REALM, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    geom_smooth(method = 'lm') +
    scale_color_brewer(palette="Set1")
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 4262 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 4262 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 4262 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20turnover%20v%20temp%20trend-6.png)<!-- -->

### Plot Richness trend vs. temperature trends

``` r
# Richness trend vs. temperature trend (across all years)
ggplot(trends, aes(temptrend_comb_allyr, Strend, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    scale_color_brewer(palette="Set1")
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/richness%20trend%20v%20temp%20trend-1.png)<!-- -->

``` r
# Richness trend vs. abs temperature trend (across all years)
ggplot(trends, aes(abs(temptrend_comb_allyr), Strend, group = REALM, color = REALM)) +
    geom_point(size = 0.2, alpha = 0.5) + 
    geom_smooth() +
    scale_color_brewer(palette="Set1")
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3132 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 3132 rows containing missing values (geom_point).

![](turnover_vs_temperature_MEmodels_files/figure-gfm/richness%20trend%20v%20temp%20trend-2.png)<!-- -->

### Compare covariates across realms

Marine are in generally warmer locations (seawater doesn’t freeze)
Marine have much lower seasonality. Marine and freshwater have some very
small masses (plankton), but much of dataset is similar to terrestrial.
Marine has a lot of slow, crawling organisms, but land has plants. Land
also has birds
(fast).

``` r
i <- trends[, !duplicated(STUDY_ID)]; sum(i)
```

    ## [1] 332

``` r
beanplot(rarefyID_y ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Latitude (degN)', ll = 0.05)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->

``` r
beanplot(tempave_comb ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature (degC)', ll = 0.05)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-2.png)<!-- -->

``` r
beanplot(seas_comb ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Seasonality (degC)', ll = 0.05)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-3.png)<!-- -->

``` r
beanplot(npp ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'NPP', ll = 0.05)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-4.png)<!-- -->

``` r
beanplot(mass_geomean ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Mass (g)', ll = 0.05)
```

    ## log="y" selected

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-5.png)<!-- -->

``` r
beanplot(speed_geomean ~ REALM, data = trends[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Speed (km/hr)', ll = 0.05, log = '')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/compare%20across%20realms-6.png)<!-- -->

### Center and scale covariates

``` r
trends[, temptrend_comb_allyr.sc := scale(temptrend_comb_allyr)]
trends[, temptrend_comb_allyr_abs.sc := scale(log(abs(temptrend_comb_allyr)))]
trends[, seas_comb.sc := scale(seas_comb)]
trends[, tempave_comb.sc := scale(tempave_comb)]
trends[, npp.sc := scale(npp)]
trends[, mass_geomean.sc := scale(log(mass_geomean))]
trends[, speed_geomean.sc := scale(log(speed_geomean+1))]

# histograms to examine
trends[, hist(temptrend_comb_allyr.sc)]
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-1.png)<!-- -->

    ## $breaks
    ##  [1] -12 -10  -8  -6  -4  -2   0   2   4   6   8  10  12  14  16
    ## 
    ## $counts
    ##  [1]     1    38    72   323  1003 22397 25251   979   216    47     5
    ## [12]     1     1     1
    ## 
    ## $density
    ##  [1] 9.933446e-06 3.774709e-04 7.152081e-04 3.208503e-03 9.963246e-03
    ##  [6] 2.224794e-01 2.508294e-01 9.724844e-03 2.145624e-03 4.668720e-04
    ## [11] 4.966723e-05 9.933446e-06 9.933446e-06 9.933446e-06
    ## 
    ## $mids
    ##  [1] -11  -9  -7  -5  -3  -1   1   3   5   7   9  11  13  15
    ## 
    ## $xname
    ## [1] "temptrend_comb_allyr.sc"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
trends[, hist(temptrend_comb_allyr_abs.sc)]
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-2.png)<!-- -->

    ## $breaks
    ##  [1] -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3
    ## 
    ## $counts
    ##  [1]     1     0     0     2    19    63   320  1254  5484 16535 19238
    ## [12]  6612   807
    ## 
    ## $density
    ##  [1] 1.986689e-05 0.000000e+00 0.000000e+00 3.973378e-05 3.774709e-04
    ##  [6] 1.251614e-03 6.357405e-03 2.491308e-02 1.089500e-01 3.284991e-01
    ## [11] 3.821993e-01 1.313599e-01 1.603258e-02
    ## 
    ## $mids
    ##  [1] -9.5 -8.5 -7.5 -6.5 -5.5 -4.5 -3.5 -2.5 -1.5 -0.5  0.5  1.5  2.5
    ## 
    ## $xname
    ## [1] "temptrend_comb_allyr_abs.sc"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
trends[, hist(seas_comb.sc)]
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-3.png)<!-- -->

    ## $breaks
    ##  [1] -2.0 -1.5 -1.0 -0.5  0.0  0.5  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5
    ## [15]  5.0  5.5  6.0
    ## 
    ## $counts
    ##  [1]  1033  6796  9996  7839 11136 10131   407   466   566  1181   603
    ## [12]   120    49     8     3     1
    ## 
    ## $density
    ##  [1] 4.104500e-02 2.700308e-01 3.971789e-01 3.114731e-01 4.424754e-01
    ##  [6] 4.025430e-01 1.617165e-02 1.851594e-02 2.248932e-02 4.692560e-02
    ## [11] 2.395947e-02 4.768054e-03 1.946955e-03 3.178703e-04 1.192014e-04
    ## [16] 3.973378e-05
    ## 
    ## $mids
    ##  [1] -1.75 -1.25 -0.75 -0.25  0.25  0.75  1.25  1.75  2.25  2.75  3.25
    ## [12]  3.75  4.25  4.75  5.25  5.75
    ## 
    ## $xname
    ## [1] "seas_comb.sc"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
trends[, hist(tempave_comb.sc)]
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-4.png)<!-- -->

    ## $breaks
    ##  [1] -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5  0.0  0.5  1.0  1.5  2.0  2.5
    ## 
    ## $counts
    ##  [1]     1     5    10  2666  5184  7788 13242  7903  4441  2070  6617
    ## [12]   408
    ## 
    ## $density
    ##  [1] 3.973378e-05 1.986689e-04 3.973378e-04 1.059303e-01 2.059799e-01
    ##  [6] 3.094467e-01 5.261548e-01 3.140161e-01 1.764577e-01 8.224893e-02
    ## [11] 2.629184e-01 1.621138e-02
    ## 
    ## $mids
    ##  [1] -3.25 -2.75 -2.25 -1.75 -1.25 -0.75 -0.25  0.25  0.75  1.25  1.75
    ## [12]  2.25
    ## 
    ## $xname
    ## [1] "tempave_comb.sc"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
trends[, hist(npp.sc)]
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-5.png)<!-- -->

    ## $breaks
    ##  [1] -1.5 -1.0 -0.5  0.0  0.5  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5  5.0
    ## [15]  5.5  6.0  6.5
    ## 
    ## $counts
    ##  [1]  3838 13798 15770  8777  4760  1936  1602   944   695   534   315
    ## [12]   161    91    55    37     1
    ## 
    ## $density
    ##  [1] 0.1439771917 0.5176126346 0.5915894512 0.3292568556 0.1785647297
    ##  [6] 0.0726263270 0.0600967851 0.0354128372 0.0260719511 0.0200322617
    ## [11] 0.0118167836 0.0060396894 0.0034137375 0.0020632479 0.0013880032
    ## [16] 0.0000375136
    ## 
    ## $mids
    ##  [1] -1.25 -0.75 -0.25  0.25  0.75  1.25  1.75  2.25  2.75  3.25  3.75
    ## [12]  4.25  4.75  5.25  5.75  6.25
    ## 
    ## $xname
    ## [1] "npp.sc"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
trends[, hist(mass_geomean.sc)]
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-6.png)<!-- -->

    ## $breaks
    ##  [1] -11 -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4
    ## 
    ## $counts
    ##  [1]     1     0    52    35     0     3   162   955   976  2426 19667
    ## [12] 19735  8622   441     3
    ## 
    ## $density
    ##  [1] 1.884020e-05 0.000000e+00 9.796903e-04 6.594069e-04 0.000000e+00
    ##  [6] 5.652059e-05 3.052112e-03 1.799239e-02 1.838803e-02 4.570632e-02
    ## [11] 3.705302e-01 3.718113e-01 1.624402e-01 8.308527e-03 5.652059e-05
    ## 
    ## $mids
    ##  [1] -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5
    ## [12]   0.5   1.5   2.5   3.5
    ## 
    ## $xname
    ## [1] "mass_geomean.sc"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

``` r
trends[, hist(speed_geomean.sc)]
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/center%20and%20scale-7.png)<!-- -->

    ## $breaks
    ##  [1] -2.0 -1.5 -1.0 -0.5  0.0  0.5  1.0  1.5  2.0  2.5  3.0  3.5  4.0
    ## 
    ## $counts
    ##  [1]  4492  8479  5032  2155 10592 14006  7044   746    68    20    11
    ## [12]     4
    ## 
    ## $density
    ##  [1] 0.1706395183 0.3220953864 0.1911527284 0.0818629034 0.4023628179
    ##  [6] 0.5320518908 0.2675834299 0.0283386199 0.0025831450 0.0007597485
    ## [11] 0.0004178617 0.0001519497
    ## 
    ## $mids
    ##  [1] -1.75 -1.25 -0.75 -0.25  0.25  0.75  1.25  1.75  2.25  2.75  3.25
    ## [12]  3.75
    ## 
    ## $xname
    ## [1] "speed_geomean.sc"
    ## 
    ## $equidist
    ## [1] TRUE
    ## 
    ## attr(,"class")
    ## [1] "histogram"

### Model choice: Jaccard turnover trend vs. temperature trend (across all years)

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

Then choose the variance structure

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

Finally, compare fixed effects among
models

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

### Validate the chosen model

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

### Plot the chosen model

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

Dissimilarity slope covariate data prep and visualization
================

Uses slope of dissimilarity vs. time from calc\_turnover.Rmd.

``` r
if(Sys.info()['nodename'] == 'annotate.sebs.rutgers.edu'){
  library('mgcv', lib.loc = '/usr/lib64/R/library') # when running on Annotate. Need to load 1.8-26, not 1.8-33.
} else {
  library('mgcv')
}
```

    ## Loading required package: nlme

    ## This is mgcv 1.8-26. For overview type 'help("mgcv-package")'.

``` r
library(data.table)
library(ggplot2)
library(beanplot) # for beanplots
library(gridExtra) # to combine ggplots together
library(grid) # to combine ggplots together
require(scales) # for custom axis scales
```

    ## Loading required package: scales

``` r
require(here)
```

    ## Loading required package: here

    ## here() starts at /local/home/malinp/climate_community_crossrealm

``` r
#knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.

signedsqrt_trans <- function() trans_new('signedsqrt', 
                                         transform = function(x) sign(x)*sqrt(abs(x)), 
                                         inverse = function(x) sign(x)*x^2)
```

# Load data

## Add covariates to BT data

## Trim data

``` r
norig <- nrow(bt); norig
```

    ## [1] 120654

``` r
# trim out timeseries with too few species
bt[Nspp <= 2 | is.na(Nspp), .N] #0
```

    ## [1] 0

``` r
bt <- bt[Nspp > 2, ]
  
# or few individuals
bt[Nave < 10, .N]
```

    ## [1] 3297

``` r
bt <- bt[Nave >= 10 | is.na(Nave), ]


nrow(bt)
```

    ## [1] 117357

``` r
nrow(bt)/norig
```

    ## [1] 0.9726739

``` r
bt[, length(unique(STUDY_ID))]
```

    ## [1] 204

``` r
bt[, length(unique(rarefyID))]
```

    ## [1] 10434

## Set up useful variables and transformations

## Write out

Only if file doesn’t yet exist

# Check variable distributions

## Response variables

``` r
bt[measure =='Jbeta' & duration_group == '3', summary(disstrend)]
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -1.00000 -0.05417  0.01868  0.01601  0.12500  0.65233

``` r
bt[measure =='Jtu' & duration_group == '3', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -1.000000 -0.095238  0.000000  0.009226  0.125000  1.000000

``` r
bt[measure =='Horn' & duration_group == '3', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -1.000000 -0.056011  0.009088  0.013092  0.145991  0.999840

``` r
bt[measure =='Jbeta' & duration_group == '5', summary(disstrend)]
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -0.18000 -0.01237  0.00966  0.01168  0.03961  0.22500

``` r
bt[measure =='Jtu' & duration_group == '5', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.285000 -0.025749  0.001921  0.008775  0.043333  0.321702

``` r
bt[measure =='Horn' & duration_group == '5', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.299190 -0.022781  0.003126  0.012460  0.044560  0.289490

``` r
bt[measure =='Jbeta' & duration_group == '10', summary(disstrend)]
```

    ##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
    ## -0.0318462  0.0002491  0.0061082  0.0080168  0.0142269  0.0958410

``` r
bt[measure =='Jtu' & duration_group == '10', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.106219 -0.003782  0.003558  0.004162  0.013034  0.097424

``` r
bt[measure =='Horn' & duration_group == '10', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.050177 -0.001275  0.005531  0.009390  0.017285  0.140489

``` r
bt[measure =='Jbeta' & duration_group == '20', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.005446  0.002775  0.003934  0.004658  0.005762  0.028083

``` r
bt[measure =='Jtu' & duration_group == '20', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.019038  0.001050  0.002657  0.002810  0.004431  0.031571

``` r
bt[measure =='Horn' & duration_group == '20', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.007307  0.002772  0.004881  0.006182  0.008285  0.033522

``` r
bt[measure =='Jbeta' & duration_group == 'All', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.366667  0.000000  0.005893  0.008879  0.015689  0.250000

``` r
bt[measure =='Jtu' & duration_group == 'All', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.420000 -0.003611  0.003987  0.006832  0.015914  0.433333

``` r
bt[measure =='Horn' & duration_group == 'All', summary(disstrend)]
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.397921 -0.001918  0.005509  0.011275  0.020494  0.383642

![](assemble_slope_covariates_files/figure-gfm/histograms%20response-1.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-2.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-3.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-4.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-5.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-6.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-7.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-8.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-9.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-10.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-11.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-12.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-13.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-14.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-15.png)<!-- -->

## Unscaled covariates

![](assemble_slope_covariates_files/figure-gfm/histograms%20unscaled-1.png)<!-- -->

## Scaled covariates

![](assemble_slope_covariates_files/figure-gfm/histograms%20scaled-1.png)<!-- -->

# Response variable SE

``` r
ggplot(bt, aes(disstrend, trendse)) + geom_point() + facet_wrap(vars(duration_group, measure), scales = 'free', ncol = 3)
```

![](assemble_slope_covariates_files/figure-gfm/response%20SE-1.png)<!-- -->

# Check correlations among variables

Pearson’s r is in the lower triangle
![](assemble_slope_covariates_files/figure-gfm/pairs-1.png)<!-- -->

# Compare covariates across realms

    ## [1] 10434

    ## log="y" selected
    ## log="y" selected

![](assemble_slope_covariates_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->

  - Marine are in generally warmer locations (seawater doesn’t freeze)
  - Marine have much lower seasonality.
  - Marine and freshwater have some very small masses (plankton), but
    much of dataset is similar to terrestrial.
  - Marine has a lot of slow, crawling organisms, but land has plants.
    Land also has birds (fast).

# Plot slope vs. duration

Use the slopes with unstandardized durations

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](assemble_slope_covariates_files/figure-gfm/plot%20slope%20vs%20duration-1.png)<!-- -->

# Plot dissimilarity vs. explanatory variables

Lines are ggplot smoother fits

## Jtu

### Realm

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Temperature trend

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### abs(Temperature trend)

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### Temperature

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Metabolic temperature

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Seasonality

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Microclimates

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### Mass

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Endotherms

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### Richness

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### Thermal bias

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

### NPP

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

### Human

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Jbeta (total)

### Realm

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

### Temperature trend

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

### abs(Temperature trend)

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### Temperature

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### Metabolic temperature

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

### Seasonality

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

### Microclimate

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

### Mass

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

### Endotherms

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

### Richness

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

### Thermal bias

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### NPP

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

### Human

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

## Horn

### Realm

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

### Temperature trend

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

### abs(Temperature trend)

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

### Temperature

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

### Metabolic temperature

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

### Seasonality

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

### Microclimates

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

### Mass

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

### Endotherms

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

### Richness

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

### Thermal bias

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

### NPP

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

### Human

![](assemble_slope_covariates_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

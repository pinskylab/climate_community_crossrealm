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

    ## This is mgcv 1.8-33. For overview type 'help("mgcv-package")'.

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

    ## here() starts at /Users/mpinsky/Documents/Rutgers/Community_and_climate/crossrealm/climate_community_crossrealm

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

    ## [1] 56530

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

    ## [1] 1191

``` r
bt <- bt[Nave >= 10 | is.na(Nave), ]


nrow(bt)
```

    ## [1] 55339

``` r
nrow(bt)/norig
```

    ## [1] 0.9789315

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

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -1.000000 -0.018167  0.008529  0.013488  0.063542  0.652327

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -1.000000 -0.032807  0.000000  0.008329  0.057876  1.000000

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.999840 -0.070581 -0.005387 -0.012287  0.022541  1.000000

![](assemble_slope_covariates_files/figure-gfm/histograms%20response-1.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-2.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-3.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-4.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-5.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-6.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-7.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-8.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-9.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-10.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-11.png)<!-- -->![](assemble_slope_covariates_files/figure-gfm/histograms%20response-12.png)<!-- -->

## Unscaled covariates

![](assemble_slope_covariates_files/figure-gfm/histograms%20unscaled-1.png)<!-- -->

## Scaled covariates

![](assemble_slope_covariates_files/figure-gfm/histograms%20scaled-1.png)<!-- -->

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

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](assemble_slope_covariates_files/figure-gfm/plot%20slope%20vs%20duration-1.png)<!-- -->

# Plot dissimilarity vs. explanatory variables

Lines are ggplot smoother fits \#\# Jtu \#\#\# Realm
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

Turnover covariate data prep
================

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
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # tell RStudio to use project root directory as the root for this notebook. Needed since we are storing code in a separate directory.

signedsqrt_trans <- function() trans_new('signedsqrt', 
                                         transform = function(x) sign(x)*sqrt(abs(x)), 
                                         inverse = function(x) sign(x)*x^2)
```

# Load data

``` r
# BioTime
load('data/biotime_blowes/bt_malin.Rdata') # loads bt_malin
bt <- data.table(bt_malin); rm(bt_malin) # rename

# BioTime data type
load('data/biotime_blowes/time_series_data_type.Rdata') # loads rarefyID_type
bt <- merge(bt, rarefyID_type, by = 'rarefyID', all.x = TRUE) # merge with biotime

# BioTime temporal turnover null model
load('data/biotime_blowes/z_scores_SAD.Rdata') # loads sim_summary

# Temperature average, trends, and seasonality
temperature <- fread('output/temperature_byrarefyID.csv.gz')

# microclimates
microclim <- fread('output/microclimates.csv.gz', drop = 1)

# NPP
npp <- fread('output/npplandocean.csv.gz')

# Body size
bs <- fread('output/mass_byrarefyid.csv.gz', drop = 1)
bs[, ':='(STUDY_ID = NULL, REALM = NULL, taxa_mod = NULL)] # remove unnecessary columns 

# Mobility
speed <- fread('output/speed_byrarefyID.csv.gz', drop = 1)
speed[, ':='(STUDY_ID = NULL, REALM = NULL, taxa_mod = NULL)] # remove unnecessary columns 

# Lifespan
lsp <- fread('output/lifespan_byrarefyID.csv.gz')

# CTI
cti <- fread('output/cti_byrarefyID.csv.gz')
    
# consumer vs. producer
consfrac <- fread('output/consfrac_byrarefyID.csv.gz')

# richness
rich <- fread('output/richness_by_rarefyID.csv.gz') # number of species

# endotherm vs. ectotherm
endofrac <- fread('output/endofrac_byrarefyID.csv.gz') # endotherm vs. ectotherm classifications

# human impact
human <- fread('output/humanimpact_by_rarefyID.csv.gz')

# %veg
veg <- as.data.table(readRDS('output/vct_by_rarefyID.rds'))
veg[, veg := (`tree cover % (mean)` + 0.5 * `non-tree veg. % (mean)`)/100] # veg index from 0 (all non-veg) to 1 (all tree). Non-tree veg counts as 0.5.
```

## Plot a turnover example

First is a long time-series, next one shows an example result of
removing the first year of self-comparison.

``` r
#bt[, .(n = .N), by = rarefyID][n > 50,]

ggplot(bt[rarefyID == '339_1085477'], aes(YEAR, Jtu_base)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    xlab('Year') + ylab('Jaccard dissimilarity')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20turnover-1.png)<!-- -->

``` r
ggplot(bt[rarefyID == '489_16754'], aes(YEAR, Jtu_base)) +
    geom_point() +
    geom_smooth(method = 'lm', se = FALSE) +
    xlab('Year') + ylab('Jaccard dissimilarity') +
  geom_smooth(data = bt[rarefyID == '489_16754' & YEAR != 2000], method = 'lm', se = FALSE, color = 'red')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20turnover-2.png)<!-- -->

## Assemble dataset of beta diversity trends (temporal turnover) and covariates

Also trim to count data

``` r
# function to calc trends, removing first year with 0
calctrendrem0 <- function(y, YEAR, nm = 'y'){
  # turn off warnings for the following
  defaultW <- getOption("warn")
  options(warn = -1)
  
  if(length(YEAR)>2){
    o <- order(YEAR)
    YEAR2 <- YEAR[o][2:length(YEAR)]
    y2 <- y[o][2:length(y)]
    
    if(sum(!is.na(y2)) >= 2){ # make sure enough values to fit a line
      mod <- lm(y2 ~ YEAR2)
      out <- list(y = coef(mod)[2], # coef for the slope
                  y_se = sqrt(diag(vcov(mod)))[2]) # SE
      names(out) <- c(nm, paste0(nm, '_se'))
      options(warn = defaultW)
      return(out)
      
    } else {
      out <- list(y = NA_real_, y_se = NA_real_)
      names(out) <- c(nm, paste0(nm, '_se'))
      
      options(warn = defaultW)
      return(out)
    }
    
  } else {
    out <- list(y = NA_real_, y_se = NA_real_)
    names(out) <- c(nm, paste0(nm, '_se'))
  
    options(warn = defaultW)
    return(out)
  }
  

}

calcfirstlast <- function(y, YEAR, nm = 'y'){ # function to get distance from last to first year
  # turn off warnings for the following
  defaultW <- getOption("warn")
  options(warn = -1)
  
  if(length(YEAR)>1){
    o <- order(YEAR)
    YEAR2 <- YEAR[o][2:length(YEAR)]
    y2 <- y[o][2:length(y)]
    
    out <- list(y = y2[length(y2)], # dissimilarity for last year
                y_se = NA_real_)
  } else {
    out <- list(y = NA_real_, y_se = NA_real_)
  }

  names(out) <- c(nm, paste0(nm, '_se'))
  options(warn = defaultW)
  return(out)
}



setkey(bt, STUDY_ID, rarefyID, YEAR)

if(file.exists('temp/trendstemp.rds')){
    trends <- readRDS('temp/trendstemp.rds')
    print('Loaded from file')
} else {
    print('Calculating from scratch')
    
    trends5 <- bt[BROAD_TYPE == 'count', calctrendrem0(Jtu_base, YEAR, 'Jtutrendrem0'), 
                  by = .(REALM, Biome, taxa_mod, STUDY_ID, 
                         rarefyID, rarefyID_x, rarefyID_y)] # calculate trend in Jaccard turnover without first year
    trends6 <- bt[BROAD_TYPE == 'count', calctrendrem0(Jbeta_base, YEAR, 'Jbetatrendrem0'),
                  by = .(rarefyID)]
    trends7 <- bt[BROAD_TYPE == 'count', calctrendrem0(1-Horn_base, YEAR, 'Horntrendrem0'),
                  by = .(rarefyID)]
    
    trends8 <- bt[BROAD_TYPE == 'count', calcfirstlast(Jtu_base, YEAR, 'Jtulast'),
                  by = .(rarefyID)]
    trends9 <- bt[BROAD_TYPE == 'count', calcfirstlast(Jbeta_base, YEAR, 'Jbetalast'),
                  by = .(rarefyID)]
    trends10 <- bt[BROAD_TYPE == 'count', calcfirstlast(1-Horn_base, YEAR, 'Hornlast'),
                   by = .(rarefyID)]
    
    
    nyrBT <-  bt[BROAD_TYPE == 'count', .(nyrBT = length(YEAR),
                     minyrBT = min(YEAR),
                     maxyrBT = max(YEAR),
                     medianyrBT = median(YEAR),
                     meanyrBT = mean(YEAR)),
                 by = .(rarefyID)] # number of years in time-series
    
    trends <- merge(trends5, trends6)
    trends <- merge(trends, trends7)
    trends <- merge(trends, trends8)
    trends <- merge(trends, trends9)
    trends <- merge(trends, trends10)
    trends <- merge(trends, nyrBT)
    
    saveRDS(trends, file = 'temp/trendstemp.rds')
    
}
```

    ## [1] "Loaded from file"

``` r
nrow(trends)
```

    ## [1] 52016

### Standardize slope values against null model simulations

``` r
# remove some duplicates
sim_summary <- as.data.table(sim_summary)
sim_summary <- sim_summary[!duplicated(sim_summary), ]

# hist
hist(sim_summary[metric == 'Jaccard (total)', sd], breaks = 200, col = 'grey')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20sds-1.png)<!-- -->

``` r
hist(sim_summary[metric == 'Nestedness', sd], breaks = 200, col = 'grey')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20sds-2.png)<!-- -->

``` r
hist(sim_summary[metric == 'Morisita-Horn', sd], breaks = 200, col = 'grey')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20sds-3.png)<!-- -->

``` r
# plot
sim_summary[metric == 'Jaccard (total)', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20sds-4.png)<!-- -->

``` r
sim_summary[metric == 'Nestedness', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20sds-5.png)<!-- -->

``` r
sim_summary[metric == 'Morisita-Horn', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20sds-6.png)<!-- -->

``` r
# vs #yrs
sim_summary[metric == 'Jaccard (total)', plot(Nyrs, sd+0.001, log = 'y', main = 'Null SD declines with #yrs')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20sds-7.png)<!-- -->

``` r
# vs #spp
sim_summary[metric == 'Jaccard (total)', plot(Nspp, sd+0.001, log = 'y')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20sds-8.png)<!-- -->

Merge null model

``` r
# merge
trends <- merge(trends, sim_summary[metric == 'Jaccard (total)', .(Jbeta_exp = expected, Jbeta_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
trends <- merge(trends, sim_summary[metric == 'Nestedness', .(Jtu_exp = expected, Jtu_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
trends <- merge(trends, sim_summary[metric == 'Morisita-Horn', .(Horn_exp = expected, Horn_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
```

Set a threshold and only standardize those timeseries with SD \> 0.001

``` r
# standardize, only where sd >0
trends[, Jtutrendz := NA_real_]
trends[, Jbetatrendz := NA_real_]
trends[, Horntrendz := NA_real_]

trends[Jtu_sd > 0.001, Jtutrendz := (Jtutrendrem0 - Jtu_exp)/Jtu_sd]
trends[Jbeta_sd > 0.001, Jbetatrendz := (Jbetatrendrem0 - Jbeta_exp)/Jbeta_sd]
trends[Horn_sd > 0.001, Horntrendz := (Horntrendrem0 - Horn_exp)/Horn_sd]

nrow(trends)
```

    ## [1] 52016

Add covariates

``` r
# add covariates
trends <- merge(trends, temperature[, .(rarefyID, tempave, tempave_metab, temptrend, seas)], all.x = TRUE, by = 'rarefyID') # temperature ave, ave metabolic, trend, and seasonality
trends <- merge(trends, microclim[, .(rarefyID, microclim = Temp_sd20km)], all.x = TRUE, by = 'rarefyID') # microclimates
trends <- merge(trends, npp, all.x = TRUE, by = 'rarefyID') # npp
trends <- merge(trends, bs[, .(rarefyID, mass_mean_weight, mass_sd_weight)], all.x = TRUE) # body size mass (g)
trends <- merge(trends, speed[, .(rarefyID, speed_mean_weight, speed_sd_weight)], all.x = TRUE) # speed (km/hr)
trends <- merge(trends, lsp[, .(rarefyID, lifespan_mean_weight, lifespan_sd_weight)], all.x = TRUE) # lifespan (yr)
trends <- merge(trends, cti[, .(rarefyID, thermal_bias)], all.x = TRUE) # thermal bias (degC)
trends <- merge(trends, consfrac[, .(rarefyID, consfrac)], all.x = TRUE) # fraction consumers
trends <- merge(trends, rich, all.x = TRUE) # species richness
trends <- merge(trends, endofrac[, .(rarefyID, endofrac)], all.x = TRUE) # endotherm vs. ectotherm
trends <- merge(trends, human[, .(rarefyID, human_bowler = atc, human_venter = hfp, human_halpern = himp)], all.x = TRUE) # human impact
trends <- merge(trends, veg[, .(rarefyID, veg = veg)], all.x = TRUE) # vegetation index
trends[REALM == 'Marine', veg := 0] # veg index is 0 at sea
```

Remove studies with only 1 species

``` r
nrow(trends)
```

    ## [1] 52016

``` r
trends <- trends[Nspp > 1, ]
nrow(trends)
```

    ## [1] 51578

#### Output a file with plots of every Jtu timeseries (a lot\!)

Not run during knitting

``` r
rids <- trends[Nspp > 1 & nyrBT > 2, sort(unique(rarefyID))]
setkey(bt, rarefyID, YEAR)
print(paste(length(rids), ' rarefyIDs'))
filenum <- 1
plotnum <- 1
for(i in 1:length(rids)){
  jtutrendrem0 <- trends[rarefyID == rids[i], Jtutrendrem0]
  nspp <- trends[rarefyID == rids[i], Nspp]
  x <- bt[rarefyID == rids[i], YEAR]
  y <- bt[rarefyID == rids[i], Jtu_base]
  
  if(length(x) > 2 & nspp > 1){
    if(plotnum %% 400 == 1){
      if(plotnum >1) dev.off()
      png(file = paste0('figures/jtu_plots/jtu_plots', formatC(filenum, width = 3, format = 'd', flag = '0'), '.png'), width = 36, height = 36, units = 'in', res = 100)
      par(mfrow=c(20,20), mai = c(0.4, 0.5, 0.5, 0.1))
      filenum <- filenum + 1
    }
    plot(x, y, main = paste('Jtu rem0:', signif(jtutrendrem0, 3), '\nNspp:', nspp, '\nrID:', rids[i]), xlab = '', ylab = 'Jtu',
         cex.main = 0.7)
    abline(lm(y ~ x))
    
    x2 <- x[2:length(x)]
    y2 <- y[2:length(y)]
    abline(lm(y2 ~ x2), col = 'red')
    
    plotnum <- plotnum + 1
  }
}
  
dev.off()
```

### Do some basic checks of the turnover calculations

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
    ## 51574:  97_2204414 Marine     Northern_European_Seas Invertebrates
    ## 51575:  97_2273557 Marine                     Arctic Invertebrates
    ## 51576:    98_88435 Marine                     Arctic       Benthos
    ## 51577:    98_93555 Marine                     Arctic       Benthos
    ## 51578:    98_97933 Marine                     Arctic       Benthos
    ##        STUDY_ID rarefyID_x rarefyID_y Jtutrendrem0 Jtutrendrem0_se
    ##     1:      100   -3.08000   51.14000   0.00297623     0.001309059
    ##     2:      101   -3.08000   51.14000   0.00000000     0.000000000
    ##     3:      108   57.96500  -65.28500           NA              NA
    ##     4:      108   59.92750  -66.29250           NA              NA
    ##     5:      108   59.97000  -66.19500           NA              NA
    ##    ---                                                            
    ## 51574:       97   33.50280   69.51390  -0.50000000             NaN
    ## 51575:       97   82.69667   71.84167           NA              NA
    ## 51576:       98 -137.25567   68.99443           NA              NA
    ## 51577:       98 -133.90000   69.65420           NA              NA
    ## 51578:       98 -132.98750   69.40090           NA              NA
    ##        Jbetatrendrem0 Jbetatrendrem0_se Horntrendrem0 Horntrendrem0_se
    ##     1:    0.002587992      0.0009975709  0.0023147105     0.0016045472
    ##     2:    0.002208307      0.0009103225  0.0001735803     0.0009087936
    ##     3:             NA                NA            NA               NA
    ##     4:             NA                NA            NA               NA
    ##     5:             NA                NA            NA               NA
    ##    ---                                                                
    ## 51574:   -0.035714286               NaN -0.0271517119              NaN
    ## 51575:             NA                NA            NA               NA
    ## 51576:             NA                NA            NA               NA
    ## 51577:             NA                NA            NA               NA
    ## 51578:             NA                NA            NA               NA
    ##          Jtulast Jtulast_se Jbetalast Jbetalast_se  Hornlast Hornlast_se
    ##     1: 0.4324324         NA 0.5000000           NA 0.3642417          NA
    ##     2: 0.0000000         NA 0.5000000           NA 0.1541960          NA
    ##     3: 1.0000000         NA 1.0000000           NA 1.0000000          NA
    ##     4: 0.8000000         NA 0.8571429           NA 0.9177154          NA
    ##     5: 0.5000000         NA 0.7142857           NA 0.7026207          NA
    ##    ---                                                                  
    ## 51574: 0.0000000         NA 0.9285714           NA 0.9456966          NA
    ## 51575: 1.0000000         NA 1.0000000           NA 1.0000000          NA
    ## 51576: 0.6666667         NA 0.7500000           NA 0.3265231          NA
    ## 51577: 0.4000000         NA 0.8571429           NA 0.9516536          NA
    ## 51578: 0.0000000         NA 0.8750000           NA 0.8025743          NA
    ##        nyrBT minyrBT maxyrBT medianyrBT meanyrBT     Jbeta_exp
    ##     1:    31    1981    2011     1996.0 1996.000  4.991765e-05
    ##     2:    31    1981    2011     1996.0 1996.000 -5.224876e-05
    ##     3:     2    1985    1999     1992.0 1992.000           NaN
    ##     4:     2    1985    1993     1989.0 1989.000           NaN
    ##     5:     2    1985    1993     1989.0 1989.000           NaN
    ##    ---                                                        
    ## 51574:     3    1934    1952     1950.0 1945.333            NA
    ## 51575:     2    1934    1955     1944.5 1944.500            NA
    ## 51576:     2    1973    1975     1974.0 1974.000            NA
    ## 51577:     2    1973    1975     1974.0 1974.000            NA
    ## 51578:     2    1971    1973     1972.0 1972.000            NA
    ##            Jbeta_sd       Jtu_exp       Jtu_sd      Horn_exp      Horn_sd
    ##     1: 0.0008978943 -2.209763e-05 0.0006692022  3.092595e-05 0.0006016976
    ##     2: 0.0014642102 -1.270871e-05 0.0012548416 -3.276256e-05 0.0008557850
    ##     3:           NA           NaN           NA           NaN           NA
    ##     4:           NA           NaN           NA           NaN           NA
    ##     5:           NA           NaN           NA           NaN           NA
    ##    ---                                                                   
    ## 51574:           NA            NA           NA            NA           NA
    ## 51575:           NA            NA           NA            NA           NA
    ## 51576:           NA            NA           NA            NA           NA
    ## 51577:           NA            NA           NA            NA           NA
    ## 51578:           NA            NA           NA            NA           NA
    ##         Jtutrendz Jbetatrendz Horntrendz    tempave tempave_metab
    ##     1:         NA          NA         NA 12.0513496    12.0513496
    ##     2: 0.01012774    1.543874         NA 12.0513496    12.0513496
    ##     3:         NA          NA         NA -1.3200413    40.0000000
    ##     4:         NA          NA         NA -1.2716919    40.0000000
    ##     5:         NA          NA         NA -1.2716919    40.0000000
    ##    ---                                                           
    ## 51574:         NA          NA         NA  4.8661794     4.8661794
    ## 51575:         NA          NA         NA         NA            NA
    ## 51576:         NA          NA         NA         NA            NA
    ## 51577:         NA          NA         NA -0.5274576    -0.5274576
    ## 51578:         NA          NA         NA         NA            NA
    ##           temptrend      seas  microclim       npp mass_mean_weight
    ##     1:  0.041129329 3.0760983 0.23603834 1685.9467     4.576354e+03
    ##     2:  0.041129329 3.0760983 0.23603834 1685.9467     1.546288e+01
    ##     3: -0.006945833 0.4808712 0.03647966  126.3856     1.011728e+03
    ##     4: -0.004529737 0.5261527 0.01186496  152.3814     3.802604e+02
    ##     5: -0.004529737 0.5261527 0.01315645  150.4591     2.559320e+03
    ##    ---                                                             
    ## 51574:  0.009485958 2.1870667 0.30366775 1362.4417     6.900539e-02
    ## 51575:           NA        NA 0.15973145  261.5264     5.060336e-03
    ## 51576:           NA        NA 0.11977866 2694.4831     4.445541e-03
    ## 51577: -0.441934864 1.4292877 0.06630448 2420.8457     1.000115e-02
    ## 51578:           NA        NA 0.13839841  268.2280     1.649153e-02
    ##        mass_sd_weight speed_mean_weight speed_sd_weight
    ##     1:   1.608100e+04       12.02066251      9.88345276
    ##     2:   6.575757e+01        5.70189420      6.03351600
    ##     3:   6.624470e+02      126.52879895      7.02709129
    ##     4:   3.606983e+02      106.39887140     18.54357255
    ##     5:   2.480082e+03       73.39031844     57.37846293
    ##    ---                                                 
    ## 51574:   2.215162e-01        0.20260052      0.31886895
    ## 51575:   1.328297e-02        0.06805861      0.08900776
    ## 51576:   2.557748e-03        0.14273492      0.01572383
    ## 51577:   1.149557e-02        0.21783476      0.41993705
    ## 51578:   1.909184e-02        0.04479078      0.08106260
    ##        lifespan_mean_weight lifespan_sd_weight thermal_bias  consfrac Nspp
    ##     1:            24.368790         13.8561561   -0.1166378 1.0000000   83
    ##     2:             6.855411          5.5954324   -1.3322161 1.0000000   15
    ##     3:             2.322789          0.2962608    5.8210833 1.0000000    4
    ##     4:             1.835437          0.3587385    5.5271034 1.0000000    7
    ##     5:             2.628751          0.8313121    3.0517163 1.0000000    7
    ##    ---                                                                    
    ## 51574:             3.580285          2.4189608    4.0842456 1.0000000   42
    ## 51575:                   NA                 NA           NA 0.9793644   37
    ## 51576:                   NA                 NA           NA 0.4236728    9
    ## 51577:             6.481516          1.4691047    3.8125779 0.9652155   23
    ## 51578:                   NA                 NA           NA 0.9925757   17
    ##        endofrac human_bowler human_venter human_halpern veg
    ##     1:        0            6    33.002857            NA   0
    ##     2:        0            6    33.002857            NA   0
    ##     3:        1            0           NA     6.8752594   0
    ##     4:        1            0           NA     1.0418338   0
    ##     5:        1            0           NA     1.0914448   0
    ##    ---                                                     
    ## 51574:        0            0           NA    10.4499359   0
    ## 51575:        0            0           NA     0.6958275   0
    ## 51576:        0            0           NA     0.4557864   0
    ## 51577:        0            1           NA     0.1919821   0
    ## 51578:        0            1     4.757507     9.3018417   0

``` r
summary(trends$Jtutrendrem0)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##  -1.000  -0.005   0.000   0.007   0.020   1.000   13390

``` r
summary(trends$Jbetatrendrem0)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##  -1.000  -0.004   0.003   0.007   0.017   1.000   13390

``` r
summary(trends$Horntrendrem0)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##  -1.000  -0.008   0.002   0.009   0.023   1.000   13390

``` r
summary(trends$Jtutrendz)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ## -44.240  -0.281   0.039   0.323   0.884  53.387   18197

``` r
summary(trends$Jbetatrendz)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ## -31.591  -0.146   0.077   0.264   0.552  84.690   16894

``` r
summary(trends$Horntrendz)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ## -74.446  -0.310   0.051   0.516   1.073 104.775   16994

``` r
summary(trends$Jtulast)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.2000  0.5000  0.4958  0.8000  1.0000

``` r
summary(trends$Jbetalast)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.5000  0.6667  0.6757  0.8571  1.0000

``` r
summary(trends$Hornlast)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.2039  0.6378  0.5791  0.9816  1.0000

#### Histograms of temporal change

Standardized slopes have very large and small values

``` r
x <- trends[, hist(Jtutrendrem0)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-1.png)<!-- -->

``` r
x <- trends[, hist(Jbetatrendrem0)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-2.png)<!-- -->

``` r
x <- trends[, hist(Horntrendrem0)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-3.png)<!-- -->

``` r
x <- trends[, hist(Jtutrendz, breaks = 100)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-4.png)<!-- -->

``` r
x <- trends[, hist(Jbetatrendz, breaks = 100)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-5.png)<!-- -->

``` r
x <- trends[, hist(Horntrendz, breaks = 100)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-6.png)<!-- -->

``` r
x <- trends[, hist(Jtulast)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-7.png)<!-- -->

``` r
x <- trends[, hist(Jbetalast)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-8.png)<!-- -->

``` r
x <- trends[, hist(Hornlast)]
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20of%20change-9.png)<!-- -->

#### Turnover calculations are correlated, though less so for Horn

``` r
# are turnover calculations correlated?
ggplot(trends, aes(Jbetatrendrem0, Jtutrendrem0)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 13390 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-1.png)<!-- -->

``` r
ggplot(trends, aes(Jbetatrendrem0, Horntrendrem0)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 13390 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-2.png)<!-- -->

``` r
ggplot(trends, aes(Jbetatrendrem0, Jbetatrendz)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 16894 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 16894 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-3.png)<!-- -->

``` r
ggplot(trends, aes(Jtutrendrem0, Jtutrendz)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 18197 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 18197 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-4.png)<!-- -->

``` r
ggplot(trends, aes(Horntrendrem0, Horntrendz)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 16994 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 16994 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-5.png)<!-- -->

``` r
ggplot(trends, aes(Jbetatrendz, Jtutrendz)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 18197 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 18197 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-6.png)<!-- -->

``` r
ggplot(trends, aes(Jbetatrendz, Horntrendz)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 16994 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 16994 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-7.png)<!-- -->

``` r
ggplot(trends, aes(Jbetatrendrem0, Jbetalast)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 13390 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-8.png)<!-- -->

``` r
ggplot(trends, aes(Jtutrendrem0, Jtulast)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 13390 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-9.png)<!-- -->

``` r
ggplot(trends, aes(Horntrendrem0, Hornlast)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).
    
    ## Warning: Removed 13390 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-10.png)<!-- -->

``` r
ggplot(trends, aes(Jbetalast, Jtulast)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-11.png)<!-- -->

``` r
ggplot(trends, aes(Jbetalast, Hornlast)) +
    geom_point(alpha = 0.3) +
    geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](assemble_turnover_covariates_files/figure-gfm/basic%20pairwise%20graphs%20of%20turnover%20metrics-12.png)<!-- -->
Standardization has the effect of up-weighting some of the smaller
trends, as expected Strange that Jbetatrendz and Jtutrendz are not
really correlated \*last correlated to \*trendrem0. Also constrained to
a parallelogram Jtulast is lower diagonal of Jbetalast, which makes
sense Hornlast correlated but with lots of scatter to Jbetalast

#### Change compared to number of years in time-series

``` r
# number of year
trends[, summary(nyrBT)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   2.000   2.000   4.000   5.538   7.000  97.000

``` r
x <- trends[, hist(nyrBT)]
```

![](assemble_turnover_covariates_files/figure-gfm/change%20vs.%20num%20years-1.png)<!-- -->

``` r
par(mfrow=c(1,3))
trends[, plot(nyrBT, Jtutrendrem0, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Jbetatrendrem0, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Horntrendrem0, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](assemble_turnover_covariates_files/figure-gfm/change%20vs.%20num%20years-2.png)<!-- -->

``` r
trends[, plot(nyrBT, Jtutrendz, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Jbetatrendz, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Horntrendz, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](assemble_turnover_covariates_files/figure-gfm/change%20vs.%20num%20years-3.png)<!-- -->

``` r
trends[, plot(nyrBT, Jtulast, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Jbetalast, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Hornlast, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](assemble_turnover_covariates_files/figure-gfm/change%20vs.%20num%20years-4.png)<!-- -->

#### Change compared to number of species in time-series

``` r
# number of species
trends[, summary(Nspp)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     2.0     6.0    11.0    19.1    24.0  1427.0

``` r
x <- trends[, hist(Nspp)]
```

![](assemble_turnover_covariates_files/figure-gfm/change%20vs.%20number%20of%20species-1.png)<!-- -->

``` r
par(mfrow=c(1,3))
trends[, plot(Nspp, Jtutrendrem0, log = 'x', col = '#00000033')]; abline(h=0)
```

    ## NULL

``` r
trends[, plot(Nspp, Jbetatrendrem0, log = 'x', col = '#00000033')]; abline(h=0)
```

    ## NULL

``` r
trends[, plot(Nspp, Horntrendrem0, log = 'x', col = '#00000033')]; abline(h=0)
```

    ## NULL

![](assemble_turnover_covariates_files/figure-gfm/change%20vs.%20number%20of%20species-2.png)<!-- -->

``` r
trends[, plot(Nspp, Jtutrendz, log = 'x', col = '#00000033')]; abline(h=0)
```

    ## NULL

``` r
trends[, plot(Nspp, Jbetatrendz, log = 'x', col = '#00000033')]; abline(h=0)
```

    ## NULL

``` r
trends[, plot(Nspp, Horntrendz, log = 'x', col = '#00000033')]; abline(h=0)
```

    ## NULL

![](assemble_turnover_covariates_files/figure-gfm/change%20vs.%20number%20of%20species-3.png)<!-- -->

``` r
trends[, plot(Nspp, Jtulast, log = 'x', col = '#00000033')]
```

    ## NULL

``` r
trends[, plot(Nspp, Jbetalast, log = 'x', col = '#00000033')]
```

    ## NULL

``` r
trends[, plot(Nspp, Hornlast, log = 'x', col = '#00000033')]
```

![](assemble_turnover_covariates_files/figure-gfm/change%20vs.%20number%20of%20species-4.png)<!-- -->

    ## NULL

#### Average change compared to \#years, \#species in time-series

``` r
# number of years
ggplot(trends, aes(nyrBT, Jtutrendrem0, color = 'Jtu trend')) +
  geom_smooth() +
  geom_smooth(aes(y = Jbetatrendrem0, color = 'Jbeta trend')) +
  geom_smooth(aes(y = Horntrendrem0, color = 'Horn trend')) +
  scale_x_log10() +
  labs(y = 'Slope') +
  geom_abline(intercept = 0, slope = 0)
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).

![](assemble_turnover_covariates_files/figure-gfm/ave%20change%20vs.%20num%20years-1.png)<!-- -->

``` r
ggplot(trends, aes(nyrBT, Jtutrendz, color = 'Jtu z')) +
  geom_smooth() +
  geom_smooth(aes(y = Jbetatrendz, color = 'Jbeta z')) +
  geom_smooth(aes(y = Horntrendz, color = 'Horn z')) +
  scale_x_log10() +
  labs(y = 'Dissimilarity') +
  geom_abline(intercept = 0, slope = 0)
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 18197 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 16894 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 16994 rows containing non-finite values (stat_smooth).

![](assemble_turnover_covariates_files/figure-gfm/ave%20change%20vs.%20num%20years-2.png)<!-- -->

``` r
ggplot(trends, aes(nyrBT, Jtulast, color = 'Jtu last')) +
  geom_smooth() +
  geom_smooth(aes(y = Jbetalast, color = 'Jbeta last')) +
  geom_smooth(aes(y = Hornlast, color = 'Horn last')) +
  scale_x_log10() +
  labs(y = 'Dissimilarity') +
  geom_abline(intercept = 0, slope = 0)
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](assemble_turnover_covariates_files/figure-gfm/ave%20change%20vs.%20num%20years-3.png)<!-- -->

``` r
# number of species
ggplot(trends, aes(Nspp, Jtutrendrem0, color = 'Jtu trend')) +
  geom_smooth() +
  geom_smooth(aes(y = Jbetatrendrem0, color = 'Jbeta trend')) +
  geom_smooth(aes(y = Horntrendrem0, color = 'Horn trend')) +
  scale_x_log10() +
  labs(y = 'Slope') +
  geom_abline(intercept = 0, slope = 0)
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 13390 rows containing non-finite values (stat_smooth).

![](assemble_turnover_covariates_files/figure-gfm/ave%20change%20vs.%20num%20years-4.png)<!-- -->

``` r
ggplot(trends, aes(Nspp, Jtutrendz, color = 'Jtu z')) +
  geom_smooth() +
  geom_smooth(aes(y = Jbetatrendz, color = 'Jbeta z')) +
  geom_smooth(aes(y = Horntrendz, color = 'Horn z')) +
  scale_x_log10() +
  labs(y = 'Slope') +
  geom_abline(intercept = 0, slope = 0)
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 18197 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 16894 rows containing non-finite values (stat_smooth).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 16994 rows containing non-finite values (stat_smooth).

![](assemble_turnover_covariates_files/figure-gfm/ave%20change%20vs.%20num%20years-5.png)<!-- -->

``` r
ggplot(trends, aes(Nspp, Jtulast, color = 'Jtu last')) +
  geom_smooth() +
  geom_smooth(aes(y = Jbetalast, color = 'Jbeta last')) +
  geom_smooth(aes(y = Hornlast, color = 'Horn last')) +
  scale_x_log10() +
  labs(y = 'Dissimilarity') +
  geom_abline(intercept = 0, slope = 0)
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](assemble_turnover_covariates_files/figure-gfm/ave%20change%20vs.%20num%20years-6.png)<!-- -->

## Compare covariates across realms

``` r
i <- trends[, !duplicated(rarefyID)]; sum(i)
```

    ## [1] 51578

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
beanplot(veg ~ REALM, data = trends[i & REALM !='Marine',], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'NPP', ll = 0.05)
```

![](assemble_turnover_covariates_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->

Marine are in generally warmer locations (seawater doesn’t freeze)
Marine have much lower seasonality. Marine and freshwater have some very
small masses (plankton), but much of dataset is similar to terrestrial.
Marine has a lot of slow, crawling organisms, but land has plants. Land
also has birds (fast).

## Plot turnover vs. temperature

![](assemble_turnover_covariates_files/figure-gfm/plot%20turnover%20vs%20temp%20trend-1.png)<!-- -->

### Time-series length and temperature trend?

``` r
ggplot(trends, aes(temptrend, nyrBT)) +
  geom_point() +
  geom_smooth() +
  scale_y_log10()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 2778 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 2778 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/ts%20length%20and%20temp%20trend-1.png)<!-- -->

## Plot turnover vs. explanatory variables

Lines are ggplot smoother fits.

    ## Warning: Computation failed in `stat_smooth()`:
    ## x has insufficient unique values to support 5 knots: reduce k.

![](assemble_turnover_covariates_files/figure-gfm/plot%20turnover%20v%20explanatory%20vars-1.png)<!-- -->

Strong trends with temperature change, but trends are pretty symmetric
around no trend in temperature, which implies warming or cooling drives
similar degree of community turnover. Some indication of less turnover
for larger organisms (mass) Higher turnover on land with higher
seasonality? More turnover for shorter-lived organisms? No really clear
differences among realms.

### Write out

``` r
write.csv(trends, gzfile('output/turnover_w_covariates.csv.gz'), row.names = FALSE)
```

### Useful variables

``` r
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

### Do the variables look ok?

#### Unscaled

``` r
# histograms to examine
cexmain = 0.6
par(mfrow = c(5,4))
invisible(trends[, hist(minyrBT, main = 'Start year', cex.main = cexmain)])
invisible(trends[, hist(maxyrBT - minyrBT, main = 'Duration (years)', cex.main = cexmain)])
invisible(trends[, hist(nyrBT, main = 'Number of sampled years', cex.main = cexmain)])
invisible(trends[, hist(mass_mean_weight, main = 'Mass (g)', cex.main = cexmain)])
invisible(trends[, hist(speed_mean_weight, main = 'Speed (km/hr)', cex.main = cexmain)])
invisible(trends[, hist(lifespan_mean_weight, main = 'Lifespan (yr)', cex.main = cexmain)])
invisible(trends[, hist(tempave_metab, main = 'Metabolic temperature (°C)', cex.main = cexmain)])
invisible(trends[, hist(consfrac, main = 'Consumers (fraction)', cex.main = cexmain)])
invisible(trends[, hist(endofrac, main = 'Endotherms (fraction)', cex.main = cexmain)])
invisible(trends[, hist(tempave, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(trends[, hist(temptrend, main = 'Temperature trend (°C/yr)', cex.main = cexmain)])
invisible(trends[, hist(seas, main = 'Seasonality (°C)', cex.main = cexmain)])
invisible(trends[, hist(microclim, main = 'Microclimates (°C)', cex.main = cexmain)])
invisible(trends[, hist(Nspp, main = 'Species richness', cex.main = cexmain)])
invisible(trends[, hist(thermal_bias, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(trends[, hist(npp, main = 'Net primary productivity', cex.main = cexmain)])
invisible(trends[, hist(veg, main = 'Vegetation index', cex.main = cexmain)])
invisible(trends[, hist(human_bowler, main = 'Human impact score (Bowler)', cex.main = cexmain)])
invisible(trends[, hist(human_venter, main = 'Human impact score (Venter)', cex.main = cexmain)])
invisible(trends[, hist(human_halpern, main = 'Human impact score (Halpern)', cex.main = cexmain)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20unscaled-1.png)<!-- -->

#### Scaled

``` r
# histograms to examine
cexmain = 0.6
par(mfrow = c(5,4))
invisible(trends[, hist(tempave.sc, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(trends[, hist(tempave_metab.sc, main = 'Metabolic temperature (°C)', cex.main = cexmain)])
invisible(trends[, hist(seas.sc, main = 'Seasonality (°C)', cex.main = cexmain)])
invisible(trends[, hist(microclim.sc, main = 'log Microclimates (°C)', cex.main = cexmain)])
invisible(trends[, hist(temptrend.sc, main = 'Temperature trend (°C/yr)', cex.main = cexmain)])
invisible(trends[, hist(temptrend_abs.sc, main = 'abs(Temperature trend) (°C/yr)', cex.main = cexmain)])
invisible(trends[, hist(mass.sc, main = 'log Mass (g)', cex.main = cexmain)])
invisible(trends[, hist(speed.sc, main = 'log Speed (km/hr)', cex.main = cexmain)])
invisible(trends[, hist(lifespan.sc, main = 'log Lifespan (yr)', cex.main = cexmain)])
invisible(trends[, hist(consumerfrac.sc, main = 'Consumers (fraction)', cex.main = cexmain)])
invisible(trends[, hist(endothermfrac.sc, main = 'Endotherms (fraction)', cex.main = cexmain)])
invisible(trends[, hist(nspp.sc, main = 'log Species richness', cex.main = cexmain)])
invisible(trends[, hist(thermal_bias.sc, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(trends[, hist(npp.sc, main = 'log Net primary productivity', cex.main = cexmain)])
invisible(trends[, hist(veg.sc, main = 'log Vegetation index', cex.main = cexmain)])
invisible(trends[, hist(human_bowler.sc, main = 'log Human impact score (Bowler)', cex.main = cexmain)])
invisible(trends[, hist(human_footprint.sc, main = 'log Human impact score (Venter & Halpern)', cex.main = cexmain)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20scaled-1.png)<!-- -->

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
pairs(formula = ~ tempave.sc + tempave_metab.sc + seas.sc + microclim.sc + temptrend.sc + temptrend_abs.sc + mass.sc + speed.sc + lifespan.sc + consumerfrac.sc + endothermfrac.sc + nspp.sc + thermal_bias.sc + npp.sc + veg.sc + human_bowler.sc + human_footprint.sc, data = trends, gap = 1/10, cex = 0.2, col = '#00000022', lower.panel = panel.cor)
```

![](assemble_turnover_covariates_files/figure-gfm/pairs-1.png)<!-- -->

Mass and lifespan look tightly correlated, but r only 0.56…?
Tempave\_metab and lifespan don’t look tightly correlated, but r= -0.81
Tempave\_metab and speed don’t look tightly correlated, but r= -0.83
Lifespan and speed don’t look tightly correlated, but r = 0.73

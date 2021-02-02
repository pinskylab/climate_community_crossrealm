Turnover vs. temperature plots
================

# Prep

<details>

<summary>Click to expand code</summary>

``` r
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
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

# set realm order
trends[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = FALSE)]

# set up sign of temperature change
trends[, tsign := factor(sign(tempchange))]

# realm that combined Terrestrial and Freshwater, for interacting with human impact
trends[, REALM2 := REALM]
levels(trends$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")

# group Marine invertebrates/plants in with All
trends[, taxa_mod2 := taxa_mod]
trends[taxa_mod == 'Marine invertebrates/plants', taxa_mod2 := 'All']

# calculate duration
trends[, duration := year2 - year1]
```

</details>

# Plot dissimilarity vs. temperature change

## All temporal scales

![](turnover_vs_temperature_plots_files/figure-gfm/plot%20diss%20vs%20temp%20change-1.png)<!-- -->

## Only 1 year differences

![](turnover_vs_temperature_plots_files/figure-gfm/plot%20diss%20vs%20temp%20change%201%20yr-1.png)<!-- -->

## Only 10 year differences

![](turnover_vs_temperature_plots_files/figure-gfm/plot%20diss%20vs%20temp%20change%2010%20yr-1.png)<!-- -->

## Only 20 year differences

![](turnover_vs_temperature_plots_files/figure-gfm/plot%20diss%20vs%20temp%20change%2020%20yr-1.png)<!-- -->

# Focus on marine

## All temporal scales

![](turnover_vs_temperature_plots_files/figure-gfm/plot%20diss%20vs%20temp%20change%20marine-1.png)<!-- -->

## Only 1 year differences

![](turnover_vs_temperature_plots_files/figure-gfm/plot%20diss%20vs%20temp%20change%201%20yr%20marine-1.png)<!-- -->

## Only 10 year differences

![](turnover_vs_temperature_plots_files/figure-gfm/plot%20diss%20vs%20temp%20change%2010%20yr%20marine-1.png)<!-- -->

# Temperature trend and time-series length

``` r
ggplot(trends, aes(year2 - year1, tempchange, color = REALM, group = REALM)) +
  geom_point(alpha = 0.2) +
  geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 15988 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 15988 rows containing missing values (geom_point).

![](turnover_vs_temperature_plots_files/figure-gfm/ts%20length%20and%20temp%20trend-1.png)<!-- -->

# Summary plots and stats

## Average community dissimilarity

``` r
trends[abs(tempchange) >= 1, tempchangetext1 := 'Change >=1']
trends[abs(tempchange) <= 0.5, tempchangetext1 := 'Stable <=0.5']
trends[tempchange <= -1, tempchangetext2 := 'Cool <= -1']
trends[tempchange >= 1, tempchangetext2 := 'Warm >= 1']

# reshape to long format
measurenms <- c('Jbeta', 'Jtu', 'Jne', 'Horn')
idnms <- setdiff(names(trends), measurenms)
trends2 <- melt(trends, id.vars = idnms, measure.vars = measurenms)

# changing vs. stable
trendsum1 <- trends2[!is.na(tempchangetext1), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = tempchangetext1, type = variable, REALM = REALM)] # turnover per year for locations changing temperature
trendsum1[, duration := NA_integer_]


# warming vs. cooling
trendsum2 <- trends2[!is.na(tempchangetext2), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = tempchangetext2, type = variable, REALM = REALM)] # inc. direction
trendsum2[, duration := NA_integer_]

# 1, 10, or 20 year intervals
trendsum3 <- trends2[!is.na(tempchangetext1) & duration %in% c(1, 10, 20), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = tempchangetext1, type = variable, duration, REALM = REALM)] # inc. time interval
setorder(trendsum3, type, duration, text)

# combine
trendsum4 <- rbind(trendsum1, trendsum2, trendsum3)
setorder(trendsum4, type, text, duration)

write.csv(trendsum4, file = 'output/trendsummary.csv', row.names = FALSE)

trendsum4
```

    ##              text  type       REALM       ave           se      n duration
    ##   1:   Change >=1 Jbeta      Marine 0.7013735 0.0006182196 145870       NA
    ##   2:   Change >=1 Jbeta Terrestrial 0.3449956 0.0005597050  78188       NA
    ##   3:   Change >=1 Jbeta  Freshwater 0.5320614 0.0061905188   2247       NA
    ##   4:   Change >=1 Jbeta      Marine 0.6086006 0.0037583737   5267        1
    ##   5:   Change >=1 Jbeta Terrestrial 0.3473984 0.0025307024   5195        1
    ##  ---                                                                      
    ## 116: Stable <=0.5  Horn Terrestrial 0.2104779 0.0033216777   1771       20
    ## 117: Stable <=0.5  Horn  Freshwater 0.5203639 0.0498359061     56       20
    ## 118:    Warm >= 1  Horn      Marine 0.6540410 0.0010744683 100409       NA
    ## 119:    Warm >= 1  Horn Terrestrial 0.1909579 0.0007059340  55767       NA
    ## 120:    Warm >= 1  Horn  Freshwater 0.5054501 0.0088474184   1564       NA

## Plots of turnover rates binned by warming rates

![](turnover_vs_temperature_plots_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-1.png)<!-- -->![](turnover_vs_temperature_plots_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-2.png)<!-- -->![](turnover_vs_temperature_plots_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-3.png)<!-- -->

# Plot dissimilarity binned by tempchange vs duration

## Jaccard turnover dissimilarity

Very slow to render (10 min?)

``` r
ggplot(trends[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change')
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20temperature%20trend%20vs.%20duration-1.png)<!-- -->

## By realm

``` r
terr <- trends[REALM == 'Terrestrial', ]
p1 <- ggplot(terr[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = terr[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = terr[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = terr[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change', title = 'Terrestrial')

mar <- trends[REALM == 'Marine', ]
p2 <- ggplot(mar[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = mar[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = mar[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = mar[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change', title = 'Marine')

fre <- trends[REALM == 'Freshwater', ]
p3 <- ggplot(fre[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = fre[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = fre[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = fre[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change', title = 'Freshwater')

grid.arrange(p1, p2, p3, ncol = 3)
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20temperature%20trend%20vs.%20duration%20by%20realm-1.png)<!-- -->

## By body size

``` r
small <- trends[mass_mean_weight < 1, ]
p1 <- ggplot(small[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = small[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = small[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = small[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change', title = '<1 g')

med <- trends[mass_mean_weight < 10000 & mass_mean_weight >= 1, ]
p2 <- ggplot(med[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = med[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = med[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = med[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change', title = '1g to 10kg')

lg <- trends[mass_mean_weight >= 10000, ]
p3 <- ggplot(lg[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = lg[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = lg[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = lg[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change', title = '>10kg')

grid.arrange(p1, p2, p3, ncol = 3)
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20temperature%20trend%20vs.%20duration%20by%20body%20size-1.png)<!-- -->

## By endo/ectotherm

``` r
ecto <- trends[endofrac <= 0.5, ]
p1 <- ggplot(ecto[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = ecto[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = ecto[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = ecto[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change', title = 'Mostly ectotherms')

endo <- trends[endofrac > 0.5, ]
p2 <- ggplot(endo[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(aes(color = '<=0.05')) +
  geom_smooth(data = endo[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], aes(color = '0.05-0.1')) +
  geom_smooth(data = endo[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], aes(color = '0.1-0.5')) +
  geom_smooth(data = endo[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change', title = 'Mostly endotherms')

grid.arrange(p1, p2, ncol = 2)
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20temperature%20trend%20vs.%20duration%20by%20endofrac-1.png)<!-- -->

# Plot dissimilarity binned by duration vs tempchange

## Jaccard turnover dissimilarity

``` r
ggplot(trends[duration == 1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = trends[duration ==5,], aes(color = '5')) +
  geom_smooth(data = trends[duration==10,], aes(color = '10')) +
  geom_smooth(data = trends[duration==15,], aes(color = '15')) +
  geom_smooth(data = trends[duration>15,], aes(color = '>15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration')
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20tempchange%20by%20duration-1.png)<!-- -->

## By realm

``` r
terr <- trends[REALM == 'Terrestrial', ]
p1 <- ggplot(terr[duration==1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = terr[duration==5,], aes(color = '5')) +
  geom_smooth(data = terr[duration==10,], aes(color = '10')) +
  geom_smooth(data = terr[duration==15,], aes(color = '15')) +
  geom_smooth(data = terr[duration>15,], aes(color = '>15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration', title = 'Terrestrial')

mar <- trends[REALM == 'Marine', ]
p2 <- ggplot(mar[duration==1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = mar[duration==5,], aes(color = '5')) +
  geom_smooth(data = mar[duration==10,], aes(color = '10')) +
  geom_smooth(data = mar[duration==15,], aes(color = '15')) +
  geom_smooth(data = mar[duration>15,], aes(color = '>15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration', title = 'Marine')

fre <- trends[REALM == 'Freshwater', ]
p3 <- ggplot(fre[duration==1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = fre[duration==5,], aes(color = '5')) +
  geom_smooth(data = fre[duration==10,], aes(color = '10')) +
  geom_smooth(data = fre[duration==15,], aes(color = '15')) +
  geom_smooth(data = fre[duration>15,], aes(color = '>15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration', title = 'Freshwater')

grid.arrange(p1, p2, p3, ncol = 3)
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20tempchange%20by%20duration%20by%20realm-1.png)<!-- -->

## By body size

``` r
small <- trends[mass_mean_weight < 1, ]
p1 <- ggplot(small[duration==1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = small[duration==5,], aes(color = '5')) +
  geom_smooth(data = small[duration==10,], aes(color = '10')) +
  geom_smooth(data = small[duration==15,], aes(color = '15')) +
  geom_smooth(data = small[duration>15,], aes(color = '>15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration', title = '<1 g')

med <- trends[mass_mean_weight < 10000 & mass_mean_weight >= 1, ]
p2 <- ggplot(med[duration==1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = med[duration==5,], aes(color = '5')) +
  geom_smooth(data = med[duration==10,], aes(color = '10')) +
  geom_smooth(data = med[duration==15,], aes(color = '15')) +
  geom_smooth(data = med[duration>15,], aes(color = '>15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration', title = '1g to 10kg')

lg <- trends[mass_mean_weight >= 10000, ]
p3 <- ggplot(lg[duration==1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = lg[duration==5,], aes(color = '5')) +
  geom_smooth(data = lg[duration==10,], aes(color = '10')) +
  geom_smooth(data = lg[duration==15,], aes(color = '15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration', title = '>10kg')

grid.arrange(p1, p2, p3, ncol = 3)
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20temperature%20trend%20by%20duration%20by%20body%20size-1.png)<!-- -->

## By endo/ectotherm

``` r
ecto <- trends[endofrac <= 0.5, ]
p1 <- ggplot(ecto[duration==1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = ecto[duration==5,], aes(color = '5')) +
  geom_smooth(data = ecto[duration==10,], aes(color = '10')) +
  geom_smooth(data = ecto[duration==15,], aes(color = '15')) +
  geom_smooth(data = ecto[duration>15,], aes(color = '>15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration', title = 'Mostly ectotherms')

endo <- trends[endofrac > 0.5, ]
p2 <- ggplot(endo[duration==1, ], aes(x=abs(tempchange), y=Jtu)) +
  geom_smooth(aes(color = '1')) +
  geom_smooth(data = endo[duration==5,], aes(color = '5')) +
  geom_smooth(data = endo[duration==10,], aes(color = '10')) +
  geom_smooth(data = endo[duration==15,], aes(color = '15')) +
  geom_smooth(data = ecto[duration>15,], aes(color = '>15')) +
  labs(y = 'Jaccard dissimilarity', color = 'Duration', title = 'Mostly endotherms')

grid.arrange(p1, p2, ncol = 2)
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20temperature%20trend%20by%20duration%20by%20endofrac-1.png)<!-- -->
\#\# By realm and body size

``` r
trends[duration == 1, dur_bin := '[1] yr']
trends[duration > 1 & duration <= 5, dur_bin := '(1-5] yr']
trends[duration > 5 & duration <= 10, dur_bin := '(5-10] yr']
trends[duration > 10 & duration <= 20, dur_bin := '(10-20] yr']
trends[duration > 20, dur_bin := '>20 yr']
trends[mass_mean_weight <= 1, mass_bin := '<=1 g']
trends[mass_mean_weight > 1 & mass_mean_weight <= 10000, mass_bin := '(1 g-10 kg]']
trends[mass_mean_weight > 10000, mass_bin := '>10 kg']
ggplot(trends[!is.na(mass_bin), ], aes(x=abs(tempchange), y=Jtu, color = dur_bin, group = dur_bin)) +
  geom_smooth() +
  facet_grid(mass_bin~REALM, scales = 'free') +
  labs(y = 'Jaccard dissimilarity', title = '')
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20tempchange%20by%20duration%20mass-1.png)<!-- -->

## By realm, body size AND endo/ecto

``` r
trends[duration == 1, dur_bin := '[1] yr']
trends[duration > 1 & duration <= 5, dur_bin := '(1-5] yr']
trends[duration > 5 & duration <= 10, dur_bin := '(5-10] yr']
trends[duration > 10 & duration <= 20, dur_bin := '(10-20] yr']
trends[duration > 20, dur_bin := '>20 yr']
trends[mass_mean_weight <= 1, mass_bin := '<=1 g']
trends[mass_mean_weight > 1 & mass_mean_weight <= 10000, mass_bin := '(1 g-10 kg]']
trends[mass_mean_weight > 10000, mass_bin := '>10 kg']
trends[endofrac > 0.5, endoTF := 'endo']
trends[endofrac <= 0.5, endoTF := 'ecto']
ggplot(trends[!is.na(mass_bin), ], aes(x=abs(tempchange), y=Jtu, color = dur_bin, group = dur_bin)) +
  geom_smooth() +
  facet_grid(mass_bin + endoTF~REALM, scales = 'free') +
  labs(y = 'Jaccard dissimilarity', title = '')
```

![](turnover_vs_temperature_plots_files/figure-gfm/diss%20vs.%20tempchange%20by%20duration%20mass%20endofrac-1.png)<!-- -->

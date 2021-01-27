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
trends[abs(tempchange) >= 0.5, tempchangetext1 := 'Change >=0.5']
trends[abs(tempchange) <= 0.05, tempchangetext1 := 'Stable <=0.05']
trends[tempchange <= -0.5, tempchangetext2 := 'Cool <= -0.5']
trends[tempchange >= 0.5, tempchangetext2 := 'Warm >= 0.5']

# reshape to long format
measurenms <- c('Jbeta', 'Jtu', 'Jne', 'Horn')
idnms <- setdiff(names(trends), measurenms)
trends2 <- melt(trends, id.vars = idnms, measure.vars = measurenms)

# changing vs. stable
trendsum1 <- trends2[!is.na(tempchangetext1), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = tempchangetext1, type = variable)] # turnover per year for locations changing temperature
trendsum1[, duration := NA_integer_]


# warming vs. cooling
trendsum2 <- trends2[!is.na(tempchangetext2), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = tempchangetext2, type = variable)] # inc. direction
trendsum2[, duration := NA_integer_]

# 1, 10, or 20 year intervals
trendsum3 <- trends2[!is.na(tempchangetext1) & duration %in% c(1, 10, 20), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = tempchangetext1, type = variable, duration)] # inc. time interval
setorder(trendsum3, type, duration, text)

# combine
trendsum4 <- rbind(trendsum1, trendsum2, trendsum3)
setorder(trendsum4, type, text, duration)

write.csv(trendsum4, file = 'output/trendsummary.csv', row.names = FALSE)

trendsum4
```

    ##              text  type       ave           se      n duration
    ##  1:  Change >=0.5 Jbeta 0.6002095 0.0003610318 589418       NA
    ##  2:  Change >=0.5 Jbeta 0.5359723 0.0014581158  39169        1
    ##  3:  Change >=0.5 Jbeta 0.5723561 0.0018278849  22423       10
    ##  4:  Change >=0.5 Jbeta 0.6402557 0.0027132290   9958       20
    ##  5:  Cool <= -0.5 Jbeta 0.6066283 0.0006551779 188968       NA
    ##  6: Stable <=0.05 Jbeta 0.6051975 0.0009785149  83940       NA
    ##  7: Stable <=0.05 Jbeta 0.5625294 0.0027096435  11618        1
    ##  8: Stable <=0.05 Jbeta 0.6399932 0.0046113191   3774       10
    ##  9: Stable <=0.05 Jbeta 0.5934181 0.0105802353    698       20
    ## 10:   Warm >= 0.5 Jbeta 0.5971806 0.0004321189 400450       NA
    ## 11:  Change >=0.5   Jtu 0.4459941 0.0004281782 589418       NA
    ## 12:  Change >=0.5   Jtu 0.3797385 0.0016539601  39169        1
    ## 13:  Change >=0.5   Jtu 0.4161877 0.0021343203  22423       10
    ## 14:  Change >=0.5   Jtu 0.4865784 0.0032528815   9958       20
    ## 15:  Cool <= -0.5   Jtu 0.4662617 0.0007791891 188968       NA
    ## 16: Stable <=0.05   Jtu 0.4262683 0.0012346728  83940       NA
    ## 17: Stable <=0.05   Jtu 0.3731232 0.0033285498  11618        1
    ## 18: Stable <=0.05   Jtu 0.4690068 0.0060822619   3774       10
    ## 19: Stable <=0.05   Jtu 0.4478041 0.0123191172    698       20
    ## 20:   Warm >= 0.5   Jtu 0.4364300 0.0005111580 400450       NA
    ## 21:  Change >=0.5   Jne 0.1468606 0.0002732904 589418       NA
    ## 22:  Change >=0.5   Jne 0.1470520 0.0010593737  39169        1
    ## 23:  Change >=0.5   Jne 0.1500010 0.0014011557  22423       10
    ## 24:  Change >=0.5   Jne 0.1453703 0.0021878906   9958       20
    ## 25:  Cool <= -0.5   Jne 0.1335554 0.0004590099 188968       NA
    ## 26: Stable <=0.05   Jne 0.1725461 0.0007931707  83940       NA
    ## 27: Stable <=0.05   Jne 0.1836436 0.0021939665  11618        1
    ## 28: Stable <=0.05   Jne 0.1671645 0.0036860889   3774       10
    ## 29: Stable <=0.05   Jne 0.1360153 0.0075265148    698       20
    ## 30:   Warm >= 0.5   Jne 0.1531392 0.0003385034 400450       NA
    ## 31:  Change >=0.5  Horn 0.5179450 0.0004985121 573102       NA
    ## 32:  Change >=0.5  Horn 0.4294316 0.0019581135  38120        1
    ## 33:  Change >=0.5  Horn 0.4816639 0.0025324227  21867       10
    ## 34:  Change >=0.5  Horn 0.5768760 0.0037623253   9634       20
    ## 35:  Cool <= -0.5  Horn 0.5261015 0.0008946141 183776       NA
    ## 36: Stable <=0.05  Horn 0.5075724 0.0013575413  81011       NA
    ## 37: Stable <=0.05  Horn 0.4341489 0.0036774869  11137        1
    ## 38: Stable <=0.05  Horn 0.5463031 0.0064589171   3634       10
    ## 39: Stable <=0.05  Horn 0.5497366 0.0146317342    647       20
    ## 40:   Warm >= 0.5  Horn 0.5140948 0.0006000470 389326       NA
    ##              text  type       ave           se      n duration

## Plots of turnover rates binned by warming rates

``` r
p1 <- ggplot(trends[!is.na(tempchangetext1), ], aes(tempchangetext1, Horn)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = 'grey') +
  labs(x = '', y = 'Horn dissimilarity', tag = 'A', title = 'Rate of temperature change') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10)) +
  geom_abline(intercept = 0, slope = 0)

p2 <- ggplot(trends[!is.na(tempchangetext2), ], aes(tempchangetext2, Horn)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = 'grey') +
  labs(x = '', y = 'Horn dissimilarity', tag = 'B', title = 'Rate of temperature change') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10)) +
  geom_abline(intercept = 0, slope = 0)

p3 <- ggplot(trendsum4[type %in% c('Jtu', 'Jbeta', 'Horn') & is.na(duration)], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width= 0.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Dissimilarity', title = 'All time intervals') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p4 <- ggplot(trendsum4[type %in% c('Jtu', 'Jbeta', 'Horn') & duration==1, ], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width= 0.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Dissimilarity', title = '1-year interval') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p5 <- ggplot(trendsum4[type %in% c('Jtu', 'Jbeta', 'Horn') & duration==10, ], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width= 0.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Dissimilarity', title = '10-year interval') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p1
```

![](turnover_vs_temperature_plots_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-1.png)<!-- -->

``` r
p2
```

![](turnover_vs_temperature_plots_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-2.png)<!-- -->

``` r
p3
```

![](turnover_vs_temperature_plots_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-3.png)<!-- -->

``` r
p4
```

![](turnover_vs_temperature_plots_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-4.png)<!-- -->

``` r
p5
```

![](turnover_vs_temperature_plots_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-5.png)<!-- -->

# Plot turnover binned by warming rates vs duration

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

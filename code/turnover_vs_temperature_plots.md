Turnover vs. temperature plots
================

# Prep

<details>

<summary>Click to expand code</summary>

``` r
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(nlme) # for ME models
library(maps) # for map
library(gridExtra) # to combine ggplots together
library(grid) # to combine ggplots together
library(RColorBrewer)
library(MASS) # for stepAIC
library(piecewiseSEM) # for rsquared() for nlme models
```

    ## Registered S3 methods overwritten by 'lme4':
    ##   method                          from
    ##   cooks.distance.influence.merMod car 
    ##   influence.merMod                car 
    ##   dfbeta.influence.merMod         car 
    ##   dfbetas.influence.merMod        car

    ## 
    ##   This is piecewiseSEM version 2.1.0.
    ## 
    ## 
    ##   Questions or bugs can be addressed to <LefcheckJ@si.edu>.

``` r
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

    ##              text  type        ave           se      n duration
    ##  1:  Change >=0.5 Jbeta 0.52130794 0.0004472551 339366       NA
    ##  2:  Change >=0.5 Jbeta 0.46748252 0.0017228371  23728        1
    ##  3:  Change >=0.5 Jbeta 0.50503948 0.0023030548  13346       10
    ##  4:  Change >=0.5 Jbeta 0.52894590 0.0030939929   6219       20
    ##  5:  Cool <= -0.5 Jbeta 0.51878078 0.0007747764 112999       NA
    ##  6: Stable <=0.05 Jbeta 0.53683708 0.0016887871  28142       NA
    ##  7: Stable <=0.05 Jbeta 0.52019540 0.0049632196   3523        1
    ##  8: Stable <=0.05 Jbeta 0.51097127 0.0088768415   1014       10
    ##  9: Stable <=0.05 Jbeta 0.54462493 0.0161558529    362       20
    ## 10:   Warm >= 0.5 Jbeta 0.52256946 0.0005477171 226367       NA
    ## 11:  Change >=0.5   Jtu 0.39697280 0.0004853488 339366       NA
    ## 12:  Change >=0.5   Jtu 0.34374062 0.0018412668  23728        1
    ## 13:  Change >=0.5   Jtu 0.38357734 0.0025106579  13346       10
    ## 14:  Change >=0.5   Jtu 0.39976299 0.0033330986   6219       20
    ## 15:  Cool <= -0.5   Jtu 0.39646555 0.0008493252 112999       NA
    ## 16: Stable <=0.05   Jtu 0.41079301 0.0018963926  28142       NA
    ## 17: Stable <=0.05   Jtu 0.39673059 0.0055696448   3523        1
    ## 18: Stable <=0.05   Jtu 0.39575557 0.0098195132   1014       10
    ## 19: Stable <=0.05   Jtu 0.44369265 0.0177423259    362       20
    ## 20:   Warm >= 0.5   Jtu 0.39722601 0.0005913475 226367       NA
    ## 21:  Change >=0.5   Jne 0.11975196 0.0002839410 339366       NA
    ## 22:  Change >=0.5   Jne 0.12011161 0.0010987511  23728        1
    ## 23:  Change >=0.5   Jne 0.11791130 0.0014279691  13346       10
    ## 24:  Change >=0.5   Jne 0.12430301 0.0021485320   6219       20
    ## 25:  Cool <= -0.5   Jne 0.11793869 0.0004867605 112999       NA
    ## 26: Stable <=0.05   Jne 0.12240402 0.0010540973  28142       NA
    ## 27: Stable <=0.05   Jne 0.12118209 0.0029777338   3523        1
    ## 28: Stable <=0.05   Jne 0.11249637 0.0052166820   1014       10
    ## 29: Stable <=0.05   Jne 0.09560664 0.0075501827    362       20
    ## 30:   Warm >= 0.5   Jne 0.12065711 0.0003495027 226367       NA
    ## 31:  Change >=0.5  Horn 0.41598316 0.0006042510 327079       NA
    ## 32:  Change >=0.5  Horn 0.34528808 0.0022800180  22980        1
    ## 33:  Change >=0.5  Horn 0.39639168 0.0030767514  12923       10
    ## 34:  Change >=0.5  Horn 0.42733560 0.0042280176   5921       20
    ## 35:  Cool <= -0.5  Horn 0.41133237 0.0010486230 109598       NA
    ## 36: Stable <=0.05  Horn 0.44018401 0.0022642391  26629       NA
    ## 37: Stable <=0.05  Horn 0.41125764 0.0066359834   3302        1
    ## 38: Stable <=0.05  Horn 0.40920880 0.0118359968    965       10
    ## 39: Stable <=0.05  Horn 0.47199668 0.0200272914    359       20
    ## 40:   Warm >= 0.5  Horn 0.41832689 0.0007392664 217481       NA
    ##              text  type        ave           se      n duration

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

## Plots of turnover rates binned by warming rates and duration

### Jaccard turnover dissimilarity

Very slow to render (10 min?)

``` r
ggplot(trends[abs(tempchange) <= 0.05, ], aes(x=duration, y=Jtu)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(tempchange) <= 0.1 & abs(tempchange) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(tempchange) <= 0.5 & abs(tempchange) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(tempchange) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard dissimilarity', color = '°C change')
```

![](turnover_vs_temperature_plots_files/figure-gfm/slope%20vs.%20temperature%20trend%20vs.%20duration-1.png)<!-- -->

Drivers of variation in the community response to temperature change
across realms
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
trends[, tsign := factor(sign(temptrend))]

# realm that combined Terrestrial and Freshwater, for interacting with human impact
trends[, REALM2 := REALM]
levels(trends$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")

# group Marine invertebrates/plants in with All
trends[, taxa_mod2 := taxa_mod]
trends[taxa_mod == 'Marine invertebrates/plants', taxa_mod2 := 'All']

# calculate duration
trends[, duration := maxyrBT - minyrBT + 1]

# trim to data with >= 3 yrs
trends <- trends[nyrBT >= 3, ]
```

## Log-transform some variables, then center and scale.

``` r
trends[, tempave.sc := scale(tempave)]
trends[, tempave_metab.sc := scale(tempave_metab)]
trends[, seas.sc := scale(seas)]
trends[, microclim.sc := scale(log(microclim))]
trends[, temptrend.sc := scale(temptrend, center = FALSE)] # do not center
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
trends[, duration.sc := scale(log(duration))]
trends[, human_bowler.sc := scale(log(human_bowler+1)), by = REALM2] # separate scaling by realm
trends[REALM2 == 'TerrFresh', human_footprint.sc := scale(log(human_venter+1))]
trends[REALM2 == 'Marine', human_footprint.sc := scale(log(human_halpern))]
```

</details>

# Summary plots and stats

## Examine how many data points are available

### Just turnover

``` r
cat('Overall # time-series: ', nrow(trends), '\n')
```

    ## Overall # time-series:  32581

``` r
cat('# studies: ', trends[, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  261

``` r
cat('Data points: ', trends[, sum(nyrBT)], '\n')
```

    ## Data points:  225221

``` r
trends[, table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##         482       29410        2689

``` r
trends[, table(taxa_mod)]
```

    ## taxa_mod
    ##                         All                  Amphibians                     Benthos                       Birds                        Fish               Invertebrates                     Mammals Marine invertebrates/plants                       Plant                    Reptiles 
    ##                        1163                         329                        4294                        4934                       19921                        1163                         458                         104                         211                           4

``` r
trends[, table(taxa_mod, REALM)]
```

    ##                              REALM
    ## taxa_mod                      Freshwater Marine Terrestrial
    ##   All                                  0   1160           3
    ##   Amphibians                           2      0         327
    ##   Benthos                              0   4294           0
    ##   Birds                                0   2848        2086
    ##   Fish                               465  19456           0
    ##   Invertebrates                       13   1072          78
    ##   Mammals                              0    417          41
    ##   Marine invertebrates/plants          0    104           0
    ##   Plant                                1     59         151
    ##   Reptiles                             1      0           3

### With all covariates (Bowler for human)

``` r
# the cases we can compare
apply(trends[, .(Jtutrendrem0, REALM, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)], MARGIN = 2, FUN = function(x) sum(!is.na(x)))
```

    ##     Jtutrendrem0            REALM       tempave.sc tempave_metab.sc          seas.sc     microclim.sc     temptrend.sc          mass.sc         speed.sc      lifespan.sc  consumerfrac.sc endothermfrac.sc          nspp.sc  thermal_bias.sc           npp.sc           veg.sc  human_bowler.sc 
    ##            32581            32581            31307            31307            31307            32365            31307            32536            32555            31734            32581            32581            32581            30945            32498            32499            32581

``` r
i <- trends[, complete.cases(Jtutrendrem0, tempave.sc, tempave_metab.sc, seas.sc, microclim.sc, temptrend.sc, mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
cat('Overall # time-series: ', sum(i), '\n')
```

    ## Overall # time-series:  30756

``` r
cat('# studies: ', trends[i, length(unique(STUDY_ID))], '\n')
```

    ## # studies:  198

``` r
cat('Data points: ', trends[i, sum(nyrBT)], '\n')
```

    ## Data points:  210294

``` r
trends[i, table(REALM)]
```

    ## REALM
    ##  Freshwater      Marine Terrestrial 
    ##         470       28054        2232

``` r
trends[i, table(taxa_mod)]
```

    ## taxa_mod
    ##                         All                  Amphibians                     Benthos                       Birds                        Fish               Invertebrates                     Mammals Marine invertebrates/plants                       Plant                    Reptiles 
    ##                        1138                          12                        4272                        4378                       19172                        1086                         453                         104                         139                           2

``` r
trends[i, table(taxa_mod, REALM)]
```

    ##                              REALM
    ## taxa_mod                      Freshwater Marine Terrestrial
    ##   All                                  0   1136           2
    ##   Amphibians                           2      0          10
    ##   Benthos                              0   4272           0
    ##   Birds                                0   2364        2014
    ##   Fish                               459  18713           0
    ##   Invertebrates                        8   1014          64
    ##   Mammals                              0    417          36
    ##   Marine invertebrates/plants          0    104           0
    ##   Plant                                1     34         104
    ##   Reptiles                             0      0           2

## Average rates of turnover

``` r
trends[abs(temptrend) >= 0.5, temptrendtext1 := 'Changing >=0.5']
trends[abs(temptrend) <= 0.05, temptrendtext1 := 'Stable <=0.05']
trends[temptrend <= -0.5, temptrendtext2 := 'Cooling <= -0.5']
trends[temptrend >= 0.5, temptrendtext2 := 'Warming >= 0.5']
trends[abs(temptrend) >= 0.5 & abs(temptrend) < 1, temptrendtext3 := 'Slow Changing >=0.5 & <1']
trends[abs(temptrend) >= 1, temptrendtext3 := 'Fast Changing >=1']

# reshape to long format
measurenms <- c('Jtutrendrem0', 'Jtutrendrem0_se', 'Jbetatrendrem0', 'Jbetatrendrem0_se', 'Horntrendrem0', 'Horntrendrem0_se', 'Jtutrendz', 'Jbetatrendz', 'Horntrendz', 'Jtulast', 'Jbetalast', 'Hornlast', 'Jtuexp', 'Jbetaexp', 'Hornexp', 'Jtumm', 'Jbetamm', 'Hornmm')
idnms <- setdiff(names(trends), measurenms)
trends2 <- melt(trends, id.vars = idnms, measure.vars = measurenms)

# changing vs. stable
trendsum1 <- trends2[!is.na(temptrendtext1), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = temptrendtext1, type = variable)] # turnover per year for locations changing temperature


# warming vs. cooling
trendsum2 <- trends2[!is.na(temptrendtext2), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = temptrendtext2, type = variable)] # inc. direction

# slow vs. fast changing
trendsum3 <- trends2[!is.na(temptrendtext3), 
                    .(ave = mean(value, na.rm=TRUE), 
                      se = sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value))),
                      n = sum(!is.na(value))),
                    by = .(text = temptrendtext3, type = variable)] # inc. rate

# combine
trendsum4 <- rbind(trendsum1, trendsum2, trendsum3)
setorder(trendsum4, type, text)

write.csv(trendsum4, file = 'output/trendsummary.csv', row.names = FALSE)

trendsum4
```

    ##                          text         type          ave           se     n
    ##   1:           Changing >=0.5 Jtutrendrem0  0.017411084 0.0107816733   355
    ##   2:          Cooling <= -0.5 Jtutrendrem0  0.013137429 0.0123577611   247
    ##   3:        Fast Changing >=1 Jtutrendrem0  0.045305532 0.0295638382    39
    ##   4: Slow Changing >=0.5 & <1 Jtutrendrem0  0.013968415 0.0115500736   316
    ##   5:            Stable <=0.05 Jtutrendrem0  0.005066829 0.0005795107 22045
    ##  ---                                                                      
    ## 104:          Cooling <= -0.5       Hornmm  1.211630617 0.1293621221   246
    ## 105:        Fast Changing >=1       Hornmm  1.398230132 0.3094748502    39
    ## 106: Slow Changing >=0.5 & <1       Hornmm  1.062721270 0.1104267346   315
    ## 107:            Stable <=0.05       Hornmm 10.953692562 0.4532973081 21954
    ## 108:           Warming >= 0.5       Hornmm  0.844694848 0.1698202469   108

## Plots of turnover rates binned by warming rates

``` r
p1 <- ggplot(trends[!is.na(temptrendtext1), ], aes(temptrendtext1, Horntrendrem0)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = 'grey') +
  labs(x = '', y = 'Slope (Horn)', tag = 'A', title = 'Rate of temperature change') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10)) +
  geom_abline(intercept = 0, slope = 0)

p2 <- ggplot(trends[!is.na(text), ], aes(temptrendtext2, Horntrendrem0)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = 'grey') +
  labs(x = '', y = 'Slope (Horn)', tag = 'B', title = 'Rate of temperature change') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10)) +
  geom_abline(intercept = 0, slope = 0)

p3 <- ggplot(trendsum4[type %in% c('Jtutrendrem0', 'Jbetatrendrem0', 'Horntrendrem0')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Slope', title = 'Slope of dissimilarity vs. time') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  coord_cartesian(ylim = c(-0.01,0.08))

p4 <- ggplot(trendsum4[type %in% c('Jtutrendz', 'Jbetatrendz', 'Horntrendz')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Standardized slope', title = 'Null model standardized slope') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p5 <- ggplot(trendsum4[type %in% c('Jtulast', 'Jbetalast', 'Hornlast')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Dissimilarity', title = 'First to last year') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p6 <- ggplot(trendsum4[type %in% c('Jtuexp', 'Jbetaexp', 'Hornexp')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Years', title = 'Exponential half-saturation') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

p7 <- ggplot(trendsum4[type %in% c('Jtumm', 'Jbetamm', 'Hornmm')], 
             aes(text, ave, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
  labs(x = '', y = 'Years', title = 'MM half-saturation') +
  geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 2, 
             layout_matrix = rbind(c(1,2), c(3,3), c(4,4), c(5,5), c(6,6), c(7,7)),
             heights=c(unit(0.1, "npc"), unit(0.18, "npc"), 
                       unit(0.18, "npc"), unit(0.18, "npc"), 
                       unit(0.18, "npc"), unit(0.18, "npc")))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/turnover%20vs.%20temperature%20violin%20plot-1.png)<!-- -->

## Plots of turnover rates binned by warming rates and duration

### Slope of dissimilarity

``` r
p1 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Jtutrendrem0)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard turnover slope', color = '°C/year')

p2 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Jbetatrendrem0)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Jaccard total slope', color = '°C/year')

p3 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Horntrendrem0)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Horn slope', color = '°C/year')

grid.arrange(p1, p2, p3, ncol = 3)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/slope%20vs.%20temperature%20trend%20vs.%20duration-1.png)<!-- -->

### Standardized slope of dissimilarity

``` r
p1 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Jtutrendz)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Standardized slope', color = '°C/year', title = 'Jaccard turnover')

p2 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Jbetatrendz)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Standardized slope', color = '°C/year', title = 'Jaccard total')

p3 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Horntrendz)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Standardized slope', color = '°C/year', title = 'Morisita horn')

grid.arrange(p1, p2, p3, ncol = 3)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/std%20slope%20vs.%20temperature%20trend%20vs.%20duration-1.png)<!-- -->

### First-last dissimilarity

``` r
p1 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Jtulast)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Standardized slope', color = '°C/year', title = 'Jaccard turnover')

p2 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Jbetalast)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Standardized slope', color = '°C/year', title = 'Jaccard total')

p3 <- ggplot(trends[abs(temptrend) <= 0.05, ], aes(x=duration, y=Hornlast)) +
  geom_smooth(method = 'gam', aes(color = '<=0.05')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.1 & abs(temptrend) > 0.05,], method = 'loess', aes(color = '0.05-0.1')) +
  geom_smooth(data = trends[abs(temptrend) <= 0.5 & abs(temptrend) > 0.1,], method = 'loess', aes(color = '0.1-0.5')) +
  geom_smooth(data = trends[abs(temptrend) >= 0.5,], method = 'lm', aes(color = '>=0.5')) +
  scale_x_log10() + 
  labs(y = 'Standardized slope', color = '°C/year', title = 'Morisita horn')

grid.arrange(p1, p2, p3, ncol = 3)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/first-last%20dissimilarity%20vs.%20temperature%20trend%20vs.%20duration-1.png)<!-- -->

# Models

## Choose the variance structure for mixed effects models

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
fixed <- formula(Jtutrendrem0 ~ temptrend_abs.sc*REALM +
                     temptrend_abs.sc*tsign + 
                     temptrend_abs.sc*tempave_metab.sc + 
                     temptrend_abs.sc*seas.sc + 
                     temptrend_abs.sc*microclim.sc + 
                     temptrend_abs.sc*mass.sc + 
                     temptrend_abs.sc*speed.sc + 
                     temptrend_abs.sc*consumerfrac.sc +
                     temptrend_abs.sc*nspp.sc +
                     temptrend_abs.sc*thermal_bias.sc:tsign +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*veg.sc +
                     temptrend_abs.sc*duration.sc +
                     temptrend_abs.sc*human_bowler.sc:REALM2)
i <- trends[, complete.cases(Jtutrendrem0, temptrend_abs.sc, REALM, tsign, tempave_metab.sc, seas.sc, 
                             microclim.sc, mass.sc, speed.sc, consumerfrac.sc, nspp.sc,
                             thermal_bias.sc, npp.sc, veg.sc, human_bowler.sc)]
mods <- vector('list', 0)
mods[[1]] <- gls(fixed, data = trends[i,])
mods[[2]] <- gls(fixed, data = trends[i,], weights = varPower(-0.5, ~nyrBT))
mods[[3]] <- gls(fixed, data = trends[i,], weights = varPower(0.5, ~temptrend_abs.sc))

mods[[4]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2, control = lmeControl(opt = "optim"))
mods[[5]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID, control = lmeControl(opt = "optim"))
mods[[6]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID, control = lmeControl(opt = "optim"))
mods[[7]] <- lme(fixed, data = trends[i,], random = ~1|STUDY_ID/rarefyID, control = lmeControl(opt = "optim"))
mods[[8]] <- lme(fixed, data = trends[i,], random = ~1|taxa_mod2/STUDY_ID/rarefyID, control = lmeControl(opt = "optim"))

mods[[9]] <- lme(fixed, data = trends[i,], random = ~temptrend_abs.sc | taxa_mod)
mods[[10]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)) # includes overdispersion. new formula so that random slope is only for study level (not enough data to extend to rarefyID).

mods[[11]] <- lme(fixed, data = trends[i,], random = list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT))
mods[[12]] <- lme(fixed, data = trends[i,], random = list(taxa_mod2 = ~ temptrend_abs.sc, STUDY_ID = ~ 1, rarefyID = ~1), weights = varPower(-0.5, ~nyrBT))

aics <- sapply(mods, AIC)
minaics <- aics - min(aics)
minaics
which.min(aics)
```

Chooses the random slopes (temptrend\_abs) & intercepts for STUDY\_ID,
overdispersion, and variance scaled to number of years. We haven’t dealt
with potential testing on the boundary issues here yet.

## Temperature-only models

### Fit the models

``` r
randef <- list(STUDY_ID = ~ abs(temptrend), rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# *trendrem0
if(file.exists('temp/modonlyTtrendrem0.rds')){
  modonlyTtrendrem0 <- readRDS('temp/modonlyTtrendrem0.rds')
} else {
  i <- trends[, complete.cases(Jtutrendrem0, REALM, temptrend)]
  modonlyTtrendrem0 <- lme(Jtutrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrendrem0, file = 'temp/modonlyTtrendrem0.rds')
}

if(file.exists('temp/modonlyTtrendJbetarem0.rds')){
  modonlyTtrendJbetarem0 <- readRDS('temp/modonlyTtrendJbetarem0.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendrem0, REALM, temptrend)]
  modonlyTtrendJbetarem0 <- lme(Jbetatrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetarem0, file = 'temp/modonlyTtrendJbetarem0.rds')
}

if(file.exists('temp/modonlyTtrendHornrem0.rds')){
  modonlyTtrendHornrem0 <- readRDS('temp/modonlyTtrendHornrem0.rds')
} else {
  i <- trends[, complete.cases(Horntrendrem0, REALM, temptrend)]
  modonlyTtrendHornrem0 <- lme(Horntrendrem0 ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrendHornrem0, file = 'temp/modonlyTtrendHornrem0.rds')
}

# *trendz
if(file.exists('temp/modonlyTtrendJtuz.rds')){
  modonlyTtrendJtuz <- readRDS('temp/modonlyTtrendJtuz.rds')
} else {
  i <- trends[, complete.cases(Jtutrendz, REALM, temptrend)]
  modonlyTtrendJtuz <- lme(Jtutrendz ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML',
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJtuz, file = 'temp/modonlyTtrendJtuz.rds')
}

if(file.exists('temp/modonlyTtrendJbetaz.rds')){
  modonlyTtrendJbetaz <- readRDS('temp/modonlyTtrendJbetaz.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendz, REALM, temptrend)]
  modonlyTtrendJbetaz <- lme(Jbetatrendz ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetaz, file = 'temp/modonlyTtrendJbetaz.rds')
}

if(file.exists('temp/modonlyTtrendHornz.rds')){
  modonlyTtrendHornz <- readRDS('temp/modonlyTtrendHornz.rds')
} else {
  i <- trends[, complete.cases(Horntrendz, REALM, temptrend)]
  modonlyTtrendHornz <- lme(Horntrendz ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrendHornz, file = 'temp/modonlyTtrendHornz.rds')
}

# *exp
if(file.exists('temp/modonlyTtrendJtuexp.rds')){
  modonlyTtrendJtuexp <- readRDS('temp/modonlyTtrendJtuexp.rds')
} else {
  i <- trends[, complete.cases(Jtuexp, REALM, temptrend)]
  modonlyTtrendJtuexp <- lme(log(Jtuexp) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML',
                   control=lmeControl(msMaxIter = 100, maxIter = 100, opt = 'optim'))
  saveRDS(modonlyTtrendJtuexp, file = 'temp/modonlyTtrendJtuexp.rds')
}

if(file.exists('temp/modonlyTtrendJbetaexp.rds')){
  modonlyTtrendJbetaexp <- readRDS('temp/modonlyTtrendJbetaexp.rds')
} else {
  i <- trends[, complete.cases(Jbetaexp, REALM, temptrend)]
  modonlyTtrendJbetaexp <- lme(log(Jbetaexp) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetaexp, file = 'temp/modonlyTtrendJbetaexp.rds')
}

if(file.exists('temp/modonlyTtrendHornexp.rds')){
  modonlyTtrendHornexp <- readRDS('temp/modonlyTtrendHornexp.rds')
} else {
  i <- trends[, complete.cases(Hornexp, REALM, temptrend)]
  modonlyTtrendHornexp <- lme(log(Hornexp) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML',
                   control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modonlyTtrendHornexp, file = 'temp/modonlyTtrendHornexp.rds')
}


# *mm
if(file.exists('temp/modonlyTtrendJtumm.rds')){
  modonlyTtrendJtumm <- readRDS('temp/modonlyTtrendJtumm.rds')
} else {
  i <- trends[, complete.cases(Jtumm, REALM, temptrend)]
  modonlyTtrendJtumm <- lme(log(Jtumm+1) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modonlyTtrendJtumm, file = 'temp/modonlyTtrendJtumm.rds')
}

if(file.exists('temp/modonlyTtrendJbetamm.rds')){
  modonlyTtrendJbetamm <- readRDS('temp/modonlyTtrendJbetamm.rds')
} else {
  i <- trends[, complete.cases(Jbetamm, REALM, temptrend)]
  modonlyTtrendJbetamm <- lme(log(Jbetamm+1) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendJbetamm, file = 'temp/modonlyTtrendJbetamm.rds')
}

if(file.exists('temp/modonlyTtrendHornmm.rds')){
  modonlyTtrendHornmm <- readRDS('temp/modonlyTtrendHornmm.rds')
} else {
  i <- trends[, complete.cases(Hornmm, REALM, temptrend)]
  modonlyTtrendHornmm <- lme(log(Hornmm+1) ~ abs(temptrend)*REALM,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modonlyTtrendHornmm, file = 'temp/modonlyTtrendHornmm.rds')
}
```

### Summary: Slope of dissimilarity

``` r
summary(modonlyTtrendrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -94827.26 -94727.05 47425.63
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev      Corr  
    ## (Intercept)    0.007952925 (Intr)
    ## abs(temptrend) 0.164763893 -0.959
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01157533 2.163894
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.246691 
    ## Fixed effects: Jtutrendrem0 ~ abs(temptrend) * REALM 
    ##                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                      0.00060984 0.00369993 31054  0.1648237  0.8691
    ## abs(temptrend)                   0.05886907 0.08331500 31054  0.7065842  0.4798
    ## REALMMarine                      0.00543120 0.00394538   247  1.3765972  0.1699
    ## REALMTerrestrial                 0.00278420 0.00411181   247  0.6771220  0.4990
    ## abs(temptrend):REALMMarine      -0.05091734 0.08813249 31054 -0.5777362  0.5634
    ## abs(temptrend):REALMTerrestrial -0.01905325 0.09142020 31054 -0.2084140  0.8349
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.788                                
    ## REALMMarine                     -0.938  0.739                         
    ## REALMTerrestrial                -0.900  0.709  0.844                  
    ## abs(temptrend):REALMMarine       0.745 -0.945 -0.802 -0.670           
    ## abs(temptrend):REALMTerrestrial  0.718 -0.911 -0.674 -0.801  0.862    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.78048404 -0.23246011 -0.01107399  0.27719738  6.55447835 
    ## 
    ## Number of Observations: 31307
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31307

``` r
summary(modonlyTtrendJbetarem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -122330.3 -122230.1 61177.15
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev      Corr  
    ## (Intercept)    0.006543123 (Intr)
    ## abs(temptrend) 0.071122544 -0.049
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev: 0.006325163 1.133689
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.094327 
    ## Fixed effects: Jbetatrendrem0 ~ abs(temptrend) * REALM 
    ##                                       Value  Std.Error    DF   t-value p-value
    ## (Intercept)                      0.00282279 0.00275536 31054  1.024471  0.3056
    ## abs(temptrend)                   0.15794305 0.05188059 31054  3.044357  0.0023
    ## REALMMarine                      0.00336986 0.00293945   247  1.146426  0.2527
    ## REALMTerrestrial                 0.00427163 0.00302954   247  1.409993  0.1598
    ## abs(temptrend):REALMMarine      -0.13690052 0.05408802 31054 -2.531069  0.0114
    ## abs(temptrend):REALMTerrestrial -0.13349019 0.05606887 31054 -2.380826  0.0173
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.567                                
    ## REALMMarine                     -0.937  0.532                         
    ## REALMTerrestrial                -0.909  0.516  0.853                  
    ## abs(temptrend):REALMMarine       0.544 -0.959 -0.546 -0.495           
    ## abs(temptrend):REALMTerrestrial  0.525 -0.925 -0.492 -0.553  0.888    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.79146172 -0.31481245 -0.01043308  0.33402768  7.91183985 
    ## 
    ## Number of Observations: 31307
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31307

``` r
summary(modonlyTtrendHornrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC       BIC   logLik
    ##   -88924.7 -88824.48 44474.35
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev     Corr  
    ## (Intercept)    0.01323437 (Intr)
    ## abs(temptrend) 0.23503241 0.139 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01934775 2.662513
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.419505 
    ## Fixed effects: Horntrendrem0 ~ abs(temptrend) * REALM 
    ##                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                      0.00358655 0.00587222 31054  0.6107649  0.5414
    ## abs(temptrend)                   0.21157330 0.11625477 31054  1.8199107  0.0688
    ## REALMMarine                      0.00591612 0.00622487   247  0.9504006  0.3428
    ## REALMTerrestrial                 0.00557399 0.00650476   247  0.8569095  0.3923
    ## abs(temptrend):REALMMarine      -0.08375236 0.12268083 31054 -0.6826850  0.4948
    ## abs(temptrend):REALMTerrestrial -0.13582251 0.12790730 31054 -1.0618824  0.2883
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.450                                
    ## REALMMarine                     -0.943  0.424                         
    ## REALMTerrestrial                -0.903  0.406  0.852                  
    ## abs(temptrend):REALMMarine       0.426 -0.948 -0.428 -0.385           
    ## abs(temptrend):REALMTerrestrial  0.409 -0.909 -0.386 -0.448  0.861    
    ## 
    ## Standardized Within-Group Residuals:
    ##          Min           Q1          Med           Q3          Max 
    ## -5.417759130 -0.228422814 -0.009748498  0.240197159  6.061908050 
    ## 
    ## Number of Observations: 31307
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31307

### Summary: Standardized slope of dissimilarity

``` r
summary(modonlyTtrendJtuz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC   logLik
    ##   128091.8 128190.6 -64033.9
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev   Corr  
    ## (Intercept)    4.443672 (Intr)
    ## abs(temptrend) 6.329423 0.001 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    2.181317  9311.15
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -7.642962 
    ## Fixed effects: Jtutrendz ~ abs(temptrend) * REALM 
    ##                                     Value Std.Error    DF    t-value p-value
    ## (Intercept)                     -0.005489  1.660118 27691 -0.0033063  0.9974
    ## abs(temptrend)                   3.164749  4.031553 27691  0.7849951  0.4325
    ## REALMMarine                      1.053327  1.742023   171  0.6046575  0.5462
    ## REALMTerrestrial                 1.786041  1.745811   171  1.0230438  0.3077
    ## abs(temptrend):REALMMarine      -5.098876  4.221906 27691 -1.2077189  0.2272
    ## abs(temptrend):REALMTerrestrial -5.025365  4.535863 27691 -1.1079182  0.2679
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.119                                
    ## REALMMarine                     -0.953  0.113                         
    ## REALMTerrestrial                -0.951  0.113  0.906                  
    ## abs(temptrend):REALMMarine       0.113 -0.955 -0.116 -0.108           
    ## abs(temptrend):REALMTerrestrial  0.105 -0.889 -0.100 -0.136  0.849    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.334605e+01 -1.139447e-02 -4.082560e-05  6.281415e-04  1.393971e+01 
    ## 
    ## Number of Observations: 27868
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    174                  27868

``` r
summary(modonlyTtrendJbetaz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##       AIC      BIC    logLik
    ##   95583.3 95682.16 -47779.65
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)     3.183263 (Intr)
    ## abs(temptrend) 18.140923 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)    Residual
    ## StdDev:    1.198246 39081402602
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -22.02594 
    ## Fixed effects: Jbetatrendz ~ abs(temptrend) * REALM 
    ##                                     Value Std.Error    DF    t-value p-value
    ## (Intercept)                      1.048980  1.168466 27780  0.8977406  0.3693
    ## abs(temptrend)                   1.328315  8.422812 27780  0.1577044  0.8747
    ## REALMMarine                     -0.256308  1.229235   176 -0.2085098  0.8351
    ## REALMTerrestrial                 0.852803  1.236422   176  0.6897347  0.4913
    ## abs(temptrend):REALMMarine      -6.310774  8.826995 27780 -0.7149403  0.4747
    ## abs(temptrend):REALMTerrestrial -6.240139  9.276016 27780 -0.6727176  0.5011
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.205                                
    ## REALMMarine                     -0.951  0.195                         
    ## REALMTerrestrial                -0.945  0.194  0.898                  
    ## abs(temptrend):REALMMarine       0.196 -0.954 -0.197 -0.185           
    ## abs(temptrend):REALMTerrestrial  0.186 -0.908 -0.177 -0.229  0.866    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.259502e+09 -7.302550e-04 -1.106162e-07  2.961634e-06  2.843046e+09 
    ## 
    ## Number of Observations: 27962
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    179                  27962

``` r
summary(modonlyTtrendHornz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##      AIC      BIC logLik
    ##   148842 148940.8 -74409
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)     3.958082 (Intr)
    ## abs(temptrend) 37.358258 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)     Residual
    ## StdDev:    3.397087 106031655443
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -22.91116 
    ## Fixed effects: Horntrendz ~ abs(temptrend) * REALM 
    ##                                      Value Std.Error    DF    t-value p-value
    ## (Intercept)                       1.160245  1.682387 27691  0.6896423  0.4904
    ## abs(temptrend)                   14.386296 16.949480 27691  0.8487751  0.3960
    ## REALMMarine                       0.841914  1.758304   167  0.4788219  0.6327
    ## REALMTerrestrial                  1.125087  1.798922   167  0.6254231  0.5325
    ## abs(temptrend):REALMMarine      -21.204956 17.822950 27691 -1.1897556  0.2342
    ## abs(temptrend):REALMTerrestrial  -9.501499 18.649869 27691 -0.5094673  0.6104
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.230                                
    ## REALMMarine                     -0.957  0.220                         
    ## REALMTerrestrial                -0.935  0.215  0.895                  
    ## abs(temptrend):REALMMarine       0.219 -0.951 -0.225 -0.204           
    ## abs(temptrend):REALMTerrestrial  0.209 -0.909 -0.200 -0.259  0.864    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -5.338773e+09 -2.267566e-04 -1.621989e-08  3.280843e-06  2.186012e+07 
    ## 
    ## Number of Observations: 27864
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    170                  27864

### Summary: Half-saturation exponential fits

``` r
summary(modonlyTtrendJtuexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   82973.57 83069.33 -41474.78
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)    0.9875993 (Intr)
    ## abs(temptrend) 2.6544585 -0.782
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.570786 18.51434
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.780927 
    ## Fixed effects: log(Jtuexp) ~ abs(temptrend) * REALM 
    ##                                     Value Std.Error    DF    t-value p-value
    ## (Intercept)                      0.599799 0.3783257 21379  1.5854038  0.1129
    ## abs(temptrend)                  -4.663902 1.8820317 21379 -2.4781210  0.0132
    ## REALMMarine                     -0.714627 0.3996836   217 -1.7879816  0.0752
    ## REALMTerrestrial                 0.012499 0.4128602   217  0.0302745  0.9759
    ## abs(temptrend):REALMMarine       1.884189 1.9535542 21379  0.9644926  0.3348
    ## abs(temptrend):REALMTerrestrial  1.495374 2.0372538 21379  0.7340144  0.4629
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.586                                
    ## REALMMarine                     -0.947  0.555                         
    ## REALMTerrestrial                -0.916  0.537  0.867                  
    ## abs(temptrend):REALMMarine       0.565 -0.963 -0.588 -0.517           
    ## abs(temptrend):REALMTerrestrial  0.541 -0.924 -0.513 -0.590  0.890    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.5025451540 -0.0394533801 -0.0003422187  0.0655541535  1.1594445321 
    ## 
    ## Number of Observations: 21602
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    220                  21602

``` r
summary(modonlyTtrendJbetaexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   73942.52 74038.01 -36959.26
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)    0.9192722 (Intr)
    ## abs(temptrend) 3.2263373 -0.85 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.214096 4.588968
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.292111 
    ## Fixed effects: log(Jbetaexp) ~ abs(temptrend) * REALM 
    ##                                     Value Std.Error    DF   t-value p-value
    ## (Intercept)                      0.490416 0.3241854 20885  1.512765  0.1304
    ## abs(temptrend)                  -5.783240 1.8660448 20885 -3.099197  0.0019
    ## REALMMarine                     -1.134735 0.3447447   225 -3.291524  0.0012
    ## REALMTerrestrial                -0.291380 0.3563609   225 -0.817654  0.4144
    ## abs(temptrend):REALMMarine       3.519076 1.9494360 20885  1.805176  0.0711
    ## abs(temptrend):REALMTerrestrial  2.122820 2.0331743 20885  1.044092  0.2965
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.658                                
    ## REALMMarine                     -0.940  0.619                         
    ## REALMTerrestrial                -0.910  0.599  0.855                  
    ## abs(temptrend):REALMMarine       0.630 -0.957 -0.662 -0.573           
    ## abs(temptrend):REALMTerrestrial  0.604 -0.918 -0.568 -0.662  0.879    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -1.86244678 -0.16187559  0.01884547  0.24753908  1.96519825 
    ## 
    ## Number of Observations: 21116
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    228                  21116

``` r
summary(modonlyTtrendHornexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   78503.15 78598.01 -39239.57
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev   Corr  
    ## (Intercept)    1.138044 (Intr)
    ## abs(temptrend) 1.640645 -0.911
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:     1.62421 40.53926
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##    power 
    ## -3.36273 
    ## Fixed effects: log(Hornexp) ~ abs(temptrend) * REALM 
    ##                                     Value Std.Error    DF   t-value p-value
    ## (Intercept)                      0.691658 0.4147774 19818  1.667540  0.0954
    ## abs(temptrend)                  -6.261463 1.6719731 19818 -3.744954  0.0002
    ## REALMMarine                     -0.557617 0.4392395   214 -1.269506  0.2056
    ## REALMTerrestrial                 0.385835 0.4501513   214  0.857124  0.3923
    ## abs(temptrend):REALMMarine       3.835362 1.7112816 19818  2.241222  0.0250
    ## abs(temptrend):REALMTerrestrial  3.549273 1.7301246 19818  2.051455  0.0402
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.506                                
    ## REALMMarine                     -0.944  0.478                         
    ## REALMTerrestrial                -0.921  0.466  0.870                  
    ## abs(temptrend):REALMMarine       0.495 -0.977 -0.510 -0.456           
    ## abs(temptrend):REALMTerrestrial  0.489 -0.966 -0.462 -0.514  0.944    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.4311778659 -0.0328922575  0.0002377582  0.0612146060  1.6562942910 
    ## 
    ## Number of Observations: 20038
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    217                  20038

### Summary: Half-saturation MM fits

``` r
summary(modonlyTtrendJtumm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   94365.67 94465.85 -47170.84
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)    0.7999516 (Intr)
    ## abs(temptrend) 2.0246507 -0.949
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:    1.090041 0.4024477
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -6.845854 
    ## Fixed effects: log(Jtumm + 1) ~ abs(temptrend) * REALM 
    ##                                      Value Std.Error    DF   t-value p-value
    ## (Intercept)                      1.1142870 0.2716646 30961  4.101702  0.0000
    ## abs(temptrend)                  -1.5947708 1.1066878 30961 -1.441030  0.1496
    ## REALMMarine                     -0.3746800 0.2886013   247 -1.298262  0.1954
    ## REALMTerrestrial                 0.3919001 0.2936656   247  1.334512  0.1833
    ## abs(temptrend):REALMMarine       0.8822214 1.1496742 30961  0.767366  0.4429
    ## abs(temptrend):REALMTerrestrial -0.7948202 1.1781360 30961 -0.674642  0.4999
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.692                                
    ## REALMMarine                     -0.941  0.651                         
    ## REALMTerrestrial                -0.925  0.640  0.871                  
    ## abs(temptrend):REALMMarine       0.666 -0.963 -0.700 -0.616           
    ## abs(temptrend):REALMTerrestrial  0.650 -0.939 -0.612 -0.699  0.904    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -6.602942e-02 -5.940902e-06 -3.900428e-08  1.227987e-06  6.131904e-02 
    ## 
    ## Number of Observations: 31214
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31214

``` r
summary(modonlyTtrendJbetamm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   67296.72 67396.91 -33636.36
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)    0.8967893 (Intr)
    ## abs(temptrend) 1.7739631 -0.975
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:  0.03150223 0.8792083
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.1283621 
    ## Fixed effects: log(Jbetamm + 1) ~ abs(temptrend) * REALM 
    ##                                      Value Std.Error    DF   t-value p-value
    ## (Intercept)                      1.2211102 0.2414200 30981  5.058032  0.0000
    ## abs(temptrend)                  -2.0277165 0.7755426 30981 -2.614578  0.0089
    ## REALMMarine                     -0.8518779 0.2615314   247 -3.257269  0.0013
    ## REALMTerrestrial                -0.0774185 0.2600212   247 -0.297739  0.7662
    ## abs(temptrend):REALMMarine       1.7380613 0.8107604 30981  2.143742  0.0321
    ## abs(temptrend):REALMTerrestrial  0.4224251 0.8210222 30981  0.514511  0.6069
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.697                                
    ## REALMMarine                     -0.923  0.643                         
    ## REALMTerrestrial                -0.928  0.647  0.857                  
    ## abs(temptrend):REALMMarine       0.666 -0.957 -0.711 -0.619           
    ## abs(temptrend):REALMTerrestrial  0.658 -0.945 -0.607 -0.704  0.904    
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -6.2496027 -0.6390446 -0.3197549  0.3643529  6.2223284 
    ## 
    ## Number of Observations: 31234
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31234

``` r
summary(modonlyTtrendHornmm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   99714.74 99814.91 -49845.37
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev   Corr  
    ## (Intercept)    1.177000 (Intr)
    ## abs(temptrend) 1.438911 -0.953
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.126723 2.869408
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##    power 
    ## -1.39652 
    ## Fixed effects: log(Hornmm + 1) ~ abs(temptrend) * REALM 
    ##                                      Value Std.Error    DF   t-value p-value
    ## (Intercept)                      1.4431310 0.3429821 30940  4.207599  0.0000
    ## abs(temptrend)                  -2.7211693 1.0512925 30940 -2.588403  0.0096
    ## REALMMarine                     -0.6687804 0.3686631   247 -1.814069  0.0709
    ## REALMTerrestrial                 0.4909528 0.3694109   247  1.329016  0.1851
    ## abs(temptrend):REALMMarine       2.4054805 1.0798956 30940  2.227512  0.0259
    ## abs(temptrend):REALMTerrestrial  1.0415061 1.0903816 30940  0.955176  0.3395
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT a():REALMM
    ## abs(temptrend)                  -0.514                                
    ## REALMMarine                     -0.930  0.478                         
    ## REALMTerrestrial                -0.928  0.477  0.864                  
    ## abs(temptrend):REALMMarine       0.500 -0.974 -0.524 -0.465           
    ## abs(temptrend):REALMTerrestrial  0.496 -0.964 -0.461 -0.523  0.939    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -0.99832886 -0.16228129 -0.05493216  0.08821033  2.40925799 
    ## 
    ## Number of Observations: 31193
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31193

### Plot the temp-only coefficients

``` r
# make table of coefficients
coefs1 <- as.data.frame(summary(modonlyTtrendrem0)$tTable)
coefs2 <- as.data.frame(summary(modonlyTtrendJbetarem0)$tTable)
coefs3 <- as.data.frame(summary(modonlyTtrendHornrem0)$tTable)
coefs4 <- as.data.frame(summary(modonlyTtrendJtuz)$tTable)
coefs5 <- as.data.frame(summary(modonlyTtrendJbetaz)$tTable)
coefs6 <- as.data.frame(summary(modonlyTtrendHornz)$tTable)
coefs1$response <- 'Jtu'
coefs2$response <- 'Jbeta'
coefs3$response <- 'Horn'
coefs4$response <- 'Jtu'
coefs5$response <- 'Jbeta'
coefs6$response <- 'Horn'
coefs1$fit <- 'rem0'
coefs2$fit <- 'rem0'
coefs3$fit <- 'rem0'
coefs4$fit <- 'z'
coefs5$fit <- 'z'
coefs6$fit <- 'z'
rows1 <- which(grepl('temptrend', rownames(coefs1))) # extract temperature effect
cols <- c('Value', 'Std.Error', 'response', 'fit')
allcoefs <- rbind(coefs1[rows1, cols], coefs2[rows1, cols], coefs3[rows1, cols],
                  coefs4[rows1, cols], coefs5[rows1, cols], coefs6[rows1, cols])
allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to marine effects
allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to terrestrial effects

allcoefs$lCI <- allcoefs$Value - allcoefs$Std.Error # lower confidence interval
allcoefs$uCI <- allcoefs$Value + allcoefs$Std.Error
allcoefs$y <- c(3, 2, 1) + rep(c(0, -0.1, -0.2), c(3, 3, 3)) # y-values
allcoefs$realm <- rep(c('Freshwater', 'Marine', 'Terrestrial'), 6)
allcoefs$group <- paste0(allcoefs$response, allcoefs$fit)

pd <- position_dodge(width = 0.5)
ggplot(subset(allcoefs, fit == 'rem0'), aes(x = realm, y = Value, color = response, shape = fit, group = group)) +
  geom_point(position = pd) + 
  geom_errorbar(aes(ymin=lCI, ymax=uCI), width=.1, position=pd) +
  coord_flip() +
  labs(y = 'Turnover per |°C/yr|')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/modonlyT%20coefs-1.png)<!-- -->

``` r
ggplot(subset(allcoefs, fit == 'z'), aes(x = realm, y = Value, color = response, shape = fit, group = group)) +
  geom_point(position = pd) + 
  geom_errorbar(aes(ymin=lCI, ymax=uCI), width=.1, position=pd) +
  coord_flip() +
  labs(y = 'Standardized Turnover per |°C/yr|')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/modonlyT%20coefs-2.png)<!-- -->

## Temperature\&duration models

### Fit the models

``` r
randef <- list(STUDY_ID = ~ abs(temptrend), rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# *trendrem0
if(file.exists('temp/modTDtrendJturem0.rds')){
  modTDtrendJturem0 <- readRDS('temp/modTDtrendJturem0.rds')
} else {
  i <- trends[, complete.cases(Jtutrendrem0, REALM, temptrend)]
  modTDtrendJturem0 <- lme(Jtutrendrem0 ~ abs(temptrend)*REALM + duration,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modTDtrendJturem0, file = 'temp/modTDtrendJturem0.rds')
}

if(file.exists('temp/modTDtrendJbetarem0.rds')){
  modTDtrendJbetarem0 <- readRDS('temp/modTDtrendJbetarem0.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendrem0, REALM, temptrend)]
  modTDtrendJbetarem0 <- lme(Jbetatrendrem0 ~ abs(temptrend)*REALM + duration,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modTDtrendJbetarem0, file = 'temp/modTDtrendJbetarem0.rds')
}

if(file.exists('temp/modTDtrendHornrem0.rds')){
  modTDtrendHornrem0 <- readRDS('temp/modTDtrendHornrem0.rds')
} else {
  i <- trends[, complete.cases(Horntrendrem0, REALM, temptrend)]
  modTDtrendHornrem0 <- lme(Horntrendrem0 ~ abs(temptrend)*REALM + duration,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modTDtrendHornrem0, file = 'temp/modTDtrendHornrem0.rds')
}

# *trendz
if(file.exists('temp/modTDtrendJtuz.rds')){
  modTDtrendJtuz <- readRDS('temp/modTDtrendJtuz.rds')
} else {
  i <- trends[, complete.cases(Jtutrendz, REALM, temptrend)]
  modTDtrendJtuz <- lme(Jtutrendz ~ abs(temptrend)*REALM + duration,
                   random = randef, weights = varef, data = trends[i,], method = 'REML',
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modTDtrendJtuz, file = 'temp/modTDtrendJtuz.rds')
}

if(file.exists('temp/modTDtrendJbetaz.rds')){
  modTDtrendJbetaz <- readRDS('temp/modTDtrendJbetaz.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendz, REALM, temptrend)]
  modTDtrendJbetaz <- lme(Jbetatrendz ~ abs(temptrend)*REALM + duration,
                   random = randef, weights = varef, data = trends[i,], method = 'REML', 
                   control=lmeControl(msMaxIter = 100, maxIter = 100))
  saveRDS(modTDtrendJbetaz, file = 'temp/modTDtrendJbetaz.rds')
}

if(file.exists('temp/modTDtrendHornz.rds')){
  modTDtrendHornz <- readRDS('temp/modTDtrendHornz.rds')
} else {
  i <- trends[, complete.cases(Horntrendz, REALM, temptrend)]
  modTDtrendHornz <- lme(Horntrendz ~ abs(temptrend)*REALM + duration,
                   random = randef, weights = varef, data = trends[i,], method = 'REML')
  saveRDS(modTDtrendHornz, file = 'temp/modTDtrendHornz.rds')
}
```

### Summary: Slope of dissimilarity

``` r
summary(modTDtrendJturem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -94806.48 -94697.91 47416.24
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev     Corr  
    ## (Intercept)    0.00796915 (Intr)
    ## abs(temptrend) 0.16577501 -0.96 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01157376  2.16344
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.246535 
    ## Fixed effects: Jtutrendrem0 ~ abs(temptrend) * REALM + duration 
    ##                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                      0.00107885 0.00374647 31053  0.2879649  0.7734
    ## abs(temptrend)                   0.05698260 0.08353780 31053  0.6821176  0.4952
    ## REALMMarine                      0.00538249 0.00394827   247  1.3632521  0.1740
    ## REALMTerrestrial                 0.00281128 0.00411437   247  0.6832836  0.4951
    ## duration                        -0.00001951 0.00002374 31053 -0.8218123  0.4112
    ## abs(temptrend):REALMMarine      -0.05064091 0.08835293 31053 -0.5731662  0.5665
    ## abs(temptrend):REALMTerrestrial -0.02005547 0.09166191 31053 -0.2187983  0.8268
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT duratn a():REALMM
    ## abs(temptrend)                  -0.783                                       
    ## REALMMarine                     -0.929  0.740                                
    ## REALMTerrestrial                -0.887  0.709  0.843                         
    ## duration                        -0.154  0.027  0.016 -0.011                  
    ## abs(temptrend):REALMMarine       0.737 -0.945 -0.802 -0.671 -0.004           
    ## abs(temptrend):REALMTerrestrial  0.708 -0.910 -0.674 -0.802  0.016  0.861    
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -6.7848569 -0.2329439 -0.0116805  0.2765163  6.5511232 
    ## 
    ## Number of Observations: 31307
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31307

``` r
summary(modTDtrendJbetarem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##       AIC       BIC   logLik
    ##   -122317 -122208.5 61171.52
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev      Corr  
    ## (Intercept)    0.006433222 (Intr)
    ## abs(temptrend) 0.072854449 -0.048
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev: 0.006331461 1.134306
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.094866 
    ## Fixed effects: Jbetatrendrem0 ~ abs(temptrend) * REALM + duration 
    ##                                       Value  Std.Error    DF   t-value p-value
    ## (Intercept)                      0.00393307 0.00276230 31053  1.423839  0.1545
    ## abs(temptrend)                   0.15516447 0.05211986 31053  2.977070  0.0029
    ## REALMMarine                      0.00322780 0.00291733   247  1.106424  0.2696
    ## REALMTerrestrial                 0.00436785 0.00300781   247  1.452169  0.1477
    ## duration                        -0.00004732 0.00001573 31053 -3.008979  0.0026
    ## abs(temptrend):REALMMarine      -0.13829445 0.05435922 31053 -2.544085  0.0110
    ## abs(temptrend):REALMTerrestrial -0.13528328 0.05637326 31053 -2.399777  0.0164
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT duratn a():REALMM
    ## abs(temptrend)                  -0.565                                       
    ## REALMMarine                     -0.931  0.533                                
    ## REALMTerrestrial                -0.899  0.516  0.852                         
    ## duration                        -0.141  0.025  0.023 -0.007                  
    ## abs(temptrend):REALMMarine       0.538 -0.958 -0.547 -0.495  0.002           
    ## abs(temptrend):REALMTerrestrial  0.519 -0.924 -0.492 -0.554  0.006  0.886    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.79417832 -0.31740519 -0.01193816  0.33042562  7.90998791 
    ## 
    ## Number of Observations: 31307
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31307

``` r
summary(modTDtrendHornrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC      BIC   logLik
    ##   -88908.66 -88800.1 44467.33
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev     Corr  
    ## (Intercept)    0.01296688 (Intr)
    ## abs(temptrend) 0.23571800 0.136 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01934892 2.661858
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.419368 
    ## Fixed effects: Horntrendrem0 ~ abs(temptrend) * REALM + duration 
    ##                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                      0.00523692 0.00586639 31053  0.8926993  0.3720
    ## abs(temptrend)                   0.20687769 0.11621673 31053  1.7801025  0.0751
    ## REALMMarine                      0.00562360 0.00616271   247  0.9125214  0.3624
    ## REALMTerrestrial                 0.00563565 0.00644175   247  0.8748627  0.3825
    ## duration                        -0.00007331 0.00003363 31053 -2.1800405  0.0293
    ## abs(temptrend):REALMMarine      -0.08478634 0.12263560 31053 -0.6913680  0.4893
    ## abs(temptrend):REALMTerrestrial -0.13737366 0.12784924 31053 -1.0744973  0.2826
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT duratn a():REALMM
    ## abs(temptrend)                  -0.450                                       
    ## REALMMarine                     -0.938  0.426                                
    ## REALMTerrestrial                -0.894  0.407  0.852                         
    ## duration                        -0.132  0.020  0.022 -0.004                  
    ## abs(temptrend):REALMMarine       0.423 -0.947 -0.429 -0.386  0.003           
    ## abs(temptrend):REALMTerrestrial  0.405 -0.909 -0.387 -0.449  0.006  0.861    
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.42088463 -0.22950254 -0.01052715  0.23841443  6.05788241 
    ## 
    ## Number of Observations: 31307
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    250                  31307

### Summary: Standardized slope of dissimilarity

``` r
summary(modTDtrendJtuz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   128046.9 128153.9 -64010.44
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev   Corr  
    ## (Intercept)    4.448517 (Intr)
    ## abs(temptrend) 6.205057 0.002 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:     2.17854 8489.824
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -7.558455 
    ## Fixed effects: Jtutrendz ~ abs(temptrend) * REALM + duration 
    ##                                     Value Std.Error    DF   t-value p-value
    ## (Intercept)                     -0.390450  1.662013 27690 -0.234926  0.8143
    ## abs(temptrend)                   4.138669  3.982931 27690  1.039101  0.2988
    ## REALMMarine                      1.162999  1.743263   171  0.667139  0.5056
    ## REALMTerrestrial                 1.840754  1.746842   171  1.053761  0.2935
    ## duration                         0.018692  0.002472 27690  7.561076  0.0000
    ## abs(temptrend):REALMMarine      -5.067722  4.168020 27690 -1.215858  0.2240
    ## abs(temptrend):REALMTerrestrial -5.482772  4.477897 27690 -1.224408  0.2208
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT duratn a():REALMM
    ## abs(temptrend)                  -0.118                                       
    ## REALMMarine                     -0.953  0.112                                
    ## REALMTerrestrial                -0.951  0.112  0.906                         
    ## duration                        -0.031  0.033  0.008  0.004                  
    ## abs(temptrend):REALMMarine       0.112 -0.955 -0.115 -0.107  0.000           
    ## abs(temptrend):REALMTerrestrial  0.105 -0.889 -0.099 -0.134 -0.015  0.849    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.336934e+01 -8.255770e-03 -1.944539e-05  1.315311e-03  1.395269e+01 
    ## 
    ## Number of Observations: 27868
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    174                  27868

``` r
summary(modTDtrendJbetaz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   95510.45 95617.55 -47742.23
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev    Corr  
    ## (Intercept)     3.176372 (Intr)
    ## abs(temptrend) 17.200966 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)     Residual
    ## StdDev:    1.196077 273824383770
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -23.79672 
    ## Fixed effects: Jbetatrendz ~ abs(temptrend) * REALM + duration 
    ##                                     Value Std.Error    DF   t-value p-value
    ## (Intercept)                      0.766367  1.162158 27779  0.659434  0.5096
    ## abs(temptrend)                   2.012849  8.046389 27779  0.250156  0.8025
    ## REALMMarine                     -0.155590  1.222513   176 -0.127271  0.8989
    ## REALMTerrestrial                 0.866941  1.229099   176  0.705346  0.4815
    ## duration                         0.012781  0.001374 27779  9.303166  0.0000
    ## abs(temptrend):REALMMarine      -6.257063  8.431757 27779 -0.742083  0.4580
    ## abs(temptrend):REALMTerrestrial -6.075665  8.869183 27779 -0.685031  0.4933
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT duratn a():REALMM
    ## abs(temptrend)                  -0.202                                       
    ## REALMMarine                     -0.950  0.192                                
    ## REALMTerrestrial                -0.945  0.191  0.898                         
    ## duration                        -0.026  0.009  0.009  0.001                  
    ## abs(temptrend):REALMMarine       0.193 -0.954 -0.194 -0.182  0.000           
    ## abs(temptrend):REALMTerrestrial  0.183 -0.907 -0.174 -0.226  0.002  0.866    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -2.216413e+09 -4.491508e-04 -2.614992e-08  5.486808e-06  1.590771e+11 
    ## 
    ## Number of Observations: 27962
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    179                  27962

``` r
summary(modTDtrendHornz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC    BIC    logLik
    ##   148805.9 148913 -74389.96
    ## 
    ## Random effects:
    ##  Formula: ~abs(temptrend) | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                StdDev   Corr  
    ## (Intercept)     4.05154 (Intr)
    ## abs(temptrend) 35.83091 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)    Residual
    ## StdDev:    3.393784 21260928001
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##    power 
    ## -21.4449 
    ## Fixed effects: Horntrendz ~ abs(temptrend) * REALM + duration 
    ##                                      Value Std.Error    DF   t-value p-value
    ## (Intercept)                       0.709720  1.711704 27690  0.414628  0.6784
    ## abs(temptrend)                   14.952621 16.393423 27690  0.912111  0.3617
    ## REALMMarine                       0.946606  1.787718   167  0.529505  0.5972
    ## REALMTerrestrial                  1.107501  1.826982   167  0.606191  0.5452
    ## duration                          0.025270  0.003661 27690  6.902895  0.0000
    ## abs(temptrend):REALMMarine      -20.347019 17.235275 27690 -1.180545  0.2378
    ## abs(temptrend):REALMTerrestrial  -8.985469 18.053231 27690 -0.497721  0.6187
    ##  Correlation: 
    ##                                 (Intr) abs(t) REALMM REALMT duratn a():REALMM
    ## abs(temptrend)                  -0.229                                       
    ## REALMMarine                     -0.956  0.219                                
    ## REALMTerrestrial                -0.935  0.214  0.895                         
    ## duration                        -0.044  0.010  0.013  0.002                  
    ## abs(temptrend):REALMMarine       0.218 -0.951 -0.224 -0.204  0.001           
    ## abs(temptrend):REALMTerrestrial  0.208 -0.908 -0.199 -0.258  0.002  0.864    
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -4.454434e+07 -2.687993e-04 -2.168357e-08  2.872747e-06  5.509914e+08 
    ## 
    ## Number of Observations: 27864
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    170                  27864

### Plot the temp coefficients from TD models

``` r
# make table of coefficients
coefs1 <- as.data.frame(summary(modTDtrendJturem0)$tTable)
coefs2 <- as.data.frame(summary(modTDtrendJbetarem0)$tTable)
coefs3 <- as.data.frame(summary(modTDtrendHornrem0)$tTable)
coefs4 <- as.data.frame(summary(modTDtrendJtuz)$tTable)
coefs5 <- as.data.frame(summary(modTDtrendJbetaz)$tTable)
coefs6 <- as.data.frame(summary(modTDtrendHornz)$tTable)
coefs1$response <- 'Jtu'
coefs2$response <- 'Jbeta'
coefs3$response <- 'Horn'
coefs4$response <- 'Jtu'
coefs5$response <- 'Jbeta'
coefs6$response <- 'Horn'
coefs1$fit <- 'rem0'
coefs2$fit <- 'rem0'
coefs3$fit <- 'rem0'
coefs4$fit <- 'z'
coefs5$fit <- 'z'
coefs6$fit <- 'z'
rows1 <- which(grepl('temptrend', rownames(coefs1))) # extract temperature effect
cols <- c('Value', 'Std.Error', 'response', 'fit')
allcoefs <- rbind(coefs1[rows1, cols], coefs2[rows1, cols], coefs3[rows1, cols],
                  coefs4[rows1, cols], coefs5[rows1, cols], coefs6[rows1, cols])
allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMMarine', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to marine effects
allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] <- 
  allcoefs$Value[grepl('REALMTerrestrial', rownames(allcoefs))] + 
  allcoefs$Value[!grepl('REALM', rownames(allcoefs))] # add intercept to terrestrial effects

allcoefs$lCI <- allcoefs$Value - allcoefs$Std.Error # lower confidence interval
allcoefs$uCI <- allcoefs$Value + allcoefs$Std.Error
allcoefs$realm <- rep(c('Freshwater', 'Marine', 'Terrestrial'), 6)
allcoefs$group <- paste0(allcoefs$response, allcoefs$fit)

pd <- position_dodge(width = 0.5)
ggplot(subset(allcoefs, fit == 'rem0'), aes(x = realm, y = Value, color = response, shape = fit, group = group)) +
  geom_point(position = pd) + 
  geom_errorbar(aes(ymin=lCI, ymax=uCI), width=.1, position=pd) +
  coord_flip() +
  labs(y = 'Turnover per |°C/yr|')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/modTD%20coefs-1.png)<!-- -->

``` r
ggplot(subset(allcoefs, fit == 'z'), aes(x = realm, y = Value, color = response, shape = fit, group = group)) +
  geom_point(position = pd) + 
  geom_errorbar(aes(ymin=lCI, ymax=uCI), width=.1, position=pd) +
  coord_flip() +
  labs(y = 'Standardized Turnover per |°C/yr|')
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/modTD%20coefs-2.png)<!-- -->

## Full models

Try static covariates plus interactions of abs temperature trend with
each covariate:

  - realm
  - speed
  - mass
  - average metabolic temperature
  - consumer fraction
  - environmental temperature
  - seasonality
  - microclimates
  - thermal bias
  - NPP
  - vegetation
  - duration
  - human footprint

Except for thermal bias: interact with temperature trend (not abs)

### Bowler vs Venter/Halpern human impact

Bowler has lower AIC.

``` r
# using Bowler for human impact
i1 <- trends[, complete.cases(Jtutrendrem0, REALM, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc, human_footprint.sc)]

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

if(file.exists('temp/modTfullbowlerrem0.rds')){
  modTfullbowlerrem0 <- readRDS('temp/modTfullbowlerrem0.rds')
} else {

modTfullbowlerrem0 <- lme(Jtutrendrem0 ~ temptrend_abs.sc*REALM +
                     temptrend_abs.sc*tsign + 
                     temptrend_abs.sc*tempave_metab.sc + 
                     temptrend_abs.sc*seas.sc + 
                     temptrend_abs.sc*microclim.sc + 
                     temptrend_abs.sc*mass.sc + 
                     temptrend_abs.sc*speed.sc + 
                     temptrend_abs.sc*consumerfrac.sc +
                     temptrend_abs.sc*nspp.sc +
                     temptrend_abs.sc*thermal_bias.sc:tsign +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*veg.sc +
                     temptrend_abs.sc*duration.sc +
                     temptrend_abs.sc*human_bowler.sc:REALM2,
                   random = randef, weights = varef, data = trends[i1,], method = 'REML')
  saveRDS(modTfullbowlerrem0, file = 'temp/modTfullbowlerrem0.rds')
}

# using Venter/Halpern for human impact
if(file.exists('temp/modTfullfootprintrem0.rds')){
  modTfullfootprintrem0 <- readRDS('temp/modTfullfootprintrem0.rds')
} else {
modTfullfootprintrem0 <- lme(Jtutrendrem0 ~ temptrend_abs.sc*REALM + 
                     temptrend_abs.sc*tsign + 
                     temptrend_abs.sc*tempave_metab.sc + 
                     temptrend_abs.sc*seas.sc + 
                     temptrend_abs.sc*microclim.sc + 
                     temptrend_abs.sc*mass.sc + 
                     temptrend_abs.sc*speed.sc + 
                     temptrend_abs.sc*consumerfrac.sc +
                     temptrend_abs.sc*nspp.sc +
                     temptrend_abs.sc*thermal_bias.sc:tsign +
                     temptrend_abs.sc*npp.sc +
                     temptrend_abs.sc*veg.sc +
                     temptrend_abs.sc*duration.sc +
                     temptrend_abs.sc*human_footprint.sc:REALM2,
                   random = randef, weights = varef, data = trends[i1,], method = 'REML',
                   control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullfootprintrem0, file = 'temp/modTfullfootprintrem0.rds')

}
AIC(modTfullbowlerrem0, modTfullfootprintrem0)
```

    ##                       df       AIC
    ## modTfullbowlerrem0    42 -101579.8
    ## modTfullfootprintrem0 42 -101568.7

### Fit/load full models

``` r
randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

terms <- c('temptrend_abs.sc*REALM', 
           'temptrend_abs.sc*tsign',
           'temptrend_abs.sc*tempave_metab.sc',
           'temptrend_abs.sc*seas.sc',
           'temptrend_abs.sc*microclim.sc',
           'temptrend_abs.sc*mass.sc',
           'temptrend_abs.sc*speed.sc', 
           'temptrend_abs.sc*consumerfrac.sc',
           'temptrend_abs.sc*nspp.sc',
           'temptrend_abs.sc*thermal_bias.sc:tsign',
           'temptrend_abs.sc*npp.sc',
           'temptrend_abs.sc*veg.sc',
           'temptrend_abs.sc*duration.sc',
           'temptrend_abs.sc*human_bowler.sc:REALM2')

completerows <- trends[, complete.cases(REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]


# rem0
if(file.exists('temp/modTfullJturem0.rds')){
  modTfullJturem0 <- readRDS('temp/modTfullJturem0.rds')
} else {
  i <- trends[, complete.cases(Jtutrendrem0)] & completerows
  thisformula <- as.formula(paste0('Jtutrendrem0 ~ ', paste(terms, collapse = ' + ')))
  modTfullJturem0 <- eval(bquote(lme(.(thisformula),
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML')))
  saveRDS(modTfullJturem0, file = 'temp/modTfullJturem0.rds')
}

if(file.exists('temp/modTfullJbetarem0.rds')){
  modTfullJbetarem0 <- readRDS('temp/modTfullJbetarem0.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendrem0)] & completerows
  thisformula <- as.formula(paste0('Jbetatrendrem0 ~ ', paste(terms, collapse = ' + ')))
  modTfullJbetarem0 <- eval(bquote(lme(.(thisformula),
                           random = randef, weights = varef, data = trends[i,], 
                           method = 'REML')))
  saveRDS(modTfullJbetarem0, file = 'temp/modTfullJbetarem0.rds')
}

if(file.exists('temp/modTfullHornrem0.rds')){
  modTfullHornrem0 <- readRDS('temp/modTfullHornrem0.rds')
} else {
  i <- trends[, complete.cases(Horntrendrem0)] & completerows
  thisformula <- as.formula(paste0('Horntrendrem0 ~ ', paste(terms, collapse = ' + ')))
  modTfullHornrem0 <- eval(bquote(lme(.(thisformula),
                      random = randef, weights = varef, data = trends[i,], 
                      method = 'REML')))
  saveRDS(modTfullHornrem0, file = 'temp/modTfullHornrem0.rds')
}

# trendz
if(file.exists('temp/modTfullJtuz.rds')){
  modTfullJtuz <- readRDS('temp/modTfullJtuz.rds')
} else {
  i <- trends[, complete.cases(Jtutrendz)] & completerows
  thisformula <- as.formula(paste0('Jtutrendz ~ ', paste(terms, collapse = ' + ')))
  modTfullJtuz <- eval(bquote(lme(.(thisformula),
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                   control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))))
  saveRDS(modTfullJtuz, file = 'temp/modTfullJtuz.rds')
}

if(file.exists('temp/modTfullJbetaz.rds')){
  modTfullJbetaz <- readRDS('temp/modTfullJbetaz.rds')
} else {
  i <- trends[, complete.cases(Jbetatrendz)] & completerows
  thisformula <- as.formula(paste0('Jbetatrendz ~ ', paste(terms, collapse = ' + ')))
  modTfullJbetaz <- eval(bquote(lme(.(thisformula),
                           random = randef, weights = varef, data = trends[i,], 
                           method = 'REML')))
  saveRDS(modTfullJbetaz, file = 'temp/modTfullJbetaz.rds')
}

if(file.exists('temp/modTfullHornz.rds')){
  modTfullHornz <- readRDS('temp/modTfullHornz.rds')
} else {
  i <- trends[, complete.cases(Horntrendz)] & completerows
  thisformula <- as.formula(paste0('Horntrendz ~ ', paste(terms, collapse = ' + ')))
  modTfullHornz <- eval(bquote(lme(.(thisformula),
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                         control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))))
  saveRDS(modTfullHornz, file = 'temp/modTfullHornz.rds')
}

# last
if(file.exists('temp/modTfullJtulast.rds')){
  modTfullJtulast <- readRDS('temp/modTfullJtulast.rds')
} else {
  i <- trends[, complete.cases(Jtulast)] & completerows
  thisformula <- as.formula(paste0('Jtulast ~ ', paste(terms, collapse = ' + ')))
  modTfullJtulast <- eval(bquote(lme(.(thisformula),
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                   control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))))
  saveRDS(modTfullJtulast, file = 'temp/modTfullJtulast.rds')
}

if(file.exists('temp/modTfullJbetalast.rds')){
  modTfullJbetalast <- readRDS('temp/modTfullJbetalast.rds')
} else {
  i <- trends[, complete.cases(Jbetalast)] & completerows
  thisformula <- as.formula(paste0('Jbetalast ~ ', paste(terms, collapse = ' + ')))
  modTfullJbetalast <- eval(bquote(lme(.(thisformula),
                           random = randef, weights = varef, data = trends[i,], 
                           method = 'REML')))
  saveRDS(modTfullJbetalast, file = 'temp/modTfullJbetalast.rds')
}

if(file.exists('temp/modTfullHornlast.rds')){
  modTfullHornlast <- readRDS('temp/modTfullHornlast.rds')
} else {
  i <- trends[, complete.cases(Hornlast)] & completerows
  thisformula <- as.formula(paste0('Hornlast ~ ', paste(terms, collapse = ' + ')))
  modTfullHornlast <- eval(bquote(lme(.(thisformula),
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                         control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))))
  saveRDS(modTfullHornlast, file = 'temp/modTfullHornlast.rds')
}


# exp
if(file.exists('temp/modTfullJtuexp.rds')){
  modTfullJtuexp <- readRDS('temp/modTfullJtuexp.rds')
} else {
  i <- trends[, complete.cases(Jtuexp)] & completerows
  thisformula <- as.formula(paste0('log(Jtuexp) ~ ', paste(terms, collapse = ' + ')))
  modTfullJtuexp <- eval(bquote(lme(.(thisformula),
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                   control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))))
  saveRDS(modTfullJtuexp, file = 'temp/modTfullJtuexp.rds')
}

if(file.exists('temp/modTfullJbetaexp.rds')){
  modTfullJbetaexp <- readRDS('temp/modTfullJbetaexp.rds')
} else {
  i <- trends[, complete.cases(Jbetaexp)] & completerows
  thisformula <- as.formula(paste0('log(Jbetaexp) ~ ', paste(terms, collapse = ' + ')))
  modTfullJbetaexp <- eval(bquote(lme(.(thisformula),
                           random = randef, weights = varef, data = trends[i,], 
                           method = 'REML')))
  saveRDS(modTfullJbetaexp, file = 'temp/modTfullJbetaexp.rds')
}

if(file.exists('temp/modTfullHornexp.rds')){
  modTfullHornexp <- readRDS('temp/modTfullHornexp.rds')
} else {
  i <- trends[, complete.cases(Hornexp)] & completerows
  thisformula <- as.formula(paste0('log(Hornexp) ~ ', paste(terms, collapse = ' + ')))
  modTfullHornexp <- eval(bquote(lme(.(thisformula),
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                         control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))))
  saveRDS(modTfullHornexp, file = 'temp/modTfullHornexp.rds')
}

# mm
if(file.exists('temp/modTfullJtumm.rds')){
  modTfullJtumm <- readRDS('temp/modTfullJtumm.rds')
} else {
  i <- trends[, complete.cases(Jtumm)] & completerows
  thisformula <- as.formula(paste0('log(Jtumm+1) ~ ', paste(terms, collapse = ' + ')))
  modTfullJtumm <- eval(bquote(lme(.(thisformula),
                       random = randef, weights = varef, data = trends[i,], 
                       method = 'REML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))))
  saveRDS(modTfullJtumm, file = 'temp/modTfullJtumm.rds')
}

if(file.exists('temp/modTfullJbetamm.rds')){
  modTfullJbetamm <- readRDS('temp/modTfullJbetamm.rds')
} else {
  i <- trends[, complete.cases(Jbetamm)] & completerows
  thisformula <- as.formula(paste0('log(Jbetamm+1) ~ ', paste(terms, collapse = ' + ')))
  modTfullJbetamm <- eval(bquote(lme(.(thisformula),
                         random = randef, weights = varef, data = trends[i,], 
                         method = 'REML',
                         control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))))
  saveRDS(modTfullJbetamm, file = 'temp/modTfullJbetamm.rds')
}

if(file.exists('temp/modTfullHornmm.rds')){
  modTfullHornmm <- readRDS('temp/modTfullHornmm.rds')
} else {
  i <- trends[, complete.cases(Hornmm)] & completerows
  thisformula <- as.formula(paste0('log(Hornmm+1) ~ ', paste(terms, collapse = ' + ')))
  modTfullHornmm <- eval(bquote(lme(.(thisformula),
                        random = randef, weights = varef, data = trends[i,], 
                        method = 'REML',
                        control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))))
  saveRDS(modTfullHornmm, file = 'temp/modTfullHornmm.rds')
}
```

#### Summary rem0

``` r
summary(modTfullJturem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -93123.78 -92773.81 46603.89
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.009644446 (Intr)
    ## temptrend_abs.sc 0.020144635 -0.977
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01153819 2.197503
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.262967 
    ## Fixed effects: Jtutrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.002481664 0.005369119 30525  0.462211  0.6439
    ## temptrend_abs.sc                                  0.009012417 0.012343709 30525  0.730122  0.4653
    ## REALMMarine                                       0.005710929 0.005769345   195  0.989875  0.3235
    ## REALMTerrestrial                                  0.005551065 0.005132316   195  1.081591  0.2808
    ## tsign1                                           -0.001908202 0.000641857 30525 -2.972939  0.0030
    ## tempave_metab.sc                                 -0.002570908 0.000880134 30525 -2.921042  0.0035
    ## seas.sc                                          -0.000379555 0.000659735 30525 -0.575315  0.5651
    ## microclim.sc                                      0.000820134 0.000346923 30525  2.364022  0.0181
    ## mass.sc                                          -0.000910860 0.000664487 30525 -1.370771  0.1705
    ## speed.sc                                          0.001336053 0.000618452 30525  2.160320  0.0308
    ## consumerfrac.sc                                  -0.000569149 0.000333709 30525 -1.705528  0.0881
    ## nspp.sc                                          -0.000490547 0.000560685 30525 -0.874906  0.3816
    ## npp.sc                                            0.000631900 0.000490255 30525  1.288921  0.1974
    ## veg.sc                                           -0.000516012 0.000965994 30525 -0.534177  0.5932
    ## duration.sc                                       0.000115399 0.000525291 30525  0.219686  0.8261
    ## temptrend_abs.sc:REALMMarine                     -0.010752874 0.013061622 30525 -0.823242  0.4104
    ## temptrend_abs.sc:REALMTerrestrial                -0.003722700 0.012358226 30525 -0.301233  0.7632
    ## temptrend_abs.sc:tsign1                          -0.002171439 0.001941946 30525 -1.118177  0.2635
    ## temptrend_abs.sc:tempave_metab.sc                 0.002193201 0.002326819 30525  0.942575  0.3459
    ## temptrend_abs.sc:seas.sc                          0.001424990 0.001564288 30525  0.910951  0.3623
    ## temptrend_abs.sc:microclim.sc                    -0.003228339 0.000985204 30525 -3.276822  0.0011
    ## temptrend_abs.sc:mass.sc                         -0.001229783 0.001502883 30525 -0.818283  0.4132
    ## temptrend_abs.sc:speed.sc                        -0.000864735 0.001646884 30525 -0.525073  0.5995
    ## temptrend_abs.sc:consumerfrac.sc                  0.004292700 0.000761042 30525  5.640559  0.0000
    ## temptrend_abs.sc:nspp.sc                          0.000768309 0.001404359 30525  0.547089  0.5843
    ## tsign-1:thermal_bias.sc                          -0.000878946 0.000701330 30525 -1.253255  0.2101
    ## tsign1:thermal_bias.sc                           -0.000087927 0.000487335 30525 -0.180424  0.8568
    ## temptrend_abs.sc:npp.sc                           0.002863289 0.001365981 30525  2.096142  0.0361
    ## temptrend_abs.sc:veg.sc                          -0.001744900 0.002180742 30525 -0.800141  0.4236
    ## temptrend_abs.sc:duration.sc                     -0.001067058 0.001066617 30525 -1.000414  0.3171
    ## human_bowler.sc:REALM2TerrFresh                  -0.000105798 0.000736866 30525 -0.143579  0.8858
    ## human_bowler.sc:REALM2Marine                      0.001149545 0.000403707 30525  2.847475  0.0044
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.000853820 0.001451395 30525 -0.588275  0.5564
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.000554015 0.001240152 30525  0.446731  0.6551
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.001211326 0.001436734 30525 -0.843111  0.3992
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.004756799 0.001086198 30525 -4.379311  0.0000
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.766                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.950  0.723                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.693  0.538  0.641                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.095  0.074  0.007  0.016                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.051 -0.022 -0.059 -0.245  0.059                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.143  0.138  0.199 -0.074 -0.076  0.091                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.069  0.055  0.080  0.021  0.018 -0.203  0.112                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.094 -0.059 -0.057 -0.042 -0.022  0.073  0.153  0.015                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.080 -0.074 -0.044 -0.084 -0.060 -0.152 -0.076  0.070 -0.548                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.019  0.019  0.013  0.115  0.021 -0.135 -0.047  0.015  0.027 -0.124                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.046 -0.025 -0.079 -0.069  0.056 -0.211  0.008 -0.100  0.067  0.131  0.102                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.048 -0.046 -0.072  0.039  0.036  0.021 -0.155 -0.145 -0.012  0.164 -0.034 -0.140                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.558  0.380  0.579  0.021 -0.010  0.019  0.083 -0.016  0.020 -0.031 -0.008  0.018 -0.191                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.041  0.015  0.043  0.024 -0.165  0.098 -0.027 -0.027 -0.075  0.015  0.010 -0.266  0.056 -0.010                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.735 -0.950 -0.779 -0.501 -0.018  0.032 -0.185 -0.067  0.027  0.042 -0.015  0.051  0.072 -0.403 -0.027                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.520 -0.662 -0.477 -0.796 -0.017  0.191  0.071  0.008  0.030  0.074 -0.096  0.040 -0.036  0.022 -0.017  0.607                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.068 -0.171 -0.024 -0.018 -0.513 -0.052  0.023 -0.013  0.007  0.031 -0.009 -0.027 -0.037 -0.015  0.040  0.056      0.023                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.020  0.018  0.025  0.176 -0.059 -0.781 -0.077  0.136 -0.052  0.121  0.114  0.180 -0.023 -0.021 -0.039 -0.011     -0.240      0.037                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.133 -0.158 -0.170  0.079  0.057 -0.059 -0.773 -0.112 -0.102  0.048  0.024  0.008  0.131 -0.101 -0.017  0.206     -0.168     -0.019  0.065                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.042 -0.077 -0.052 -0.001 -0.028  0.143 -0.076 -0.770  0.016 -0.059  0.000  0.070  0.148  0.013  0.020  0.095     -0.024      0.047 -0.102   0.121                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.058  0.056  0.020  0.018 -0.001 -0.073 -0.092  0.008 -0.741  0.426 -0.006 -0.032 -0.013 -0.003  0.067 -0.016      0.009     -0.014  0.176   0.072            -0.040                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.071  0.090  0.040  0.080  0.058  0.144  0.044 -0.050  0.336 -0.717  0.083 -0.120 -0.099 -0.002 -0.014 -0.059     -0.113     -0.023 -0.333  -0.015             0.053            -0.571                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.024 -0.034 -0.018 -0.100 -0.009  0.107  0.033  0.014 -0.007  0.102 -0.782 -0.084  0.015  0.000 -0.041  0.024      0.115      0.038 -0.142  -0.059            -0.029            -0.026            -0.092                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.014  0.022  0.037  0.048 -0.039  0.179  0.013  0.068 -0.019 -0.119 -0.061 -0.759  0.088 -0.027  0.211 -0.056     -0.085      0.023 -0.163  -0.007            -0.050             0.074             0.108             0.124                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.065 -0.042 -0.062 -0.055 -0.121  0.181 -0.244 -0.044 -0.027 -0.025  0.000 -0.085 -0.012 -0.004  0.100  0.054      0.042     -0.020 -0.139   0.180             0.034             0.022             0.027            -0.006            0.056                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.092 -0.077 -0.107 -0.111  0.067  0.288 -0.355 -0.137 -0.060 -0.008 -0.025 -0.102 -0.071  0.005  0.097  0.082      0.100      0.023 -0.246   0.235             0.089             0.020             0.033             0.023            0.061             0.386                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.042  0.072  0.061 -0.032 -0.023 -0.012  0.122  0.142 -0.016 -0.093 -0.005  0.088 -0.753  0.147 -0.040 -0.114      0.072     -0.001  0.035  -0.197            -0.279             0.038             0.076            -0.009           -0.096             0.000  0.054                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.379 -0.517 -0.397  0.036 -0.009 -0.029 -0.115  0.009 -0.017  0.018  0.009 -0.022  0.180 -0.705  0.002  0.550     -0.060      0.055  0.019   0.155             0.020            -0.001            -0.003            -0.002            0.025             0.017 -0.015 -0.254                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.004  0.049 -0.040 -0.055  0.053  0.014 -0.034 -0.015  0.079 -0.047 -0.006  0.266 -0.031  0.023 -0.553  0.027      0.045     -0.171  0.033   0.009             0.044            -0.102             0.027             0.011           -0.365            -0.035 -0.085  0.014            -0.009                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.117  0.156  0.121 -0.041  0.020  0.036  0.049  0.078  0.049 -0.053 -0.040 -0.014 -0.116  0.220 -0.021 -0.162      0.040     -0.030 -0.028  -0.115             0.004            -0.007             0.010             0.026            0.023             0.058  0.101  0.065            -0.286            0.011                                                               
    ## human_bowler.sc:REALM2Marine                      0.048 -0.043 -0.052 -0.046 -0.010  0.110 -0.191  0.022 -0.035  0.072 -0.088 -0.055 -0.170  0.025  0.007  0.047      0.038      0.004 -0.075   0.156            -0.018             0.039            -0.070             0.103            0.052             0.114  0.161  0.126            -0.016            0.032            0.026                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.038  0.030  0.052  0.044 -0.065 -0.159  0.233  0.085  0.010  0.021  0.008  0.035 -0.032  0.003 -0.059 -0.068     -0.060      0.241  0.197  -0.268            -0.128             0.053            -0.091             0.014           -0.003            -0.515 -0.350  0.063            -0.021            0.012           -0.055      -0.090                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.056  0.081  0.064  0.074  0.001 -0.250  0.224  0.112  0.019  0.029  0.007  0.073  0.057 -0.016 -0.054 -0.079     -0.105     -0.073  0.343  -0.271            -0.225             0.052            -0.162             0.008           -0.056            -0.297 -0.731 -0.022             0.013            0.052           -0.075      -0.138       0.440                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.189 -0.282 -0.198  0.040 -0.011 -0.066 -0.136 -0.042 -0.023  0.041  0.025  0.020  0.114 -0.324  0.024  0.297     -0.066      0.057  0.079   0.184             0.053             0.016            -0.048            -0.026           -0.024            -0.043 -0.066 -0.155             0.499           -0.032           -0.650      -0.008       0.095  0.088               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.034  0.044  0.038  0.041 -0.003 -0.102  0.146 -0.006  0.047 -0.075  0.153  0.052  0.132 -0.017 -0.007 -0.049     -0.047     -0.005  0.113  -0.175            -0.034            -0.035             0.066            -0.189           -0.061            -0.076 -0.140 -0.179             0.032           -0.039           -0.013      -0.790       0.101  0.155  0.013        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -6.87855240 -0.23743170 -0.01238702  0.27081091  6.63391080 
    ## 
    ## Number of Observations: 30756
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30756

``` r
summary(modTfullJbetarem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC     BIC   logLik
    ##   -120468.9 -120119 60276.46
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.008978352 (Intr)
    ## temptrend_abs.sc 0.007179505 -0.371
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev: 0.006180149 1.145676
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.109054 
    ## Fixed effects: Jbetatrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                         Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.000526319 0.003997326 30525  0.131668  0.8952
    ## temptrend_abs.sc                                  0.025586661 0.007440079 30525  3.439031  0.0006
    ## REALMMarine                                       0.006490272 0.004274850   195  1.518245  0.1306
    ## REALMTerrestrial                                  0.012819339 0.004066917   195  3.152102  0.0019
    ## tsign1                                           -0.000976013 0.000406444 30525 -2.401346  0.0163
    ## tempave_metab.sc                                 -0.003550400 0.000586361 30525 -6.054976  0.0000
    ## seas.sc                                          -0.000434312 0.000414010 30525 -1.049037  0.2942
    ## microclim.sc                                      0.000766065 0.000217850 30525  3.516483  0.0004
    ## mass.sc                                          -0.001175952 0.000470879 30525 -2.497356  0.0125
    ## speed.sc                                          0.000760916 0.000419846 30525  1.812371  0.0699
    ## consumerfrac.sc                                   0.000050126 0.000212522 30525  0.235861  0.8135
    ## nspp.sc                                          -0.000306051 0.000358942 30525 -0.852648  0.3939
    ## npp.sc                                            0.000222679 0.000305043 30525  0.729993  0.4654
    ## veg.sc                                           -0.000099515 0.000588104 30525 -0.169213  0.8656
    ## duration.sc                                      -0.001252768 0.000344091 30525 -3.640809  0.0003
    ## temptrend_abs.sc:REALMMarine                     -0.025132990 0.007779287 30525 -3.230758  0.0012
    ## temptrend_abs.sc:REALMTerrestrial                -0.024306266 0.007531347 30525 -3.227347  0.0013
    ## temptrend_abs.sc:tsign1                          -0.001427815 0.001279393 30525 -1.116009  0.2644
    ## temptrend_abs.sc:tempave_metab.sc                 0.003612654 0.001479523 30525  2.441770  0.0146
    ## temptrend_abs.sc:seas.sc                          0.001116919 0.000977648 30525  1.142455  0.2533
    ## temptrend_abs.sc:microclim.sc                    -0.002526463 0.000631422 30525 -4.001227  0.0001
    ## temptrend_abs.sc:mass.sc                          0.000634951 0.000953284 30525  0.666067  0.5054
    ## temptrend_abs.sc:speed.sc                        -0.000610754 0.001066487 30525 -0.572679  0.5669
    ## temptrend_abs.sc:consumerfrac.sc                  0.000913361 0.000488389 30525  1.870150  0.0615
    ## temptrend_abs.sc:nspp.sc                          0.001776058 0.000884584 30525  2.007788  0.0447
    ## tsign-1:thermal_bias.sc                          -0.000360076 0.000452426 30525 -0.795879  0.4261
    ## tsign1:thermal_bias.sc                           -0.000086001 0.000318544 30525 -0.269982  0.7872
    ## temptrend_abs.sc:npp.sc                           0.001604281 0.000890181 30525  1.802197  0.0715
    ## temptrend_abs.sc:veg.sc                          -0.001788501 0.001395844 30525 -1.281305  0.2001
    ## temptrend_abs.sc:duration.sc                     -0.000207074 0.000676551 30525 -0.306074  0.7596
    ## human_bowler.sc:REALM2TerrFresh                   0.000262674 0.000445616 30525  0.589462  0.5556
    ## human_bowler.sc:REALM2Marine                      0.000241525 0.000252690 30525  0.955812  0.3392
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.001312179 0.000937816 30525 -1.399185  0.1618
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.000199820 0.000790124 30525  0.252896  0.8004
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.001044344 0.000909954 30525 -1.147689  0.2511
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.000657002 0.000694171 30525 -0.946456  0.3439
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.598                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.943  0.572                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.763  0.390  0.708                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.077  0.069  0.003  0.006                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.090 -0.015 -0.085 -0.232  0.064                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.117  0.096  0.167 -0.070 -0.080  0.058                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.070  0.067  0.081  0.034  0.018 -0.232  0.121                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.133 -0.032 -0.091 -0.053  0.000  0.069  0.110  0.023                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.033 -0.050 -0.015 -0.072 -0.077 -0.075 -0.026  0.048 -0.612                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.003  0.019  0.004  0.068  0.006 -0.109 -0.033  0.017  0.012 -0.068                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.029 -0.044 -0.065 -0.058  0.051 -0.225  0.035 -0.095  0.019  0.154  0.111                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.025 -0.014 -0.044  0.042  0.043  0.030 -0.125 -0.117 -0.025  0.173 -0.036 -0.145                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.447  0.393  0.467  0.014 -0.016  0.006  0.081 -0.010  0.020 -0.034 -0.011  0.015 -0.195                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.036  0.016  0.039  0.014 -0.150  0.129 -0.027 -0.037 -0.052  0.001  0.009 -0.257  0.049 -0.008                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.581 -0.961 -0.596 -0.367 -0.018  0.020 -0.123 -0.075  0.013  0.036 -0.015  0.068  0.049 -0.422 -0.021                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.394 -0.634 -0.363 -0.589 -0.004  0.164  0.136  0.008 -0.011  0.094 -0.066  0.081 -0.062  0.043 -0.009  0.585                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.056 -0.182 -0.022  0.001 -0.463 -0.054  0.003 -0.025  0.016  0.027 -0.004 -0.008 -0.047 -0.018  0.043  0.062      0.019                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                 0.010  0.002 -0.020  0.131 -0.091 -0.509 -0.107  0.080 -0.080  0.169  0.095  0.125 -0.024 -0.019 -0.056  0.010     -0.280      0.039                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.058 -0.126 -0.078  0.116  0.055 -0.060 -0.687 -0.108 -0.081  0.029  0.020  0.008  0.104 -0.091 -0.020  0.167     -0.245     -0.023  0.076                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.040 -0.089 -0.045 -0.017 -0.033  0.155 -0.066 -0.757  0.002 -0.040 -0.010  0.062  0.124  0.015  0.018  0.104     -0.035      0.047 -0.048   0.135                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.033  0.025  0.024 -0.016 -0.010 -0.080 -0.079  0.007 -0.612  0.392  0.014 -0.036 -0.011 -0.008  0.077  0.007      0.059     -0.008  0.218   0.062            -0.036                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.046  0.087  0.049  0.086  0.069  0.122  0.061 -0.032  0.313 -0.590  0.057 -0.105 -0.079  0.005 -0.010 -0.063     -0.139     -0.019 -0.426   0.000             0.041            -0.583                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.002 -0.029 -0.008 -0.049 -0.002  0.084  0.027  0.014 -0.008  0.059 -0.736 -0.083  0.011  0.003 -0.022  0.021      0.099      0.028 -0.119  -0.071            -0.030            -0.039            -0.091                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.041  0.050  0.052  0.077 -0.027  0.121 -0.002  0.059 -0.042 -0.106 -0.060 -0.643  0.079 -0.024  0.173 -0.088     -0.128      0.013 -0.127  -0.021            -0.048             0.089             0.101             0.110                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.069 -0.028 -0.065 -0.057 -0.122  0.200 -0.257 -0.053 -0.016 -0.007  0.004 -0.104 -0.015 -0.002  0.107  0.035      0.024     -0.023 -0.041   0.156             0.034             0.016             0.003            -0.007            0.043                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.110 -0.061 -0.116 -0.115  0.067  0.332 -0.399 -0.157 -0.024  0.001 -0.021 -0.138 -0.074  0.000  0.122  0.049      0.079      0.040 -0.126   0.201             0.084             0.011             0.016             0.016            0.049             0.391                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.012  0.048 -0.002 -0.057 -0.025  0.001  0.068  0.124 -0.001 -0.091  0.003  0.072 -0.704  0.133 -0.036 -0.098      0.105      0.016  0.059  -0.200            -0.289             0.037             0.044            -0.003           -0.080             0.008  0.078                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.304 -0.533 -0.318  0.040 -0.006 -0.028 -0.099  0.007 -0.013  0.015  0.007 -0.021  0.175 -0.716  0.007  0.576     -0.082      0.054  0.010   0.161             0.027            -0.002             0.009            -0.006            0.018             0.010 -0.022 -0.264                                                                                                  
    ## temptrend_abs.sc:duration.sc                      0.004  0.050 -0.042 -0.053  0.049  0.037 -0.033 -0.015  0.076 -0.032  0.002  0.231 -0.027  0.015 -0.467  0.034      0.058     -0.167  0.042   0.001             0.029            -0.121             0.007             0.030           -0.357            -0.023 -0.074  0.024            -0.008                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.088  0.138  0.090 -0.030  0.020  0.027  0.033  0.083  0.035 -0.044 -0.026 -0.004 -0.131  0.211 -0.021 -0.150      0.059     -0.022 -0.015  -0.113            -0.010             0.002             0.002             0.018            0.014             0.048  0.108  0.074            -0.274            0.022                                                               
    ## human_bowler.sc:REALM2Marine                      0.045 -0.037 -0.050 -0.046 -0.011  0.123 -0.198  0.022 -0.034  0.077 -0.096 -0.049 -0.185  0.027  0.010  0.042      0.037      0.010 -0.053   0.155            -0.023             0.040            -0.078             0.114            0.044             0.122  0.168  0.140            -0.020            0.042            0.028                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.017  0.011  0.017  0.035 -0.072 -0.077  0.193  0.053  0.000  0.037 -0.001  0.033 -0.030  0.008 -0.070 -0.045     -0.055      0.267  0.181  -0.276            -0.109             0.058            -0.130             0.028           -0.005            -0.431 -0.281  0.073            -0.028            0.032           -0.056      -0.077                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.026  0.065  0.019  0.055 -0.009 -0.125  0.190  0.073 -0.008  0.050 -0.003  0.056  0.057 -0.019 -0.073 -0.053     -0.108     -0.074  0.328  -0.271            -0.195             0.067            -0.205             0.026           -0.047            -0.228 -0.639 -0.015             0.013            0.067           -0.064      -0.126       0.438                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.131 -0.272 -0.140  0.050 -0.014 -0.040 -0.115 -0.053 -0.019  0.031  0.013  0.013  0.108 -0.325  0.024  0.293     -0.097      0.054  0.073   0.180             0.067             0.011            -0.041            -0.019           -0.022            -0.032 -0.063 -0.160             0.497           -0.029           -0.647      -0.009       0.100  0.089               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.029  0.044  0.035  0.037  0.000 -0.088  0.157 -0.009  0.042 -0.070  0.169  0.038  0.146 -0.019  0.008 -0.049     -0.046     -0.011  0.097  -0.178            -0.040            -0.031             0.067            -0.188           -0.050            -0.081 -0.136 -0.175             0.032           -0.042           -0.020      -0.776       0.087  0.144  0.011        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -8.85386366 -0.32283102 -0.01276817  0.32602346  7.98731227 
    ## 
    ## Number of Observations: 30756
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30756

``` r
summary(modTfullHornrem0)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -87368.55 -87018.57 43726.27
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.01717909 (Intr)
    ## temptrend_abs.sc 0.02379603 0.103 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01937619 2.727363
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.444923 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                        Value   Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.00540876 0.008414279 30525  0.642807  0.5204
    ## temptrend_abs.sc                                  0.03375995 0.015401478 30525  2.191994  0.0284
    ## REALMMarine                                       0.00596674 0.008995781   195  0.663283  0.5079
    ## REALMTerrestrial                                  0.01818339 0.008394240   195  2.166175  0.0315
    ## tsign1                                           -0.00200925 0.000797929 30525 -2.518074  0.0118
    ## tempave_metab.sc                                 -0.00793053 0.001176921 30525 -6.738371  0.0000
    ## seas.sc                                          -0.00292862 0.000872141 30525 -3.357970  0.0008
    ## microclim.sc                                     -0.00006305 0.000430068 30525 -0.146606  0.8834
    ## mass.sc                                          -0.00026278 0.000913556 30525 -0.287645  0.7736
    ## speed.sc                                          0.00054092 0.000852474 30525  0.634531  0.5257
    ## consumerfrac.sc                                   0.00077674 0.000438250 30525  1.772362  0.0763
    ## nspp.sc                                           0.00068961 0.000715905 30525  0.963272  0.3354
    ## npp.sc                                           -0.00057765 0.000621356 30525 -0.929661  0.3526
    ## veg.sc                                            0.00070219 0.001355901 30525  0.517877  0.6045
    ## duration.sc                                      -0.00219805 0.000662421 30525 -3.318210  0.0009
    ## temptrend_abs.sc:REALMMarine                     -0.01997624 0.016177705 30525 -1.234801  0.2169
    ## temptrend_abs.sc:REALMTerrestrial                -0.01207568 0.016722987 30525 -0.722101  0.4702
    ## temptrend_abs.sc:tsign1                          -0.00658239 0.002153176 30525 -3.057061  0.0022
    ## temptrend_abs.sc:tempave_metab.sc                -0.00147870 0.003129498 30525 -0.472504  0.6366
    ## temptrend_abs.sc:seas.sc                          0.00091985 0.001836646 30525  0.500832  0.6165
    ## temptrend_abs.sc:microclim.sc                     0.00109056 0.001152427 30525  0.946314  0.3440
    ## temptrend_abs.sc:mass.sc                         -0.00179471 0.001782155 30525 -1.007044  0.3139
    ## temptrend_abs.sc:speed.sc                         0.00173048 0.002043513 30525  0.846814  0.3971
    ## temptrend_abs.sc:consumerfrac.sc                  0.00562026 0.000913647 30525  6.151461  0.0000
    ## temptrend_abs.sc:nspp.sc                         -0.00002630 0.001648425 30525 -0.015953  0.9873
    ## tsign-1:thermal_bias.sc                          -0.00143407 0.000878404 30525 -1.632585  0.1026
    ## tsign1:thermal_bias.sc                           -0.00058676 0.000641648 30525 -0.914457  0.3605
    ## temptrend_abs.sc:npp.sc                           0.00454029 0.001643264 30525  2.762971  0.0057
    ## temptrend_abs.sc:veg.sc                          -0.00398901 0.002350458 30525 -1.697121  0.0897
    ## temptrend_abs.sc:duration.sc                      0.00110693 0.001230765 30525  0.899383  0.3685
    ## human_bowler.sc:REALM2TerrFresh                   0.00190747 0.001070478 30525  1.781890  0.0748
    ## human_bowler.sc:REALM2Marine                      0.00050202 0.000504744 30525  0.994610  0.3199
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.00201115 0.001643566 30525 -1.223648  0.2211
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.00117701 0.001460915 30525  0.805667  0.4204
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.00283530 0.001544340 30525 -1.835933  0.0664
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.00328081 0.001280872 30525 -2.561392  0.0104
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.467                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.952  0.447                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.749  0.329  0.696                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.064  0.055 -0.003  0.011                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.085 -0.015 -0.080 -0.230  0.026                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.095  0.043  0.151 -0.094 -0.088  0.073                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.068  0.050  0.080  0.026  0.014 -0.213  0.113                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.104 -0.010 -0.070 -0.039 -0.007  0.078  0.111  0.020                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.042 -0.053 -0.020 -0.072 -0.058 -0.086 -0.033  0.058 -0.581                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.009  0.023  0.005  0.077  0.034 -0.116 -0.047  0.014  0.001 -0.079                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.044 -0.056 -0.076 -0.076  0.048 -0.197  0.031 -0.076  0.016  0.142  0.094                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.022  0.004 -0.043  0.063  0.012 -0.007 -0.201 -0.169 -0.023  0.149 -0.036 -0.141                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.495  0.271  0.517  0.021 -0.016  0.010  0.091  0.001  0.012 -0.016 -0.003  0.014 -0.176                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.029  0.016  0.040  0.029 -0.167  0.095 -0.036 -0.029 -0.057  0.024  0.006 -0.270  0.035 -0.003                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.454 -0.951 -0.459 -0.308 -0.014  0.023 -0.067 -0.058 -0.004  0.036 -0.016  0.073  0.026 -0.292 -0.023                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.293 -0.713 -0.269 -0.474 -0.006  0.150  0.135  0.014 -0.030  0.087 -0.058  0.078 -0.064  0.037 -0.022  0.663                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.041 -0.132 -0.008  0.002 -0.486 -0.043  0.014 -0.012  0.018  0.016 -0.012 -0.009 -0.039 -0.005  0.047  0.038      0.016                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                 0.003  0.021 -0.007  0.130 -0.066 -0.460 -0.078  0.067 -0.050  0.136  0.080  0.097 -0.034 -0.012 -0.020 -0.022     -0.275      0.019                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.044 -0.075 -0.070  0.116  0.059 -0.045 -0.637 -0.106 -0.073  0.026  0.015  0.014  0.129 -0.092 -0.012  0.125     -0.226     -0.042  0.042                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.031 -0.059 -0.038 -0.006 -0.028  0.124 -0.062 -0.708  0.006 -0.043 -0.007  0.047  0.142  0.016  0.013  0.073     -0.049      0.023 -0.028   0.136                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.021  0.027  0.015 -0.022  0.009 -0.065 -0.066  0.014 -0.567  0.335  0.024 -0.030 -0.023 -0.007  0.071  0.017      0.051     -0.016  0.143   0.050            -0.051                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.034  0.086  0.033  0.076  0.043  0.112  0.048 -0.045  0.272 -0.560  0.066 -0.095 -0.054 -0.007 -0.020 -0.046     -0.113     -0.001 -0.333   0.013             0.055            -0.508                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.005 -0.033 -0.007 -0.054 -0.016  0.085  0.029  0.018  0.003  0.068 -0.718 -0.069  0.013  0.000 -0.022  0.021      0.091      0.034 -0.100  -0.060            -0.039            -0.058            -0.108                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.044  0.046  0.054  0.069 -0.019  0.103  0.000  0.046 -0.031 -0.098 -0.050 -0.621  0.072 -0.019  0.172 -0.086     -0.088      0.012 -0.139  -0.030            -0.032             0.048             0.108             0.104                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.058 -0.016 -0.056 -0.058 -0.130  0.241 -0.232 -0.062 -0.011 -0.030 -0.007 -0.091 -0.010  0.001  0.087  0.023      0.020     -0.017 -0.052   0.141             0.041             0.012             0.011             0.001            0.030                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.091 -0.034 -0.100 -0.107  0.071  0.381 -0.333 -0.167 -0.020 -0.033 -0.017 -0.113 -0.079 -0.003  0.101  0.028      0.054      0.026 -0.123   0.164             0.082             0.016             0.026             0.014            0.038             0.414                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.015  0.021 -0.004 -0.062 -0.011  0.004  0.093  0.141 -0.009 -0.069  0.007  0.067 -0.652  0.108 -0.025 -0.065      0.100      0.016  0.077  -0.226            -0.318             0.061             0.022            -0.002           -0.070            -0.013  0.048                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.284 -0.424 -0.299  0.034 -0.010 -0.027 -0.109  0.009 -0.012  0.011  0.003 -0.017  0.174 -0.610  0.000  0.462     -0.077      0.052  0.005   0.187             0.011             0.005             0.015            -0.006            0.010             0.014 -0.016 -0.263                                                                                                  
    ## temptrend_abs.sc:duration.sc                      0.005  0.044 -0.042 -0.051  0.042  0.066 -0.037 -0.010  0.064 -0.038  0.006  0.219 -0.022  0.010 -0.467  0.037      0.052     -0.129  0.013   0.010             0.023            -0.089             0.007             0.040           -0.338            -0.012 -0.064  0.031            -0.004                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.062  0.058  0.065 -0.029  0.013  0.033  0.037  0.067  0.039 -0.043 -0.024 -0.006 -0.108  0.141 -0.014 -0.067      0.052     -0.011 -0.011  -0.095             0.002            -0.002             0.008             0.015            0.020             0.047  0.085  0.049            -0.194            0.018                                                               
    ## human_bowler.sc:REALM2Marine                      0.037 -0.024 -0.044 -0.033 -0.004  0.115 -0.194 -0.002 -0.031  0.058 -0.074 -0.053 -0.150  0.017  0.009  0.030      0.021      0.005 -0.056   0.150            -0.009             0.036            -0.064             0.097            0.046             0.127  0.167  0.099            -0.010            0.042            0.018                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.014  0.014  0.019  0.035 -0.066 -0.107  0.178  0.064  0.006  0.042  0.001  0.026 -0.042  0.008 -0.052 -0.048     -0.061      0.267  0.278  -0.280            -0.133             0.034            -0.128             0.030           -0.014            -0.442 -0.278  0.102            -0.034            0.031           -0.037      -0.083                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.018  0.045  0.018  0.055 -0.012 -0.154  0.163  0.080  0.002  0.068 -0.002  0.046  0.034 -0.012 -0.051 -0.046     -0.113     -0.059  0.430  -0.243            -0.204             0.044            -0.208             0.024           -0.065            -0.226 -0.578  0.020             0.012            0.051           -0.046      -0.127       0.502                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.111 -0.211 -0.122  0.053 -0.011 -0.050 -0.129 -0.041 -0.018  0.035  0.009  0.021  0.106 -0.265  0.015  0.231     -0.106      0.046  0.094   0.228             0.044             0.003            -0.043            -0.015           -0.039            -0.038 -0.063 -0.158             0.511           -0.016           -0.528      -0.001       0.102  0.105               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.022  0.033  0.030  0.029 -0.006 -0.084  0.141  0.000  0.039 -0.060  0.153  0.040  0.102 -0.009  0.003 -0.038     -0.031     -0.002  0.104  -0.174            -0.052            -0.030             0.061            -0.181           -0.054            -0.083 -0.129 -0.130             0.021           -0.037           -0.011      -0.749       0.093  0.148  0.001        
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -5.4616815 -0.2353987 -0.0133081  0.2335245  6.1296080 
    ## 
    ## Number of Observations: 30756
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30756

``` r
rsquared(modTfullJturem0)
```

    ##       Response   family     link method     Marginal  Conditional
    ## 1 Jtutrendrem0 gaussian identity   none 1.276013e-05 0.0001376072

``` r
rsquared(modTfullJbetarem0)
```

    ##         Response   family     link method     Marginal  Conditional
    ## 1 Jbetatrendrem0 gaussian identity   none 1.899003e-05 0.0001459506

``` r
rsquared(modTfullHornrem0)
```

    ##        Response   family     link method     Marginal  Conditional
    ## 1 Horntrendrem0 gaussian identity   none 2.481326e-05 0.0001856183

#### Summary trendz

``` r
summary(modTfullJtuz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   125009.4 125354.5 -62462.72
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      3.6897259 (Intr)
    ## temptrend_abs.sc 0.7777827 -0.001
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    2.133608  1755.81
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -6.113408 
    ## Fixed effects: Jtutrendz ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF   t-value p-value
    ## (Intercept)                                       0.3199969 1.4652417 27194  0.218392  0.8271
    ## temptrend_abs.sc                                  0.3853620 0.5762759 27194  0.668711  0.5037
    ## REALMMarine                                       0.6113024 1.5422166   138  0.396379  0.6924
    ## REALMTerrestrial                                  1.7564164 1.5402026   138  1.140380  0.2561
    ## tsign1                                           -0.1488229 0.0451003 27194 -3.299817  0.0010
    ## tempave_metab.sc                                 -0.5076519 0.0817706 27194 -6.208245  0.0000
    ## seas.sc                                          -0.1509686 0.0540862 27194 -2.791261  0.0053
    ## microclim.sc                                     -0.0244059 0.0252693 27194 -0.965831  0.3341
    ## mass.sc                                          -0.0623083 0.0532879 27194 -1.169277  0.2423
    ## speed.sc                                          0.1428145 0.0547117 27194  2.610311  0.0091
    ## consumerfrac.sc                                   0.0214851 0.0450895 27194  0.476499  0.6337
    ## nspp.sc                                           0.4138425 0.0412309 27194 10.037204  0.0000
    ## npp.sc                                            0.1536723 0.0379173 27194  4.052831  0.0001
    ## veg.sc                                           -0.0625537 0.1168322 27194 -0.535415  0.5924
    ## duration.sc                                       0.1454903 0.0320363 27194  4.541424  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.5478625 0.5963821 27194 -0.918643  0.3583
    ## temptrend_abs.sc:REALMTerrestrial                -0.5530646 0.6880260 27194 -0.803843  0.4215
    ## temptrend_abs.sc:tsign1                           0.1513398 0.0604558 27194  2.503313  0.0123
    ## temptrend_abs.sc:tempave_metab.sc                 0.1025680 0.1102800 27194  0.930069  0.3523
    ## temptrend_abs.sc:seas.sc                          0.0246616 0.0531181 27194  0.464279  0.6425
    ## temptrend_abs.sc:microclim.sc                    -0.0313677 0.0345897 27194 -0.906850  0.3645
    ## temptrend_abs.sc:mass.sc                          0.0184161 0.0562016 27194  0.327678  0.7432
    ## temptrend_abs.sc:speed.sc                        -0.1129194 0.0635237 27194 -1.777594  0.0755
    ## temptrend_abs.sc:consumerfrac.sc                  0.1963727 0.0792475 27194  2.477966  0.0132
    ## temptrend_abs.sc:nspp.sc                         -0.0776294 0.0475871 27194 -1.631312  0.1028
    ## tsign-1:thermal_bias.sc                          -0.1282425 0.0478992 27194 -2.677344  0.0074
    ## tsign1:thermal_bias.sc                           -0.1040113 0.0398225 27194 -2.611874  0.0090
    ## temptrend_abs.sc:npp.sc                          -0.0354106 0.0526668 27194 -0.672351  0.5014
    ## temptrend_abs.sc:veg.sc                          -0.0192787 0.0651040 27194 -0.296121  0.7671
    ## temptrend_abs.sc:duration.sc                     -0.0152259 0.0448603 27194 -0.339408  0.7343
    ## human_bowler.sc:REALM2TerrFresh                   0.0997487 0.0963692 27194  1.035068  0.3006
    ## human_bowler.sc:REALM2Marine                      0.0376089 0.0302897 27194  1.241640  0.2144
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0569854 0.0465103 27194  1.225220  0.2205
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.0166564 0.0442543 27194 -0.376379  0.7066
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0397563 0.0400355 27194 -0.993027  0.3207
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.1156316 0.0508889 27194 -2.272235  0.0231
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.192                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.955  0.185                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.896  0.147  0.850                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.017  0.039 -0.006  0.009                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.047 -0.031 -0.038 -0.088 -0.050                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.031  0.013  0.058 -0.053 -0.106  0.077                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.034  0.042  0.040  0.009 -0.020 -0.174  0.125                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.017  0.023 -0.006  0.008  0.013  0.126  0.001  0.011                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.036 -0.054 -0.019 -0.041 -0.059  0.051  0.019  0.047 -0.504                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.011  0.022  0.010  0.035  0.011 -0.066 -0.026  0.035  0.012 -0.019                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.013 -0.037 -0.027 -0.031  0.054 -0.199 -0.003 -0.033 -0.067  0.138  0.075                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.014  0.011 -0.024  0.037  0.001 -0.077 -0.345 -0.272 -0.014  0.085 -0.039 -0.146                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.243  0.170  0.258  0.003 -0.026  0.004  0.133  0.044  0.015 -0.005  0.012  0.029 -0.164                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.001  0.004  0.012  0.025 -0.187  0.040 -0.062 -0.021 -0.013  0.008 -0.015 -0.307  0.032 -0.011                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.189 -0.954 -0.193 -0.140 -0.014  0.035 -0.038 -0.050 -0.037  0.033 -0.018  0.051  0.010 -0.191 -0.020                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.125 -0.767 -0.114 -0.197 -0.017  0.147  0.125  0.003 -0.076  0.110 -0.022  0.064 -0.065  0.037 -0.042  0.724                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.005 -0.096  0.006 -0.003 -0.510  0.015  0.007 -0.011  0.003  0.029 -0.003 -0.006 -0.035  0.032  0.093  0.036      0.033                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.021  0.065  0.018  0.068 -0.027 -0.420 -0.042  0.086 -0.032  0.042  0.050  0.065 -0.039  0.006  0.025 -0.062     -0.334      0.008                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.029 -0.042 -0.043  0.035  0.051 -0.003 -0.527 -0.129 -0.005 -0.020  0.001  0.003  0.149 -0.131 -0.001  0.085     -0.234     -0.030  0.010                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.021 -0.046 -0.024  0.000 -0.003  0.090 -0.093 -0.557  0.002 -0.022 -0.025  0.016  0.142 -0.028 -0.007  0.052     -0.064      0.039 -0.049   0.234                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                          0.007 -0.027 -0.008 -0.026 -0.002 -0.050 -0.007  0.003 -0.571  0.307  0.011  0.024 -0.004 -0.021  0.032  0.059      0.114     -0.019  0.084  -0.012            -0.058                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.013  0.093  0.008  0.044  0.035  0.040 -0.005 -0.044  0.284 -0.525  0.008 -0.082 -0.022 -0.019 -0.017 -0.051     -0.194     -0.026 -0.141   0.056             0.076            -0.440                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.005 -0.042 -0.002  0.014 -0.005  0.047  0.001 -0.020  0.013 -0.003 -0.499 -0.035  0.029 -0.006  0.007  0.037      0.041      0.002 -0.079   0.007             0.037            -0.045            -0.032                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.007  0.068  0.009  0.026 -0.003  0.078 -0.005  0.019  0.019 -0.077 -0.041 -0.511  0.070 -0.052  0.132 -0.082     -0.108     -0.022 -0.086  -0.026             0.009            -0.023             0.130             0.048                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.024 -0.019 -0.023 -0.030 -0.093  0.415 -0.151 -0.091  0.057 -0.032 -0.010 -0.092 -0.061  0.003  0.061  0.023      0.026     -0.014 -0.101   0.141             0.050            -0.015             0.026             0.010            0.027                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.039 -0.041 -0.039 -0.046  0.048  0.569 -0.192 -0.218  0.070 -0.028 -0.015 -0.089 -0.109 -0.013  0.067  0.031      0.050      0.063 -0.164   0.155             0.107            -0.012             0.039             0.015            0.031             0.503                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.000 -0.021  0.004 -0.022 -0.022  0.012  0.112  0.132 -0.004 -0.022  0.033  0.059 -0.545  0.109 -0.002 -0.015      0.094      0.049  0.132  -0.197            -0.294             0.048            -0.038            -0.050           -0.075             0.012  0.025                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.129 -0.263 -0.139  0.011  0.003 -0.009 -0.157 -0.017 -0.021 -0.008 -0.012 -0.035  0.217 -0.561 -0.016  0.310     -0.099      0.020 -0.032   0.325             0.019             0.037             0.066             0.011            0.040             0.025  0.010 -0.377                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.012  0.095 -0.004 -0.016  0.016  0.104 -0.017  0.002  0.013 -0.032  0.010  0.137 -0.008  0.012 -0.402  0.020      0.043     -0.036 -0.045   0.029            -0.015            -0.030             0.010            -0.004           -0.172             0.011 -0.048  0.009             0.023                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.025  0.007  0.030 -0.035 -0.003  0.035  0.097  0.018 -0.003  0.003 -0.002 -0.014 -0.074  0.176 -0.025 -0.010      0.046      0.002 -0.011  -0.075             0.013             0.004             0.011             0.004            0.044             0.019  0.027  0.024            -0.116            0.027                                                               
    ## human_bowler.sc:REALM2Marine                      0.015 -0.021 -0.019 -0.008 -0.008  0.100 -0.195 -0.087 -0.023  0.045 -0.026 -0.047 -0.072 -0.010 -0.022  0.034      0.023      0.035 -0.091   0.156             0.048             0.044            -0.052             0.018            0.043             0.147  0.149  0.050             0.018            0.058           -0.005                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.014  0.033  0.016  0.022 -0.052 -0.186  0.156  0.057 -0.023  0.049  0.011  0.032 -0.008  0.008 -0.022 -0.050     -0.077      0.252  0.341  -0.313            -0.117             0.049            -0.098            -0.025           -0.014            -0.522 -0.325  0.082            -0.085            0.016            0.010      -0.131                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.025  0.066  0.025  0.032  0.038 -0.240  0.151  0.102 -0.014  0.043  0.014  0.032  0.025  0.025 -0.048 -0.055     -0.130     -0.169  0.423  -0.247            -0.177             0.034            -0.115            -0.027           -0.033            -0.295 -0.558  0.055            -0.035            0.036           -0.012      -0.151       0.603                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.035 -0.114 -0.041  0.033 -0.001 -0.046 -0.163 -0.029  0.011  0.005 -0.001  0.021  0.114 -0.207  0.009  0.136     -0.137      0.028  0.066   0.346             0.039            -0.033            -0.013             0.004           -0.067            -0.018 -0.033 -0.170             0.523           -0.011           -0.579       0.018       0.021  0.062               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.013  0.035  0.014  0.014  0.026 -0.077  0.087  0.038  0.034 -0.048  0.020  0.031  0.043  0.013  0.020 -0.053     -0.049     -0.043  0.183  -0.215            -0.101            -0.068             0.064            -0.040           -0.056            -0.091 -0.089 -0.079            -0.017           -0.064            0.001      -0.641       0.168  0.192 -0.028        
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.362924e+01 -1.127015e-02 -1.727242e-05  7.828064e-03  1.425198e+01 
    ## 
    ## Number of Observations: 27368
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    141                  27368

``` r
summary(modTfullJbetaz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   91887.47 92232.67 -45901.73
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev   Corr  
    ## (Intercept)      2.841575 (Intr)
    ## temptrend_abs.sc 2.152057 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)    Residual
    ## StdDev:    1.140902 65195843421
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -22.46541 
    ## Fixed effects: Jbetatrendz ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF   t-value p-value
    ## (Intercept)                                       1.0419829 1.1029987 27282  0.944682  0.3448
    ## temptrend_abs.sc                                  0.5032103 1.1023417 27282  0.456492  0.6480
    ## REALMMarine                                      -0.3887028 1.1604275   141 -0.334965  0.7381
    ## REALMTerrestrial                                  1.4418607 1.1772037   141  1.224818  0.2227
    ## tsign1                                           -0.0244646 0.0243348 27282 -1.005335  0.3147
    ## tempave_metab.sc                                 -0.3152433 0.0474347 27282 -6.645838  0.0000
    ## seas.sc                                          -0.0939898 0.0293734 27282 -3.199822  0.0014
    ## microclim.sc                                      0.0009104 0.0136800 27282  0.066548  0.9469
    ## mass.sc                                          -0.0013881 0.0293036 27282 -0.047371  0.9622
    ## speed.sc                                          0.0440341 0.0303868 27282  1.449118  0.1473
    ## consumerfrac.sc                                   0.0325658 0.0248230 27282  1.311920  0.1896
    ## nspp.sc                                           0.2719835 0.0226205 27282 12.023779  0.0000
    ## npp.sc                                            0.0339400 0.0205831 27282  1.648924  0.0992
    ## veg.sc                                           -0.0510387 0.0620195 27282 -0.822946  0.4105
    ## duration.sc                                       0.0894864 0.0175522 27282  5.098306  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.9342878 1.1470671 27282 -0.814501  0.4154
    ## temptrend_abs.sc:REALMTerrestrial                -0.8527722 1.2197702 27282 -0.699125  0.4845
    ## temptrend_abs.sc:tsign1                           0.0377858 0.0333519 27282  1.132944  0.2572
    ## temptrend_abs.sc:tempave_metab.sc                 0.1502009 0.0799273 27282  1.879220  0.0602
    ## temptrend_abs.sc:seas.sc                         -0.0044728 0.0305493 27282 -0.146412  0.8836
    ## temptrend_abs.sc:microclim.sc                    -0.0311101 0.0193615 27282 -1.606802  0.1081
    ## temptrend_abs.sc:mass.sc                         -0.0235335 0.0324570 27282 -0.725065  0.4684
    ## temptrend_abs.sc:speed.sc                         0.0227957 0.0386969 27282  0.589085  0.5558
    ## temptrend_abs.sc:consumerfrac.sc                 -0.0214257 0.0476355 27282 -0.449784  0.6529
    ## temptrend_abs.sc:nspp.sc                          0.0021793 0.0275420 27282  0.079128  0.9369
    ## tsign-1:thermal_bias.sc                          -0.0690733 0.0262506 27282 -2.631307  0.0085
    ## tsign1:thermal_bias.sc                           -0.0165024 0.0219485 27282 -0.751869  0.4521
    ## temptrend_abs.sc:npp.sc                           0.0326444 0.0297971 27282  1.095557  0.2733
    ## temptrend_abs.sc:veg.sc                          -0.0567352 0.0361775 27282 -1.568247  0.1168
    ## temptrend_abs.sc:duration.sc                      0.0193911 0.0247685 27282  0.782896  0.4337
    ## human_bowler.sc:REALM2TerrFresh                   0.0742050 0.0500691 27282  1.482051  0.1383
    ## human_bowler.sc:REALM2Marine                     -0.0088858 0.0164354 27282 -0.540654  0.5888
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0239330 0.0278760 27282  0.858552  0.3906
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.0175372 0.0267885 27282 -0.654654  0.5127
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0523487 0.0220535 27282 -2.373711  0.0176
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0009293 0.0283414 27282 -0.032789  0.9738
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.275                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.953  0.262                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.909  0.250  0.864                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.013  0.011 -0.005  0.006                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.045 -0.020 -0.037 -0.075 -0.052                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.034  0.009  0.053 -0.022 -0.108  0.055                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.027  0.016  0.031  0.007 -0.021 -0.184  0.138                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.026  0.003 -0.016 -0.007  0.010  0.116 -0.006  0.009                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.028 -0.020 -0.015 -0.033 -0.062  0.112  0.022  0.039 -0.483                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.007  0.006  0.007  0.026  0.011 -0.053 -0.019  0.035  0.006  0.001                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.002 -0.008 -0.013 -0.014  0.059 -0.205 -0.001 -0.029 -0.075  0.126  0.077                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.017  0.001 -0.025  0.020  0.003 -0.068 -0.342 -0.271 -0.009  0.085 -0.042 -0.149                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.178  0.054  0.189  0.012 -0.025  0.003  0.128  0.054  0.016 -0.004  0.010  0.033 -0.173                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.001  0.001  0.010  0.015 -0.188  0.038 -0.062 -0.023  0.002  0.007 -0.020 -0.313  0.036 -0.010                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.266 -0.960 -0.265 -0.240 -0.003  0.020 -0.019 -0.019 -0.009  0.011 -0.005  0.014  0.007 -0.061 -0.006                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.241 -0.897 -0.229 -0.291 -0.006  0.065  0.028 -0.002 -0.020  0.041 -0.012  0.016 -0.017  0.004 -0.009  0.861                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.005 -0.028  0.004 -0.003 -0.509  0.018  0.012 -0.008  0.004  0.031 -0.002 -0.012 -0.040  0.030  0.092  0.009      0.011                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.020  0.029  0.019  0.053 -0.013 -0.517 -0.013  0.099 -0.019 -0.063  0.027  0.075 -0.038  0.002  0.024 -0.030     -0.113      0.005                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.029 -0.020 -0.038  0.014  0.052  0.031 -0.539 -0.146  0.002 -0.010  0.002  0.003  0.155 -0.135  0.002  0.036     -0.056     -0.039 -0.050                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.016 -0.015 -0.018  0.001 -0.002  0.087 -0.104 -0.564  0.009 -0.022 -0.026  0.010  0.151 -0.033 -0.001  0.017     -0.020      0.032 -0.036   0.248                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.005 -0.001  0.001 -0.007  0.003 -0.055  0.004  0.008 -0.584  0.273  0.009  0.037 -0.013 -0.018  0.007  0.012      0.025     -0.017  0.070  -0.030            -0.070                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.013  0.034  0.006  0.031  0.039 -0.063 -0.004 -0.034  0.236 -0.551 -0.016 -0.058 -0.028 -0.024 -0.020 -0.015     -0.064     -0.024  0.076   0.036             0.077            -0.340                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.003 -0.010 -0.005  0.003 -0.008  0.027 -0.002 -0.018  0.014 -0.023 -0.522 -0.037  0.030 -0.009  0.019  0.011      0.028      0.001 -0.011   0.004             0.040            -0.041             0.023                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                          0.004  0.012 -0.001  0.008 -0.011  0.102 -0.006  0.013  0.034 -0.060 -0.045 -0.530  0.077 -0.056  0.150 -0.019     -0.021     -0.008 -0.112  -0.026             0.019            -0.047             0.085             0.049                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.025 -0.012 -0.023 -0.029 -0.096  0.437 -0.164 -0.103  0.059 -0.004 -0.009 -0.102 -0.056  0.002  0.064  0.014      0.018     -0.013 -0.164   0.157             0.058            -0.023            -0.014             0.008            0.039                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.036 -0.020 -0.036 -0.043  0.046  0.590 -0.209 -0.229  0.069  0.009 -0.012 -0.101 -0.101 -0.012  0.067  0.017      0.029      0.061 -0.234   0.173             0.111            -0.017            -0.013             0.011            0.048             0.521                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.006  0.000  0.009 -0.012 -0.025  0.006  0.116  0.137 -0.009 -0.029  0.035  0.065 -0.558  0.117 -0.005 -0.011      0.024      0.057  0.114  -0.206            -0.305             0.057            -0.021            -0.051           -0.086             0.008  0.021                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.095 -0.080 -0.103  0.002  0.004 -0.010 -0.163 -0.026 -0.024 -0.011 -0.013 -0.034  0.227 -0.560 -0.018  0.096     -0.023      0.019 -0.025   0.332             0.031             0.040             0.071             0.015            0.038             0.027  0.010 -0.383                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.005  0.023 -0.007 -0.012  0.020  0.101 -0.015  0.003 -0.002 -0.029  0.014  0.150 -0.013  0.010 -0.415  0.008      0.011     -0.037 -0.045   0.027            -0.020            -0.003             0.013            -0.013           -0.192             0.010 -0.044  0.014             0.027                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.007 -0.001  0.011 -0.026 -0.001  0.032  0.086  0.031 -0.003  0.002 -0.003 -0.015 -0.078  0.151 -0.025  0.000      0.010      0.001 -0.011  -0.077             0.010             0.006             0.009             0.001            0.048             0.019  0.029  0.023            -0.110            0.024                                                               
    ## human_bowler.sc:REALM2Marine                      0.014 -0.008 -0.016 -0.008 -0.007  0.101 -0.193 -0.090 -0.020  0.051 -0.022 -0.050 -0.066 -0.009 -0.024  0.011      0.007      0.030 -0.081   0.155             0.054             0.038            -0.061             0.009            0.045             0.149  0.150  0.036             0.021            0.060           -0.003                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.016  0.019  0.018  0.024 -0.045 -0.286  0.167  0.082 -0.020 -0.011  0.006  0.042 -0.017  0.008 -0.021 -0.027     -0.043      0.236  0.480  -0.333            -0.122             0.054             0.015            -0.002           -0.035            -0.531 -0.364  0.089            -0.086            0.013            0.009      -0.131                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.023  0.028  0.023  0.033  0.037 -0.342  0.159  0.122 -0.011 -0.019  0.008  0.049  0.016  0.020 -0.042 -0.028     -0.060     -0.148  0.557  -0.266            -0.171             0.040             0.005            -0.001           -0.064            -0.327 -0.573  0.058            -0.033            0.029           -0.013      -0.145       0.666                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.024 -0.036 -0.028  0.021 -0.001 -0.050 -0.167 -0.039  0.014  0.004  0.000  0.026  0.122 -0.202  0.009  0.043     -0.038      0.025  0.070   0.358             0.051            -0.040            -0.007             0.006           -0.074            -0.016 -0.034 -0.177             0.530           -0.010           -0.566       0.020       0.020  0.066               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.010  0.012  0.011  0.009  0.024 -0.089  0.090  0.044  0.032 -0.061  0.014  0.036  0.031  0.014  0.021 -0.015     -0.014     -0.036  0.174  -0.216            -0.111            -0.058             0.082            -0.023           -0.061            -0.098 -0.096 -0.053            -0.024           -0.064            0.001      -0.649       0.177  0.194 -0.032        
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.685874e+10 -3.666735e-04 -8.949927e-09  1.491834e-04  8.696827e+07 
    ## 
    ## Number of Observations: 27459
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    144                  27459

``` r
summary(modTfullHornz)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC   logLik
    ##   145200.6 145545.7 -72558.3
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev   Corr  
    ## (Intercept)      3.091432 (Intr)
    ## temptrend_abs.sc 4.013309 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)     Residual
    ## StdDev:    3.317881 9.261969e+14
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -31.07832 
    ## Fixed effects: Horntrendz ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                      Value Std.Error    DF   t-value p-value
    ## (Intercept)                                       1.166602 1.5202988 27198  0.767351  0.4429
    ## temptrend_abs.sc                                  2.704641 2.0183161 27198  1.340048  0.1802
    ## REALMMarine                                       0.312947 1.5973687   135  0.195914  0.8450
    ## REALMTerrestrial                                  2.572130 1.5758759   135  1.632191  0.1050
    ## tsign1                                            0.030679 0.0654226 27198  0.468932  0.6391
    ## tempave_metab.sc                                 -0.864349 0.1199135 27198 -7.208107  0.0000
    ## seas.sc                                          -0.285819 0.0764780 27198 -3.737270  0.0002
    ## microclim.sc                                      0.008775 0.0368636 27198  0.238039  0.8119
    ## mass.sc                                          -0.012565 0.0747057 27198 -0.168197  0.8664
    ## speed.sc                                          0.033678 0.0783876 27198  0.429637  0.6675
    ## consumerfrac.sc                                   0.272184 0.0671965 27198  4.050572  0.0001
    ## nspp.sc                                           0.559956 0.0589826 27198  9.493583  0.0000
    ## npp.sc                                            0.008923 0.0547109 27198  0.163093  0.8704
    ## veg.sc                                           -0.183982 0.1692921 27198 -1.086773  0.2771
    ## duration.sc                                       0.109692 0.0449683 27198  2.439307  0.0147
    ## temptrend_abs.sc:REALMMarine                     -3.257722 2.1066393 27198 -1.546407  0.1220
    ## temptrend_abs.sc:REALMTerrestrial                -1.732626 2.2477153 27198 -0.770839  0.4408
    ## temptrend_abs.sc:tsign1                          -0.092820 0.0800466 27198 -1.159571  0.2462
    ## temptrend_abs.sc:tempave_metab.sc                 0.042168 0.1870595 27198  0.225427  0.8216
    ## temptrend_abs.sc:seas.sc                         -0.000484 0.0701070 27198 -0.006898  0.9945
    ## temptrend_abs.sc:microclim.sc                    -0.002262 0.0457535 27198 -0.049442  0.9606
    ## temptrend_abs.sc:mass.sc                         -0.192086 0.0765253 27198 -2.510095  0.0121
    ## temptrend_abs.sc:speed.sc                         0.143217 0.0897832 27198  1.595139  0.1107
    ## temptrend_abs.sc:consumerfrac.sc                 -0.039956 0.1136541 27198 -0.351557  0.7252
    ## temptrend_abs.sc:nspp.sc                         -0.035299 0.0644717 27198 -0.547518  0.5840
    ## tsign-1:thermal_bias.sc                          -0.201247 0.0681957 27198 -2.951027  0.0032
    ## tsign1:thermal_bias.sc                           -0.088816 0.0567513 27198 -1.564996  0.1176
    ## temptrend_abs.sc:npp.sc                           0.138924 0.0704182 27198  1.972843  0.0485
    ## temptrend_abs.sc:veg.sc                          -0.155755 0.0851464 27198 -1.829265  0.0674
    ## temptrend_abs.sc:duration.sc                      0.080285 0.0617585 27198  1.299989  0.1936
    ## human_bowler.sc:REALM2TerrFresh                   0.102229 0.1382630 27198  0.739383  0.4597
    ## human_bowler.sc:REALM2Marine                     -0.053400 0.0441538 27198 -1.209408  0.2265
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.043521 0.0652948 27198 -0.666536  0.5051
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.051770 0.0631383 27198 -0.819948  0.4123
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.116546 0.0523174 27198 -2.227664  0.0259
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine     0.117681 0.0689859 27198  1.705875  0.0880
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.296                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.962  0.283                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.851  0.264  0.808                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.024  0.014 -0.008  0.014                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.065 -0.021 -0.056 -0.130 -0.052                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.047  0.004  0.082 -0.069 -0.099  0.067                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.050  0.020  0.058  0.013 -0.021 -0.168  0.121                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.019  0.013 -0.006  0.012  0.008  0.119  0.005  0.010                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.051 -0.024 -0.027 -0.059 -0.057  0.048  0.003  0.049 -0.477                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.019  0.009  0.017  0.061  0.020 -0.071 -0.040  0.034  0.023 -0.038                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.022 -0.015 -0.039 -0.044  0.054 -0.189 -0.008 -0.030 -0.061  0.132  0.068                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.023  0.006 -0.037  0.050  0.000 -0.055 -0.349 -0.287 -0.014  0.092 -0.031 -0.145                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.349  0.072  0.371  0.008 -0.025  0.007  0.142  0.048  0.011  0.001  0.016  0.031 -0.172                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.001  0.000  0.017  0.034 -0.185  0.026 -0.060 -0.016 -0.010  0.011 -0.012 -0.305  0.038 -0.019                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.286 -0.956 -0.289 -0.252 -0.004  0.023 -0.015 -0.023 -0.018  0.014 -0.007  0.021  0.003 -0.080 -0.007                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.243 -0.887 -0.230 -0.319 -0.009  0.086  0.043 -0.004 -0.034  0.054 -0.017  0.023 -0.025  0.012 -0.012  0.848                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.007 -0.036  0.009 -0.008 -0.519  0.025  0.010 -0.008 -0.001  0.034 -0.008 -0.012 -0.033  0.037  0.097  0.012      0.019                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.028  0.031  0.026  0.085 -0.009 -0.508 -0.010  0.095 -0.020 -0.036  0.038  0.068 -0.043  0.005  0.023 -0.034     -0.138     -0.011                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.048 -0.023 -0.066  0.038  0.048  0.024 -0.528 -0.140 -0.004 -0.006  0.008  0.007  0.155 -0.141 -0.010  0.042     -0.075     -0.030 -0.048                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.032 -0.019 -0.036  0.001  0.002  0.081 -0.100 -0.553  0.006 -0.024 -0.023  0.016  0.156 -0.039 -0.014  0.021     -0.025      0.035 -0.042   0.257                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                          0.003 -0.006 -0.008 -0.021  0.001 -0.056 -0.003  0.004 -0.584  0.271  0.003  0.033 -0.007 -0.018  0.013  0.019      0.032     -0.013  0.071  -0.024            -0.063                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.022  0.042  0.009  0.049  0.040 -0.030  0.006 -0.036  0.243 -0.544 -0.003 -0.056 -0.031 -0.023 -0.025 -0.019     -0.076     -0.038  0.050   0.033             0.073            -0.349                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.008 -0.015 -0.008 -0.009 -0.012  0.038  0.005 -0.017  0.017 -0.017 -0.503 -0.034  0.024 -0.009  0.013  0.015      0.049      0.009 -0.042   0.004             0.032            -0.047             0.020                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                          0.004  0.015  0.001  0.020 -0.009  0.093 -0.003  0.019  0.033 -0.057 -0.043 -0.516  0.068 -0.060  0.136 -0.023     -0.027     -0.010 -0.106  -0.026             0.004            -0.052             0.073             0.057                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.033 -0.012 -0.032 -0.046 -0.082  0.408 -0.145 -0.092  0.058 -0.036 -0.012 -0.080 -0.052  0.005  0.048  0.013      0.021     -0.017 -0.142   0.154             0.051            -0.024             0.002             0.010            0.032                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.056 -0.024 -0.057 -0.072  0.044  0.561 -0.182 -0.219  0.067 -0.035 -0.018 -0.078 -0.095 -0.010  0.053  0.020      0.036      0.069 -0.216   0.172             0.112            -0.016             0.007             0.015            0.041             0.491                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.010 -0.001  0.016 -0.024 -0.022 -0.005  0.117  0.144 -0.006 -0.030  0.029  0.059 -0.547  0.113 -0.007 -0.014      0.029      0.046  0.129  -0.205            -0.308             0.047            -0.020            -0.046           -0.068             0.007  0.012                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.187 -0.101 -0.203  0.011  0.008 -0.006 -0.172 -0.028 -0.026 -0.013 -0.014 -0.035  0.227 -0.568 -0.015  0.120     -0.035      0.010 -0.032   0.354             0.037             0.048             0.074             0.016            0.043             0.026  0.019 -0.393                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.015  0.037 -0.007 -0.021  0.016  0.100 -0.020  0.000  0.003 -0.036  0.010  0.136 -0.014  0.017 -0.407  0.008      0.014     -0.024 -0.036   0.036            -0.017            -0.007             0.018            -0.008           -0.165             0.014 -0.046  0.010             0.021                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.034 -0.006  0.042 -0.056  0.000  0.033  0.098  0.021 -0.004  0.002 -0.008 -0.015 -0.077  0.196 -0.020  0.003      0.024      0.003 -0.008  -0.078             0.012             0.005             0.010             0.007            0.048             0.019  0.031  0.025            -0.125            0.020                                                               
    ## human_bowler.sc:REALM2Marine                      0.020 -0.008 -0.025 -0.007 -0.009  0.093 -0.199 -0.094 -0.024  0.044 -0.023 -0.051 -0.058 -0.011 -0.029  0.012      0.007      0.037 -0.074   0.152             0.060             0.045            -0.057             0.011            0.044             0.150  0.150  0.029             0.025            0.057           -0.005                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.020  0.018  0.024  0.034 -0.049 -0.266  0.161  0.071 -0.021  0.009  0.009  0.034 -0.016  0.007 -0.015 -0.028     -0.047      0.232  0.448  -0.336            -0.114             0.062            -0.008            -0.016           -0.035            -0.536 -0.357  0.086            -0.086            0.008            0.011      -0.134                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.038  0.033  0.039  0.050  0.047 -0.323  0.154  0.116 -0.008  0.000  0.014  0.040  0.011  0.028 -0.042 -0.032     -0.070     -0.180  0.526  -0.274            -0.168             0.040            -0.009            -0.021           -0.061            -0.314 -0.575  0.069            -0.051            0.028           -0.013      -0.148       0.659                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.050 -0.041 -0.059  0.046  0.000 -0.048 -0.168 -0.031  0.016  0.004  0.004  0.024  0.117 -0.212  0.006  0.049     -0.052      0.018  0.058   0.357             0.044            -0.040            -0.010             0.001           -0.071            -0.016 -0.031 -0.171             0.521           -0.008           -0.599       0.019       0.015  0.052               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.016  0.013  0.017  0.013  0.030 -0.081  0.084  0.047  0.036 -0.055  0.016  0.034  0.023  0.016  0.025 -0.019     -0.016     -0.050  0.169  -0.215            -0.117            -0.073             0.081            -0.025           -0.059            -0.095 -0.094 -0.043            -0.032           -0.063            0.002      -0.632       0.188  0.208 -0.030        
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.308624e+18 -2.794205e-03  0.000000e+00  1.982341e-04  8.063733e+19 
    ## 
    ## Number of Observations: 27369
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    138                  27369

``` r
rsquared(modTfullJtuz)
```

    ##    Response   family     link method     Marginal Conditional
    ## 1 Jtutrendz gaussian identity   none 1.772871e-07 6.25997e-06

``` r
rsquared(modTfullJbetaz)
```

    ##      Response   family     link method     Marginal  Conditional
    ## 1 Jbetatrendz gaussian identity   none 7.690367e-23 3.334853e-21

``` r
rsquared(modTfullHornz)
```

    ##     Response   family     link method     Marginal  Conditional
    ## 1 Horntrendz gaussian identity   none 1.849542e-30 4.400296e-29

#### Summary last

``` r
summary(modTfullJtulast)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   4054.981 4404.953 -1985.491
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.18708646 (Intr)
    ## temptrend_abs.sc 0.06635838 -0.382
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 0.0004019716 0.3392482
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.1637998 
    ## Fixed effects: Jtulast ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.3876638 0.07142010 30525   5.427936  0.0000
    ## temptrend_abs.sc                                  0.0086627 0.04969309 30525   0.174324  0.8616
    ## REALMMarine                                       0.0888679 0.07627144   195   1.165153  0.2454
    ## REALMTerrestrial                                 -0.1842583 0.07247114   195  -2.542506  0.0118
    ## tsign1                                           -0.0563815 0.00481553 30525 -11.708267  0.0000
    ## tempave_metab.sc                                  0.0517887 0.00799029 30525   6.481451  0.0000
    ## seas.sc                                           0.0221193 0.00503827 30525   4.390247  0.0000
    ## microclim.sc                                      0.0049208 0.00261373 30525   1.882672  0.0598
    ## mass.sc                                          -0.0138706 0.00518054 30525  -2.677450  0.0074
    ## speed.sc                                         -0.0173844 0.00542357 30525  -3.205348  0.0014
    ## consumerfrac.sc                                   0.0083317 0.00329292 30525   2.530181  0.0114
    ## nspp.sc                                           0.0726261 0.00412330 30525  17.613592  0.0000
    ## npp.sc                                            0.0077751 0.00379842 30525   2.046944  0.0407
    ## veg.sc                                           -0.0032246 0.00968396 30525  -0.332979  0.7392
    ## duration.sc                                       0.0339713 0.00323098 30525  10.514245  0.0000
    ## temptrend_abs.sc:REALMMarine                      0.0015758 0.05133145 30525   0.030699  0.9755
    ## temptrend_abs.sc:REALMTerrestrial                 0.0112985 0.05734980 30525   0.197010  0.8438
    ## temptrend_abs.sc:tsign1                           0.0270470 0.00597280 30525   4.528357  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                -0.0056422 0.01004077 30525  -0.561928  0.5742
    ## temptrend_abs.sc:seas.sc                          0.0005748 0.00496672 30525   0.115730  0.9079
    ## temptrend_abs.sc:microclim.sc                    -0.0027948 0.00338365 30525  -0.825975  0.4088
    ## temptrend_abs.sc:mass.sc                          0.0062388 0.00503800 30525   1.238345  0.2156
    ## temptrend_abs.sc:speed.sc                        -0.0007282 0.00598067 30525  -0.121761  0.9031
    ## temptrend_abs.sc:consumerfrac.sc                  0.0046806 0.00412248 30525   1.135379  0.2562
    ## temptrend_abs.sc:nspp.sc                         -0.0122923 0.00459616 30525  -2.674468  0.0075
    ## tsign-1:thermal_bias.sc                          -0.0160099 0.00504722 30525  -3.172016  0.0015
    ## tsign1:thermal_bias.sc                           -0.0168764 0.00398773 30525  -4.232072  0.0000
    ## temptrend_abs.sc:npp.sc                          -0.0044615 0.00505965 30525  -0.881775  0.3779
    ## temptrend_abs.sc:veg.sc                           0.0039608 0.00618271 30525   0.640626  0.5218
    ## temptrend_abs.sc:duration.sc                     -0.0003253 0.00437374 30525  -0.074383  0.9407
    ## human_bowler.sc:REALM2TerrFresh                  -0.0198129 0.00788688 30525  -2.512130  0.0120
    ## human_bowler.sc:REALM2Marine                     -0.0045040 0.00310654 30525  -1.449839  0.1471
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0041820 0.00459597 30525   0.909922  0.3629
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.0018265 0.00437041 30525  -0.417930  0.6760
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.0045586 0.00376562 30525   1.210593  0.2261
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0070298 0.00467967 30525  -1.502211  0.1331
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.442                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.946  0.416                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.809  0.362  0.753                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.045  0.043 -0.005  0.018                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.093 -0.037 -0.078 -0.169 -0.034                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.096  0.041  0.144 -0.043 -0.090  0.035                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.058  0.046  0.066  0.006 -0.015 -0.162  0.118                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.084 -0.001 -0.048 -0.042  0.018  0.119  0.026  0.003                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.050 -0.052 -0.015 -0.044 -0.060 -0.026  0.026  0.058 -0.524                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.010  0.022  0.007  0.068  0.046 -0.091 -0.057  0.015 -0.022 -0.055                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.009 -0.033 -0.034 -0.043  0.037 -0.188  0.024 -0.036 -0.048  0.139  0.081                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.052 -0.005 -0.068  0.041 -0.005 -0.025 -0.305 -0.272 -0.019  0.101 -0.024 -0.143                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.419  0.183  0.439  0.017 -0.017  0.002  0.122  0.042  0.010 -0.008 -0.005  0.022 -0.190                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.002  0.006  0.026  0.036 -0.197  0.035 -0.053 -0.011 -0.028  0.008 -0.007 -0.320  0.035 -0.006                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.431 -0.951 -0.443 -0.345 -0.014  0.045 -0.074 -0.053 -0.018  0.028 -0.016  0.049  0.028 -0.204 -0.025                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.312 -0.755 -0.288 -0.450 -0.015  0.148  0.088  0.004 -0.045  0.086 -0.043  0.064 -0.049  0.025 -0.031  0.714                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.021 -0.109  0.002 -0.011 -0.520  0.002  0.007 -0.009 -0.005  0.033 -0.012 -0.005 -0.024  0.026  0.099  0.038      0.036                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.031  0.061  0.025  0.114 -0.029 -0.456 -0.028  0.085 -0.048  0.063  0.070  0.068 -0.047  0.010  0.013 -0.057     -0.326      0.014                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.051 -0.081 -0.071  0.065  0.046  0.006 -0.532 -0.122 -0.015 -0.014  0.018 -0.003  0.136 -0.122 -0.016  0.122     -0.220     -0.028  0.007                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.031 -0.053 -0.035  0.004 -0.002  0.085 -0.089 -0.559  0.003 -0.028 -0.017  0.020  0.146 -0.023 -0.015  0.058     -0.058      0.038 -0.056   0.219                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.010 -0.014  0.000 -0.022 -0.010 -0.065 -0.014  0.007 -0.568  0.327  0.061  0.024 -0.010 -0.017  0.044  0.056      0.093     -0.007  0.120  -0.006            -0.053                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.020  0.091  0.011  0.057  0.040  0.063 -0.003 -0.045  0.297 -0.536  0.023 -0.079 -0.024 -0.012 -0.013 -0.054     -0.159     -0.037 -0.231   0.039             0.068            -0.507                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.000 -0.027 -0.002 -0.021 -0.011  0.060  0.028 -0.007  0.037  0.013 -0.632 -0.040  0.006  0.008 -0.013  0.026      0.070      0.015 -0.114  -0.032            -0.004            -0.117            -0.052                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.011  0.071  0.015  0.045  0.002  0.070 -0.015  0.029  0.017 -0.079 -0.039 -0.516  0.062 -0.043  0.131 -0.086     -0.131     -0.035 -0.064  -0.017            -0.009            -0.013             0.127             0.081                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.051 -0.026 -0.048 -0.057 -0.085  0.360 -0.170 -0.080  0.044 -0.060 -0.015 -0.082 -0.033  0.001  0.049  0.032      0.026     -0.020 -0.102   0.146             0.040            -0.017             0.029             0.010            0.028                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.081 -0.058 -0.083 -0.091  0.052  0.507 -0.244 -0.216  0.057 -0.074 -0.019 -0.092 -0.061 -0.014  0.060  0.048      0.054      0.063 -0.171   0.174             0.105            -0.018             0.045             0.014            0.032             0.450                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.002 -0.005  0.008 -0.038 -0.012 -0.006  0.101  0.136 -0.011 -0.028  0.006  0.055 -0.548  0.113 -0.009 -0.035      0.091      0.021  0.138  -0.196            -0.299             0.070            -0.033             0.003           -0.051             0.007  0.010                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.199 -0.304 -0.211  0.022 -0.001 -0.005 -0.154 -0.011 -0.016 -0.002  0.013 -0.032  0.218 -0.514 -0.017  0.353     -0.098      0.028 -0.038   0.332             0.010             0.026             0.048            -0.027            0.036             0.026  0.017 -0.382                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.019  0.114 -0.013 -0.024  0.010  0.101 -0.046  0.006  0.030 -0.033  0.029  0.136 -0.001  0.012 -0.400  0.019      0.038     -0.052 -0.031   0.030            -0.013            -0.053            -0.003             0.026           -0.176             0.013 -0.054 -0.004             0.014                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.066  0.019  0.068 -0.016  0.011  0.028  0.037  0.071  0.016 -0.018 -0.020 -0.007 -0.097  0.136 -0.011 -0.024      0.023     -0.003 -0.004  -0.063            -0.001             0.000             0.011             0.010            0.032             0.031  0.055  0.036            -0.110            0.001                                                               
    ## human_bowler.sc:REALM2Marine                      0.032 -0.032 -0.034 -0.018 -0.015  0.092 -0.195 -0.073 -0.027  0.041 -0.019 -0.048 -0.091  0.002 -0.014  0.038      0.028      0.033 -0.073   0.150             0.040             0.053            -0.061             0.050            0.032             0.144  0.151  0.052             0.013            0.032            0.005                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.021  0.037  0.025  0.031 -0.052 -0.175  0.168  0.051 -0.025  0.057  0.015  0.034 -0.010  0.011 -0.024 -0.057     -0.068      0.246  0.294  -0.313            -0.101             0.072            -0.131            -0.018           -0.020            -0.536 -0.329  0.067            -0.080            0.003           -0.002      -0.130                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.042  0.076  0.041  0.051  0.049 -0.217  0.167  0.097 -0.017  0.053  0.016  0.038  0.011  0.025 -0.054 -0.062     -0.120     -0.188  0.371  -0.254            -0.164             0.061            -0.144            -0.024           -0.028            -0.282 -0.570  0.058            -0.040            0.039           -0.028      -0.147       0.591                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.066 -0.152 -0.072  0.039 -0.006 -0.038 -0.133 -0.045 -0.001  0.018  0.014  0.016  0.117 -0.182 -0.001  0.174     -0.122      0.029  0.055   0.336             0.035            -0.020            -0.026            -0.019           -0.054            -0.023 -0.040 -0.177             0.539           -0.002           -0.524       0.012       0.038  0.067               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.022  0.052  0.023  0.019  0.029 -0.064  0.102  0.039  0.041 -0.057  0.068  0.025  0.043  0.010 -0.005 -0.057     -0.053     -0.040  0.148  -0.215            -0.089            -0.089             0.093            -0.139           -0.033            -0.094 -0.100 -0.076            -0.016           -0.004            0.003      -0.632       0.176  0.201 -0.027        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -4.09949639 -0.69495289  0.06893221  0.70470990  3.34093227 
    ## 
    ## Number of Observations: 30756
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30756

``` r
summary(modTfullJbetalast)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##         AIC       BIC   logLik
    ##   -21274.32 -20924.35 10679.16
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.1911675 (Intr)
    ## temptrend_abs.sc 0.0644066 -0.353
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 0.0001567023 0.2253046
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##    power 
    ## -0.16687 
    ## Fixed effects: Jbetalast ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.6480870 0.06097241 30525  10.629183  0.0000
    ## temptrend_abs.sc                                  0.0300696 0.04170541 30525   0.721000  0.4709
    ## REALMMarine                                      -0.0172125 0.06547538   195  -0.262885  0.7929
    ## REALMTerrestrial                                 -0.2197198 0.06392942   195  -3.436912  0.0007
    ## tsign1                                           -0.0350154 0.00319123 30525 -10.972411  0.0000
    ## tempave_metab.sc                                  0.0430106 0.00565710 30525   7.602944  0.0000
    ## seas.sc                                           0.0032365 0.00340297 30525   0.951082  0.3416
    ## microclim.sc                                      0.0058100 0.00173620 30525   3.346396  0.0008
    ## mass.sc                                          -0.0044232 0.00360042 30525  -1.228525  0.2193
    ## speed.sc                                         -0.0165685 0.00375062 30525  -4.417523  0.0000
    ## consumerfrac.sc                                   0.0004154 0.00222265 30525   0.186905  0.8517
    ## nspp.sc                                           0.0487814 0.00278605 30525  17.509188  0.0000
    ## npp.sc                                            0.0026629 0.00253874 30525   1.048906  0.2942
    ## veg.sc                                           -0.0157968 0.00648937 30525  -2.434262  0.0149
    ## duration.sc                                       0.0234522 0.00216123 30525  10.851335  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.0324635 0.04325238 30525  -0.750561  0.4529
    ## temptrend_abs.sc:REALMTerrestrial                -0.0443928 0.04797908 30525  -0.925254  0.3548
    ## temptrend_abs.sc:tsign1                           0.0184937 0.00397435 30525   4.653263  0.0000
    ## temptrend_abs.sc:tempave_metab.sc                 0.0056697 0.00751639 30525   0.754315  0.4507
    ## temptrend_abs.sc:seas.sc                          0.0023506 0.00338776 30525   0.693863  0.4878
    ## temptrend_abs.sc:microclim.sc                    -0.0026298 0.00226958 30525  -1.158713  0.2466
    ## temptrend_abs.sc:mass.sc                          0.0059629 0.00350948 30525   1.699074  0.0893
    ## temptrend_abs.sc:speed.sc                         0.0002207 0.00415737 30525   0.053094  0.9577
    ## temptrend_abs.sc:consumerfrac.sc                 -0.0037849 0.00280944 30525  -1.347220  0.1779
    ## temptrend_abs.sc:nspp.sc                         -0.0059433 0.00314063 30525  -1.892385  0.0584
    ## tsign-1:thermal_bias.sc                          -0.0107610 0.00338966 30525  -3.174667  0.0015
    ## tsign1:thermal_bias.sc                           -0.0139546 0.00271182 30525  -5.145830  0.0000
    ## temptrend_abs.sc:npp.sc                           0.0026227 0.00342015 30525   0.766843  0.4432
    ## temptrend_abs.sc:veg.sc                          -0.0007133 0.00412852 30525  -0.172767  0.8628
    ## temptrend_abs.sc:duration.sc                      0.0034056 0.00291968 30525   1.166413  0.2435
    ## human_bowler.sc:REALM2TerrFresh                  -0.0151605 0.00526753 30525  -2.878102  0.0040
    ## human_bowler.sc:REALM2Marine                     -0.0049220 0.00206445 30525  -2.384148  0.0171
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0053799 0.00314054 30525   1.713063  0.0867
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0014120 0.00299765 30525   0.471034  0.6376
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.0029331 0.00251115 30525   1.168027  0.2428
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0007023 0.00314053 30525  -0.223627  0.8230
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.431                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.936  0.402                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.849  0.368  0.787                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.035  0.033 -0.005  0.013                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.087 -0.040 -0.069 -0.137 -0.038                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.083  0.036  0.122 -0.022 -0.094  0.033                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.047  0.039  0.053  0.006 -0.016 -0.169  0.121                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.081 -0.004 -0.047 -0.045  0.020  0.116  0.010  0.003                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.041 -0.045 -0.011 -0.034 -0.065  0.031  0.040  0.051 -0.520                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.003  0.017  0.004  0.049  0.043 -0.076 -0.047  0.015 -0.026 -0.032                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.000 -0.025 -0.022 -0.028  0.039 -0.198  0.025 -0.032 -0.062  0.134  0.082                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.045 -0.005 -0.057  0.027 -0.003 -0.029 -0.310 -0.271 -0.016  0.097 -0.027 -0.144                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.326  0.145  0.340  0.013 -0.018 -0.002  0.122  0.045  0.012 -0.007 -0.006  0.022 -0.188                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.000  0.005  0.019  0.026 -0.195  0.037 -0.054 -0.012 -0.015  0.003 -0.011 -0.329  0.038 -0.006                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.417 -0.953 -0.427 -0.351 -0.010  0.043 -0.062 -0.045 -0.015  0.024 -0.014  0.039  0.024 -0.163 -0.021                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.332 -0.804 -0.306 -0.439 -0.013  0.139  0.066  0.001 -0.038  0.081 -0.033  0.048 -0.037  0.020 -0.026  0.763                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.017 -0.083  0.003 -0.009 -0.519  0.007  0.011 -0.007 -0.006  0.034 -0.010 -0.009 -0.026  0.026  0.098  0.025      0.028                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.028  0.067  0.023  0.100 -0.024 -0.468 -0.020  0.090 -0.033  0.026  0.062  0.071 -0.043  0.010  0.016 -0.063     -0.302      0.008                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.042 -0.071 -0.058  0.046  0.049  0.017 -0.533 -0.129 -0.006 -0.017  0.012 -0.003  0.139 -0.124 -0.016  0.107     -0.167     -0.034 -0.012                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.025 -0.043 -0.028  0.006  0.000  0.086 -0.094 -0.561  0.005 -0.024 -0.015  0.018  0.150 -0.026 -0.014  0.048     -0.051      0.032 -0.051   0.229                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.011 -0.004  0.001 -0.017 -0.010 -0.055 -0.005  0.007 -0.568  0.315  0.057  0.030 -0.011 -0.018  0.033  0.045      0.075     -0.006  0.086  -0.015            -0.058                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.017  0.083  0.008  0.052  0.041  0.028 -0.009 -0.043  0.284 -0.537  0.016 -0.073 -0.025 -0.013 -0.011 -0.045     -0.150     -0.036 -0.145   0.043             0.068            -0.467                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.002 -0.022 -0.001 -0.012 -0.010  0.054  0.023 -0.006  0.035  0.007 -0.637 -0.042  0.008  0.007 -0.009  0.023      0.057      0.013 -0.098  -0.027            -0.005            -0.105            -0.048                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.007  0.053  0.011  0.031 -0.002  0.079 -0.014  0.027  0.024 -0.073 -0.039 -0.520  0.062 -0.045  0.138 -0.070     -0.096     -0.026 -0.075  -0.021            -0.005            -0.024             0.110             0.079                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.049 -0.026 -0.044 -0.049 -0.087  0.386 -0.170 -0.088  0.050 -0.039 -0.014 -0.092 -0.035 -0.001  0.052  0.030      0.029     -0.019 -0.118   0.151             0.046            -0.019             0.018             0.011            0.033                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.075 -0.052 -0.073 -0.078  0.049  0.538 -0.242 -0.223  0.065 -0.044 -0.015 -0.106 -0.063 -0.017  0.064  0.043      0.054      0.063 -0.192   0.178             0.110            -0.019             0.030             0.015            0.041             0.469                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.001 -0.002  0.006 -0.029 -0.015 -0.003  0.103  0.139 -0.011 -0.028  0.008  0.056 -0.550  0.115 -0.011 -0.032      0.071      0.026  0.129  -0.203            -0.308             0.065            -0.027             0.001           -0.053             0.007  0.011                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.154 -0.241 -0.164  0.017  0.001 -0.004 -0.156 -0.015 -0.018 -0.005  0.011 -0.031  0.220 -0.515 -0.017  0.283     -0.077      0.025 -0.038   0.339             0.017             0.031             0.054            -0.024            0.034             0.027  0.019 -0.386                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.015  0.089 -0.011 -0.018  0.010  0.101 -0.041  0.004  0.022 -0.028  0.030  0.141 -0.005  0.010 -0.405  0.018      0.030     -0.050 -0.034   0.031            -0.015            -0.041            -0.002             0.025           -0.183             0.016 -0.049 -0.001             0.017                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.058  0.015  0.059 -0.005  0.009  0.027  0.044  0.070  0.009 -0.013 -0.014 -0.006 -0.099  0.136 -0.013 -0.019      0.018     -0.002 -0.005  -0.066            -0.001             0.003             0.008             0.008            0.035             0.028  0.051  0.036            -0.109            0.004                                                               
    ## human_bowler.sc:REALM2Marine                      0.027 -0.027 -0.028 -0.016 -0.015  0.094 -0.192 -0.075 -0.024  0.045 -0.016 -0.049 -0.085  0.001 -0.013  0.031      0.024      0.033 -0.071   0.146             0.043             0.048            -0.066             0.046            0.034             0.146  0.153  0.044             0.015            0.032            0.004                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.021  0.043  0.024  0.032 -0.051 -0.206  0.172  0.061 -0.025  0.041  0.018  0.038 -0.012  0.012 -0.022 -0.060     -0.081      0.243  0.346  -0.323            -0.107             0.066            -0.100            -0.021           -0.024            -0.535 -0.341  0.072            -0.085            0.001           -0.001      -0.131                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.038  0.071  0.036  0.049  0.047 -0.252  0.169  0.105 -0.016  0.034  0.018  0.047  0.010  0.026 -0.051 -0.062     -0.124     -0.179  0.425  -0.261            -0.165             0.050            -0.110            -0.026           -0.041            -0.294 -0.572  0.062            -0.045            0.035           -0.027      -0.147       0.614                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.053 -0.123 -0.058  0.029 -0.004 -0.039 -0.139 -0.046  0.005  0.014  0.010  0.018  0.119 -0.182 -0.001  0.141     -0.099      0.026  0.058   0.344             0.041            -0.028            -0.021            -0.016           -0.059            -0.020 -0.037 -0.179             0.540           -0.002           -0.524       0.014       0.032  0.064               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.018  0.045  0.017  0.015  0.029 -0.066  0.099  0.041  0.039 -0.061  0.065  0.028  0.036  0.011 -0.006 -0.047     -0.047     -0.039  0.146  -0.209            -0.095            -0.079             0.101            -0.132           -0.037            -0.095 -0.100 -0.060            -0.019           -0.001            0.003      -0.634       0.180  0.203 -0.028        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.90186096 -0.58375041  0.06240915  0.67203220  4.35067517 
    ## 
    ## Number of Observations: 30756
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30756

``` r
summary(modTfullHornlast)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   14979.44 15329.41 -7447.718
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.24177490 (Intr)
    ## temptrend_abs.sc 0.08912252 -0.537
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 0.0005354858 0.3336718
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.0521168 
    ## Fixed effects: Hornlast ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.4839726 0.09281608 30525   5.214318  0.0000
    ## temptrend_abs.sc                                  0.0466911 0.06102181 30525   0.765155  0.4442
    ## REALMMarine                                      -0.0024070 0.09897047   195  -0.024320  0.9806
    ## REALMTerrestrial                                 -0.2415477 0.09420311   195  -2.564116  0.0111
    ## tsign1                                           -0.0314983 0.00569624 30525  -5.529669  0.0000
    ## tempave_metab.sc                                  0.0653347 0.00976111 30525   6.693365  0.0000
    ## seas.sc                                           0.0220347 0.00606842 30525   3.631042  0.0003
    ## microclim.sc                                      0.0190821 0.00314842 30525   6.060836  0.0000
    ## mass.sc                                          -0.0002872 0.00614698 30525  -0.046714  0.9627
    ## speed.sc                                         -0.0795207 0.00653604 30525 -12.166482  0.0000
    ## consumerfrac.sc                                   0.0175392 0.00405791 30525   4.322226  0.0000
    ## nspp.sc                                           0.0575305 0.00492769 30525  11.674945  0.0000
    ## npp.sc                                            0.0173340 0.00456738 30525   3.795182  0.0001
    ## veg.sc                                           -0.0258411 0.01227618 30525  -2.104980  0.0353
    ## duration.sc                                       0.0379497 0.00378329 30525  10.030868  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.0122677 0.06311399 30525  -0.194374  0.8459
    ## temptrend_abs.sc:REALMTerrestrial                -0.0822457 0.07004721 30525  -1.174147  0.2403
    ## temptrend_abs.sc:tsign1                           0.0138542 0.00679426 30525   2.039108  0.0414
    ## temptrend_abs.sc:tempave_metab.sc                 0.0027632 0.01182581 30525   0.233662  0.8152
    ## temptrend_abs.sc:seas.sc                          0.0118294 0.00566446 30525   2.088356  0.0368
    ## temptrend_abs.sc:microclim.sc                     0.0000812 0.00385254 30525   0.021083  0.9832
    ## temptrend_abs.sc:mass.sc                          0.0081140 0.00578361 30525   1.402925  0.1606
    ## temptrend_abs.sc:speed.sc                         0.0210096 0.00687490 30525   3.055988  0.0022
    ## temptrend_abs.sc:consumerfrac.sc                 -0.0060846 0.00481988 30525  -1.262407  0.2068
    ## temptrend_abs.sc:nspp.sc                         -0.0083740 0.00525471 30525  -1.593623  0.1110
    ## tsign-1:thermal_bias.sc                          -0.0166753 0.00595206 30525  -2.801592  0.0051
    ## tsign1:thermal_bias.sc                           -0.0117748 0.00479459 30525  -2.455857  0.0141
    ## temptrend_abs.sc:npp.sc                           0.0044904 0.00578173 30525   0.776654  0.4374
    ## temptrend_abs.sc:veg.sc                           0.0106639 0.00708277 30525   1.505611  0.1322
    ## temptrend_abs.sc:duration.sc                      0.0058785 0.00507171 30525   1.159068  0.2464
    ## human_bowler.sc:REALM2TerrFresh                   0.0009301 0.00999028 30525   0.093096  0.9258
    ## human_bowler.sc:REALM2Marine                     -0.0116930 0.00372528 30525  -3.138840  0.0017
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.0060117 0.00526746 30525  -1.141298  0.2538
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.0111413 0.00503106 30525  -2.214504  0.0268
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.0060363 0.00433105 30525   1.393737  0.1634
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine     0.0047608 0.00536996 30525   0.886563  0.3753
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.507                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.948  0.477                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.817  0.428  0.762                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.041  0.040 -0.004  0.018                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.089 -0.044 -0.073 -0.157 -0.041                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.094  0.049  0.139 -0.038 -0.089  0.042                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.055  0.045  0.063  0.005 -0.016 -0.157  0.116                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.078 -0.010 -0.044 -0.043  0.018  0.121  0.020  0.002                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.049 -0.050 -0.015 -0.038 -0.059 -0.013  0.026  0.059 -0.510                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.009  0.020  0.007  0.067  0.047 -0.088 -0.056  0.014 -0.021 -0.054                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.007 -0.025 -0.030 -0.040  0.037 -0.184  0.021 -0.030 -0.050  0.135  0.077                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.052 -0.012 -0.069  0.038 -0.005 -0.033 -0.321 -0.283 -0.018  0.095 -0.020 -0.142                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.409  0.178  0.429  0.014 -0.017  0.003  0.129  0.044  0.009 -0.005 -0.003  0.025 -0.187                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.001  0.003  0.022  0.035 -0.194  0.027 -0.054 -0.011 -0.024  0.010 -0.006 -0.315  0.037 -0.010                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.493 -0.951 -0.510 -0.409 -0.013  0.048 -0.082 -0.051 -0.011  0.024 -0.015  0.041  0.034 -0.199 -0.023                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.374 -0.774 -0.347 -0.511 -0.016  0.149  0.072  0.003 -0.033  0.077 -0.043  0.054 -0.043  0.021 -0.030  0.733                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.018 -0.098  0.004 -0.011 -0.528  0.009  0.009 -0.008 -0.006  0.032 -0.013 -0.006 -0.023  0.029  0.100  0.033      0.035                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.034  0.065  0.028  0.112 -0.022 -0.484 -0.022  0.088 -0.046  0.047  0.067  0.071 -0.043  0.011  0.016 -0.060     -0.317      0.010                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.054 -0.082 -0.073  0.057  0.047  0.011 -0.538 -0.124 -0.012 -0.014  0.018 -0.004  0.146 -0.127 -0.016  0.121     -0.197     -0.029 -0.003                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.032 -0.051 -0.035  0.005  0.001  0.083 -0.092 -0.560  0.004 -0.028 -0.016  0.018  0.154 -0.028 -0.016  0.056     -0.055      0.036 -0.059   0.224                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.012 -0.009  0.000 -0.018 -0.010 -0.065 -0.012  0.006 -0.581  0.322  0.059  0.027 -0.010 -0.016  0.040  0.051      0.083     -0.006  0.110  -0.008            -0.052                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.020  0.087  0.009  0.053  0.040  0.050 -0.006 -0.046  0.294 -0.549  0.021 -0.074 -0.023 -0.014 -0.015 -0.049     -0.151     -0.038 -0.206   0.040             0.068            -0.495                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.000 -0.025 -0.003 -0.020 -0.012  0.058  0.028 -0.007  0.038  0.013 -0.634 -0.039  0.006  0.006 -0.012  0.025      0.067      0.014 -0.111  -0.031            -0.003            -0.117            -0.050                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.007  0.064  0.010  0.040  0.002  0.072 -0.014  0.027  0.021 -0.075 -0.039 -0.522  0.062 -0.048  0.133 -0.078     -0.119     -0.034 -0.065  -0.016            -0.008            -0.018             0.121             0.079                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.048 -0.029 -0.045 -0.054 -0.081  0.377 -0.158 -0.082  0.049 -0.059 -0.017 -0.080 -0.036  0.002  0.043  0.034      0.030     -0.020 -0.123   0.148             0.042            -0.020             0.027             0.011            0.029                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.075 -0.060 -0.076 -0.086  0.048  0.523 -0.217 -0.213  0.063 -0.073 -0.019 -0.088 -0.066 -0.012  0.052  0.051      0.058      0.065 -0.197   0.174             0.109            -0.021             0.044             0.015            0.035             0.460                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.006 -0.002  0.012 -0.035 -0.012 -0.006  0.110  0.143 -0.010 -0.026  0.006  0.055 -0.552  0.113 -0.010 -0.035      0.083      0.021  0.136  -0.200            -0.303             0.068            -0.032             0.003           -0.050             0.007  0.008                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.203 -0.285 -0.216  0.019  0.000 -0.004 -0.161 -0.017 -0.018 -0.005  0.012 -0.033  0.223 -0.533 -0.016  0.331     -0.089      0.025 -0.040   0.338             0.016             0.028             0.050            -0.027            0.039             0.026  0.019 -0.384                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.019  0.111 -0.011 -0.022  0.009  0.099 -0.043  0.006  0.026 -0.031  0.028  0.133 -0.003  0.014 -0.403  0.017      0.035     -0.047 -0.033   0.030            -0.015            -0.049            -0.003             0.023           -0.165             0.014 -0.051 -0.003             0.014                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.071  0.025  0.074 -0.016  0.009  0.028  0.050  0.061  0.013 -0.014 -0.017 -0.008 -0.092  0.151 -0.011 -0.029      0.018     -0.002 -0.005  -0.066             0.001             0.001             0.010             0.010            0.036             0.027  0.048  0.035            -0.113            0.002                                                               
    ## human_bowler.sc:REALM2Marine                      0.031 -0.032 -0.032 -0.015 -0.014  0.086 -0.198 -0.079 -0.027  0.039 -0.016 -0.050 -0.079 -0.002 -0.017  0.038      0.026      0.033 -0.071   0.151             0.045             0.052            -0.060             0.046            0.033             0.146  0.151  0.044             0.017            0.031            0.002                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.022  0.040  0.025  0.032 -0.052 -0.197  0.165  0.054 -0.027  0.053  0.016  0.035 -0.010  0.011 -0.020 -0.058     -0.071      0.242  0.310  -0.317            -0.102             0.072            -0.124            -0.021           -0.021            -0.549 -0.342  0.068            -0.081            0.001            0.002      -0.132                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.043  0.077  0.042  0.051  0.051 -0.243  0.162  0.102 -0.018  0.049  0.017  0.040  0.010  0.027 -0.051 -0.064     -0.120     -0.190  0.387  -0.258            -0.165             0.059            -0.135            -0.026           -0.032            -0.297 -0.583  0.061            -0.045            0.036           -0.024      -0.148       0.600                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.069 -0.142 -0.074  0.037 -0.005 -0.038 -0.139 -0.043  0.001  0.016  0.013  0.017  0.118 -0.187 -0.001  0.162     -0.111      0.026  0.053   0.336             0.035            -0.022            -0.025            -0.019           -0.056            -0.022 -0.038 -0.175             0.530           -0.002           -0.552       0.014       0.032  0.062               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.022  0.051  0.022  0.017  0.031 -0.064  0.101  0.045  0.040 -0.056  0.061  0.028  0.035  0.013 -0.004 -0.055     -0.050     -0.041  0.147  -0.215            -0.094            -0.089             0.094            -0.132           -0.034            -0.097 -0.100 -0.066            -0.020           -0.003            0.004      -0.631       0.179  0.203 -0.028        
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -3.1758666 -0.7868590  0.0732762  0.7828231  2.6922253 
    ## 
    ## Number of Observations: 30756
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30756

``` r
rsquared(modTfullJtulast)
```

    ##   Response   family     link method   Marginal Conditional
    ## 1  Jtulast gaussian identity   none 0.05449033   0.2942132

``` r
rsquared(modTfullJbetalast)
```

    ##    Response   family     link method   Marginal Conditional
    ## 1 Jbetalast gaussian identity   none 0.05631666    0.474508

``` r
rsquared(modTfullHornlast)
```

    ##   Response   family     link method  Marginal Conditional
    ## 1 Hornlast gaussian identity   none 0.0690544   0.4149697

#### Summary exp

``` r
summary(modTfullJtuexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   80514.95 80849.29 -40215.48
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.824904663 (Intr)
    ## temptrend_abs.sc 0.006059369 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:     1.53867 148.5071
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -4.702845 
    ## Fixed effects: log(Jtuexp) ~ temptrend_abs.sc * REALM + temptrend_abs.sc * tsign +      temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.2179012 0.4127713 21001   0.527898  0.5976
    ## temptrend_abs.sc                                  0.4381675 0.2177658 21001   2.012104  0.0442
    ## REALMMarine                                      -0.0079716 0.4395651   171  -0.018135  0.9856
    ## REALMTerrestrial                                  1.3662814 0.3885639   171   3.516234  0.0006
    ## tsign1                                           -0.1056573 0.0351987 21001  -3.001734  0.0027
    ## tempave_metab.sc                                 -0.3217130 0.0506670 21001  -6.349553  0.0000
    ## seas.sc                                          -0.1627161 0.0349263 21001  -4.658849  0.0000
    ## microclim.sc                                     -0.0999135 0.0186730 21001  -5.350704  0.0000
    ## mass.sc                                           0.0352586 0.0336777 21001   1.046940  0.2951
    ## speed.sc                                          0.0113253 0.0364660 21001   0.310571  0.7561
    ## consumerfrac.sc                                   0.0615543 0.0223711 21001   2.751516  0.0059
    ## nspp.sc                                          -0.4618508 0.0278539 21001 -16.581179  0.0000
    ## npp.sc                                           -0.0035704 0.0267945 21001  -0.133250  0.8940
    ## veg.sc                                            0.0200687 0.0745542 21001   0.269182  0.7878
    ## duration.sc                                       0.5723267 0.0220193 21001  25.992058  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.3144930 0.2171823 21001  -1.448060  0.1476
    ## temptrend_abs.sc:REALMTerrestrial                -0.2367383 0.2202724 21001  -1.074753  0.2825
    ## temptrend_abs.sc:tsign1                           0.0680726 0.0400978 21001   1.697663  0.0896
    ## temptrend_abs.sc:tempave_metab.sc                 0.0120160 0.0337546 21001   0.355980  0.7219
    ## temptrend_abs.sc:seas.sc                          0.0238713 0.0288718 21001   0.826803  0.4084
    ## temptrend_abs.sc:microclim.sc                     0.0488741 0.0194391 21001   2.514211  0.0119
    ## temptrend_abs.sc:mass.sc                          0.0202724 0.0285294 21001   0.710581  0.4774
    ## temptrend_abs.sc:speed.sc                        -0.0108812 0.0336145 21001  -0.323706  0.7462
    ## temptrend_abs.sc:consumerfrac.sc                 -0.0594020 0.0240029 21001  -2.474780  0.0133
    ## temptrend_abs.sc:nspp.sc                          0.0772743 0.0249605 21001   3.095863  0.0020
    ## tsign-1:thermal_bias.sc                          -0.1227645 0.0353982 21001  -3.468096  0.0005
    ## tsign1:thermal_bias.sc                           -0.0309613 0.0273078 21001  -1.133789  0.2569
    ## temptrend_abs.sc:npp.sc                           0.0097879 0.0276696 21001   0.353743  0.7235
    ## temptrend_abs.sc:veg.sc                          -0.0877197 0.0423139 21001  -2.073071  0.0382
    ## temptrend_abs.sc:duration.sc                      0.1048629 0.0292614 21001   3.583664  0.0003
    ## human_bowler.sc:REALM2TerrFresh                   0.1234357 0.0593481 21001   2.079858  0.0376
    ## human_bowler.sc:REALM2Marine                      0.0456744 0.0221125 21001   2.065543  0.0389
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0508621 0.0254689 21001   1.997030  0.0458
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0254364 0.0244166 21001   1.041766  0.2975
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0565149 0.0256205 21001  -2.205845  0.0274
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0499449 0.0284749 21001  -1.753997  0.0794
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.346                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.961  0.332                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.736  0.201  0.685                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.061  0.069  0.001  0.031                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.072  0.001 -0.067 -0.178 -0.045                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.076  0.016  0.132 -0.101 -0.067  0.055                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.069  0.085  0.081 -0.006 -0.004 -0.127  0.126                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.068  0.015 -0.034 -0.028  0.005  0.147  0.070  0.021                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.060 -0.041 -0.027 -0.051 -0.048 -0.163 -0.030  0.061 -0.537                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.020  0.030  0.008  0.097  0.085 -0.131 -0.100 -0.001 -0.019 -0.061                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.049 -0.069 -0.070 -0.102  0.065 -0.152  0.001 -0.024 -0.022  0.153  0.122                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.047  0.008 -0.070  0.094 -0.027 -0.024 -0.350 -0.262 -0.009  0.093  0.013 -0.159                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.547  0.292  0.573  0.000 -0.010  0.007  0.107  0.036 -0.003 -0.003  0.000  0.013 -0.179                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.004  0.006  0.031  0.048 -0.194  0.030 -0.051 -0.027 -0.076  0.040 -0.030 -0.277  0.028  0.004                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.354 -0.962 -0.355 -0.198 -0.028  0.029 -0.047 -0.087 -0.019  0.033 -0.013  0.070  0.010 -0.329 -0.026                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.156 -0.634 -0.138 -0.288 -0.025  0.047  0.172  0.038 -0.068  0.090 -0.028  0.167 -0.161  0.064 -0.039  0.609                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.037 -0.189 -0.009 -0.019 -0.514  0.010 -0.017 -0.030  0.025  0.016 -0.072  0.003  0.005  0.014  0.088  0.082      0.058                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.002  0.011 -0.009  0.042 -0.048 -0.258 -0.107  0.036 -0.208  0.267  0.133 -0.001  0.017 -0.005  0.006  0.023     -0.144      0.017                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.043 -0.076 -0.064  0.100  0.022 -0.006 -0.464 -0.147 -0.048  0.025  0.050  0.010  0.163 -0.103 -0.018  0.112     -0.386      0.015  0.091                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.040 -0.128 -0.048  0.018 -0.026  0.083 -0.132 -0.491 -0.031 -0.034  0.004  0.003  0.049 -0.020 -0.014  0.125     -0.134      0.097 -0.027   0.339                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                          0.010 -0.035 -0.012 -0.033  0.025 -0.119 -0.073 -0.023 -0.511  0.348  0.061  0.010 -0.026 -0.010  0.099  0.053      0.137     -0.095  0.468   0.060             0.000                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.018  0.075  0.028  0.047  0.028  0.152  0.071 -0.024  0.350 -0.474 -0.046 -0.092 -0.018 -0.002 -0.033 -0.092     -0.199     -0.014 -0.619  -0.034             0.056            -0.693                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.008 -0.016  0.014 -0.016 -0.063  0.095  0.084  0.002  0.043 -0.052 -0.542 -0.121 -0.049  0.014  0.041  0.009      0.039      0.119 -0.258  -0.086            -0.025            -0.156             0.104                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.033  0.139  0.029  0.095 -0.032  0.000 -0.001 -0.002 -0.038 -0.071 -0.123 -0.429  0.115 -0.035  0.038 -0.110     -0.323     -0.059  0.068  -0.032             0.047             0.089             0.125             0.252                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.033 -0.006 -0.035 -0.059 -0.051  0.295 -0.144 -0.064  0.018 -0.093 -0.008 -0.063 -0.029  0.003  0.058  0.015      0.012     -0.033  0.043   0.118             0.016             0.028             0.022            -0.001            0.040                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.067 -0.076 -0.074 -0.102  0.012  0.443 -0.200 -0.189  0.018 -0.107 -0.016 -0.021 -0.076  0.001  0.051  0.058      0.043      0.123  0.031   0.151             0.078             0.036             0.010            -0.003           -0.077             0.413                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.007 -0.019 -0.003 -0.103  0.032  0.022  0.128  0.030 -0.014 -0.022 -0.060  0.078 -0.482  0.107 -0.032 -0.009      0.320     -0.052 -0.032  -0.262            -0.103             0.086            -0.018             0.121           -0.118             0.019  0.052                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.284 -0.503 -0.295  0.019 -0.012 -0.004 -0.089 -0.035  0.016 -0.009  0.015 -0.031  0.191 -0.545 -0.010  0.565     -0.149      0.041  0.011   0.239             0.061            -0.024             0.033            -0.048            0.045             0.013  0.009 -0.391                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.034  0.220 -0.003 -0.025  0.001  0.116 -0.041  0.024  0.074 -0.060  0.046  0.070  0.005  0.011 -0.328 -0.011      0.035     -0.056 -0.009   0.010            -0.015            -0.137             0.009             0.011           -0.058            -0.016 -0.065  0.011            -0.013                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.071 -0.011  0.075 -0.036  0.006  0.020 -0.007  0.124  0.022 -0.011 -0.027  0.001 -0.110  0.143 -0.006  0.009      0.036      0.007  0.014   0.008            -0.121             0.016            -0.003             0.020            0.004             0.040  0.041  0.099            -0.070            0.014                                                               
    ## human_bowler.sc:REALM2Marine                      0.030 -0.034 -0.034 -0.020 -0.005  0.082 -0.200 -0.100 -0.021  0.020 -0.059 -0.030 -0.131  0.012 -0.027  0.040      0.042      0.025 -0.043   0.165             0.046             0.043            -0.059             0.119            0.017             0.143  0.161  0.171            -0.039            0.032            0.030                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.001  0.015  0.003  0.018 -0.050 -0.005  0.100  0.005  0.041  0.035 -0.019  0.005 -0.009  0.004 -0.075 -0.042     -0.038      0.242 -0.057  -0.250            -0.076            -0.082            -0.043             0.046            0.041            -0.524 -0.247  0.000            -0.042            0.060           -0.007      -0.145                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.038  0.148  0.032  0.048  0.130 -0.027  0.124  0.036  0.048  0.013  0.013 -0.072  0.021  0.011 -0.082 -0.111     -0.131     -0.356  0.006  -0.247            -0.120            -0.071            -0.019            -0.012            0.224            -0.225 -0.530 -0.033            -0.035            0.074            0.001      -0.177       0.442                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.085 -0.221 -0.091  0.033 -0.011  0.000 -0.062 -0.124  0.010 -0.001  0.021 -0.003  0.135 -0.193 -0.003  0.248     -0.144      0.032  0.001   0.240             0.209            -0.058             0.006            -0.043           -0.013            -0.030 -0.026 -0.278             0.517           -0.021           -0.509      -0.030       0.049  0.019               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.025  0.067  0.026  0.029  0.019 -0.050  0.110  0.052  0.018 -0.027  0.139 -0.004  0.144 -0.009  0.017 -0.060     -0.077     -0.027  0.105  -0.265            -0.094            -0.054             0.090            -0.289           -0.013            -0.101 -0.120 -0.323             0.078           -0.011           -0.039      -0.590       0.218  0.270  0.042        
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.4027424041 -0.0121954760  0.0000114472  0.0337219642  1.1557977502 
    ## 
    ## Number of Observations: 21208
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    174                  21208

``` r
summary(modTfullJbetaexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   71205.05 71538.44 -35560.53
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.8554396 (Intr)
    ## temptrend_abs.sc 0.2439011 -0.825
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:       1.188 5.389029
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.496251 
    ## Fixed effects: log(Jbetaexp) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF    t-value p-value
    ## (Intercept)                                      -0.1876228 0.3753932 20517  -0.499803  0.6172
    ## temptrend_abs.sc                                  0.2734133 0.2289477 20517   1.194217  0.2324
    ## REALMMarine                                      -0.1211093 0.4007652   181  -0.302195  0.7629
    ## REALMTerrestrial                                  1.2997199 0.3676258   181   3.535443  0.0005
    ## tsign1                                           -0.1278899 0.0303541 20517  -4.213263  0.0000
    ## tempave_metab.sc                                 -0.3207036 0.0472450 20517  -6.788093  0.0000
    ## seas.sc                                          -0.1780568 0.0312132 20517  -5.704543  0.0000
    ## microclim.sc                                     -0.0551856 0.0165067 20517  -3.343214  0.0008
    ## mass.sc                                           0.0181275 0.0307105 20517   0.590270  0.5550
    ## speed.sc                                         -0.0500452 0.0327535 20517  -1.527933  0.1265
    ## consumerfrac.sc                                  -0.0248085 0.0194944 20517  -1.272599  0.2032
    ## nspp.sc                                          -0.5971659 0.0253019 20517 -23.601597  0.0000
    ## npp.sc                                            0.0010995 0.0238128 20517   0.046172  0.9632
    ## veg.sc                                            0.0338301 0.0606210 20517   0.558058  0.5768
    ## duration.sc                                       0.5431752 0.0202847 20517  26.777526  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.1391537 0.2333756 20517  -0.596265  0.5510
    ## temptrend_abs.sc:REALMTerrestrial                -0.5181060 0.2518857 20517  -2.056909  0.0397
    ## temptrend_abs.sc:tsign1                           0.0372215 0.0372560 20517   0.999074  0.3178
    ## temptrend_abs.sc:tempave_metab.sc                 0.1009891 0.0455399 20517   2.217595  0.0266
    ## temptrend_abs.sc:seas.sc                          0.0681084 0.0286186 20517   2.379862  0.0173
    ## temptrend_abs.sc:microclim.sc                     0.0055250 0.0202739 20517   0.272516  0.7852
    ## temptrend_abs.sc:mass.sc                         -0.0207843 0.0278039 20517  -0.747529  0.4548
    ## temptrend_abs.sc:speed.sc                         0.0361310 0.0332937 20517   1.085218  0.2778
    ## temptrend_abs.sc:consumerfrac.sc                  0.0423555 0.0236469 20517   1.791166  0.0733
    ## temptrend_abs.sc:nspp.sc                          0.2002581 0.0256763 20517   7.799336  0.0000
    ## tsign-1:thermal_bias.sc                          -0.0767892 0.0315679 20517  -2.432511  0.0150
    ## tsign1:thermal_bias.sc                            0.0116906 0.0245905 20517   0.475412  0.6345
    ## temptrend_abs.sc:npp.sc                           0.0415487 0.0289380 20517   1.435782  0.1511
    ## temptrend_abs.sc:veg.sc                          -0.0424486 0.0387222 20517  -1.096236  0.2730
    ## temptrend_abs.sc:duration.sc                      0.0962572 0.0267534 20517   3.597943  0.0003
    ## human_bowler.sc:REALM2TerrFresh                   0.0866782 0.0486010 20517   1.783466  0.0745
    ## human_bowler.sc:REALM2Marine                      0.0410067 0.0194164 20517   2.111958  0.0347
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0053437 0.0270927 20517   0.197238  0.8436
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0005054 0.0257431 20517   0.019633  0.9843
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0295203 0.0234567 20517  -1.258501  0.2082
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0536310 0.0273730 20517  -1.959266  0.0501
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.561                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.953  0.527                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.765  0.439  0.712                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.057  0.066 -0.002  0.026                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.086 -0.029 -0.076 -0.178 -0.045                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.108  0.060  0.165 -0.061 -0.091  0.028                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.070  0.066  0.080  0.005 -0.003 -0.152  0.109                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.084 -0.013 -0.045 -0.045  0.015  0.117  0.047  0.004                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.057 -0.060 -0.019 -0.042 -0.060 -0.088  0.006  0.066 -0.530                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.013  0.025  0.006  0.084  0.055 -0.104 -0.070  0.009 -0.021 -0.064                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.021 -0.037 -0.050 -0.066  0.048 -0.173  0.022 -0.046 -0.038  0.142  0.086                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.066 -0.021 -0.087  0.055 -0.010 -0.018 -0.309 -0.273 -0.014  0.099 -0.018 -0.146                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.496  0.250  0.520  0.018 -0.016  0.013  0.128  0.047  0.008 -0.007 -0.005  0.018 -0.204                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.002  0.001  0.030  0.043 -0.187  0.046 -0.052 -0.010 -0.046  0.028 -0.001 -0.307  0.028 -0.001                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.554 -0.952 -0.571 -0.425 -0.027  0.046 -0.107 -0.072 -0.007  0.033 -0.011  0.056  0.050 -0.280 -0.026                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.385 -0.668 -0.355 -0.565 -0.024  0.134  0.108  0.008 -0.031  0.078 -0.056  0.092 -0.069  0.029 -0.034  0.633                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.031 -0.169 -0.005 -0.012 -0.522  0.001  0.002 -0.019 -0.002  0.037 -0.026  0.001 -0.014  0.021  0.095  0.076      0.044                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.033  0.033  0.025  0.108 -0.040 -0.471 -0.041  0.087 -0.094  0.122  0.104  0.069 -0.051  0.010  0.005 -0.027     -0.290      0.046                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.053 -0.087 -0.074  0.089  0.043  0.000 -0.541 -0.109 -0.031 -0.001  0.031 -0.010  0.140 -0.107 -0.012  0.131     -0.313     -0.008  0.048                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.043 -0.084 -0.045  0.000 -0.017  0.087 -0.079 -0.556  0.002 -0.036 -0.021  0.024  0.129 -0.028 -0.027  0.087     -0.057      0.063 -0.080   0.198                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.013 -0.031  0.001 -0.022 -0.008 -0.085 -0.030  0.006 -0.587  0.362  0.068  0.025 -0.014 -0.014  0.064  0.062      0.115     -0.019  0.237   0.017            -0.051                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.025  0.093  0.017  0.055  0.048  0.093  0.010 -0.044  0.334 -0.552 -0.002 -0.085 -0.021 -0.007 -0.030 -0.077     -0.164     -0.057 -0.389   0.013             0.069            -0.626                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.002 -0.021  0.000 -0.023 -0.016  0.068  0.036 -0.012  0.034 -0.002 -0.643 -0.053 -0.002  0.008 -0.016  0.016      0.072      0.041 -0.175  -0.048             0.008            -0.132             0.015                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.021  0.099  0.023  0.070 -0.001  0.055 -0.021  0.031  0.012 -0.089 -0.058 -0.525  0.074 -0.038  0.123 -0.104     -0.218     -0.052 -0.039   0.006            -0.012             0.000             0.152             0.118                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.043 -0.017 -0.044 -0.055 -0.065  0.329 -0.175 -0.071  0.031 -0.081 -0.009 -0.064 -0.032  0.011  0.058  0.027      0.025     -0.024 -0.090   0.140             0.025            -0.014             0.038             0.004            0.024                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.077 -0.067 -0.083 -0.094  0.032  0.472 -0.243 -0.192  0.038 -0.097 -0.015 -0.070 -0.057  0.000  0.064  0.056      0.060      0.082 -0.151   0.166             0.080            -0.015             0.047             0.013            0.015             0.435                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.006  0.003  0.012 -0.054  0.000 -0.009  0.104  0.125 -0.010 -0.028 -0.007  0.058 -0.545  0.117 -0.005 -0.045      0.141     -0.008  0.151  -0.186            -0.253             0.079            -0.042             0.028           -0.053             0.016  0.018                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.244 -0.429 -0.257  0.023 -0.007 -0.011 -0.138 -0.020 -0.014  0.004  0.017 -0.031  0.219 -0.524 -0.014  0.485     -0.121      0.052 -0.026   0.286             0.018             0.013             0.026            -0.036            0.045             0.000 -0.001 -0.387                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.028  0.174 -0.013 -0.024  0.005  0.102 -0.058  0.008  0.041 -0.045  0.030  0.128  0.017  0.013 -0.383  0.004      0.040     -0.063 -0.025   0.025             0.008            -0.081            -0.002             0.029           -0.153             0.006 -0.059 -0.013             0.000                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.077  0.045  0.081 -0.020  0.012  0.023  0.025  0.080  0.028 -0.025 -0.027 -0.005 -0.098  0.147 -0.011 -0.052      0.028     -0.004 -0.007  -0.062            -0.003            -0.009             0.016             0.017            0.019             0.034  0.055  0.041            -0.124            0.003                                                               
    ## human_bowler.sc:REALM2Marine                      0.033 -0.036 -0.035 -0.017 -0.003  0.086 -0.197 -0.075 -0.029  0.028 -0.016 -0.049 -0.102  0.006 -0.017  0.044      0.035      0.032 -0.075   0.166             0.030             0.061            -0.052             0.057            0.027             0.139  0.155  0.084            -0.002            0.035            0.006                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.008 -0.003  0.012  0.025 -0.055 -0.124  0.160  0.021 -0.022  0.072  0.007  0.028  0.001 -0.008 -0.043 -0.022     -0.040      0.249  0.158  -0.280            -0.054             0.072            -0.157            -0.001           -0.012            -0.544 -0.312  0.030            -0.019            0.015           -0.010      -0.134                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.041  0.074  0.039  0.045  0.073 -0.155  0.157  0.064 -0.012  0.060  0.014  0.020  0.014  0.016 -0.068 -0.050     -0.098     -0.225  0.209  -0.231            -0.124             0.072            -0.153            -0.021            0.007            -0.274 -0.578  0.041            -0.006            0.054           -0.036      -0.155       0.532                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.083 -0.216 -0.088  0.042 -0.008 -0.030 -0.113 -0.053 -0.009  0.024  0.020  0.004  0.118 -0.192  0.007  0.242     -0.144      0.045  0.043   0.295             0.042            -0.007            -0.031            -0.028           -0.019            -0.043 -0.049 -0.187             0.543           -0.020           -0.533       0.002       0.086  0.087               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.024  0.060  0.025  0.020  0.020 -0.058  0.114  0.040  0.046 -0.045  0.083  0.020  0.073  0.003 -0.004 -0.063     -0.066     -0.048  0.146  -0.243            -0.067            -0.107             0.084            -0.182           -0.021            -0.096 -0.109 -0.154             0.016           -0.013            0.003      -0.625       0.184  0.217 -0.009        
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -1.85753440 -0.13226179  0.01924062  0.22922469  1.72827486 
    ## 
    ## Number of Observations: 20734
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    184                  20734

``` r
summary(modTfullHornexp)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##       AIC      BIC    logLik
    ##   75913.1 76244.25 -37914.55
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      1.0139013 (Intr)
    ## temptrend_abs.sc 0.1488556 -0.682
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:    1.579336 43.33736
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -3.497225 
    ## Fixed effects: log(Hornexp) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.5196239 0.4640640 19448   1.119725  0.2628
    ## temptrend_abs.sc                                  0.4175405 0.2464312 19448   1.694349  0.0902
    ## REALMMarine                                       0.1210043 0.4949882   174   0.244459  0.8072
    ## REALMTerrestrial                                  1.8357639 0.4491982   174   4.086757  0.0001
    ## tsign1                                           -0.1033734 0.0380505 19448  -2.716745  0.0066
    ## tempave_metab.sc                                 -0.6684625 0.0577799 19448 -11.569127  0.0000
    ## seas.sc                                          -0.1728180 0.0389170 19448  -4.440685  0.0000
    ## microclim.sc                                     -0.1322821 0.0207934 19448  -6.361747  0.0000
    ## mass.sc                                          -0.0220703 0.0385218 19448  -0.572931  0.5667
    ## speed.sc                                          0.2181842 0.0410568 19448   5.314198  0.0000
    ## consumerfrac.sc                                  -0.0310141 0.0249688 19448  -1.242113  0.2142
    ## nspp.sc                                          -0.4921007 0.0317646 19448 -15.492091  0.0000
    ## npp.sc                                           -0.0712807 0.0297957 19448  -2.392311  0.0168
    ## veg.sc                                            0.0592727 0.0764798 19448   0.775012  0.4383
    ## duration.sc                                       0.5512291 0.0244927 19448  22.505824  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.3683345 0.2471952 19448  -1.490055  0.1362
    ## temptrend_abs.sc:REALMTerrestrial                -0.0798951 0.2663271 19448  -0.299989  0.7642
    ## temptrend_abs.sc:tsign1                          -0.0169873 0.0446288 19448  -0.380636  0.7035
    ## temptrend_abs.sc:tempave_metab.sc                 0.1079240 0.0468410 19448   2.304049  0.0212
    ## temptrend_abs.sc:seas.sc                         -0.0236867 0.0332040 19448  -0.713368  0.4756
    ## temptrend_abs.sc:microclim.sc                     0.0509543 0.0239826 19448   2.124632  0.0336
    ## temptrend_abs.sc:mass.sc                         -0.0118959 0.0334069 19448  -0.356091  0.7218
    ## temptrend_abs.sc:speed.sc                        -0.0342199 0.0389074 19448  -0.879521  0.3791
    ## temptrend_abs.sc:consumerfrac.sc                  0.0773801 0.0292050 19448   2.649553  0.0081
    ## temptrend_abs.sc:nspp.sc                          0.1713642 0.0294151 19448   5.825717  0.0000
    ## tsign-1:thermal_bias.sc                          -0.0614796 0.0389757 19448  -1.577384  0.1147
    ## tsign1:thermal_bias.sc                           -0.0728357 0.0308736 19448  -2.359153  0.0183
    ## temptrend_abs.sc:npp.sc                           0.0956295 0.0335763 19448   2.848127  0.0044
    ## temptrend_abs.sc:veg.sc                          -0.1821673 0.0454695 19448  -4.006362  0.0001
    ## temptrend_abs.sc:duration.sc                      0.0872082 0.0330025 19448   2.642474  0.0082
    ## human_bowler.sc:REALM2TerrFresh                   0.1772166 0.0620915 19448   2.854118  0.0043
    ## human_bowler.sc:REALM2Marine                      0.0089024 0.0244112 19448   0.364685  0.7154
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.0034447 0.0310794 19448  -0.110835  0.9117
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0158569 0.0297935 19448   0.532227  0.5946
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.1114842 0.0283045 19448  -3.938746  0.0001
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0671970 0.0326048 19448  -2.060954  0.0393
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.438                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.955  0.414                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.766  0.309  0.713                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.058  0.074 -0.002  0.029                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.079 -0.011 -0.072 -0.175 -0.055                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.087  0.051  0.144 -0.087 -0.080  0.026                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.067  0.072  0.078  0.003  0.000 -0.140  0.113                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.065  0.014 -0.031 -0.030  0.009  0.149  0.054  0.008                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.061 -0.060 -0.026 -0.042 -0.051 -0.140 -0.002  0.071 -0.527                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.016  0.031  0.009  0.084  0.062 -0.121 -0.061  0.012 -0.015 -0.069                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.030 -0.048 -0.058 -0.075  0.046 -0.155  0.010 -0.039 -0.030  0.145  0.085                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.052 -0.011 -0.074  0.074 -0.029 -0.022 -0.327 -0.278 -0.017  0.100 -0.017 -0.142                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.510  0.262  0.533  0.020 -0.014  0.009  0.120  0.037  0.007 -0.013 -0.003  0.013 -0.189                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.001 -0.005  0.030  0.046 -0.181  0.042 -0.053 -0.023 -0.054  0.022 -0.006 -0.310  0.029  0.002                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.440 -0.953 -0.450 -0.302 -0.028  0.033 -0.096 -0.078 -0.029  0.042 -0.019  0.062  0.039 -0.299 -0.022                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.259 -0.617 -0.234 -0.425 -0.027  0.096  0.153  0.011 -0.062  0.090 -0.048  0.125 -0.105  0.058 -0.046  0.582                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.036 -0.205 -0.007 -0.015 -0.523  0.011 -0.012 -0.028  0.009  0.029 -0.037  0.003  0.003  0.020  0.102  0.096      0.059                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.020  0.029  0.012  0.079 -0.041 -0.371 -0.057  0.075 -0.143  0.188  0.110  0.029 -0.030  0.004  0.002 -0.009     -0.255      0.033                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.056 -0.130 -0.080  0.099  0.033 -0.005 -0.510 -0.114 -0.031 -0.001  0.019 -0.002  0.144 -0.133 -0.004  0.177     -0.370      0.018  0.053                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.040 -0.104 -0.045 -0.004 -0.026  0.088 -0.083 -0.542 -0.009 -0.042 -0.017  0.024  0.102 -0.014 -0.007  0.103     -0.062      0.080 -0.080   0.221                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                          0.002 -0.047 -0.007 -0.037  0.008 -0.108 -0.039 -0.003 -0.558  0.365  0.056  0.009 -0.007 -0.013  0.080  0.075      0.143     -0.057  0.319   0.024            -0.028                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.020  0.089  0.020  0.057  0.035  0.130  0.025 -0.043  0.355 -0.522  0.003 -0.077 -0.029 -0.004 -0.021 -0.089     -0.191     -0.034 -0.470   0.013             0.078            -0.678                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.002 -0.029  0.001 -0.014 -0.022  0.068  0.035 -0.010  0.025 -0.005 -0.612 -0.056 -0.007  0.005 -0.002  0.027      0.063      0.052 -0.183  -0.030            -0.001            -0.115             0.001                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.022  0.098  0.021  0.081 -0.004  0.024 -0.014  0.026 -0.012 -0.072 -0.065 -0.483  0.081 -0.045  0.114 -0.093     -0.267     -0.060  0.024   0.009            -0.009             0.043             0.110             0.132                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.033 -0.007 -0.034 -0.058 -0.066  0.318 -0.158 -0.060  0.045 -0.104 -0.016 -0.069 -0.034  0.014  0.057  0.015      0.017     -0.016 -0.044   0.119             0.022            -0.014             0.047             0.004            0.033                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.073 -0.072 -0.080 -0.095  0.027  0.462 -0.231 -0.192  0.052 -0.126 -0.029 -0.060 -0.054 -0.007  0.066  0.056      0.042      0.099 -0.087   0.160             0.089            -0.019             0.061             0.013           -0.006             0.421                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.002  0.010  0.008 -0.068  0.021  0.001  0.112  0.101 -0.006 -0.024 -0.010  0.058 -0.516  0.114 -0.025 -0.053      0.201     -0.030  0.106  -0.203            -0.205             0.071            -0.037             0.035           -0.074             0.010  0.016                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.243 -0.451 -0.256  0.033 -0.007 -0.013 -0.152 -0.006 -0.010  0.007  0.012 -0.037  0.214 -0.520 -0.006  0.519     -0.181      0.050 -0.007   0.355            -0.005             0.010             0.019            -0.028            0.075             0.004  0.012 -0.409                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.035  0.222 -0.004 -0.021  0.010  0.097 -0.046  0.018  0.048 -0.049  0.025  0.116  0.004  0.012 -0.370 -0.017      0.038     -0.073  0.003  -0.002            -0.025            -0.098            -0.007             0.027           -0.138            -0.005 -0.068  0.017            -0.013                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.078  0.032  0.082 -0.031  0.005  0.028  0.036  0.080  0.020 -0.022 -0.026  0.005 -0.105  0.157 -0.009 -0.041      0.041      0.002 -0.004  -0.058            -0.017            -0.009             0.018             0.018            0.022             0.041  0.056  0.049            -0.123            0.000                                                               
    ## human_bowler.sc:REALM2Marine                      0.031 -0.043 -0.034 -0.013 -0.008  0.076 -0.203 -0.083 -0.035  0.024 -0.022 -0.046 -0.105  0.004 -0.015  0.050      0.033      0.040 -0.056   0.165             0.028             0.062            -0.061             0.065            0.034             0.144  0.156  0.115            -0.004            0.024            0.006                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.003 -0.005  0.002  0.017 -0.048 -0.074  0.136  0.014 -0.015  0.076  0.003  0.030 -0.001 -0.008 -0.052 -0.019     -0.024      0.226  0.093  -0.264            -0.060             0.040            -0.136             0.006           -0.021            -0.525 -0.279  0.028            -0.033            0.038           -0.013      -0.127                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.040  0.104  0.037  0.040  0.092 -0.111  0.146  0.070 -0.011  0.069  0.016 -0.002  0.011  0.021 -0.080 -0.073     -0.099     -0.271  0.154  -0.235            -0.143             0.063            -0.147            -0.016            0.059            -0.247 -0.556  0.031            -0.022            0.069           -0.031      -0.155       0.499                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.078 -0.225 -0.085  0.049 -0.003 -0.026 -0.122 -0.056 -0.010  0.022  0.017 -0.001  0.127 -0.195  0.006  0.260     -0.177      0.039  0.034   0.332             0.054             0.000            -0.030            -0.023           -0.012            -0.043 -0.041 -0.219             0.557           -0.019           -0.537       0.003       0.075  0.069               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.027  0.078  0.029  0.021  0.024 -0.046  0.114  0.042  0.039 -0.041  0.096  0.014  0.096  0.006  0.000 -0.079     -0.067     -0.057  0.110  -0.248            -0.063            -0.103             0.108            -0.199           -0.031            -0.088 -0.105 -0.222             0.022           -0.002           -0.001      -0.607       0.169  0.212 -0.011        
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -1.4694105498 -0.0264616172  0.0002205077  0.0558747989  1.6792999453 
    ## 
    ## Number of Observations: 19658
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    177                  19658

``` r
rsquared(modTfullJtuexp)
```

    ##   Response   family     link method    Marginal  Conditional
    ## 1   Jtuexp gaussian identity   none 2.10169e-05 0.0001591993

``` r
rsquared(modTfullJbetaexp)
```

    ##   Response   family     link method   Marginal Conditional
    ## 1 Jbetaexp gaussian identity   none 0.01318546  0.08277463

``` r
rsquared(modTfullHornexp)
```

    ##   Response   family     link method     Marginal Conditional
    ## 1  Hornexp gaussian identity   none 0.0002989413 0.002182572

#### Summary mm

``` r
summary(modTfullJtumm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   91574.83 91924.67 -45745.41
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.6868309 (Intr)
    ## temptrend_abs.sc 0.1719101 -0.973
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:    1.067187 0.2791371
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -6.549746 
    ## Fixed effects: log(Jtumm + 1) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.9239176 0.29418448 30432   3.140606  0.0017
    ## temptrend_abs.sc                                  0.1382240 0.13927871 30432   0.992427  0.3210
    ## REALMMarine                                       0.0201198 0.31284521   195   0.064312  0.9488
    ## REALMTerrestrial                                  1.3628582 0.29028258   195   4.694936  0.0000
    ## tsign1                                           -0.0078633 0.01966722 30432  -0.399816  0.6893
    ## tempave_metab.sc                                 -0.3437056 0.03177561 30432 -10.816650  0.0000
    ## seas.sc                                          -0.1390088 0.02070902 30432  -6.712474  0.0000
    ## microclim.sc                                     -0.0240373 0.01083890 30432  -2.217691  0.0266
    ## mass.sc                                           0.0390995 0.02051577 30432   1.905829  0.0567
    ## speed.sc                                         -0.0258566 0.02193738 30432  -1.178655  0.2385
    ## consumerfrac.sc                                   0.0129393 0.01377460 30432   0.939358  0.3476
    ## nspp.sc                                          -0.3769968 0.01671560 30432 -22.553588  0.0000
    ## npp.sc                                            0.0247409 0.01554833 30432   1.591227  0.1116
    ## veg.sc                                           -0.0364737 0.04347056 30432  -0.839043  0.4015
    ## duration.sc                                       0.2672078 0.01267444 30432  21.082416  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.0893644 0.14131296 30432  -0.632387  0.5271
    ## temptrend_abs.sc:REALMTerrestrial                -0.3832214 0.14832038 30432  -2.583741  0.0098
    ## temptrend_abs.sc:tsign1                          -0.0008952 0.02240801 30432  -0.039952  0.9681
    ## temptrend_abs.sc:tempave_metab.sc                 0.0552764 0.02265810 30432   2.439586  0.0147
    ## temptrend_abs.sc:seas.sc                          0.0356042 0.01698093 30432   2.096719  0.0360
    ## temptrend_abs.sc:microclim.sc                     0.0085779 0.01167868 30432   0.734489  0.4627
    ## temptrend_abs.sc:mass.sc                         -0.0171547 0.01700950 30432  -1.008535  0.3132
    ## temptrend_abs.sc:speed.sc                         0.0206043 0.01971980 30432   1.044855  0.2961
    ## temptrend_abs.sc:consumerfrac.sc                 -0.0004379 0.01371439 30432  -0.031931  0.9745
    ## temptrend_abs.sc:nspp.sc                          0.0802689 0.01477945 30432   5.431115  0.0000
    ## tsign-1:thermal_bias.sc                          -0.0674955 0.02028837 30432  -3.326809  0.0009
    ## tsign1:thermal_bias.sc                           -0.0008059 0.01624664 30432  -0.049607  0.9604
    ## temptrend_abs.sc:npp.sc                          -0.0043218 0.01675424 30432  -0.257953  0.7964
    ## temptrend_abs.sc:veg.sc                          -0.0012097 0.02333030 30432  -0.051850  0.9586
    ## temptrend_abs.sc:duration.sc                      0.0197257 0.01679616 30432   1.174419  0.2402
    ## human_bowler.sc:REALM2TerrFresh                   0.0333704 0.03556767 30432   0.938224  0.3481
    ## human_bowler.sc:REALM2Marine                      0.0001579 0.01271139 30432   0.012425  0.9901
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0138149 0.01559689 30432   0.885747  0.3758
    ## temptrend_abs.sc:tsign1:thermal_bias.sc           0.0047364 0.01481990 30432   0.319597  0.7493
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0001052 0.01451059 30432  -0.007248  0.9942
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0097784 0.01649897 30432  -0.592665  0.5534
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.627                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.954  0.592                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.796  0.505  0.745                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.046  0.061 -0.002  0.022                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.082 -0.033 -0.070 -0.151 -0.045                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.093  0.070  0.140 -0.055 -0.082  0.052                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.058  0.065  0.066  0.004 -0.008 -0.147  0.108                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.072 -0.015 -0.038 -0.042  0.013  0.139  0.036  0.004                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.053 -0.055 -0.019 -0.034 -0.056 -0.061  0.011  0.064 -0.518                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.015  0.025  0.011  0.078  0.056 -0.102 -0.062  0.015 -0.015 -0.063                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.015 -0.030 -0.037 -0.056  0.035 -0.172  0.019 -0.028 -0.041  0.139  0.087                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.055 -0.023 -0.072  0.053 -0.013 -0.040 -0.334 -0.279 -0.019  0.093 -0.012 -0.146                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.459  0.265  0.481  0.011 -0.014  0.010  0.127  0.040  0.008 -0.008  0.001  0.021 -0.184                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.002 -0.003  0.023  0.042 -0.191  0.023 -0.062 -0.011 -0.042  0.022 -0.015 -0.305  0.033 -0.012                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.622 -0.954 -0.644 -0.494 -0.019  0.049 -0.116 -0.070 -0.007  0.030 -0.014  0.045  0.048 -0.296 -0.026                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.460 -0.689 -0.429 -0.625 -0.030  0.110  0.107  0.009 -0.021  0.067 -0.049  0.099 -0.089  0.037 -0.044  0.660                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.022 -0.156 -0.001 -0.008 -0.532  0.007 -0.005 -0.019  0.011  0.027 -0.032  0.003 -0.010  0.027  0.102  0.063      0.042                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.033  0.030  0.024  0.082 -0.024 -0.507 -0.063  0.085 -0.150  0.152  0.112  0.065 -0.008  0.000  0.012 -0.016     -0.203      0.030                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.055 -0.094 -0.071  0.074  0.041 -0.005 -0.540 -0.108 -0.037  0.004  0.021 -0.010  0.162 -0.121  0.001  0.134     -0.323      0.002  0.053                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.037 -0.091 -0.037 -0.004 -0.021  0.089 -0.083 -0.551 -0.010 -0.034 -0.019  0.018  0.114 -0.025 -0.022  0.088     -0.052      0.074 -0.078   0.209                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.010 -0.034 -0.002 -0.015  0.005 -0.100 -0.039 -0.003 -0.594  0.367  0.052  0.023 -0.008 -0.018  0.073  0.060      0.115     -0.061  0.351   0.042            -0.018                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.022  0.077  0.016  0.044  0.043  0.093  0.017 -0.042  0.354 -0.559 -0.008 -0.083 -0.026 -0.003 -0.031 -0.076     -0.165     -0.038 -0.479  -0.001             0.065            -0.686                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                  0.003 -0.021 -0.007 -0.020 -0.030  0.075  0.032 -0.012  0.032 -0.005 -0.617 -0.069 -0.009  0.002  0.014  0.022      0.045      0.067 -0.212  -0.031             0.007            -0.114             0.048                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.013  0.100  0.013  0.063  0.006  0.048 -0.013  0.020  0.007 -0.083 -0.071 -0.522  0.085 -0.043  0.117 -0.089     -0.256     -0.060  0.013  -0.004             0.004             0.027             0.142             0.169                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.042 -0.023 -0.040 -0.055 -0.083  0.362 -0.142 -0.074  0.046 -0.075 -0.016 -0.073 -0.035  0.007  0.046  0.032      0.031     -0.015 -0.116   0.125             0.038            -0.011             0.036             0.008            0.026                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.068 -0.073 -0.072 -0.085  0.043  0.502 -0.199 -0.198  0.061 -0.094 -0.022 -0.071 -0.071 -0.004  0.047  0.064      0.054      0.080 -0.177   0.160             0.091            -0.015             0.048             0.018            0.004             0.441                                                                                                                
    ## temptrend_abs.sc:npp.sc                          -0.005 -0.013  0.009 -0.055  0.007  0.009  0.121  0.107 -0.005 -0.024 -0.013  0.060 -0.538  0.109 -0.008 -0.024      0.203     -0.013  0.067  -0.213            -0.193             0.065            -0.027             0.048           -0.075             0.010  0.027                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.237 -0.419 -0.250  0.022  0.002 -0.011 -0.149 -0.016 -0.014  0.001  0.011 -0.037  0.214 -0.550 -0.004  0.475     -0.143      0.028 -0.002   0.309             0.009             0.011             0.020            -0.028            0.060             0.009  0.011 -0.390                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.024  0.187 -0.009 -0.022  0.004  0.103 -0.040  0.004  0.043 -0.042  0.039  0.112  0.007  0.015 -0.372  0.002      0.040     -0.056 -0.027   0.014            -0.001            -0.099             0.006             0.003           -0.118             0.009 -0.050 -0.019             0.000                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.075  0.049  0.078 -0.029  0.010  0.027  0.042  0.051  0.013 -0.014 -0.021 -0.005 -0.084  0.155 -0.005 -0.056      0.026     -0.005 -0.009  -0.058             0.017             0.001             0.011             0.014            0.028             0.027  0.043  0.030            -0.117           -0.005                                                               
    ## human_bowler.sc:REALM2Marine                      0.031 -0.040 -0.034 -0.012 -0.009  0.082 -0.209 -0.079 -0.022  0.025 -0.025 -0.045 -0.101  0.001 -0.019  0.048      0.030      0.028 -0.070   0.178             0.030             0.046            -0.051             0.084            0.024             0.145  0.152  0.113            -0.004            0.035           -0.001                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.013  0.009  0.014  0.023 -0.040 -0.132  0.135  0.025 -0.015  0.061  0.006  0.029  0.000  0.000 -0.042 -0.035     -0.032      0.236  0.100  -0.264            -0.068             0.025            -0.129             0.003            0.011            -0.569 -0.308  0.021            -0.042            0.016            0.002      -0.139                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.039  0.090  0.037  0.040  0.081 -0.171  0.142  0.077 -0.015  0.057  0.020  0.010  0.021  0.017 -0.068 -0.069     -0.091     -0.263  0.162  -0.232            -0.141             0.041            -0.125            -0.036            0.079            -0.281 -0.591  0.015            -0.023            0.041           -0.018      -0.158       0.493                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.072 -0.186 -0.076  0.039 -0.001 -0.025 -0.115 -0.044  0.000  0.014  0.015  0.003  0.118 -0.182  0.004  0.210     -0.141      0.033  0.025   0.282             0.038            -0.023            -0.017            -0.023           -0.020            -0.036 -0.029 -0.200             0.499           -0.007           -0.588       0.006       0.053  0.042               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.025  0.063  0.027  0.014  0.021 -0.053  0.118  0.038  0.026 -0.035  0.086  0.017  0.093  0.007  0.002 -0.066     -0.054     -0.033  0.127  -0.272            -0.057            -0.075             0.082            -0.236           -0.020            -0.098 -0.105 -0.239             0.029           -0.018            0.009      -0.613       0.198  0.228 -0.010        
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -2.487398e-02 -6.346473e-06 -4.382351e-08  1.861011e-06  2.006661e-02 
    ## 
    ## Number of Observations: 30663
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30663

``` r
summary(modTfullJbetamm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##        AIC      BIC    logLik
    ##   64704.11 65053.98 -32310.05
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev    Corr  
    ## (Intercept)      0.7623677 (Intr)
    ## temptrend_abs.sc 0.1985801 -0.991
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept)  Residual
    ## StdDev:   0.3703729 0.8132841
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##      power 
    ## -0.1989195 
    ## Fixed effects: log(Jbetamm + 1) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value  Std.Error    DF    t-value p-value
    ## (Intercept)                                       0.7494244 0.24450129 30454   3.065114  0.0022
    ## temptrend_abs.sc                                  0.0175761 0.10493144 30454   0.167501  0.8670
    ## REALMMarine                                      -0.2576908 0.26224215   195  -0.982644  0.3270
    ## REALMTerrestrial                                  0.8679334 0.25369149   195   3.421216  0.0008
    ## tsign1                                           -0.0393899 0.01286322 30454  -3.062214  0.0022
    ## tempave_metab.sc                                 -0.3375189 0.02212634 30454 -15.254169  0.0000
    ## seas.sc                                          -0.0937220 0.01359616 30454  -6.893271  0.0000
    ## microclim.sc                                     -0.0233472 0.00696381 30454  -3.352656  0.0008
    ## mass.sc                                           0.0032908 0.01433759 30454   0.229526  0.8185
    ## speed.sc                                         -0.0451568 0.01497826 30454  -3.014825  0.0026
    ## consumerfrac.sc                                  -0.0178069 0.00880022 30454  -2.023459  0.0430
    ## nspp.sc                                          -0.2524030 0.01111210 30454 -22.714247  0.0000
    ## npp.sc                                            0.0080889 0.01009398 30454   0.801357  0.4229
    ## veg.sc                                            0.0146100 0.02686187 30454   0.543895  0.5865
    ## duration.sc                                       0.1763543 0.00849750 30454  20.753669  0.0000
    ## temptrend_abs.sc:REALMMarine                      0.0194624 0.10757739 30454   0.180915  0.8564
    ## temptrend_abs.sc:REALMTerrestrial                -0.2662536 0.11301977 30454  -2.355814  0.0185
    ## temptrend_abs.sc:tsign1                           0.0141858 0.01537099 30454   0.922894  0.3561
    ## temptrend_abs.sc:tempave_metab.sc                 0.0770503 0.01564641 30454   4.924467  0.0000
    ## temptrend_abs.sc:seas.sc                          0.0218798 0.01175680 30454   1.861033  0.0627
    ## temptrend_abs.sc:microclim.sc                    -0.0012926 0.00803691 30454  -0.160839  0.8722
    ## temptrend_abs.sc:mass.sc                         -0.0040027 0.01173697 30454  -0.341031  0.7331
    ## temptrend_abs.sc:speed.sc                         0.0253071 0.01366665 30454   1.851741  0.0641
    ## temptrend_abs.sc:consumerfrac.sc                  0.0139787 0.00930010 30454   1.503070  0.1328
    ## temptrend_abs.sc:nspp.sc                          0.0655740 0.01016116 30454   6.453396  0.0000
    ## tsign-1:thermal_bias.sc                          -0.0703259 0.01360417 30454  -5.169437  0.0000
    ## tsign1:thermal_bias.sc                           -0.0206176 0.01084826 30454  -1.900546  0.0574
    ## temptrend_abs.sc:npp.sc                           0.0058216 0.01149041 30454   0.506647  0.6124
    ## temptrend_abs.sc:veg.sc                          -0.0076368 0.01589240 30454  -0.480533  0.6309
    ## temptrend_abs.sc:duration.sc                      0.0016858 0.01124752 30454   0.149882  0.8809
    ## human_bowler.sc:REALM2TerrFresh                   0.0144183 0.02176449 30454   0.662471  0.5077
    ## human_bowler.sc:REALM2Marine                     -0.0041463 0.00820370 30454  -0.505418  0.6133
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc          0.0138646 0.01072001 30454   1.293336  0.1959
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.0016634 0.01011832 30454  -0.164395  0.8694
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0052916 0.00976545 30454  -0.541864  0.5879
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0053453 0.01125497 30454  -0.474932  0.6348
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.695                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.937  0.647                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.849  0.589  0.788                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.036  0.054 -0.003  0.015                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.086 -0.042 -0.067 -0.118 -0.045                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.088  0.071  0.125 -0.021 -0.089  0.036                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.046  0.058  0.050  0.005 -0.008 -0.162  0.116                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.086 -0.034 -0.050 -0.056  0.017  0.135  0.016  0.004                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.042 -0.050 -0.011 -0.021 -0.065  0.016  0.035  0.053 -0.525                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.005  0.018  0.006  0.052  0.050 -0.083 -0.050  0.015 -0.021 -0.027                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                          -0.002 -0.017 -0.019 -0.030  0.040 -0.192  0.026 -0.030 -0.063  0.139  0.095                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.050 -0.023 -0.062  0.031 -0.011 -0.037 -0.322 -0.265 -0.016  0.094 -0.019 -0.148                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.336  0.226  0.350  0.012 -0.018  0.004  0.122  0.042  0.011 -0.010 -0.005  0.020 -0.186                                                                                                                                                                                                                                                                                     
    ## duration.sc                                       0.002 -0.001  0.019  0.027 -0.192  0.038 -0.061 -0.013 -0.023  0.011 -0.022 -0.327  0.031 -0.007                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.677 -0.950 -0.714 -0.571 -0.016  0.052 -0.113 -0.061  0.010  0.025 -0.010  0.034  0.046 -0.251 -0.024                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.565 -0.732 -0.525 -0.696 -0.026  0.097  0.085  0.009 -0.006  0.059 -0.039  0.081 -0.075  0.030 -0.037  0.698                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.019 -0.144 -0.004 -0.007 -0.522  0.006 -0.004 -0.019  0.010  0.028 -0.032 -0.001 -0.012  0.025  0.100  0.058      0.039                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.036  0.034  0.025  0.065 -0.024 -0.513 -0.057  0.093 -0.141  0.115  0.111  0.080 -0.008  0.001  0.003 -0.020     -0.183      0.027                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.039 -0.079 -0.048  0.057  0.046  0.009 -0.535 -0.113 -0.023 -0.004  0.017 -0.014  0.154 -0.114  0.002  0.112     -0.300     -0.002  0.047                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.028 -0.081 -0.025 -0.003 -0.022  0.099 -0.083 -0.552 -0.009 -0.028 -0.016  0.019  0.104 -0.021 -0.022  0.077     -0.048      0.072 -0.084   0.209                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                         -0.018 -0.022  0.006 -0.009  0.001 -0.093 -0.026  0.000 -0.590  0.370  0.060  0.035 -0.013 -0.016  0.062  0.047      0.108     -0.060  0.346   0.025            -0.023                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.017  0.070  0.012  0.037  0.048  0.063  0.006 -0.039  0.348 -0.553 -0.025 -0.089 -0.025 -0.002 -0.024 -0.068     -0.158     -0.036 -0.466   0.008             0.064            -0.687                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.002 -0.015 -0.005 -0.007 -0.027  0.067  0.025 -0.009  0.032 -0.020 -0.622 -0.073 -0.008  0.006  0.017  0.017      0.034      0.072 -0.214  -0.035             0.000            -0.120             0.066                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.010  0.091  0.010  0.047  0.002  0.057 -0.018  0.021  0.015 -0.086 -0.076 -0.519  0.086 -0.040  0.127 -0.083     -0.233     -0.056  0.006  -0.003             0.003             0.023             0.150             0.176                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.046 -0.029 -0.041 -0.046 -0.090  0.386 -0.157 -0.084  0.050 -0.043 -0.013 -0.091 -0.035  0.004  0.057  0.036      0.030     -0.011 -0.131   0.130             0.044            -0.012             0.024             0.007            0.032                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.073 -0.077 -0.071 -0.072  0.041  0.532 -0.232 -0.212  0.068 -0.048 -0.016 -0.095 -0.069 -0.011  0.063  0.067      0.052      0.079 -0.197   0.174             0.098            -0.015             0.031             0.016            0.012             0.463                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.001 -0.016  0.002 -0.044  0.007  0.011  0.113  0.099 -0.008 -0.023 -0.014  0.059 -0.535  0.112 -0.005 -0.016      0.189     -0.013  0.063  -0.213            -0.189             0.071            -0.029             0.053           -0.078             0.011  0.031                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.165 -0.375 -0.172  0.018  0.003 -0.009 -0.142 -0.012 -0.011  0.001  0.013 -0.035  0.210 -0.529 -0.006  0.420     -0.130      0.031 -0.003   0.302             0.005             0.002             0.023            -0.031            0.056             0.010  0.014 -0.391                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.019  0.163 -0.009 -0.016  0.004  0.103 -0.047  0.003  0.035 -0.038  0.044  0.122  0.012  0.012 -0.370  0.003      0.037     -0.062 -0.027   0.016             0.003            -0.097             0.003             0.002           -0.134             0.012 -0.049 -0.022             0.002                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.068  0.058  0.069 -0.004  0.012  0.024  0.036  0.065  0.010 -0.013 -0.016 -0.003 -0.090  0.145 -0.012 -0.063      0.013     -0.007 -0.011  -0.065             0.011             0.000             0.010             0.012            0.028             0.025  0.049  0.032            -0.118           -0.003                                                               
    ## human_bowler.sc:REALM2Marine                      0.028 -0.037 -0.030 -0.013 -0.010  0.091 -0.202 -0.075 -0.017  0.031 -0.026 -0.043 -0.109  0.003 -0.014  0.045      0.030      0.029 -0.072   0.173             0.025             0.046            -0.055             0.090            0.022             0.145  0.155  0.118            -0.007            0.039            0.001                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.016  0.012  0.014  0.018 -0.038 -0.145  0.143  0.030 -0.020  0.045  0.005  0.038  0.000  0.002 -0.046 -0.034     -0.028      0.236  0.107  -0.268            -0.073             0.028            -0.123             0.007            0.008            -0.562 -0.315  0.020            -0.044            0.015            0.000      -0.138                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.038  0.089  0.034  0.033  0.080 -0.183  0.157  0.082 -0.018  0.037  0.017  0.020  0.023  0.020 -0.077 -0.068     -0.084     -0.262  0.170  -0.242            -0.149             0.039            -0.118            -0.034            0.081            -0.290 -0.586  0.011            -0.025            0.046           -0.032      -0.160       0.497                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.056 -0.179 -0.057  0.025 -0.001 -0.024 -0.117 -0.052  0.000  0.014  0.013  0.003  0.120 -0.181  0.010  0.199     -0.131      0.035  0.029   0.299             0.049            -0.022            -0.018            -0.026           -0.022            -0.034 -0.033 -0.205             0.519           -0.014           -0.553       0.004       0.056  0.056               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.019  0.056  0.022  0.010  0.020 -0.056  0.117  0.033  0.026 -0.039  0.093  0.016  0.098  0.003  0.001 -0.059     -0.048     -0.033  0.124  -0.266            -0.051            -0.075             0.084            -0.249           -0.019            -0.098 -0.106 -0.243             0.034           -0.022            0.009      -0.615       0.195  0.229 -0.009        
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -4.8088209 -0.5419471 -0.2390022  0.3379120  5.0650986 
    ## 
    ## Number of Observations: 30685
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30685

``` r
summary(modTfullHornmm)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i, ] 
    ##       AIC      BIC    logLik
    ##   96950.5 97300.31 -48433.25
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.889280364 (Intr)
    ## temptrend_abs.sc 0.007292481 0.002 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:       1.096 2.659822
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.291161 
    ## Fixed effects: log(Hornmm + 1) ~ temptrend_abs.sc * REALM + temptrend_abs.sc *      tsign + temptrend_abs.sc * tempave_metab.sc + temptrend_abs.sc *      seas.sc + temptrend_abs.sc * microclim.sc + temptrend_abs.sc *      mass.sc + temptrend_abs.sc * speed.sc + temptrend_abs.sc *      consumerfrac.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc:REALM2 
    ##                                                       Value Std.Error    DF    t-value p-value
    ## (Intercept)                                       1.2781750 0.3303234 30411   3.869466  0.0001
    ## temptrend_abs.sc                                  0.0174368 0.1386050 30411   0.125802  0.8999
    ## REALMMarine                                      -0.2866093 0.3529802   195  -0.811970  0.4178
    ## REALMTerrestrial                                  1.4770173 0.3280320   195   4.502662  0.0000
    ## tsign1                                           -0.0370111 0.0213953 30411  -1.729867  0.0837
    ## tempave_metab.sc                                 -0.6216271 0.0334620 30411 -18.577079  0.0000
    ## seas.sc                                          -0.2155078 0.0220007 30411  -9.795495  0.0000
    ## microclim.sc                                     -0.0704297 0.0114165 30411  -6.169127  0.0000
    ## mass.sc                                          -0.0360424 0.0224308 30411  -1.606830  0.1081
    ## speed.sc                                          0.1629223 0.0236880 30411   6.877852  0.0000
    ## consumerfrac.sc                                  -0.0631736 0.0141504 30411  -4.464432  0.0000
    ## nspp.sc                                          -0.2233165 0.0178444 30411 -12.514654  0.0000
    ## npp.sc                                           -0.0421054 0.0165484 30411  -2.544382  0.0110
    ## veg.sc                                            0.0052142 0.0469732 30411   0.111003  0.9116
    ## duration.sc                                       0.1879722 0.0137411 30411  13.679522  0.0000
    ## temptrend_abs.sc:REALMMarine                     -0.0442971 0.1385629 30411  -0.319689  0.7492
    ## temptrend_abs.sc:REALMTerrestrial                -0.0210141 0.1428863 30411  -0.147068  0.8831
    ## temptrend_abs.sc:tsign1                           0.0024881 0.0255172 30411   0.097506  0.9223
    ## temptrend_abs.sc:tempave_metab.sc                 0.0638417 0.0214784 30411   2.972362  0.0030
    ## temptrend_abs.sc:seas.sc                          0.0065882 0.0185142 30411   0.355848  0.7220
    ## temptrend_abs.sc:microclim.sc                     0.0265762 0.0122523 30411   2.169079  0.0301
    ## temptrend_abs.sc:mass.sc                          0.0116003 0.0186032 30411   0.623564  0.5329
    ## temptrend_abs.sc:speed.sc                        -0.0334077 0.0219500 30411  -1.521986  0.1280
    ## temptrend_abs.sc:consumerfrac.sc                  0.0182946 0.0144096 30411   1.269606  0.2042
    ## temptrend_abs.sc:nspp.sc                          0.0465911 0.0159566 30411   2.919873  0.0035
    ## tsign-1:thermal_bias.sc                          -0.0470859 0.0219808 30411  -2.142137  0.0322
    ## tsign1:thermal_bias.sc                           -0.0305485 0.0174862 30411  -1.747006  0.0806
    ## temptrend_abs.sc:npp.sc                           0.0237002 0.0179642 30411   1.319300  0.1871
    ## temptrend_abs.sc:veg.sc                          -0.0616969 0.0261855 30411  -2.356147  0.0185
    ## temptrend_abs.sc:duration.sc                     -0.0169735 0.0183909 30411  -0.922932  0.3561
    ## human_bowler.sc:REALM2TerrFresh                   0.1726123 0.0374020 30411   4.615059  0.0000
    ## human_bowler.sc:REALM2Marine                     -0.0009122 0.0135935 30411  -0.067106  0.9465
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.0154754 0.0164796 30411  -0.939068  0.3477
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.0060095 0.0155670 30411  -0.386044  0.6995
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh -0.0655236 0.0155523 30411  -4.213115  0.0000
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.0256999 0.0183057 30411  -1.403927  0.1604
    ##  Correlation: 
    ##                                                  (Intr) tmpt_. REALMM REALMT tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s t_.:REALMM t_.:REALMT tm_.:1 tm_.:_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. h_.:REALM2T h_.:REALM2M t_.:-1 t_.:1: t_.:_.:REALM2T
    ## temptrend_abs.sc                                 -0.267                                                                                                                                                                                                                                                                                                                                                                         
    ## REALMMarine                                      -0.948  0.257                                                                                                                                                                                                                                                                                                                                                                  
    ## REALMTerrestrial                                 -0.810  0.150  0.754                                                                                                                                                                                                                                                                                                                                                           
    ## tsign1                                           -0.046  0.063 -0.002  0.019                                                                                                                                                                                                                                                                                                                                                    
    ## tempave_metab.sc                                  0.083 -0.002 -0.069 -0.137 -0.044                                                                                                                                                                                                                                                                                                                                             
    ## seas.sc                                          -0.087  0.011  0.132 -0.051 -0.075  0.047                                                                                                                                                                                                                                                                                                                                      
    ## microclim.sc                                     -0.058  0.087  0.067 -0.002 -0.008 -0.150  0.128                                                                                                                                                                                                                                                                                                                               
    ## mass.sc                                           0.077  0.020 -0.044 -0.042  0.009  0.151  0.041  0.012                                                                                                                                                                                                                                                                                                                        
    ## speed.sc                                          0.048 -0.038 -0.019 -0.035 -0.053 -0.066  0.003  0.059 -0.539                                                                                                                                                                                                                                                                                                                 
    ## consumerfrac.sc                                  -0.008  0.025  0.001  0.064  0.078 -0.109 -0.081  0.003 -0.020 -0.026                                                                                                                                                                                                                                                                                                          
    ## nspp.sc                                           0.016 -0.063 -0.039 -0.062  0.055 -0.166  0.015 -0.022 -0.039  0.146  0.129                                                                                                                                                                                                                                                                                                   
    ## npp.sc                                            0.051  0.012 -0.066  0.060 -0.022 -0.041 -0.348 -0.250 -0.008  0.088  0.006 -0.167                                                                                                                                                                                                                                                                                            
    ## veg.sc                                           -0.437  0.285  0.456  0.012 -0.007  0.005  0.112  0.044  0.005 -0.009 -0.006  0.014 -0.179                                                                                                                                                                                                                                                                                     
    ## duration.sc                                      -0.001  0.015  0.022  0.038 -0.196  0.040 -0.066 -0.025 -0.051  0.020 -0.046 -0.294  0.037 -0.005                                                                                                                                                                                                                                                                              
    ## temptrend_abs.sc:REALMMarine                      0.274 -0.963 -0.275 -0.149 -0.020  0.033 -0.039 -0.090 -0.025  0.034 -0.008  0.065  0.007 -0.320 -0.030                                                                                                                                                                                                                                                                       
    ## temptrend_abs.sc:REALMTerrestrial                 0.119 -0.648 -0.106 -0.211 -0.023  0.043  0.171  0.037 -0.077  0.091 -0.021  0.161 -0.165  0.057 -0.051  0.624                                                                                                                                                                                                                                                                
    ## temptrend_abs.sc:tsign1                           0.024 -0.172 -0.001 -0.009 -0.508  0.006 -0.013 -0.028  0.027  0.010 -0.067 -0.001  0.000  0.018  0.090  0.064      0.048                                                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tempave_metab.sc                -0.003  0.017 -0.006  0.023 -0.044 -0.230 -0.093  0.043 -0.199  0.268  0.131 -0.002  0.023 -0.005 -0.007  0.017     -0.138      0.013                                                                                                                                                                                                                                          
    ## temptrend_abs.sc:seas.sc                          0.036 -0.066 -0.054  0.068  0.032  0.007 -0.453 -0.149 -0.038  0.012  0.041  0.000  0.161 -0.104 -0.005  0.102     -0.394     -0.007  0.074                                                                                                                                                                                                                                   
    ## temptrend_abs.sc:microclim.sc                     0.038 -0.135 -0.045  0.010 -0.016  0.091 -0.133 -0.492 -0.027 -0.029  0.000  0.002  0.045 -0.031 -0.008  0.134     -0.132      0.082 -0.034   0.348                                                                                                                                                                                                                           
    ## temptrend_abs.sc:mass.sc                          0.003 -0.038 -0.004 -0.031  0.023 -0.091 -0.059 -0.020 -0.502  0.349  0.051  0.002 -0.027 -0.012  0.089  0.056      0.149     -0.094  0.446   0.049             0.000                                                                                                                                                                                                         
    ## temptrend_abs.sc:speed.sc                        -0.009  0.070  0.018  0.044  0.024  0.134  0.057 -0.024  0.343 -0.466 -0.051 -0.087 -0.016 -0.002 -0.022 -0.090     -0.196     -0.002 -0.628  -0.018             0.049            -0.676                                                                                                                                                                                       
    ## temptrend_abs.sc:consumerfrac.sc                 -0.008 -0.016  0.013 -0.003 -0.068  0.085  0.074  0.007  0.033 -0.054 -0.541 -0.126 -0.051  0.016  0.056  0.008      0.044      0.125 -0.246  -0.086            -0.022            -0.110             0.087                                                                                                                                                                     
    ## temptrend_abs.sc:nspp.sc                         -0.022  0.122  0.020  0.074 -0.038  0.004 -0.002  0.005 -0.038 -0.067 -0.130 -0.426  0.119 -0.036  0.063 -0.099     -0.318     -0.037  0.063  -0.020             0.041             0.089             0.125             0.272                                                                                                                                                   
    ## tsign-1:thermal_bias.sc                           0.042 -0.004 -0.039 -0.054 -0.079  0.347 -0.149 -0.058  0.032 -0.068 -0.016 -0.092 -0.029 -0.002  0.072  0.014      0.008     -0.026  0.037   0.117             0.014             0.029             0.025             0.009            0.050                                                                                                                                  
    ## tsign1:thermal_bias.sc                            0.068 -0.063 -0.069 -0.084  0.027  0.490 -0.206 -0.192  0.049 -0.076 -0.015 -0.054 -0.076 -0.001  0.060  0.049      0.031      0.107  0.021   0.149             0.082             0.030             0.022             0.005           -0.059             0.443                                                                                                                
    ## temptrend_abs.sc:npp.sc                           0.004 -0.023 -0.002 -0.073  0.026  0.026  0.119  0.024 -0.018 -0.016 -0.056  0.083 -0.480  0.107 -0.030 -0.008      0.330     -0.041 -0.046  -0.264            -0.097             0.079            -0.013             0.124           -0.136             0.019  0.050                                                                                                         
    ## temptrend_abs.sc:veg.sc                           0.221 -0.481 -0.228  0.012 -0.012  0.000 -0.087 -0.039  0.011 -0.005  0.020 -0.034  0.200 -0.536 -0.011  0.542     -0.155      0.036  0.015   0.244             0.071            -0.018             0.032            -0.056            0.065             0.017  0.008 -0.414                                                                                                  
    ## temptrend_abs.sc:duration.sc                     -0.022  0.205 -0.005 -0.022  0.003  0.108 -0.034  0.015  0.062 -0.043  0.058  0.084  0.004  0.010 -0.326 -0.005      0.045     -0.052  0.010   0.002            -0.014            -0.123            -0.004            -0.013           -0.087            -0.015 -0.057  0.009            -0.007                                                                                
    ## human_bowler.sc:REALM2TerrFresh                  -0.070 -0.015  0.073 -0.010  0.002  0.017  0.006  0.115  0.002 -0.003 -0.017  0.003 -0.112  0.144 -0.013  0.014      0.027      0.012  0.017   0.017            -0.108             0.026            -0.008             0.015            0.007             0.040  0.039  0.089            -0.053            0.014                                                               
    ## human_bowler.sc:REALM2Marine                      0.029 -0.029 -0.031 -0.019 -0.012  0.087 -0.205 -0.083 -0.017  0.022 -0.061 -0.031 -0.135  0.010 -0.025  0.035      0.040      0.024 -0.035   0.156             0.035             0.038            -0.060             0.141            0.017             0.148  0.161  0.167            -0.043            0.036            0.028                                              
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc         -0.002  0.010  0.004  0.018 -0.047 -0.025  0.102 -0.009  0.033  0.036 -0.004  0.021 -0.013  0.006 -0.084 -0.037     -0.033      0.257 -0.040  -0.250            -0.064            -0.079            -0.059             0.025            0.018            -0.506 -0.246 -0.007            -0.042            0.074           -0.013      -0.145                                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc          -0.029  0.125  0.023  0.034  0.121 -0.049  0.123  0.030  0.035  0.020  0.022 -0.051  0.019  0.013 -0.088 -0.090     -0.112     -0.336  0.024  -0.237            -0.119            -0.066            -0.044            -0.036            0.191            -0.232 -0.518 -0.035            -0.029            0.071            0.004      -0.176       0.444                      
    ## temptrend_abs.sc:human_bowler.sc:REALM2TerrFresh  0.065 -0.189 -0.070  0.019 -0.003  0.001 -0.059 -0.118  0.017 -0.002  0.020 -0.004  0.136 -0.179 -0.001  0.214     -0.143      0.019  0.000   0.224             0.195            -0.064             0.006            -0.046           -0.003            -0.035 -0.023 -0.276             0.481           -0.022           -0.510      -0.034       0.062  0.020               
    ## temptrend_abs.sc:human_bowler.sc:REALM2Marine    -0.019  0.059  0.020  0.021  0.022 -0.047  0.109  0.039  0.015 -0.026  0.153 -0.005  0.140 -0.009  0.023 -0.053     -0.076     -0.023  0.090  -0.245            -0.078            -0.044             0.092            -0.333           -0.014            -0.102 -0.114 -0.315             0.086           -0.019           -0.037      -0.587       0.221  0.268  0.048        
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -0.9459229 -0.1775189 -0.0596009  0.1054110  2.5803660 
    ## 
    ## Number of Observations: 30642
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                    198                  30642

``` r
rsquared(modTfullJtumm)
```

    ##   Response   family     link method   Marginal Conditional
    ## 1    Jtumm gaussian identity   none 0.08893255   0.9586234

``` r
rsquared(modTfullJbetamm)
```

    ##   Response   family     link method   Marginal Conditional
    ## 1  Jbetamm gaussian identity   none 0.06897058   0.5652316

``` r
rsquared(modTfullHornmm)
```

    ##   Response   family     link method   Marginal Conditional
    ## 1   Hornmm gaussian identity   none 0.02360949   0.2381358

### Plots from the full models

#### Plot the coefficients from rem0

``` r
coefs1 <- summary(modTfullJturem0)$tTable
coefs2 <- summary(modTfullJbetarem0)$tTable
coefs3 <- summary(modTfullHornrem0)$tTable

varstoplot <- unique(c(rownames(coefs1), rownames(coefs2), rownames(coefs3)))
varstoplot <- varstoplot[which(!grepl('Intercept', varstoplot) | grepl(':', varstoplot))] # vars to plot

rows1_1 <- which(rownames(coefs1) %in% varstoplot) # rows in coefs
rows1_2 <- which(rownames(coefs2) %in% varstoplot)
rows1_3 <- which(rownames(coefs3) %in% varstoplot)
xlims <- range(c(coefs1[rows1_1,1] - coefs1[rows1_1,2], coefs1[rows1_1,1] + coefs1[rows1_1,2], 
                  coefs2[rows1_2,1] - coefs2[rows1_2,2], coefs2[rows1_2,1] + coefs2[rows1_2,2], 
                  coefs3[rows1_3,1] - coefs3[rows1_3,2], coefs3[rows1_3,1] + coefs3[rows1_3,2]))


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
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[1], length(varstoplot) + 1 - i + offs[1]), col = cols[1])
  }
  if(varstoplot[i] %in% rownames(coefs2)){
    x = coefs2[rownames(coefs2) == varstoplot[i], 1]
    se = coefs2[rownames(coefs2) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[2], pch = pchs[2], col = cols[2])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[2], length(varstoplot) + 1 - i + offs[2]), col = cols[2])
  }
  if(varstoplot[i] %in% rownames(coefs3)){
    x = coefs3[rownames(coefs3) == varstoplot[i], 1]
    se = coefs3[rownames(coefs3) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[3], pch = pchs[3], col = cols[3])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[3], length(varstoplot) + 1 - i + offs[3]), col = cols[3])
  }
}
legend('topleft', col = cols, pch = 16, lwd = 1, legend = c('Jtu', 'Jbeta', 'Horn'), cex = 0.5)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20fullTmods-1.png)<!-- -->

#### Plot the main effects

``` r
# set up the variables to plot
# if variable is logged before scaling (see 'center and scale' above), then need to mark it here and express the limits on a log10 scale (even though log transforming is log)
vars <- data.frame(vars = c('temptrend_abs', 'tempave_metab', 'seas', 'microclim', 'mass', 'speed', 
                            'consumerfrac', 'nspp', 'thermal_bias', 'npp', 'veg', 'duration', 
                            'human_bowler', 'human_bowler'),
           min =      c(0,  0,   0.1, -2,  0,   0,   0,   0.3, -10, 1.9, 0,   0.5, 0,   0), 
           max =      c(0.2,  30,  16,  0.8, 8,   2,   1,   2.6, 10,  3.7, 0.3, 2,   1,   1),
           log =      c(F,  F,   F,   T,   T,   T,   F,   T,   F,   T,   T,   T,   T,   T),
           len =      c(100,  100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
           discrete = c(F,  F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F),
           plus =     c(0,  0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   1,   1), # what to add before log-scaling
           REALM = c(rep('Terrestrial', 13), 'Marine'),
           REALM2 = c(rep('TerrFresh', 13), 'Marine'),
           stringsAsFactors = FALSE)
baseall <- trends[, .(type = 'all', 
                      temptrend = -0.0001,
                      temptrend_abs.sc = 0.0001,
                      tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                      seas.sc = mean(seas.sc, na.rm=TRUE), 
                      microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                      speed.sc = mean(speed.sc, na.rm=TRUE), 
                      mass.sc = mean(mass.sc, na.rm=TRUE), 
                      nspp.sc = 0, 
                      thermal_bias.sc = mean(thermal_bias.sc, na.rm=TRUE), 
                      npp.sc = mean(npp.sc, na.rm=TRUE), 
                      human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                      veg.sc = mean(veg.sc, na.rm=TRUE), 
                      consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
baseterr <- trends[REALM == 'Terrestrial', 
                   .(type = 'Terrestrial', 
                     temptrend = -0.0001,
                     temptrend_abs.sc = 0.0001,
                     tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                     seas.sc = mean(seas.sc, na.rm=TRUE), 
                     microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                     speed.sc = mean(speed.sc, na.rm=TRUE), 
                     mass.sc = mean(mass.sc, na.rm=TRUE), 
                     nspp.sc = 0, 
                     thermal_bias.sc = 0, 
                     npp.sc = mean(npp.sc, na.rm=TRUE), 
                     human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                     veg.sc = mean(veg.sc, na.rm=TRUE), 
                     consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
basemar <- trends[REALM == 'Marine', 
                  .(type = 'Marine',
                    temptrend = -0.0001,
                    temptrend_abs.sc = 0.0001,
                    tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                    seas.sc = mean(seas.sc, na.rm=TRUE), 
                    microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                    speed.sc = mean(speed.sc, na.rm=TRUE), 
                    mass.sc = mean(mass.sc, na.rm=TRUE), 
                    nspp.sc = 0, 
                    thermal_bias.sc = 0, 
                    npp.sc = mean(npp.sc, na.rm=TRUE), 
                    human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                    veg.sc = mean(veg.sc, na.rm=TRUE), 
                    consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
basetab <- rbind(baseall, baseterr, basemar)
basetab[, ':='(duration.sc = 0, nyrBT = 20, STUDY_ID = 127L, rarefyID = '127_514668')]

# make the data frames for each interaction to plot                
for(j in 1:nrow(vars)){
  # set up the main effects
  if(vars$log[j]){
    thisdat <- data.frame(new = 10^seq(vars$min[j], vars$max[j], length.out = vars$len[j]),
                          var = vars$vars[j], stringsAsFactors = FALSE)
  } 
  if(!vars$log[j]){
    thisdat <- data.frame(new = seq(vars$min[j], vars$max[j], length.out = vars$len[j]),
                          var = vars$vars[j], stringsAsFactors = FALSE)
  }
  names(thisdat) <- c(vars$vars[j], 'var')

  # scale the variable
  cent <- attr(trends[[paste0(vars$vars[j], '.sc')]], 'scaled:center')
  scl <- attr(trends[[paste0(vars$vars[j], '.sc')]], 'scaled:scale')
  if(is.null(cent)) cent <- 0
  if(!is.null(cent) & !is.null(scl)){
    if(vars$log[j]) thisdat[[paste0(vars$var[j], '.sc')]] <- (log(thisdat[[vars$vars[j]]] + vars$plus[j]) - cent)/scl
    if(!vars$log[j]) thisdat[[paste0(vars$var[j], '.sc')]] <- (thisdat[[vars$var[j]]] - cent)/scl
  }

  # merge with the rest of the columns
  # use realm-specific averages for human impacts
  if(vars$vars[j] != 'tsign') colnamestouse <- setdiff(colnames(basetab), paste0(vars$vars[j], '.sc'))
  if(vars$vars[j] == 'tsign') colnamestouse <- setdiff(colnames(basetab), vars$var[j])
  if(vars$vars[j] != 'human_bowler'){
    thisdat <- cbind(thisdat, basetab[type == 'all', ..colnamestouse])
  }
  if(vars$vars[j] == 'human_bowler' & vars$REALM[j] == 'Terrestrial'){
    thisdat <- cbind(thisdat, basetab[type == 'Terrestrial', ..colnamestouse])
  }
  if(vars$vars[j] == 'human_bowler' & vars$REALM[j] == 'Marine'){
    thisdat <- cbind(thisdat, basetab[type == 'Marine', ..colnamestouse])
  }
  
  # add realm
  thisdat$REALM <- vars$REALM[j]
  thisdat$REALM2 <- vars$REALM2[j]
  
  # merge with the previous iterations
  if(j == 1) newdat <- thisdat
  if(j > 1){
    colstoadd <- setdiff(colnames(thisdat), colnames(newdat))
    for(toadd in colstoadd){
      newdat[[toadd]] <- NA
    }
    
    colstoadd2 <- setdiff(colnames(newdat), colnames(thisdat))
    for(toadd in colstoadd2){
      thisdat[[toadd]] <- NA
    }
    
    newdat <- rbind(newdat, thisdat)
  } 
}

# character so that new levels can be added
newdat$REALM <- as.character(newdat$REALM)
newdat$REALM2 <- as.character(newdat$REALM2)

# add extra rows so that all factor levels are represented (for predict.lme to work)
newdat <- rbind(newdat[1:6, ], newdat)
newdat$REALM[1:6] <- c('Marine', 'Marine', 'Freshwater', 'Freshwater', 'Terrestrial', 'Terrestrial')
newdat$REALM2[1:6] <- c('Marine', 'Marine', 'TerrFresh', 'TerrFresh', 'TerrFresh', 'TerrFresh')
newdat$temptrend_abs.sc[1:6] <- rep(0.0001, 6)
newdat$temptrend[1:6] <- rep(c(0.0001, -0.0001), 3)
newdat$var[1:6] <- 'test'

# trim to at least some temperature change (so that tsign is -1 or 1)
newdat <- newdat[newdat$temptrend_abs.sc != 0,]

# set up tsign
newdat$tsign <- factor(sign(newdat$temptrend))


# make predictions
newdat$predsJtu <- predict(object = modTfullJturem0, newdata = newdat, level = 0)
newdat$predsJbeta <- predict(object = modTfullJbetarem0, newdata = newdat, level = 0)
newdat$predsHorn <- predict(object = modTfullHornrem0, newdata = newdat, level = 0)

#compute standard error for predictions
# from https://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit
DesignmatJtu <- model.matrix(eval(eval(modTfullJturem0$call$fixed)[-2]), newdat)
DesignmatJbeta <- model.matrix(eval(eval(modTfullJbetarem0$call$fixed)[-2]), newdat)
DesignmatHorn <- model.matrix(eval(eval(modTfullHornrem0$call$fixed)[-2]), newdat)

predvarJtu <- diag(DesignmatJtu %*% modTfullJturem0$varFix %*% t(DesignmatJtu))
predvarJbeta <- diag(DesignmatJbeta %*% modTfullJbetarem0$varFix %*% t(DesignmatJbeta))
predvarHorn <- diag(DesignmatHorn %*% modTfullHornrem0$varFix %*% t(DesignmatHorn))

newdat$SE_Jtu <- sqrt(predvarJtu) 
newdat$SE_Jbeta <- sqrt(predvarJbeta) 
newdat$SE_Horn <- sqrt(predvarHorn) 

# prep the plots
varplots <- vector('list', nrow(vars))
for(j in 1:length(varplots)){
  subs <- newdat$var == vars$vars[j] # which rows of newdat
  xvar <- vars$vars[j]
  title <- vars$vars[j]
  if(vars$vars[j] %in% c('human_bowler')){
    subs <- newdat$var == vars$vars[j] & newdat$REALM2 == vars$REALM2[j]
    title <- paste0('human:', vars$REALM2[j])
  } 

  se <- 1
  thisplot <- ggplot(newdat[subs, ], 
                     aes_string(x = xvar, y = 'predsJtu')) +
    geom_line() +
    geom_ribbon(aes(ymin = predsJtu - se*SE_Jtu, ymax = predsJtu + se*SE_Jtu), alpha = 0.5, fill = "grey") +
    geom_line(aes(y = predsJbeta), color = 'red') +
    geom_ribbon(aes(ymin = predsJbeta - se*SE_Jbeta, ymax = predsJbeta + se*SE_Jbeta), alpha = 0.5, fill = "red") +
    geom_line(aes(y = predsHorn), color = 'blue') +
    geom_ribbon(aes(ymin = predsHorn - se*SE_Horn, ymax = predsHorn + se*SE_Horn), alpha = 0.5, fill = "blue") +
    #coord_cartesian(ylim = c(0, 0.4)) +
    theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm')) +
    labs(title = title) 
  varplots[[j]] <- thisplot
  if(vars$log[j] & !vars$discrete[j]){
    varplots[[j]] <- thisplot + scale_x_log10()
  }
}

grid.arrange(grobs = varplots, ncol = 3)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/main%20effect%20plots-1.png)<!-- -->

``` r
# write out the interactions
write.csv(newdat, file = 'output/maineffects.csv')
```

Grey: Jaccard turnover Red: Jaccard total Blue: Horn

#### Plot interactions (Jaccard turnover)

``` r
# set up the interactions to plot
# if variable is logged before scaling (see 'center and scale' above), then need to mark it here and express the limits on a log10 scale (even though log transforming is log)
ints <- data.frame(vars = c('tsign', 'tempave_metab', 'seas', 'microclim', 'mass', 'speed', 
                            'consumerfrac', 'nspp', 'thermal_bias', 'npp', 'veg', 'duration', 
                            'human_bowler', 'human_bowler'),
           min =      c(1,  0,   0.1, -2,  0,   0,   0,   0.3, -10, 1.9, 0,   0.5, 0,   0), 
           max =      c(2,  30,  16,  0.8, 8,   2,   1,   2.6, 10,  3.7, 0.3, 2,   1,   1),
           log =      c(F,  F,   F,   T,   T,   T,   F,   T,   F,   T,   T,   T,   T,   T),
           len =      c(2,  100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
           discrete = c(T,  F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   F),
           plus =     c(0,  0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,   1,   1), # what to add before log-scaling
           REALM = c(rep('Terrestrial', 12), 'Terrestrial', 'Marine'),
           REALM2 = c(rep('TerrFresh', 13), 'Marine'),
           stringsAsFactors = FALSE)
baseall <- trends[, .(type = 'all', tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                      seas.sc = mean(seas.sc, na.rm=TRUE), 
                      microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                      speed.sc = mean(speed.sc, na.rm=TRUE), 
                      mass.sc = mean(mass.sc, na.rm=TRUE), 
                      nspp.sc = 0, 
                      thermal_bias.sc = mean(thermal_bias.sc, na.rm=TRUE), 
                      npp.sc = mean(npp.sc, na.rm=TRUE), 
                      human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                      veg.sc = mean(veg.sc, na.rm=TRUE), 
                      consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
baseterr <- trends[REALM == 'Terrestrial', 
                   .(type = 'Terrestrial', 
                     tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                     seas.sc = mean(seas.sc, na.rm=TRUE), 
                     microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                     speed.sc = mean(speed.sc, na.rm=TRUE), 
                     mass.sc = mean(mass.sc, na.rm=TRUE), 
                     nspp.sc = 0, 
                     thermal_bias.sc = 0, 
                     npp.sc = mean(npp.sc, na.rm=TRUE), 
                     human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                     veg.sc = mean(veg.sc, na.rm=TRUE), 
                     consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
basemar <- trends[REALM == 'Marine', 
                  .(type = 'Marine',
                    tempave_metab.sc = mean(tempave_metab.sc, na.rm=TRUE), 
                    seas.sc = mean(seas.sc, na.rm=TRUE), 
                    microclim.sc = mean(microclim.sc, na.rm=TRUE), 
                    speed.sc = mean(speed.sc, na.rm=TRUE), 
                    mass.sc = mean(mass.sc, na.rm=TRUE), 
                    nspp.sc = 0, 
                    thermal_bias.sc = 0, 
                    npp.sc = mean(npp.sc, na.rm=TRUE), 
                    human_bowler.sc = mean(human_bowler.sc, na.rm=TRUE), 
                    veg.sc = mean(veg.sc, na.rm=TRUE), 
                    consumerfrac.sc = mean(consumerfrac.sc, na.rm=TRUE))]
basetab <- rbind(baseall, baseterr, basemar)
basetab[, ':='(duration.sc = 0, nyrBT = 20, STUDY_ID = 127L, rarefyID = '127_514668')]

# make the data frames for each interaction to plot                
for(j in 1:nrow(ints)){
  # set up a grid of temperature trends and the interacting variable
  if(ints$log[j]) intvars <- list(temptrend = seq(-0.2, 0.2, length.out = 100), 
                                  new = 10^seq(ints$min[j], ints$max[j], length.out = ints$len[j]),
                                   var = ints$vars[j])
  if(!ints$log[j]) intvars <- list(temptrend = seq(-0.2, 0.2, length.out = 100), 
                                   new = seq(ints$min[j], ints$max[j], length.out = ints$len[j]),
                                   var = ints$vars[j])
  names(intvars) <- c('temptrend', ints$vars[j], 'var')
  thisdat <- expand.grid(intvars)
  
  # scale the interacting variable
  cent <- attr(trends[[paste0(ints$var[j], '.sc')]], 'scaled:center')
  scl <- attr(trends[[paste0(ints$var[j], '.sc')]], 'scaled:scale')
  if(!is.null(cent) & !is.null(scl)){
    if(ints$log[j]) thisdat[[paste0(ints$var[j], '.sc')]] <- (log(thisdat[[ints$vars[j]]] + ints$plus[j]) - cent)/scl
    if(!ints$log[j]) thisdat[[paste0(ints$var[j], '.sc')]] <- (thisdat[[ints$var[j]]] - cent)/scl
  }

  # merge with the rest of the columns
  # use realm-specific averages for human impacts
  if(ints$vars[j] != 'tsign') colnamestouse <- setdiff(colnames(basetab), paste0(ints$var[j], '.sc'))
  if(ints$vars[j] == 'tsign') colnamestouse <- setdiff(colnames(basetab), ints$var[j])
  if(ints$vars[j] != 'human_bowler'){
    thisdat <- cbind(thisdat, basetab[type == 'all', ..colnamestouse])
  }
  if(ints$vars[j] == 'human_bowler' & ints$REALM[j] == 'Terrestrial'){
    thisdat <- cbind(thisdat, basetab[type == 'Terrestrial', ..colnamestouse])
  }
  if(ints$vars[j] == 'human_bowler' & ints$REALM[j] == 'Marine'){
    thisdat <- cbind(thisdat, basetab[type == 'Marine', ..colnamestouse])
  }
  
  # add realm
  thisdat$REALM <- ints$REALM[j]
  thisdat$REALM2 <- ints$REALM2[j]
  
  # merge with the previous iterations
  if(j == 1) newdat <- thisdat
  if(j > 1){
    colstoadd <- setdiff(colnames(thisdat), colnames(newdat))
    for(toadd in colstoadd){
      newdat[[toadd]] <- NA
    }
    
    colstoadd2 <- setdiff(colnames(newdat), colnames(thisdat))
    for(toadd in colstoadd2){
      thisdat[[toadd]] <- NA
    }
    
    newdat <- rbind(newdat, thisdat)
  } 
}

# character so that new levels can be added
newdat$REALM <- as.character(newdat$REALM)
newdat$REALM2 <- as.character(newdat$REALM2)

# add extra rows so that all factor levels are represented (for predict.lme to work)
newdat <- rbind(newdat[1:6, ], newdat)
newdat$REALM[1:6] <- c('Marine', 'Marine', 'Freshwater', 'Freshwater', 'Terrestrial', 'Terrestrial')
newdat$REALM2[1:6] <- c('Marine', 'Marine', 'TerrFresh', 'TerrFresh', 'TerrFresh', 'TerrFresh')
newdat$temptrend[1:6] <- c(-1, 1, -1, 1, -1, 1)

# trim to at least some temperature change (so that tsign is -1 or 1)
newdat <- newdat[newdat$temptrend != 0,]

# scale the temperature vars
newdat$temptrend.sc <- newdat$temptrend/attr(trends$temptrend.sc, 'scaled:scale') 
newdat$temptrend_abs <- abs(newdat$temptrend)
newdat$temptrend_abs.sc <- (newdat$temptrend_abs)/attr(trends$temptrend_abs.sc, 'scaled:scale')
newdat$tsign <- factor(sign(newdat$temptrend))

# make predictions
newdat$preds <- predict(object = modTfullJturem0, newdata = newdat, level = 0)

#remove the extra rows
newdat <- newdat[5:nrow(newdat), ]

# prep the plots
intplots <- vector('list', nrow(ints))
for(j in 1:length(intplots)){
  subs <- newdat$var == ints$vars[j] & newdat$temptrend > 0 # select warming side
  xvar <- 'temptrend_abs'
  title <- ints$vars[j]
  if(ints$vars[j] %in% c('tsign')){
    subs <- newdat$var == ints$vars[j]
  } 
  if(ints$vars[j] %in% c('thermal_bias')){
    subs <- newdat$var == ints$vars[j]
    xvar <- 'temptrend'
  } 
  if(ints$vars[j] %in% c('human_bowler')){
    subs <- newdat$var == ints$vars[j] & newdat$temptrend > 0 & newdat$REALM2 == ints$REALM2[j]
    title <- paste0('human:', ints$REALM2[j])
  } 

  thisplot <- ggplot(newdat[subs, ], 
                     aes_string(x = xvar, y = 'preds', 
                                group = ints$vars[j], 
                                color = ints$vars[j])) +
    geom_line() +
    coord_cartesian(ylim = c(-0.2, 0.4)) +
    theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm')) +
    labs(title = title)
  if(ints$log[j] & !ints$discrete[j]){
    intplots[[j]] <- thisplot + scale_color_distiller(palette = "YlGnBu", trans = 'log')
  }
  if(!ints$log[j] & !ints$discrete[j]){
    intplots[[j]] <- thisplot + scale_color_distiller(palette = "YlGnBu", trans = 'identity')
  }
  if(ints$discrete[j]){
    intplots[[j]] <- thisplot + scale_color_brewer(palette = "Dark2")
  }
}

#grid.arrange(grobs = intplots, '+', theme(plot.margin = unit(c(0,0,0,0), 'cm'))), ncol=2)
#do.call('grid.arrange', c(intplots, ncol = 2))
grid.arrange(grobs = intplots, ncol = 3)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/interaction%20plots%20modTfullJturem0-1.png)<!-- -->

``` r
# write out the interactions
write.csv(newdat, file = 'temp/interactions.csv')
```

#### Plot residuals against each predictor (Jaccard turnover)

``` r
resids <- resid(modTfullJturem0)
preds <- getData(modTfullJturem0)
col = '#00000033'
cex = 0.5
par(mfrow = c(5,4))
boxplot(resids ~ preds$REALM, cex = cex, col = col)
plot(preds$temptrend_abs.sc, resids, cex = cex, col = col)
plot(preds$tsign, resids, cex = cex, col = col)
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
plot(preds$veg.sc, resids, cex = cex, col = col)
plot(preds$human_bowler.sc, resids, cex = cex, col = col)
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/resids%20modTfull1-1.png)<!-- -->

### Remove each term from the full model

``` r
AICnas <- function(x){
  if(class(x) == 'NULL'){
    return(NA)
  } else {
    return(AIC(x))
  }
}

# set up the response variables to run
torun = data.frame(var = c('Jtutrendrem0', 'Jbetatrendrem0', 'Horntrendrem0', 
                           'Jtutrendz', 'Jbetatrendz', 'Horntrendz',
                           'Jtuexp', 'Jbetaexp', 'Hornexp',
                           'Jtumm', 'Jbetamm', 'Hornmm'), 
                   modvar = c('Jtutrendrem0', 'Jbetatrendrem0', 'Horntrendrem0', 
                           'Jtutrendz', 'Jbetatrendz', 'Horntrendz',
                           'log(Jtuexp)', 'log(Jbetaexp)', 'log(Hornexp)',
                           'log(Jtumm+1)', 'log(Jbetamm+1)', 'log(Hornmm+1)'),
                   run = TRUE, aics = NA)

# check if any already exist. If so, don't re-run that model.
if(file.exists('output/aics_from_full.csv')){
  aicsfromfull <- read.csv('output/aics_from_full.csv')
  
  for(i in 1:nrow(torun)){
    nm <- paste0('dAIC_', torun$var[i])
    if(nm %in% colnames(aicsfromfull)){ # check for column name
      if(!is.na(aicsfromfull[[nm]][1])){ # check that full model is not NA
        torun$run[i] <- FALSE
      }
    }
  }
}

randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

terms <- c('temptrend_abs.sc*REALM', 
           'temptrend_abs.sc*tsign',
           'temptrend_abs.sc*tempave_metab.sc',
           'temptrend_abs.sc*seas.sc',
           'temptrend_abs.sc*microclim.sc',
           'temptrend_abs.sc*mass.sc',
           'temptrend_abs.sc*speed.sc', 
           'temptrend_abs.sc*consumerfrac.sc',
           'temptrend_abs.sc*nspp.sc',
           'temptrend_abs.sc*thermal_bias.sc:tsign',
           'temptrend_abs.sc*npp.sc',
           'temptrend_abs.sc*veg.sc',
           'temptrend_abs.sc*duration.sc',
           'temptrend_abs.sc*human_bowler.sc:REALM2')

completerows <- trends[, complete.cases(REALM, tempave_metab.sc, 
                                seas.sc, microclim.sc, 
                                temptrend_abs.sc, mass.sc, speed.sc, 
                                consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                veg.sc, duration.sc, human_bowler.sc)]

if(any(torun$run)){
  for(k in which(torun$run)){

    i <- !is.na(trends[[as.character(torun$var[k])]]) & completerows
    
    modTdrops <- vector('list', length(terms)+2)
    names(modTdrops) <- c('full', '-temptrend_abs.sc', paste0('-', terms))
    
    # fit full model with ML for model comparison
    modTdrops[[1]] <- lme(formula(paste0(torun$modvar[k], ' ~ ', paste(terms, collapse = ' + '))),
                          random = randef, weights = varef, data = trends[i,], method = 'ML',
                          control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))
    
    # w/out temptrend
    modTdrops[[2]] <- lme(formula(paste0(torun$modvar[k], ' ~ ', paste(gsub('temptrend_abs.sc\\*', '', terms), 
                                                                       collapse = ' + '))),
                          random = list(STUDY_ID = ~ 1, rarefyID = ~1), weights = varef, data = trends[i,], method = 'ML',
                          control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))
    
    for(j in 1:length(terms)){
      print(j)
      tryCatch({
        modTdrops[[j+2]] <- lme(formula(paste0(torun$modvar[k], ' ~ ', paste(terms[-j], collapse = ' + '))),
                                random = randef, weights = varef, data = trends[i,], method = 'ML')
        
      }, error = function(e){
        print('going to optim (Jtu)')
        tryCatch({
          modTdrops[[j+2]] <- lme(formula(paste0(torun$modvar[k], ' ~ ', paste(terms[-j], collapse = ' + '))),
                                  random = randef, weights = varef, data = trends[i,], method = 'ML',
                                  control = lmeControl(opt = 'optim'))
          
        }, error = function(e){
          print('going to more iters (Jtu)') 
          tryCatch({
            modTdrops[[j+2]] <- lme(formula(paste0(torun$modvar[k], ' ~ ', paste(terms[-j], collapse = ' + '))),
                                    random = randef, weights = varef, data = trends[i,], method = 'ML',
                                    control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
            
          }, error= function(e){
            print('giving up on this one')
            modTdrops[[j+2]] <- NA
          })
        })
      })
    }
    
    torun$aics[k] <- list(sapply(modTdrops, AICnas))
  }
}



# if there was anything new
if(any(torun$run)){
  if(!exists('aicsfromfull')){
    aicsfromfull <- data.frame(mod = names(torun$aics[which(torun$run)[1]][[1]]))
  }
  
  # subtract full from each model AIC. Negative means term removal is supported. Positive means full is the better model.
  for(j in which(torun$run)){
    nm <- paste0('dAIC_', torun$var[j])
    aics <- torun$aics[j][[1]]
    aicsfromfull[[nm]] <- aics - aics[1]
  }

  
  # write out
  write.csv(aicsfromfull, file = 'output/aics_from_full.csv', row.names = FALSE)
}

aicsfromfull
```

    ##                                         mod dAIC_Jtutrendrem0 dAIC_Jbetatrendrem0 dAIC_Horntrendrem0 dAIC_Jtutrendz dAIC_Jbetatrendz dAIC_Horntrendz dAIC_Jtuexp dAIC_Jbetaexp dAIC_Hornexp  dAIC_Jtumm dAIC_Jbetamm dAIC_Hornmm
    ## 1                                      full          0.000000           0.0000000          0.0000000     0.00000000        0.0000000       0.0000000   0.0000000     0.0000000    0.0000000   0.0000000      0.00000   0.0000000
    ## 2                         -temptrend_abs.sc         62.200728          80.3280362        112.1176194    29.88523622      305.1488520     210.5158602  29.6469380    72.7609667   56.8277616  35.1372318     35.94037  76.2525218
    ## 3                   -temptrend_abs.sc*REALM                NA          13.5408327          2.2425092    -4.45663630        2.8891238              NA  22.2246952            NA   38.6483666  40.7924410           NA  58.7190131
    ## 4                   -temptrend_abs.sc*tsign        -23.396632           8.3097160         26.6812823     7.95574584       -2.4582723              NA   5.0503891    15.6737018    7.6595869  -3.8070037           NA  -0.2236897
    ## 5        -temptrend_abs.sc*tempave_metab.sc        -28.992763          31.1465717                 NA    37.65341409       43.2422703              NA  38.1202881    42.9415560           NA 127.3682873           NA 341.3841201
    ## 6                 -temptrend_abs.sc*seas.sc                NA          -2.5851040         11.6469985     5.23666626       11.1566418      15.5011308  20.0987205            NA   28.1472412  46.1012779           NA 112.1358726
    ## 7            -temptrend_abs.sc*microclim.sc          3.098449          12.3565788         -2.5048568    -0.02753739       -0.3926990              NA  24.8226107    11.1933986           NA   1.5393081           NA  35.0506541
    ## 8                 -temptrend_abs.sc*mass.sc        -32.437228           3.5624654         -1.7724568    -2.41821426       -3.0954303              NA  -0.8262754    -3.8188936   -3.0600878  -0.6556442           NA  -1.5110113
    ## 9                -temptrend_abs.sc*speed.sc                NA          -0.1755050         -1.5559022     3.11149019        0.9130635       0.8625964  -3.8627649            NA   32.6971458  -1.7694147           NA  47.0579398
    ## 10        -temptrend_abs.sc*consumerfrac.sc                NA           4.8686487        112.5776151     6.12974472       -2.2050354              NA   4.7509945            NA           NA  -2.4577371           NA  17.6927302
    ## 11                -temptrend_abs.sc*nspp.sc        -41.268152           0.8339248         -2.4826457   113.16458351      197.8833760              NA 288.1708849   570.7120619  233.3251121 556.7857883           NA 159.3119398
    ## 12  -temptrend_abs.sc*thermal_bias.sc:tsign        -41.169431          -2.3162423          1.9716230     7.82665526               NA       8.3961926   4.2824601     3.0093791   -0.7313611   7.5126486           NA   2.8970160
    ## 13                 -temptrend_abs.sc*npp.sc                NA           7.0497869          5.0155832    15.69389286        4.6186491       2.1463637  -3.8738460    -1.0186107    1.9520696  -0.7595978           NA   2.4049496
    ## 14                 -temptrend_abs.sc*veg.sc                NA           0.1194515         -0.6426642    -3.16679252        2.6966695       6.1104134   1.3718145    -2.8664130   13.4485193  -2.9724213           NA   3.3707273
    ## 15            -temptrend_abs.sc*duration.sc                NA          14.2669868          7.6351796    19.44480088       32.2482413       8.2984360          NA   905.0947284  631.3708169 530.1876356           NA 196.6247988
    ## 16 -temptrend_abs.sc*human_bowler.sc:REALM2        -23.675058          -5.4155448          4.7771692    -1.41942433       -1.7762558              NA   2.3598431     0.9384098   15.8953935  -4.4610292           NA  20.3709221

#### Plot deltaAICs for all models

``` r
# transform for a plot
aicsfromfulllong <- reshape(aicsfromfull, direction = 'long',
                            varying = c('dAIC_Jtutrendrem0', 'dAIC_Jbetatrendrem0', 'dAIC_Horntrendrem0',
                                        'dAIC_Jtutrendz', 'dAIC_Jbetatrendz', 'dAIC_Horntrendz',
                                        'dAIC_Jtuexp', 'dAIC_Jbetaexp', 'dAIC_Hornexp',
                                        'dAIC_Jtumm', 'dAIC_Jbetamm', 'dAIC_Hornmm'),
                            v.names = 'dAIC',
                            idvar = 'mod',
                            timevar = 'type',
                            times = c('Jtutrendrem0', 'Jbetatrendrem0', 'Horntrendrem0',
                                      'Jtutrendz', 'Jbetatrendz', 'Horntrendz',
                                      'Jtuexp', 'Jbetaexp', 'Hornexp',
                                      'Jtumm', 'Jbetamm', 'Hornmm'))
aicsfromfulllong$response <- gsub('trendrem0|trendz|exp|mm', '', aicsfromfulllong$type)
aicsfromfulllong$fit <- gsub('Jtu|Jbeta|Horn', '', aicsfromfulllong$type)

# plot
ggplot(aicsfromfulllong, aes(y=dAIC, x = mod, shape = response, color = fit, group = fit)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_y_continuous(trans = signedsqrttrans) +
  coord_flip()
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20dAICs-1.png)<!-- -->

Light grey is for Jaccard turnover, dark grey is for Jaccard total,
black is for Morisita-Horn. Clear that removing temperature trend makes
the model quite a bit worse and has the biggest effect.

## Simplify the full models

This takes a couple days on a laptop to run if temp/ files not
available. Not run when knitting

``` r
randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

terms <- ' ~ temptrend_abs.sc*REALM + temptrend_abs.sc*tsign + temptrend_abs.sc*tempave_metab.sc + temptrend_abs.sc*seas.sc + temptrend_abs.sc*microclim.sc + temptrend_abs.sc*mass.sc + temptrend_abs.sc*speed.sc + temptrend_abs.sc*consumerfrac.sc + temptrend_abs.sc*nspp.sc + temptrend_abs.sc*thermal_bias.sc:tsign + temptrend_abs.sc*npp.sc + temptrend_abs.sc*veg.sc + temptrend_abs.sc*duration.sc + temptrend_abs.sc*human_bowler.sc:REALM2'

completerows <- trends[, complete.cases(REALM, tempave_metab.sc, 
                                        seas.sc, microclim.sc, 
                                        temptrend_abs.sc, mass.sc, speed.sc, 
                                        consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                                        veg.sc, duration.sc, human_bowler.sc)]

# rem0
if(file.exists('temp/modTsimpJturem0.rds')){
  modTsimpJturem0 <- readRDS('temp/modTsimpJturem0.rds')
} else {
  i <- !is.na(trends$Jtutrendrem0) & completerows
  
  modTfullJturem0ML <- lme(formula(paste0('Jtutrendrem0', terms)),
                           random = randef, weights = varef, data = trends[i,], method = 'ML',
                           control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))
  modTsimpJturem0 <- stepAIC(modTfullJturem0ML, direction = 'backward')
  saveRDS(modTsimpJturem0, file = 'temp/modTsimpJturem0.rds')
}

if(file.exists('temp/modTsimpJbetarem0.rds')){
  modTsimpJbetarem0 <- readRDS('temp/modTsimpJbetarem0.rds')
} else {
  i <- !is.na(trends$Jbetatrendrem0) & completerows
  modTfullJbetarem0ML <- lme(formula(paste0('Jbetatrendrem0', terms)),
                       random = randef, weights = varef, data = trends[i,], method = 'ML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))
  modTsimpJbetarem0 <- stepAIC(modTfullJbetarem0ML, direction = 'backward')
  saveRDS(modTsimpJbetarem0, file = 'temp/modTsimpJbetarem0.rds')
}

if(file.exists('temp/modTsimpHornrem0.rds')){
  modTsimpHornrem0 <- readRDS('temp/modTsimpHornrem0.rds')
} else {
  i <- !is.na(trends$Horntrendrem0) & completerows
  modTfullHornrem0ML <- lme(formula(paste0('Horntrendrem0', terms)),
                            random = randef, weights = varef, data = trends[i,], method = 'ML',
                            control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))
  modTsimpHornrem0 <- stepAIC(modTfullHornrem0ML, direction = 'backward')
  saveRDS(modTsimpHornrem0, file = 'temp/modTsimpHornrem0.rds')
}

# z
if(file.exists('temp/modTsimpJtuz.rds')){
  modTsimpJtuz <- readRDS('temp/modTsimpJtuz.rds')
} else {
  i <- !is.na(trends$Jtutrendz) & completerows
  
  modTfullJtuzML <- lme(formula(paste0('Jtutrendz', terms)),
                           random = randef, weights = varef, data = trends[i,], method = 'ML',
                           control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))
  modTsimpJtuz <- stepAIC(modTfullJtuzML, direction = 'backward')
  saveRDS(modTsimpJtuz, file = 'temp/modTsimpJtuz.rds')
}

if(file.exists('temp/modTsimpJbetaz.rds')){
  modTsimpJbetaz <- readRDS('temp/modTsimpJbetaz.rds')
} else {
  i <- !is.na(trends$Jbetatrendz) & completerows
  modTfullJbetazML <- lme(formula(paste0('Jbetatrendz', terms)),
                       random = randef, weights = varef, data = trends[i,], method = 'ML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  modTsimpJbetaz <- stepAIC(modTfullJbetazML, direction = 'backward')
  saveRDS(modTsimpJbetaz, file = 'temp/modTsimpJbetaz.rds')
}

if(file.exists('temp/modTsimpHornz.rds')){
  modTsimpHornz <- readRDS('temp/modTsimpHornz.rds')
} else {
  i <- !is.na(trends$Horntrendz) & completerows
  modTfullHornzML <- lme(formula(paste0('Horntrendz', terms)),
                            random = randef, weights = varef, data = trends[i,], method = 'ML',
                            control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500, opt = 'optim'))
  modTsimpHornz <- stepAIC(modTfullHornzML, direction = 'backward')
  saveRDS(modTsimpHornz, file = 'temp/modTsimpHornz.rds')
}

summary(modTsimpJturem0)
summary(modTsimpJbetarem0)
summary(modTsimpHornrem0)
```

## Make realm-specific models

``` r
i1 <- trends[, REALM == 'Terrestrial' & complete.cases(Horntrendrem0, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc)]
i2 <- trends[, REALM == 'Freshwater' & complete.cases(Horntrendrem0, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             nspp.sc, thermal_bias.sc, npp.sc, 
                             veg.sc, duration.sc, human_bowler.sc)] # no consumerfrac
i3 <- trends[, REALM == 'Marine' & complete.cases(Horntrendrem0, tempave_metab.sc, seas.sc, microclim.sc, 
                             temptrend_abs.sc, mass.sc, speed.sc, 
                             consumerfrac.sc, nspp.sc, thermal_bias.sc, npp.sc, 
                             duration.sc, human_bowler.sc)] # no veg

print(paste('Terrestrial', sum(i1)))
```

    ## [1] "Terrestrial 2232"

``` r
print(paste('Freshwater', sum(i2)))
```

    ## [1] "Freshwater 470"

``` r
print(paste('Marine', sum(i3)))
```

    ## [1] "Marine 28054"

``` r
randef <- list(STUDY_ID = ~ temptrend_abs.sc, rarefyID = ~1)
varef <- varPower(-0.5, ~nyrBT)

# land
if(file.exists('temp/modTfullHornTerr.rds')){
  modTfullHornTerr <- readRDS('temp/modTfullHornTerr.rds')
} else {
  modTfullHornTerr <- lme(Horntrendrem0 ~ 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*consumerfrac.sc +
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*veg.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc,
                       random = randef, weights = varef, data = trends[i1,], method = 'REML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullHornTerr, file = 'temp/modTfullHornTerr.rds')
}

# freshwater
if(file.exists('temp/modTfullHornFresh.rds')){
  modTfullHornFresh <- readRDS('temp/modTfullHornFresh.rds')
} else {
  modTfullHornFresh <- lme(Horntrendrem0 ~ 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*veg.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc,
                       random = randef, weights = varef, data = trends[i2,], method = 'REML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullHornFresh, file = 'temp/modTfullHornFresh.rds')
}

# marine
if(file.exists('temp/modTfullHornMar.rds')){
  modTfullHornMar <- readRDS('temp/modTfullHornMar.rds')
} else {
  modTfullHornMar <- lme(Horntrendrem0 ~ 
                         temptrend_abs.sc*tsign +
                         temptrend_abs.sc*tempave_metab.sc + 
                         temptrend_abs.sc*seas.sc + 
                         temptrend_abs.sc*microclim.sc + 
                         temptrend_abs.sc*mass.sc + 
                         temptrend_abs.sc*speed.sc + 
                         temptrend_abs.sc*consumerfrac.sc +
                         temptrend_abs.sc*nspp.sc +
                         temptrend_abs.sc*thermal_bias.sc:tsign +
                         temptrend_abs.sc*npp.sc +
                         temptrend_abs.sc*duration.sc +
                         temptrend_abs.sc*human_bowler.sc,
                       random = randef, weights = varef, data = trends[i3,], method = 'REML',
                       control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, msMaxEval = 500))
  saveRDS(modTfullHornMar, file = 'temp/modTfullHornMar.rds')
}

summary(modTfullHornTerr)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i1, ] 
    ##         AIC       BIC   logLik
    ##   -7679.822 -7474.726 3875.911
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev      Corr  
    ## (Intercept)      0.005729249 (Intr)
    ## temptrend_abs.sc 0.021903364 0.711 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##          (Intercept)  Residual
    ## StdDev: 6.238717e-09 0.9040518
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -1.617461 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * tsign + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * veg.sc + temptrend_abs.sc * duration.sc +      temptrend_abs.sc * human_bowler.sc 
    ##                                                 Value   Std.Error   DF    t-value p-value
    ## (Intercept)                               0.010726721 0.005609111 2113  1.9123746  0.0560
    ## temptrend_abs.sc                          0.005241922 0.017430573 2113  0.3007314  0.7636
    ## tsign1                                    0.003710803 0.003340948 2113  1.1107034  0.2668
    ## tempave_metab.sc                         -0.005381245 0.003090136 2113 -1.7414270  0.0818
    ## seas.sc                                   0.000569580 0.000808835 2113  0.7041982  0.4814
    ## microclim.sc                             -0.000502527 0.000574763 2113 -0.8743207  0.3820
    ## mass.sc                                  -0.000791663 0.001747304 2113 -0.4530768  0.6505
    ## speed.sc                                  0.003527914 0.002328604 2113  1.5150340  0.1299
    ## consumerfrac.sc                          -0.000184310 0.001682919 2113 -0.1095180  0.9128
    ## nspp.sc                                  -0.001359577 0.001044352 2113 -1.3018377  0.1931
    ## npp.sc                                    0.000417629 0.001154073 2113  0.3618737  0.7175
    ## veg.sc                                    0.000641897 0.000776532 2113  0.8266203  0.4085
    ## duration.sc                              -0.001973706 0.001157250 2113 -1.7055141  0.0882
    ## human_bowler.sc                           0.000825245 0.000486397 2113  1.6966495  0.0899
    ## temptrend_abs.sc:tsign1                   0.001411494 0.006887970 2113  0.2049216  0.8377
    ## temptrend_abs.sc:tempave_metab.sc         0.005774600 0.009449869 2113  0.6110772  0.5412
    ## temptrend_abs.sc:seas.sc                  0.001648212 0.002348589 2113  0.7017882  0.4829
    ## temptrend_abs.sc:microclim.sc             0.000171024 0.001720580 2113  0.0993991  0.9208
    ## temptrend_abs.sc:mass.sc                  0.001782059 0.004450382 2113  0.4004283  0.6889
    ## temptrend_abs.sc:speed.sc                 0.001555549 0.007670522 2113  0.2027957  0.8393
    ## temptrend_abs.sc:consumerfrac.sc         -0.005547075 0.005697791 2113 -0.9735483  0.3304
    ## temptrend_abs.sc:nspp.sc                  0.004214386 0.002552153 2113  1.6513062  0.0988
    ## tsign-1:thermal_bias.sc                   0.002360737 0.002794819 2113  0.8446834  0.3984
    ## tsign1:thermal_bias.sc                    0.000348047 0.000879702 2113  0.3956424  0.6924
    ## temptrend_abs.sc:npp.sc                  -0.002322102 0.003625000 2113 -0.6405798  0.5219
    ## temptrend_abs.sc:veg.sc                  -0.001261633 0.002477730 2113 -0.5091891  0.6107
    ## temptrend_abs.sc:duration.sc             -0.001367813 0.001974040 2113 -0.6929005  0.4884
    ## temptrend_abs.sc:human_bowler.sc         -0.001215046 0.001287672 2113 -0.9435995  0.3455
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.002383313 0.003543596 2113 -0.6725691  0.5013
    ## temptrend_abs.sc:tsign1:thermal_bias.sc  -0.006133752 0.002758227 2113 -2.2238026  0.0263
    ##  Correlation: 
    ##                                          (Intr) tmpt_. tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc veg.sc drtn.s hmn_b. tm_.:1 tmptrnd_bs.sc:t_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. tmptrnd_bs.sc:h_. t_.:-1
    ## temptrend_abs.sc                         -0.624                                                                                                                                                                                                                                                                                                                          
    ## tsign1                                   -0.167 -0.068                                                                                                                                                                                                                                                                                                                   
    ## tempave_metab.sc                         -0.478  0.318 -0.324                                                                                                                                                                                                                                                                                                            
    ## seas.sc                                  -0.251  0.268 -0.155  0.091                                                                                                                                                                                                                                                                                                     
    ## microclim.sc                             -0.183  0.230  0.077 -0.003  0.485                                                                                                                                                                                                                                                                                              
    ## mass.sc                                   0.601 -0.266 -0.060 -0.145 -0.040  0.016                                                                                                                                                                                                                                                                                       
    ## speed.sc                                  0.034  0.100  0.013 -0.607 -0.053 -0.048 -0.182                                                                                                                                                                                                                                                                                
    ## consumerfrac.sc                           0.403 -0.389  0.195  0.010 -0.068  0.073  0.482 -0.697                                                                                                                                                                                                                                                                         
    ## nspp.sc                                   0.066 -0.037 -0.057 -0.185 -0.077 -0.122  0.122  0.053  0.031                                                                                                                                                                                                                                                                  
    ## npp.sc                                   -0.018 -0.041  0.106 -0.044  0.205  0.154 -0.048  0.020 -0.004 -0.075                                                                                                                                                                                                                                                           
    ## veg.sc                                   -0.240  0.291 -0.085 -0.062 -0.064 -0.101  0.020  0.183 -0.163 -0.001 -0.601                                                                                                                                                                                                                                                    
    ## duration.sc                              -0.139  0.063 -0.201  0.078  0.050 -0.029 -0.028 -0.126  0.181  0.166 -0.005 -0.076                                                                                                                                                                                                                                             
    ## human_bowler.sc                          -0.255  0.261  0.037  0.057 -0.036  0.182 -0.122  0.018 -0.146  0.045 -0.227  0.324 -0.129                                                                                                                                                                                                                                      
    ## temptrend_abs.sc:tsign1                   0.161 -0.251 -0.452  0.100 -0.096 -0.212  0.041  0.050 -0.115 -0.014 -0.112  0.013  0.126 -0.119                                                                                                                                                                                                                               
    ## temptrend_abs.sc:tempave_metab.sc         0.288 -0.534  0.208 -0.654 -0.055 -0.016 -0.013  0.415 -0.030  0.149  0.031  0.052 -0.100 -0.022 -0.088                                                                                                                                                                                                                        
    ## temptrend_abs.sc:seas.sc                  0.229 -0.375  0.085 -0.043 -0.758 -0.410  0.031  0.016  0.073  0.037 -0.099 -0.007 -0.008 -0.105  0.138  0.034                                                                                                                                                                                                                 
    ## temptrend_abs.sc:microclim.sc             0.173 -0.277 -0.080  0.005 -0.391 -0.854 -0.016  0.049 -0.071  0.051 -0.094  0.047  0.027 -0.169  0.250  0.001             0.492                                                                                                                                                                                               
    ## temptrend_abs.sc:mass.sc                 -0.360  0.384  0.046 -0.055 -0.021 -0.025 -0.634  0.344 -0.482 -0.045  0.004  0.008  0.029  0.092  0.039  0.099            -0.001             0.027                                                                                                                                                                             
    ## temptrend_abs.sc:speed.sc                 0.092 -0.043 -0.033  0.423  0.070  0.075  0.268 -0.736  0.595 -0.017 -0.007 -0.127  0.093 -0.007 -0.063 -0.575            -0.027            -0.074            -0.513                                                                                                                                                           
    ## temptrend_abs.sc:consumerfrac.sc         -0.312  0.454 -0.100 -0.060  0.000 -0.069 -0.373  0.598 -0.795 -0.011  0.004  0.090 -0.149  0.061  0.097  0.011            -0.041             0.080             0.586            -0.720                                                                                                                                         
    ## temptrend_abs.sc:nspp.sc                 -0.011 -0.074  0.032  0.113 -0.025  0.009 -0.062 -0.007 -0.032 -0.561  0.018 -0.012 -0.131  0.034  0.013 -0.047             0.026            -0.006             0.027            -0.031             0.030                                                                                                                       
    ## tsign-1:thermal_bias.sc                  -0.053 -0.120  0.654 -0.300 -0.282 -0.080 -0.039  0.031  0.156 -0.034  0.094 -0.083  0.004  0.041 -0.188  0.188             0.205             0.056             0.050            -0.062            -0.070            0.043                                                                                                      
    ## tsign1:thermal_bias.sc                    0.209 -0.215  0.001 -0.078 -0.669 -0.404  0.052  0.074 -0.014 -0.061  0.144 -0.111  0.017  0.040  0.306  0.026             0.517             0.316             0.026            -0.052             0.035            0.090             0.260                                                                                    
    ## temptrend_abs.sc:npp.sc                  -0.013  0.072 -0.105  0.053 -0.134 -0.082  0.022 -0.022 -0.004  0.005 -0.895  0.560  0.003  0.192  0.141 -0.043             0.070             0.078            -0.014             0.028            -0.027           -0.032            -0.097 -0.145                                                                             
    ## temptrend_abs.sc:veg.sc                   0.250 -0.320  0.067  0.034  0.025  0.046  0.016 -0.131  0.133 -0.016  0.554 -0.872  0.068 -0.339 -0.010 -0.075             0.034            -0.026            -0.019             0.153            -0.108           -0.019             0.076  0.123 -0.652                                                                      
    ## temptrend_abs.sc:duration.sc             -0.086  0.171  0.006 -0.093  0.028  0.051 -0.008  0.072 -0.069  0.166 -0.028  0.068 -0.033  0.110 -0.226  0.087            -0.072            -0.061             0.034            -0.059             0.054           -0.269            -0.048 -0.143  0.020            -0.042                                                    
    ## temptrend_abs.sc:human_bowler.sc          0.268 -0.337 -0.020 -0.062 -0.058 -0.182  0.103 -0.005  0.115  0.015  0.177 -0.375  0.126 -0.811  0.145  0.052             0.094             0.219            -0.084            -0.021            -0.067           -0.034            -0.026 -0.006 -0.224             0.449           -0.113                                   
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc  0.022  0.019 -0.381  0.147  0.248  0.033  0.016  0.020 -0.118 -0.011 -0.247  0.150  0.003 -0.086  0.626 -0.128            -0.337            -0.037             0.047            -0.030             0.096           -0.038            -0.442 -0.280  0.284            -0.168           -0.022            0.133                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc  -0.155  0.291  0.066  0.029  0.466  0.314 -0.015 -0.030  0.014  0.049 -0.183  0.146 -0.068  0.038 -0.369 -0.020            -0.639            -0.367            -0.012             0.024            -0.007           -0.120            -0.187 -0.839  0.186            -0.166            0.169           -0.003             0.334
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.21851197 -0.50451880 -0.05139768  0.43121416  6.16758380 
    ## 
    ## Number of Observations: 2232
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                     90                   2232

``` r
summary(modTfullHornFresh)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i2, ] 
    ##         AIC       BIC   logLik
    ##   -699.8191 -560.7146 383.9096
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev       Corr  
    ## (Intercept)      6.612208e-07 (Intr)
    ## temptrend_abs.sc 6.649456e-06 0     
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.01627269 1.814654
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##  power 
    ## -2.126 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * tsign + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * nspp.sc + temptrend_abs.sc *      thermal_bias.sc:tsign + temptrend_abs.sc * npp.sc + temptrend_abs.sc *      veg.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc 
    ##                                                Value  Std.Error  DF   t-value p-value
    ## (Intercept)                               0.02445726 0.02729062 425  0.896178  0.3707
    ## temptrend_abs.sc                         -0.10743909 0.08401182 425 -1.278857  0.2016
    ## tsign1                                    0.00372375 0.01113036 425  0.334558  0.7381
    ## tempave_metab.sc                         -0.01904753 0.02125145 425 -0.896293  0.3706
    ## seas.sc                                   0.00020744 0.00688607 425  0.030124  0.9760
    ## microclim.sc                             -0.00426453 0.00518764 425 -0.822055  0.4115
    ## mass.sc                                   0.00459537 0.00430297 425  1.067952  0.2861
    ## speed.sc                                 -0.00467024 0.00677339 425 -0.689499  0.4909
    ## nspp.sc                                  -0.00527132 0.00491548 425 -1.072393  0.2842
    ## npp.sc                                    0.00997940 0.00687300 425  1.451972  0.1472
    ## veg.sc                                   -0.01245506 0.00665398 425 -1.871820  0.0619
    ## duration.sc                               0.00105354 0.00558516 425  0.188633  0.8505
    ## human_bowler.sc                           0.00478392 0.00614358 425  0.778686  0.4366
    ## temptrend_abs.sc:tsign1                  -0.02463559 0.02712632 425 -0.908180  0.3643
    ## temptrend_abs.sc:tempave_metab.sc         0.07119119 0.06306649 425  1.128827  0.2596
    ## temptrend_abs.sc:seas.sc                  0.00930454 0.02144083 425  0.433964  0.6645
    ## temptrend_abs.sc:microclim.sc             0.01204772 0.01880174 425  0.640777  0.5220
    ## temptrend_abs.sc:mass.sc                 -0.02376320 0.01418380 425 -1.675376  0.0946
    ## temptrend_abs.sc:speed.sc                 0.02854621 0.02152794 425  1.326007  0.1855
    ## temptrend_abs.sc:nspp.sc                  0.01869024 0.01526823 425  1.224126  0.2216
    ## tsign-1:thermal_bias.sc                   0.00628044 0.01773296 425  0.354168  0.7234
    ## tsign1:thermal_bias.sc                    0.00475907 0.01073456 425  0.443341  0.6577
    ## temptrend_abs.sc:npp.sc                  -0.04307431 0.02280246 425 -1.889020  0.0596
    ## temptrend_abs.sc:veg.sc                   0.06875263 0.03163765 425  2.173127  0.0303
    ## temptrend_abs.sc:duration.sc             -0.03667602 0.01149088 425 -3.191751  0.0015
    ## temptrend_abs.sc:human_bowler.sc         -0.01600478 0.02089076 425 -0.766118  0.4440
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc  0.04400304 0.04511074 425  0.975445  0.3299
    ## temptrend_abs.sc:tsign1:thermal_bias.sc  -0.01479546 0.03900917 425 -0.379282  0.7047
    ##  Correlation: 
    ##                                          (Intr) tmpt_. tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc nspp.s npp.sc veg.sc drtn.s hmn_b. tm_.:1 tmptrnd_bs.sc:t_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:v. tmptrnd_bs.sc:d. tmptrnd_bs.sc:h_. t_.:-1
    ## temptrend_abs.sc                         -0.608                                                                                                                                                                                                                                                                                                  
    ## tsign1                                   -0.301  0.075                                                                                                                                                                                                                                                                                           
    ## tempave_metab.sc                          0.313 -0.075  0.053                                                                                                                                                                                                                                                                                    
    ## seas.sc                                  -0.094  0.045 -0.007  0.381                                                                                                                                                                                                                                                                             
    ## microclim.sc                             -0.205  0.100  0.083  0.327  0.105                                                                                                                                                                                                                                                                      
    ## mass.sc                                  -0.021  0.042  0.041 -0.183  0.216 -0.363                                                                                                                                                                                                                                                               
    ## speed.sc                                  0.175 -0.125 -0.083  0.120  0.217  0.045 -0.654                                                                                                                                                                                                                                                        
    ## nspp.sc                                   0.286 -0.103 -0.121 -0.335 -0.558 -0.066 -0.125  0.015                                                                                                                                                                                                                                                 
    ## npp.sc                                   -0.205 -0.067 -0.023 -0.077  0.463  0.180  0.038  0.120 -0.079                                                                                                                                                                                                                                          
    ## veg.sc                                   -0.574  0.614 -0.072  0.040  0.029 -0.144 -0.065  0.071 -0.247 -0.349                                                                                                                                                                                                                                   
    ## duration.sc                              -0.236  0.207 -0.173  0.087 -0.047  0.094 -0.044 -0.075 -0.217  0.048  0.161                                                                                                                                                                                                                            
    ## human_bowler.sc                           0.091  0.049 -0.043 -0.141 -0.324 -0.036  0.242 -0.186  0.060 -0.128  0.034  0.039                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tsign1                   0.150 -0.334 -0.514  0.021  0.013  0.012 -0.158  0.108  0.011 -0.017  0.017  0.073 -0.144                                                                                                                                                                                                              
    ## temptrend_abs.sc:tempave_metab.sc        -0.026 -0.049  0.029 -0.653 -0.541 -0.310  0.137 -0.151  0.317 -0.196 -0.092 -0.029  0.325 -0.084                                                                                                                                                                                                       
    ## temptrend_abs.sc:seas.sc                  0.055 -0.055  0.060 -0.430 -0.763 -0.184 -0.004 -0.254  0.335 -0.361 -0.075  0.039  0.356 -0.059  0.843                                                                                                                                                                                                
    ## temptrend_abs.sc:microclim.sc             0.082 -0.056 -0.006 -0.283 -0.201 -0.650  0.183 -0.017  0.011 -0.113  0.094 -0.008  0.133 -0.048  0.571             0.452                                                                                                                                                                              
    ## temptrend_abs.sc:mass.sc                 -0.049  0.014 -0.094  0.144  0.065  0.232 -0.655  0.426 -0.012  0.017  0.084  0.049 -0.380  0.226 -0.433            -0.379            -0.364                                                                                                                                                            
    ## temptrend_abs.sc:speed.sc                -0.087  0.144  0.112 -0.103 -0.277 -0.023  0.418 -0.649  0.006 -0.076 -0.083  0.053  0.274 -0.187  0.275             0.470             0.109            -0.707                                                                                                                                          
    ## temptrend_abs.sc:nspp.sc                 -0.203  0.126  0.099  0.198  0.256  0.033 -0.025 -0.032 -0.683  0.017  0.153  0.128 -0.147  0.021 -0.373            -0.273            -0.016             0.209            -0.008                                                                                                                        
    ## tsign-1:thermal_bias.sc                   0.444 -0.184 -0.464  0.187 -0.010  0.053 -0.066  0.043  0.156  0.130 -0.251 -0.033  0.065  0.185 -0.118            -0.076            -0.040             0.067            -0.054            -0.053                                                                                                      
    ## tsign1:thermal_bias.sc                    0.174 -0.010  0.106  0.604  0.376  0.104  0.165 -0.026 -0.226  0.114  0.076 -0.010  0.243 -0.201 -0.417            -0.366            -0.075            -0.176             0.055             0.079             0.157                                                                                    
    ## temptrend_abs.sc:npp.sc                  -0.095  0.368  0.010 -0.123 -0.337 -0.076 -0.061 -0.055  0.041 -0.668  0.397  0.023  0.062  0.110  0.360             0.429             0.208             0.025            -0.016             0.041            -0.145 -0.235                                                                             
    ## temptrend_abs.sc:veg.sc                   0.403 -0.845  0.069 -0.028 -0.059  0.059  0.028 -0.048  0.095  0.255 -0.717 -0.149 -0.008  0.019  0.094             0.091            -0.132            -0.130             0.122            -0.149             0.109 -0.034 -0.603                                                                      
    ## temptrend_abs.sc:duration.sc              0.078 -0.033  0.003  0.076  0.215  0.043  0.022  0.122  0.057  0.062 -0.098 -0.411 -0.065 -0.105 -0.063            -0.271            -0.138             0.042            -0.195            -0.219             0.086  0.116 -0.103             0.124                                                    
    ## temptrend_abs.sc:human_bowler.sc         -0.057  0.005 -0.034  0.228  0.333  0.101 -0.273  0.240 -0.083  0.077  0.080 -0.005 -0.700  0.191 -0.485            -0.523            -0.180             0.576            -0.452             0.042            -0.025 -0.068 -0.008            -0.143            0.052                                   
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.231  0.376  0.211 -0.220 -0.214 -0.060  0.117 -0.113  0.074 -0.246  0.210  0.049  0.126 -0.412  0.426             0.407             0.184            -0.290             0.189            -0.175            -0.506 -0.169  0.439            -0.271           -0.101           -0.181                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc  -0.090 -0.051 -0.086 -0.411 -0.435 -0.038 -0.182 -0.013  0.205 -0.203 -0.050  0.121 -0.047  0.244  0.615             0.590             0.150             0.059            -0.028            -0.213            -0.124 -0.756  0.405             0.060           -0.118           -0.062             0.314
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -4.1334763 -0.2822923 -0.0376483  0.2876199  4.8717935 
    ## 
    ## Number of Observations: 470
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                     18                    470

``` r
summary(modTfullHornMar)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: trends[i3, ] 
    ##         AIC       BIC   logLik
    ##   -79420.71 -79140.52 39744.36
    ## 
    ## Random effects:
    ##  Formula: ~temptrend_abs.sc | STUDY_ID
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##                  StdDev     Corr  
    ## (Intercept)      0.02420393 (Intr)
    ## temptrend_abs.sc 0.02638454 0.003 
    ## 
    ##  Formula: ~1 | rarefyID %in% STUDY_ID
    ##         (Intercept) Residual
    ## StdDev:  0.02099175 3.128558
    ## 
    ## Variance function:
    ##  Structure: Power of variance covariate
    ##  Formula: ~nyrBT 
    ##  Parameter estimates:
    ##     power 
    ## -2.555666 
    ## Fixed effects: Horntrendrem0 ~ temptrend_abs.sc * tsign + temptrend_abs.sc *      tempave_metab.sc + temptrend_abs.sc * seas.sc + temptrend_abs.sc *      microclim.sc + temptrend_abs.sc * mass.sc + temptrend_abs.sc *      speed.sc + temptrend_abs.sc * consumerfrac.sc + temptrend_abs.sc *      nspp.sc + temptrend_abs.sc * thermal_bias.sc:tsign + temptrend_abs.sc *      npp.sc + temptrend_abs.sc * duration.sc + temptrend_abs.sc *      human_bowler.sc 
    ##                                                 Value   Std.Error    DF   t-value p-value
    ## (Intercept)                               0.009782583 0.003535617 27937  2.766867  0.0057
    ## temptrend_abs.sc                          0.016369615 0.005525760 27937  2.962419  0.0031
    ## tsign1                                   -0.002951738 0.000897864 27937 -3.287510  0.0010
    ## tempave_metab.sc                         -0.011724422 0.001462964 27937 -8.014155  0.0000
    ## seas.sc                                  -0.005400659 0.001182908 27937 -4.565579  0.0000
    ## microclim.sc                             -0.000177716 0.000491809 27937 -0.361351  0.7178
    ## mass.sc                                   0.000607872 0.001129822 27937  0.538024  0.5906
    ## speed.sc                                 -0.001024927 0.001049094 27937 -0.976964  0.3286
    ## consumerfrac.sc                           0.000646884 0.000472995 27937  1.367633  0.1714
    ## nspp.sc                                   0.001410573 0.000845167 27937  1.668988  0.0951
    ## npp.sc                                   -0.000009999 0.000736324 27937 -0.013579  0.9892
    ## duration.sc                              -0.002183733 0.000743509 27937 -2.937065  0.0033
    ## human_bowler.sc                           0.000430216 0.000526376 27937  0.817316  0.4138
    ## temptrend_abs.sc:tsign1                  -0.005175693 0.002478651 27937 -2.088108  0.0368
    ## temptrend_abs.sc:tempave_metab.sc        -0.000104767 0.003811069 27937 -0.027490  0.9781
    ## temptrend_abs.sc:seas.sc                  0.003014583 0.002697391 27937  1.117592  0.2638
    ## temptrend_abs.sc:microclim.sc             0.001500967 0.001401801 27937  1.070742  0.2843
    ## temptrend_abs.sc:mass.sc                 -0.004239176 0.002094004 27937 -2.024435  0.0429
    ## temptrend_abs.sc:speed.sc                 0.004519428 0.002360550 27937  1.914566  0.0556
    ## temptrend_abs.sc:consumerfrac.sc          0.005996845 0.000962045 27937  6.233435  0.0000
    ## temptrend_abs.sc:nspp.sc                 -0.001707107 0.002042941 27937 -0.835612  0.4034
    ## tsign-1:thermal_bias.sc                  -0.001660932 0.001013973 27937 -1.638044  0.1014
    ## tsign1:thermal_bias.sc                   -0.002146300 0.000761943 27937 -2.816876  0.0049
    ## temptrend_abs.sc:npp.sc                   0.004756856 0.002013668 27937  2.362284  0.0182
    ## temptrend_abs.sc:duration.sc              0.001987202 0.001432726 27937  1.387008  0.1655
    ## temptrend_abs.sc:human_bowler.sc         -0.003252771 0.001319338 27937 -2.465458  0.0137
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.004204777 0.002132315 27937 -1.971930  0.0486
    ## temptrend_abs.sc:tsign1:thermal_bias.sc   0.002206323 0.001784090 27937  1.236666  0.2162
    ##  Correlation: 
    ##                                          (Intr) tmpt_. tsign1 tmpv_. ses.sc mcrcl. mss.sc spd.sc cnsmr. nspp.s npp.sc drtn.s hmn_b. tm_.:1 tmptrnd_bs.sc:t_. tmptrnd_bs.sc:ss. tmptrnd_bs.sc:mc. tmptrnd_bs.sc:ms. tmptrnd_bs.sc:sp. tmptrnd_bs.sc:c. tmptrnd_bs.sc:ns. t-1:_. ts1:_. tmptrnd_bs.sc:np. tmptrnd_bs.sc:d. tmptrnd_bs.sc:h_. t_.:-1
    ## temptrend_abs.sc                         -0.277                                                                                                                                                                                                                                                                                                  
    ## tsign1                                   -0.187  0.142                                                                                                                                                                                                                                                                                           
    ## tempave_metab.sc                          0.048  0.017  0.019                                                                                                                                                                                                                                                                                    
    ## seas.sc                                   0.195 -0.115 -0.090  0.197                                                                                                                                                                                                                                                                             
    ## microclim.sc                              0.049 -0.050  0.004 -0.236  0.068                                                                                                                                                                                                                                                                      
    ## mass.sc                                   0.083 -0.051  0.002  0.010  0.020  0.034                                                                                                                                                                                                                                                               
    ## speed.sc                                  0.066 -0.045 -0.064  0.071  0.055  0.022 -0.640                                                                                                                                                                                                                                                        
    ## consumerfrac.sc                          -0.003  0.019  0.033 -0.067 -0.007  0.007 -0.060  0.021                                                                                                                                                                                                                                                 
    ## nspp.sc                                  -0.106  0.077  0.060 -0.188  0.082 -0.076 -0.022  0.163  0.120                                                                                                                                                                                                                                          
    ## npp.sc                                   -0.069  0.088 -0.004 -0.092 -0.367 -0.187 -0.011  0.109 -0.040 -0.163                                                                                                                                                                                                                                   
    ## duration.sc                               0.054 -0.044 -0.159  0.082 -0.058 -0.019 -0.029 -0.006 -0.021 -0.319  0.041                                                                                                                                                                                                                            
    ## human_bowler.sc                          -0.032  0.035  0.007  0.116 -0.195  0.006 -0.031  0.058 -0.065 -0.064 -0.119  0.013                                                                                                                                                                                                                     
    ## temptrend_abs.sc:tsign1                   0.093 -0.295 -0.524 -0.059 -0.002 -0.001  0.023  0.012 -0.019 -0.026 -0.003  0.035 -0.013                                                                                                                                                                                                              
    ## temptrend_abs.sc:tempave_metab.sc        -0.014  0.025 -0.086 -0.428 -0.147  0.065  0.000  0.049  0.042  0.058  0.038  0.015 -0.063  0.082                                                                                                                                                                                                       
    ## temptrend_abs.sc:seas.sc                 -0.091  0.190  0.051 -0.130 -0.694 -0.053 -0.009 -0.029 -0.022 -0.042  0.280  0.012  0.154 -0.005  0.178                                                                                                                                                                                                
    ## temptrend_abs.sc:microclim.sc            -0.032  0.070 -0.013  0.153 -0.010 -0.734 -0.005 -0.022 -0.001  0.058  0.146 -0.002 -0.020 -0.010 -0.039             0.058                                                                                                                                                                              
    ## temptrend_abs.sc:mass.sc                 -0.010  0.130  0.009 -0.018  0.000  0.020 -0.567  0.365  0.073 -0.009 -0.037  0.047  0.037 -0.028  0.042            -0.024            -0.062                                                                                                                                                            
    ## temptrend_abs.sc:speed.sc                -0.008  0.117  0.042  0.023 -0.020 -0.043  0.297 -0.564 -0.001 -0.131 -0.033  0.018 -0.061  0.002 -0.161             0.073             0.060            -0.546                                                                                                                                          
    ## temptrend_abs.sc:consumerfrac.sc         -0.012 -0.041 -0.021  0.054  0.003  0.025  0.046 -0.006 -0.720 -0.102  0.018 -0.001  0.087  0.051 -0.049            -0.025            -0.046            -0.126            -0.018                                                                                                                        
    ## temptrend_abs.sc:nspp.sc                  0.036 -0.138 -0.027  0.080 -0.050  0.054 -0.009 -0.130 -0.076 -0.664  0.105  0.216  0.050  0.030 -0.098             0.029            -0.043             0.028             0.153             0.139                                                                                                      
    ## tsign-1:thermal_bias.sc                   0.031  0.001 -0.231  0.285 -0.155 -0.055 -0.014  0.014 -0.009 -0.114 -0.030  0.101  0.135  0.064 -0.060             0.102             0.032             0.017            -0.008             0.002            0.032                                                                                     
    ## tsign1:thermal_bias.sc                   -0.020  0.002  0.138  0.485 -0.167 -0.163 -0.024  0.022 -0.005 -0.116 -0.154  0.113  0.178 -0.073 -0.168             0.063             0.087             0.021             0.007             0.002            0.025             0.374                                                                   
    ## temptrend_abs.sc:npp.sc                   0.028 -0.130  0.012  0.051  0.231  0.170 -0.014 -0.047  0.020  0.093 -0.663 -0.022  0.073 -0.039 -0.047            -0.413            -0.332             0.083             0.000            -0.010           -0.107            -0.005  0.102                                                            
    ## temptrend_abs.sc:duration.sc             -0.113  0.263  0.052  0.093  0.007 -0.044  0.036  0.006  0.034  0.266 -0.027 -0.514  0.043 -0.121 -0.025            -0.029             0.065            -0.074            -0.026             0.018           -0.373            -0.025 -0.049  0.032                                                     
    ## temptrend_abs.sc:human_bowler.sc          0.032 -0.035 -0.023 -0.084  0.145  0.000  0.037 -0.052  0.145  0.047  0.081  0.002 -0.749  0.024  0.119            -0.171            -0.047            -0.028             0.054            -0.170           -0.057            -0.087 -0.141 -0.099            -0.038                                   
    ## temptrend_abs.sc:tsign-1:thermal_bias.sc -0.006 -0.026  0.047 -0.107  0.116  0.057  0.010  0.013 -0.001  0.007 -0.003 -0.022 -0.081  0.051  0.321            -0.149            -0.092            -0.003            -0.054             0.034            0.010            -0.468 -0.188  0.030             0.033            0.096                  
    ## temptrend_abs.sc:tsign1:thermal_bias.sc   0.015 -0.029 -0.108 -0.194  0.066  0.082  0.014  0.041 -0.018  0.014  0.098 -0.036 -0.137  0.124  0.535            -0.056            -0.198             0.003            -0.135             0.053           -0.023            -0.166 -0.565 -0.095             0.016            0.171             0.382
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -5.51920044 -0.21923803 -0.01263607  0.21747792  6.06647854 
    ## 
    ## Number of Observations: 28054
    ## Number of Groups: 
    ##               STUDY_ID rarefyID %in% STUDY_ID 
    ##                     90                  28054

### Plot the realm-specific coefficients

Also uses the full models across all realms

``` r
coefs1 <- summary(modTfullHornrem0)$tTable
coefs2 <- summary(modTfullHornTerr)$tTable
coefs3 <- summary(modTfullHornFresh)$tTable
coefs4 <- summary(modTfullHornMar)$tTable

varstoplot <- unique(c(rownames(coefs1), rownames(coefs2), rownames(coefs3), rownames(coefs4)))
varstoplot <- varstoplot[which(!grepl('Intercept', varstoplot) | grepl(':', varstoplot))] # vars to plot

rows1_1 <- which(rownames(coefs1) %in% varstoplot) # rows in coefs
rows1_2 <- which(rownames(coefs2) %in% varstoplot)
rows1_3 <- which(rownames(coefs3) %in% varstoplot)
rows1_4 <- which(rownames(coefs4) %in% varstoplot)
xlims <- range(c(coefs1[rows1_1,1] - coefs1[rows1_1,2], coefs1[rows1_1,1] + coefs1[rows1_1,2], 
                 coefs2[rows1_2,1] - coefs2[rows1_2,2], coefs2[rows1_2,1] + coefs2[rows1_2,2], 
                 coefs3[rows1_3,1] - coefs3[rows1_3,2], coefs3[rows1_3,1] + coefs3[rows1_3,2],
                 coefs4[rows1_4,1] - coefs4[rows1_4,2], coefs4[rows1_4,1] + coefs4[rows1_4,2]))


cols <- brewer.pal(4, 'Dark2') # for full, terr, fresh, mar
pchs <- c(1, 16, 16, 16)
offs <- c(0.1, 0, -0.1, -0.2) # offset vertically for each model


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
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[1], length(varstoplot) + 1 - i + offs[1]), col = cols[1])
  }
  if(varstoplot[i] %in% rownames(coefs2)){
    x = coefs2[rownames(coefs2) == varstoplot[i], 1]
    se = coefs2[rownames(coefs2) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[2], pch = pchs[2], col = cols[2])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[2], length(varstoplot) + 1 - i + offs[2]), col = cols[2])
  }
  if(varstoplot[i] %in% rownames(coefs3)){
    x = coefs3[rownames(coefs3) == varstoplot[i], 1]
    se = coefs3[rownames(coefs3) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[3], pch = pchs[3], col = cols[3])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[3], length(varstoplot) + 1 - i + offs[3]), col = cols[3])
  }
  if(varstoplot[i] %in% rownames(coefs4)){
    x = coefs4[rownames(coefs4) == varstoplot[i], 1]
    se = coefs4[rownames(coefs4) == varstoplot[i], 2]
    points(x, length(varstoplot) + 1 - i + offs[4], pch = pchs[4], col = cols[4])
    lines(x = c(x-se, x+se), y = c(length(varstoplot) + 1 - i + offs[4], length(varstoplot) + 1 - i + offs[4]), col = cols[4])
  }
}
legend('bottomleft', col = cols, pch = pchs, lwd = 1, legend = c('All', 'Terestrial', 'Freshwater', 'Marine'))
```

![](turnover_vs_temperature_MEmodels_files/figure-gfm/plot%20realm%20mods-1.png)<!-- -->

\[End text in hopes this helps the last figure show up when knitted\]

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
# Load BioTime trends calculated by calc_turnover.Rmd
trends <- readRDS('temp/trendstemp2.rds')

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

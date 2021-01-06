Turnover covariate data prep and visualization
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
# biotime taxa category and other info
load('data/biotime_blowes/bt_malin.Rdata') # load bt_malin
btinfo <- data.table(bt_malin); rm(bt_malin)
btinfo <- btinfo[!duplicated(rarefyID), .(rarefyID, rarefyID_x, rarefyID_y, Biome, taxa_mod, REALM, STUDY_ID)]

# biotime community dissimilarity data
load('data/biotime_blowes/all_pairs_beta.RData') # load rarefied_beta_medians
bt <- data.table(rarefied_beta_medians); rm(rarefied_beta_medians)
bt[, year1 := as.numeric(year1)] # not sure why it gets read in as character
bt[, year2 := as.numeric(year2)]
bt[, Horn := 1- Hornsim] # convert similarity to dissimilarity
bt[, Hornsim := NULL]
bt <- merge(bt, btinfo[, .(rarefyID, rarefyID_x, rarefyID_y, REALM, STUDY_ID, Biome, taxa_mod)]) # add lat/lon and realm

# Temperature average, changes, and seasonality
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

## Add covariates to BT data

``` r
# add covariates
bt <- merge(bt, temperature, all.x = TRUE, by = c('rarefyID', 'year1', 'year2')) # temperature ave, ave metabolic, change, and seasonality
bt <- merge(bt, microclim[, .(rarefyID, microclim = Temp_sd20km)], all.x = TRUE, by = 'rarefyID') # microclimates
bt <- merge(bt, npp, all.x = TRUE, by = 'rarefyID') # npp
bt <- merge(bt, bs[, .(rarefyID, mass_mean_weight, mass_sd_weight)], all.x = TRUE) # body size mass (g)
bt <- merge(bt, speed[, .(rarefyID, speed_mean_weight, speed_sd_weight)], all.x = TRUE) # speed (km/hr)
bt <- merge(bt, lsp[, .(rarefyID, lifespan_mean_weight, lifespan_sd_weight)], all.x = TRUE) # lifespan (yr)
bt <- merge(bt, cti[, .(rarefyID, thermal_bias)], all.x = TRUE) # thermal bias (degC)
bt <- merge(bt, consfrac[, .(rarefyID, consfrac)], all.x = TRUE) # fraction consumers
bt <- merge(bt, rich, all.x = TRUE) # species richness
bt <- merge(bt, endofrac[, .(rarefyID, endofrac)], all.x = TRUE) # endotherm vs. ectotherm
bt <- merge(bt, human[, .(rarefyID, human_bowler = atc, human_venter = hfp, human_halpern = himp)], all.x = TRUE) # human impact
bt <- merge(bt, veg[, .(rarefyID, veg = veg)], all.x = TRUE) # vegetation index
bt[REALM == 'Marine', veg := 0] # veg index is 0 at sea
```

## Write out

``` r
write.csv(bt, gzfile('output/turnover_w_covariates.csv.gz'), row.names = FALSE)
```

# Compare covariates across realms

``` r
i <- bt[, !duplicated(rarefyID)]; sum(i)
```

    ## [1] 53467

``` r
par(mfrow=c(5,3))
beanplot(rarefyID_y ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Latitude (degN)', ll = 0.05)
beanplot(tempave ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature (degC)', ll = 0.05)
beanplot(tempave_metab ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Metabolic Temperature (degC)', ll = 0.05, bw = 'nrd0') # nrd0 bandwidth to calculation gap
beanplot(seas ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Seasonality (degC)', ll = 0.05)
beanplot(microclim ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Microclimates (degC)', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(tempchange ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature change (degC)', ll = 0.05)
beanplot(mass_mean_weight ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Mass (g)', ll = 0.05, log = 'y')
beanplot(speed_mean_weight +1 ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Speed (km/hr)', ll = 0.05, log = 'y')
beanplot(lifespan_mean_weight ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Lifespan (yr)', ll = 0.05, log = 'y')
#beanplot(consfrac ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Consumers (fraction)', ll = 0.05, log = '') # too sparse
#beanplot(endofrac ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Endotherms (fraction)', ll = 0.05, log = '') # too sparse
beanplot(Nspp ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Number of species', ll = 0.05, log = 'y')
beanplot(thermal_bias ~ REALM, data = bt[i & !is.na(thermal_bias),], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Thermal bias (degC)', ll = 0.05)
beanplot(npp ~ REALM, data = bt[i,], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'NPP', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(veg ~ REALM, data = bt[i & REALM !='Marine',], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'NPP', ll = 0.05)
```

![](assemble_turnover_covariates_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->

Marine are in generally warmer locations (seawater doesn’t freeze)
Marine have much lower seasonality. Marine and freshwater have some very
small masses (plankton), but much of dataset is similar to terrestrial.
Marine has a lot of slow, crawling organisms, but land has plants. Land
also has birds (fast).

# Frequency of time differences

``` r
invisible(bt[, hist(year2 - year1, breaks = seq(0.5, 120.5, by = 1))])
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20time%20differences-1.png)<!-- -->

# Plot dissimilarity vs. time

![](assemble_turnover_covariates_files/figure-gfm/plot%20diss%20vs%20time-1.png)<!-- -->

# Plot dissimilarity vs. temperature change

## All temporal scales

![](assemble_turnover_covariates_files/figure-gfm/plot%20diss%20vs%20temp%20change-1.png)<!-- -->

## Only 1 year differences

![](assemble_turnover_covariates_files/figure-gfm/plot%20diss%20vs%20temp%20change%201%20yr-1.png)<!-- -->

## Only 10 year differences

![](assemble_turnover_covariates_files/figure-gfm/plot%20diss%20vs%20temp%20change%2010%20yr-1.png)<!-- -->

## Only 20 year differences

![](assemble_turnover_covariates_files/figure-gfm/plot%20diss%20vs%20temp%20change%2020%20yr-1.png)<!-- -->

# Temperature trend and time-series length

``` r
ggplot(bt, aes(year2 - year1, tempchange)) +
  geom_point(alpha = 0.2) +
  geom_smooth()
```

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 713895 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 713895 rows containing missing values (geom_point).

![](assemble_turnover_covariates_files/figure-gfm/ts%20length%20and%20temp%20trend-1.png)<!-- -->
\# Make a data.table summarized by rarefyID

``` r
btave <- bt[, .(Jtu = mean(Jtu), minyrBT = min(year1), maxyrBT = max(year2), nyrBT = length(Jtu) + 1,
                REALM = unique(REALM), tempchange = mean(tempchange), tempave = mean(tempave),
                tempave_metab = mean(tempave_metab), seas = mean(seas), microclim = mean(microclim),
                mass_mean_weight = mean(mass_mean_weight), speed_mean_weight = mean(speed_mean_weight),
                lifespan_mean_weight = mean(lifespan_mean_weight), consfrac = mean(consfrac),
                endofrac = mean(endofrac), Nspp = mean(Nspp), thermal_bias = mean(thermal_bias),
                npp = mean(npp), veg = mean(veg), human_bowler = mean(human_bowler),
                human_venter = mean(human_venter), human_halpern = mean(human_halpern)),  by = rarefyID]
```

# Plot dissimilarity vs. explanatory variables

Lines are ggplot smoother fits Just Jtu for now, averaged within
rarefyID
![](assemble_turnover_covariates_files/figure-gfm/plot%20diss%20v%20explanatory%20vars-1.png)<!-- -->

Strong trends with temperature change, but trends are pretty symmetric
around no trend in temperature, which implies warming or cooling drives
similar degree of community turnover. Some indication of less turnover
for larger organisms (mass) Higher turnover on land with higher
seasonality? More turnover for shorter-lived organisms? No really clear
differences among realms.

# A bit more prep for visualizing covariate distributions

## Add Useful variables

``` r
# realm that combined Terrestrial and Freshwater, for interacting with human impact
btave[, REALM2 := REALM]
levels(btave$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")
```

## Log-transform some variables, then center and scale.

``` r
btave[, tempave.sc := scale(tempave)]
btave[, tempave_metab.sc := scale(tempave_metab)]
btave[, seas.sc := scale(seas)]
btave[, microclim.sc := scale(log(microclim))]
bt[, tempchange.sc := scale(tempchange, center = FALSE)] # all the raw data
bt[, tempchange_abs.sc := scale(abs(tempchange), center = FALSE)] # do not center, so that 0 is still 0 temperature change
btave[, tempchange.sc := scale(tempchange, center = FALSE)] # averages within studies
btave[, tempchange_abs.sc := scale(abs(tempchange), center = FALSE)]
btave[, mass.sc := scale(log(mass_mean_weight))]
btave[, speed.sc := scale(log(speed_mean_weight+1))]
btave[, lifespan.sc := scale(log(lifespan_mean_weight))]
btave[, consumerfrac.sc := scale(consfrac)]
btave[, endothermfrac.sc := scale(endofrac)]
btave[, nspp.sc := scale(log(Nspp))]
btave[, thermal_bias.sc := scale(thermal_bias)]
btave[, npp.sc := scale(log(npp))]
btave[, veg.sc := scale(log(veg+1))]
btave[, human_bowler.sc := scale(log(human_bowler+1)), by = REALM2] # separate scaling by realm
btave[REALM2 == 'TerrFresh', human_footprint.sc := scale(log(human_venter+1))]
btave[REALM2 == 'Marine', human_footprint.sc := scale(log(human_halpern))]
```

# Do the variables look ok?

## Unscaled

``` r
# histograms to examine
cexmain = 0.6
par(mfrow = c(5,4))
invisible(btave[, hist(minyrBT, main = 'Start year', cex.main = cexmain)])
invisible(btave[, hist(maxyrBT - minyrBT, main = 'Duration (years)', cex.main = cexmain)])
invisible(btave[, hist(nyrBT, main = 'Number of sampled years', cex.main = cexmain)])
invisible(btave[, hist(mass_mean_weight, main = 'Mass (g)', cex.main = cexmain)])
invisible(btave[, hist(speed_mean_weight, main = 'Speed (km/hr)', cex.main = cexmain)])
invisible(btave[, hist(lifespan_mean_weight, main = 'Lifespan (yr)', cex.main = cexmain)])
invisible(btave[, hist(tempave_metab, main = 'Metabolic temperature (°C)', cex.main = cexmain)])
invisible(btave[, hist(consfrac, main = 'Consumers (fraction)', cex.main = cexmain)])
invisible(btave[, hist(endofrac, main = 'Endotherms (fraction)', cex.main = cexmain)])
invisible(btave[, hist(tempave, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(bt[, hist(tempchange, main = 'Temperature trend (°C/yr)', cex.main = cexmain)]) # all the raw data
invisible(btave[, hist(seas, main = 'Seasonality (°C)', cex.main = cexmain)])
invisible(btave[, hist(microclim, main = 'Microclimates (°C)', cex.main = cexmain)])
invisible(btave[, hist(Nspp, main = 'Species richness', cex.main = cexmain)])
invisible(btave[, hist(thermal_bias, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(btave[, hist(npp, main = 'Net primary productivity', cex.main = cexmain)])
invisible(btave[, hist(veg, main = 'Vegetation index', cex.main = cexmain)])
invisible(btave[, hist(human_bowler, main = 'Human impact score (Bowler)', cex.main = cexmain)])
invisible(btave[, hist(human_venter, main = 'Human impact score (Venter)', cex.main = cexmain)])
invisible(btave[, hist(human_halpern, main = 'Human impact score (Halpern)', cex.main = cexmain)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20unscaled-1.png)<!-- -->

## Scaled

``` r
# histograms to examine
cexmain = 0.6
par(mfrow = c(5,4))
invisible(btave[, hist(tempave.sc, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(btave[, hist(tempave_metab.sc, main = 'Metabolic temperature (°C)', cex.main = cexmain)])
invisible(btave[, hist(seas.sc, main = 'Seasonality (°C)', cex.main = cexmain)])
invisible(btave[, hist(microclim.sc, main = 'log Microclimates (°C)', cex.main = cexmain)])
invisible(bt[, hist(tempchange.sc, main = 'Temperature trend (°C/yr)', cex.main = cexmain)])
invisible(bt[, hist(tempchange_abs.sc, main = 'abs(Temperature trend) (°C/yr)', cex.main = cexmain)])
invisible(btave[, hist(mass.sc, main = 'log Mass (g)', cex.main = cexmain)])
invisible(btave[, hist(speed.sc, main = 'log Speed (km/hr)', cex.main = cexmain)])
invisible(btave[, hist(lifespan.sc, main = 'log Lifespan (yr)', cex.main = cexmain)])
invisible(btave[, hist(consumerfrac.sc, main = 'Consumers (fraction)', cex.main = cexmain)])
invisible(btave[, hist(endothermfrac.sc, main = 'Endotherms (fraction)', cex.main = cexmain)])
invisible(btave[, hist(nspp.sc, main = 'log Species richness', cex.main = cexmain)])
invisible(btave[, hist(thermal_bias.sc, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(btave[, hist(npp.sc, main = 'log Net primary productivity', cex.main = cexmain)])
invisible(btave[, hist(veg.sc, main = 'log Vegetation index', cex.main = cexmain)])
invisible(btave[, hist(human_bowler.sc, main = 'log Human impact score (Bowler)', cex.main = cexmain)])
invisible(btave[, hist(human_footprint.sc, main = 'log Human impact score (Venter & Halpern)', cex.main = cexmain)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20scaled-1.png)<!-- -->

# Check correlations among variables. Pearson’s r is on the lower diagonal.

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
pairs(formula = ~ tempave.sc + tempave_metab.sc + seas.sc + microclim.sc + tempchange.sc + tempchange_abs.sc + mass.sc + speed.sc + lifespan.sc + consumerfrac.sc + endothermfrac.sc + nspp.sc + thermal_bias.sc + npp.sc + veg.sc + human_bowler.sc + human_footprint.sc, data = btave, gap = 1/10, cex = 0.2, col = '#00000022', lower.panel = panel.cor)
```

![](assemble_turnover_covariates_files/figure-gfm/pairs-1.png)<!-- -->

Mass and lifespan look tightly correlated, but r only 0.56…?
Tempave\_metab and lifespan don’t look tightly correlated, but r= -0.81
Tempave\_metab and speed don’t look tightly correlated, but r= -0.83
Lifespan and speed don’t look tightly correlated, but r = 0.73

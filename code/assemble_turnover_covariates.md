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

``` r
# biotime taxa category and other info
load(here::here('data', 'biotime_blowes', 'bt_malin.Rdata')) # load bt_malin
btinfo <- data.table(bt_malin); rm(bt_malin)
btinfo <- btinfo[!duplicated(rarefyID), .(rarefyID, rarefyID_x, rarefyID_y, Biome, taxa_mod, REALM, STUDY_ID)]

# biotime community dissimilarity data
load(here::here('data', 'biotime_blowes', 'all_pairs_beta.Rdata')) # load rarefied_beta_medians
bt <- data.table(rarefied_beta_medians); rm(rarefied_beta_medians)
bt[, year1 := as.numeric(year1)] # not sure why it gets read in as character
bt[, year2 := as.numeric(year2)]
bt[, Horn := 1- Hornsim] # convert similarity to dissimilarity
bt[, Hornsim := NULL]
bt <- merge(bt, btinfo[, .(rarefyID, rarefyID_x, rarefyID_y, REALM, STUDY_ID, Biome, taxa_mod)]) # add lat/lon and realm

# Temperature average, changes, and seasonality
temperature <- fread(here('output', 'temperature_byrarefyID.csv.gz'))

# microclimates
microclim <- fread(here('output', 'microclimates.csv.gz'), drop = 1)

# NPP
npp <- fread(here('output', 'npplandocean.csv.gz'))

# Body size
bs <- fread(here('output', 'mass_byrarefyID.csv.gz'), drop = 1)
bs[, ':='(STUDY_ID = NULL, REALM = NULL, taxa_mod = NULL)] # remove unnecessary columns 

# Mobility
speed <- fread(here('output', 'speed_byrarefyID.csv.gz'), drop = 1)
speed[, ':='(STUDY_ID = NULL, REALM = NULL, taxa_mod = NULL)] # remove unnecessary columns 

# Lifespan
lsp <- fread(here('output', 'lifespan_byrarefyID.csv.gz'))

# CTI
cti <- fread(here('output', 'cti_byrarefyID.csv.gz'))
    
# consumer vs. producer
consfrac <- fread(here('output', 'consfrac_byrarefyID.csv.gz'))

# richness
rich <- fread(here('output', 'richness_by_rarefyID.csv.gz')) # number of species

# endotherm vs. ectotherm
endofrac <- fread(here('output', 'endofrac_byrarefyID.csv.gz')) # endotherm vs. ectotherm classifications

# human impact
human <- fread(here('output', 'humanimpact_by_rarefyID.csv.gz'))

# %veg
veg <- as.data.table(readRDS(here('output', 'vct_by_rarefyID.rds')))
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

Only if file doesn’t yet exist

``` r
if(!file.exists(here('output', 'turnover_w_covariates.csv.gz'))){
  write.csv(bt, gzfile(here('output', 'turnover_w_covariates.csv.gz')), row.names = FALSE)
}
```

# A bit more prep for visualizing covariate distributions

## Make a data.table summarized by rarefyID

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

## Add useful variables

``` r
# realm that combined Terrestrial and Freshwater, for interacting with human impact
btave[, REALM2 := REALM]
levels(btave$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")
```

## Try a logit-transform of the response variables

``` r
logit = function(x) return(log(x/(1-x)))

# have to adjust variables away from 0 and 1
jtumin = bt[Jtu > 0, min(Jtu)]
jtumax = bt[Jtu < 1, max(Jtu)]
bt[Jtu > 0 & Jtu < 1, Jtulogit := logit(Jtu)]
bt[Jtu == 0, Jtulogit := logit(jtumin)]
bt[Jtu == 1, Jtulogit := logit(jtumax)]

jbetamin = bt[Jbeta > 0, min(Jbeta)]
jbetamax = bt[Jbeta < 1, max(Jbeta)]
bt[Jbeta > 0 & Jbeta < 1, Jbetalogit := logit(Jbeta)]
bt[Jbeta == 0, Jbetalogit := logit(jbetamin)]
bt[Jbeta == 1, Jbetalogit := logit(jbetamax)]

hornmin = bt[Horn > 0, min(Horn)]
hornmax = bt[Horn < 1, max(Horn)]
bt[Horn > 0 & Horn < 1, Hornlogit := logit(Horn)]
bt[Horn == 0, Hornlogit := logit(hornmin)]
bt[Horn == 1, Hornlogit := logit(hornmax)]

bt[, summary(Jtulogit)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## -4.1821 -1.5656 -0.4055 -0.3728  0.6931  3.9120

``` r
bt[, summary(Jbetalogit)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## -5.5922 -0.4855  0.3365  1.0419  1.8718  5.7301

``` r
bt[, summary(Hornlogit)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##  -27.46   -1.85   -0.12    2.90    3.38   19.55   41879

## Log-transform some explanantory variables, then center and scale.

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

# Check variable distributions

## Response variables

``` r
bt[, summary(Jbeta)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.3810  0.5833  0.6043  0.8667  1.0000

``` r
bt[, summary(Jtu)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.1728  0.4000  0.4381  0.6667  1.0000

``` r
bt[, summary(Horn)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##    0.00    0.14    0.47    0.52    0.97    1.00   41879

``` r
# fraction 0 or 1
bt[, sum(Jbeta == 0)/.N]
```

    ## [1] 0.02257946

``` r
bt[, sum(Jtu == 0)/.N]
```

    ## [1] 0.1925047

``` r
bt[!is.na(Horn), sum(Horn == 0)/.N]
```

    ## [1] 0.007767746

``` r
bt[, sum(Jbeta == 1)/.N]
```

    ## [1] 0.1709619

``` r
bt[, sum(Jtu == 1)/.N]
```

    ## [1] 0.1709619

``` r
bt[!is.na(Horn), sum(Horn == 1)/.N]
```

    ## [1] 0.1747303

``` r
# histograms
invisible(bt[, hist(Jbeta, main = 'Jaccard total', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-1.png)<!-- -->

``` r
invisible(bt[, hist(Jtu, main = 'Jaccard turnover', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-2.png)<!-- -->

``` r
invisible(bt[, hist(Horn, main = 'Morisita-Horn', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-3.png)<!-- -->

``` r
invisible(bt[, hist(Jbetalogit, main = 'logit(Jaccard total)', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-4.png)<!-- -->

``` r
invisible(bt[, hist(Jtulogit, main = 'logit(Jaccard turnover)', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-5.png)<!-- -->

``` r
invisible(bt[, hist(Hornlogit, main = 'logit(Morisita-Horn)', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-6.png)<!-- -->

## Unscaled covariates

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

## Scaled covariates

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

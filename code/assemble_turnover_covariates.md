Pairwise dissimilarity covariate data prep and visualization
================

Trims based on data quality, separates the data points to use for
3/5/10/20 year comparisons, and adds covariate data.

``` r
if(Sys.info()[[4]]=='annotate.sebs.rutgers.edu'){ # for gam smoother
library('mgcv', lib.loc = '/usr/lib64/R/library') # when running on Annotate. Need to load 1.8-26, not 1.8-33.
} else {
  library(mgcv)
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
source(here('code', 'util.R')) # simple utility functions
```

# Load data

``` r
# biotime taxa category and other info
load(here::here('data', 'biotime_blowes', 'bt_malin.Rdata')) # load bt_malin
btinfo <- data.table(bt_malin); rm(bt_malin)
btinfo[, Nave := mean(N), by = rarefyID] # average number of individuals in the timeseries
btinfo <- btinfo[!duplicated(rarefyID), .(rarefyID, rarefyID_x, rarefyID_y, Biome, taxa_mod, Nave, REALM, STUDY_ID)]

# biotime community dissimilarity data (all pairs)
load(here::here('data', 'biotime_blowes', 'all_pairs_beta.Rdata')) # load rarefied_beta_medians
bt <- data.table(rarefied_beta_medians); rm(rarefied_beta_medians)
bt[, year1 := as.numeric(year1)] # not sure why it gets read in as character
bt[, year2 := as.numeric(year2)]
bt[, Horn := 1- Hornsim] # convert similarity to dissimilarity
bt[, Hornsim := NULL]
bt <- merge(bt, btinfo) # add lat/lon and realm

# Temperature average, changes
temperature <- fread(here('output', 'temperature_byrarefyID.csv.gz')) # from assemble_temp.Rmd

# microclimates
microclim <- fread(here('output', 'microclimates.csv.gz'), drop = 1) # from assemble_microclimates.R

# CTI
cti <- fread(here('output', 'cti_byrarefyID.csv.gz')) # from assemble_temperature_bias.R
    
# richness
rich <- fread(here('output', 'richness_by_rarefyID.csv.gz')) # number of species. from extract_richness.R

# human impact
human <- fread(here('output', 'humanimpact_by_rarefyID.csv.gz'))
```

## Add covariates to BT data

``` r
# add covariates
bt <- merge(bt, temperature, all.x = TRUE, by = c('rarefyID', 'year1', 'year2')) # temperature ave, change
bt <- merge(bt, microclim[, .(rarefyID, microclim = Temp_sd20km)], all.x = TRUE, by = 'rarefyID') # microclimates
bt <- merge(bt, cti[, .(rarefyID, thermal_bias)], all.x = TRUE) # thermal bias (degC)
bt <- merge(bt, rich, all.x = TRUE) # species richness
bt <- merge(bt, human[, .(rarefyID, human_bowler = atc)], all.x = TRUE) # human impact
```

## Trim data

At least 5 species and at least 10 individuals on average.

``` r
cat('original number of rows'); norig <- nrow(bt); norig
```

    ## original number of rows

    ## [1] 1304150

``` r
# trim out timeseries with too few species
cat('rows with <5 spp'); bt[Nspp < 5 | is.na(Nspp), .N]
```

    ## rows with <5 spp

    ## [1] 35975

``` r
bt <- bt[Nspp >= 5, ]
  
# or few individuals
cat('rows with <10 indivs'); bt[Nave < 10, .N]
```

    ## rows with <10 indivs

    ## [1] 31899

``` r
bt <- bt[Nave >= 10 | is.na(Nave), ]


cat('rows left'); nrow(bt)
```

    ## rows left

    ## [1] 1236276

``` r
cat('fraction kept'); nrow(bt)/norig
```

    ## fraction kept

    ## [1] 0.9479554

``` r
cat('number of studies'); bt[, length(unique(STUDY_ID))]
```

    ## number of studies

    ## [1] 319

``` r
cat('number of timeseries'); bt[, length(unique(rarefyID))]
```

    ## number of timeseries

    ## [1] 42255

### Choose the appropriate temperature value for each study and duration

Average temperature for the full duration. Thiel-Sen temperature slope.

``` r
btall <- bt # legacy from when we divided dataset into different ts durations (3, 5, 10, 20 years)
btall[, tempave := picklongest(year1, year2, tempave), by = rarefyID]
btall[, tempchange := median(tempchange/(year2 - year1)), by = rarefyID] # Thiel-Sen estimator of the slope: median of all pairwise differences
```

## Set up useful variables and transformations

``` r
# set realm order
btall[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = FALSE)]

# realm that combines Terrestrial and Freshwater, for interacting with human impact
btall[, REALM2 := REALM]
levels(btall$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")


# duration
btall[, duration := year2 - year1]

# group Marine invertebrates/plants in with All
btall[, taxa_mod2 := taxa_mod]
btall[taxa_mod == 'Marine invertebrates/plants', taxa_mod2 := 'All']


btall[, tsign := signneg11(tempchange)]

#######################
## Transformations

### Adjust response away from 0-1
btall[, ':='(Jtu.sc = transform01(Jtu), Jne.sc = transform01(Jne), 
           Jbeta.sc = transform01(Jbeta), Horn.sc = transform01(Horn))]


### Log-transform some variables, then center and scale. 
btall[, ':='(tempave.sc = scale(tempave),
           microclim.sc = scale(log(microclim)),
           tempchange.sc = scale(tempchange, center = TRUE), 
           tempchange_abs.sc = scale(abs(tempchange), center = TRUE),
           nspp.sc = scale(log(Nspp)),
           thermal_bias.sc = scale(thermal_bias),
           duration.sc = scale(duration),
           durationlog.sc = scale(log(duration)))]
btall[, human_bowler.sc := scale(log(human_bowler+1)), by = REALM2] # separate scaling by realm


## save the centering and scaling
scalingall <- data.frame(center = t(t(sapply(btall, attr, 'scaled:center'))),
                 scale = t(t(sapply(btall, attr, 'scaled:scale'))))
scalingall$var <- rownames(scalingall)
scalingall <- scalingall[!vapply(scalingall$center, is.null, TRUE) | !vapply(scalingall$scale, is.null, TRUE),]
scalingall$center[vapply(scalingall$center, is.null, TRUE)] <- NA # turn null to NA
scalingall$center <- vapply(scalingall$center, paste, collapse = ", ", character(1L)) # flatten list to character. not sure why it's a list.
scalingall$scale <- vapply(scalingall$scale, paste, collapse = ", ", character(1L))
```

## Dataset sizes

``` r
btall[, .N]
```

    ## [1] 1236276

``` r
btall[, length(unique(STUDY_ID)), by = REALM]
```

    ##          REALM  V1
    ## 1:      Marine 130
    ## 2: Terrestrial 168
    ## 3:  Freshwater  21

``` r
btall[, length(unique(rarefyID)), by = REALM]
```

    ##          REALM    V1
    ## 1:      Marine 38451
    ## 2: Terrestrial  3159
    ## 3:  Freshwater   645

## Write out

Only if file doesn’t yet exist

``` r
if(!file.exists(here('output', 'turnover_w_covariates.csv.gz'))){
  write.csv(btall[,.(STUDY_ID, rarefyID, rarefyID_x, rarefyID_y, REALM, Biome, taxa_mod, taxa_mod2, year1, year2, duration, 
                   Jtu.sc, Jne.sc, Jbeta.sc, Horn.sc, 
                   tempave.sc, 
                   microclim.sc, tempchange, tsign, tempchange.sc, tempchange_abs.sc,
                   nspp.sc, thermal_bias.sc,
                   duration.sc, durationlog.sc, human_bowler.sc, REALM2)], 
            gzfile(here('output', 'turnover_w_covariates.csv.gz')), row.names = FALSE)
  write.csv(scalingall, here('output','turnover_w_covariates_scaling.csv'), row.names = FALSE)
}
```

# A bit more prep for visualizing covariate distributions

## Make a data.table summarized by rarefyID

``` r
btaveall <- btall[, .(Jtuslope = lmNAcoef(Jtu, duration), Jbetaslope = lmNAcoef(Jbeta, duration), Hornslope = lmNAcoef(Horn, duration), 
                    minyrBT = min(year1), maxyrBT = max(year2), nyrBT = length(Jtu) + 1,
                    REALM = unique(REALM), REALM2 = unique(REALM2), rarefyID_y = mean(rarefyID_y),
                    tempchange = mean(tempchange), tempave = mean(tempave),
                    microclim = mean(microclim),
                    Nspp = mean(Nspp), thermal_bias = mean(thermal_bias),
                    human_bowler = mean(human_bowler),
                    tempave.sc = mean(tempave.sc),
                    tempchange.sc = mean(tempchange.sc), tempchange_abs.sc = mean(tempchange_abs.sc),
                    microclim.sc = mean(microclim.sc),
                    nspp.sc = mean(nspp.sc), thermal_bias.sc = mean(thermal_bias.sc),
                    human_bowler.sc = mean(human_bowler.sc)),
                by = rarefyID]
```

# Check variable distributions

## Response variables

``` r
btall[, summary(Jbeta)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.3750  0.5714  0.5970  0.8333  1.0000

``` r
btall[, summary(Jtu)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.1818  0.3810  0.4283  0.6667  1.0000

``` r
btall[, summary(Horn)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##    0.00    0.14    0.45    0.51    0.95    1.00   41507

``` r
# fraction 0 or 1
btall[, sum(Jbeta == 0)/.N]
```

    ## [1] 0.01621078

``` r
btall[, sum(Jtu == 0)/.N]
```

    ## [1] 0.1831282

``` r
btall[!is.na(Horn), sum(Horn == 0)/.N]
```

    ## [1] 0.002628123

``` r
btall[, sum(Jbeta == 1)/.N]
```

    ## [1] 0.1486165

``` r
btall[, sum(Jtu == 1)/.N]
```

    ## [1] 0.1486165

``` r
btall[!is.na(Horn), sum(Horn == 1)/.N]
```

    ## [1] 0.1517942

``` r
# histograms
invisible(btall[, hist(Jbeta, main = 'Jaccard total', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-1.png)<!-- -->

``` r
invisible(btall[, hist(Jtu, main = 'Jaccard turnover', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-2.png)<!-- -->

``` r
invisible(btall[, hist(Horn, main = 'Morisita-Horn', breaks = 80)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20response-3.png)<!-- -->

## Unscaled covariates

``` r
# histograms to examine
cexmain = 0.6
par(mfrow = c(5,4))
invisible(btaveall[, hist(minyrBT, main = 'Start year', cex.main = cexmain)])
invisible(btaveall[, hist(maxyrBT - minyrBT, main = 'Duration (years)', cex.main = cexmain)])
invisible(btaveall[, hist(nyrBT, main = 'Number of sampled years', cex.main = cexmain)])
invisible(btaveall[, hist(tempave, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(tempchange, main = 'Temperature trend (°C/yr)', cex.main = cexmain)]) # all the raw data
invisible(btaveall[, hist(microclim, main = 'Microclimates (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(Nspp, main = 'Species richness', cex.main = cexmain)])
invisible(btaveall[, hist(thermal_bias, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(human_bowler, main = 'Human impact score (Bowler)', cex.main = cexmain)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20unscaled-1.png)<!-- -->

## Scaled covariates

``` r
# histograms to examine
cexmain = 0.6
par(mfrow = c(5,4))
invisible(btaveall[, hist(tempave.sc, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(microclim.sc, main = 'log Microclimates (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(tempchange.sc, main = 'Temperature trend (°C/yr)', cex.main = cexmain)])
invisible(btaveall[, hist(tempchange_abs.sc, main = 'abs(Temperature trend) (°C/yr)', cex.main = cexmain)])
invisible(btaveall[, hist(nspp.sc, main = 'log Species richness', cex.main = cexmain)])
invisible(btaveall[, hist(thermal_bias.sc, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(human_bowler.sc, main = 'log Human impact score (Bowler)', cex.main = cexmain)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20scaled-1.png)<!-- -->

# Check correlations among variables

Pearson’s r is in the lower triangle

``` r
pairs(formula = ~ tempave.sc + microclim.sc + tempchange.sc + tempchange_abs.sc + nspp.sc + thermal_bias.sc + human_bowler.sc, data = btaveall, gap = 1/10, cex = 0.2, col = '#00000022', lower.panel = panel.cor)
```

![](assemble_turnover_covariates_files/figure-gfm/pairs-1.png)<!-- -->

# Compare covariates across realms

``` r
par(mfrow=c(5,3))
beanplot(rarefyID_y ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Latitude (degN)', ll = 0.05)
beanplot(tempave ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature (degC)', ll = 0.05)
beanplot(microclim ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Microclimates (degC)', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(tempchange ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature change (degC)', ll = 0.05)
beanplot(Nspp ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Number of species', ll = 0.05, log = 'y')
beanplot(thermal_bias ~ REALM, data = btaveall[!is.na(thermal_bias),], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Thermal bias (degC)', ll = 0.05)
```

![](assemble_turnover_covariates_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->

  - Marine are in generally warmer locations (seawater doesn’t freeze)
  - Marine has a lot of slow, crawling organisms, but land has plants.
    Land also has birds (fast).

# Frequency of time differences

``` r
invisible(btall[, hist(year2 - year1, breaks = seq(0.5, max(year2 - year1)+0.5, by = 1))])
```

![](assemble_turnover_covariates_files/figure-gfm/plot%20time%20differences-1.png)<!-- -->

# Plot dissimilarity vs. time

## Jaccard turnover

![](assemble_turnover_covariates_files/figure-gfm/plot%20diss%20vs%20time-1.png)<!-- -->

# Plot slope vs. explanatory variables

## Each variable individually

Lines are ggplot smoother fits Just Jtu for now, slope within rarefyID
for btaveall
![](assemble_turnover_covariates_files/figure-gfm/plot%20diss%20v%20explanatory%20vars-1.png)<!-- -->

  - Strong trends with temperature change, but trends are pretty
    symmetric around no trend in temperature, which implies warming or
    cooling drives similar degree of community turnover.

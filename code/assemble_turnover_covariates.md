Pairwise dissimilarity covariate data prep and visualization
================

Trims based on data quality, separates the data points to use for
3/5/10/20 year comparisons, and adds covariate data.

``` r
library('mgcv', lib.loc = '/usr/lib64/R/library') # when running on Annotate. Need to load 1.8-26, not 1.8-33.
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
picklongest <- function(year1, year2, y){
  dy <- year2 - year1
  keep <- which.max(dy)
  return(y[keep])
}

btall <- bt # legacy from when we divided dataset into different ts durations (3, 5, 10, 20 years)
btall[, tempave := picklongest(year1, year2, tempave), by = rarefyID]
btall[, tempave_metab := picklongest(year1, year2, tempave_metab), by = rarefyID]
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

# sign of temperature change
signneg11 <- function(x){ # assign 0 a sign of 1 so that there are only 2 levels
  out <- sign(x)
  out[out == 0] <- 1
  return(out)
}
btall[, tsign := signneg11(tempchange)]

#######################
## Transformations

### Adjust response away from 0-1
# transformation for 2 categories. Eq. 1 in Douma & Weedon 2019 MEE
transform01 <- function(x) (x * (length(x) - 1) + 0.5) / (length(x))

btall[, ':='(Jtu.sc = transform01(Jtu), Jne.sc = transform01(Jne), 
           Jbeta.sc = transform01(Jbeta), Horn.sc = transform01(Horn))]


### Log-transform some variables, then center and scale. 
btall[, ':='(tempave.sc = scale(tempave),
           tempave_metab.sc = scale(tempave_metab),
           seas.sc = scale(seas),
           microclim.sc = scale(log(microclim)),
           tempchange.sc = scale(tempchange, center = TRUE), 
           tempchange_abs.sc = scale(abs(tempchange), center = TRUE),
           mass.sc = scale(log(mass_mean_weight)),
           speed.sc = scale(log(speed_mean_weight+1)),
           lifespan.sc = scale(log(lifespan_mean_weight)),
           consumerfrac.sc = scale(consfrac),
           endothermfrac.sc = scale(endofrac),
           nspp.sc = scale(log(Nspp)),
           thermal_bias.sc = scale(thermal_bias),
           npp.sc = scale(log(npp)),
           veg.sc = scale(log(veg+1)),
           duration.sc = scale(duration),
           durationlog.sc = scale(log(duration)))]
btall[, human_bowler.sc := scale(log(human_bowler+1)), by = REALM2] # separate scaling by realm
btall[REALM2 == 'TerrFresh', human_footprint.sc := scale(log(human_venter+1))]
btall[REALM2 == 'Marine', human_footprint.sc := scale(log(human_halpern))]


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
    ## 1: Terrestrial 168
    ## 2:      Marine 130
    ## 3:  Freshwater  21

``` r
btall[, length(unique(rarefyID)), by = REALM]
```

    ##          REALM    V1
    ## 1: Terrestrial  3159
    ## 2:      Marine 38451
    ## 3:  Freshwater   645

## Write out

Only if file doesn’t yet exist

``` r
if(!file.exists(here('output', 'turnover_w_covariates.csv.gz'))){
  write.csv(btall[,.(STUDY_ID, rarefyID, rarefyID_x, rarefyID_y, REALM, Biome, taxa_mod, taxa_mod2, year1, year2, duration, 
                   Jtu.sc, Jne.sc, Jbeta.sc, Horn.sc, 
                   tempave.sc, tempave_metab.sc, seas.sc, 
                   microclim.sc, tempchange, tsign, tempchange.sc, tempchange_abs.sc,
                   mass.sc, speed.sc, lifespan.sc, consumerfrac.sc, endothermfrac.sc, 
                   nspp.sc, thermal_bias.sc, npp.sc, veg.sc,
                   duration.sc, durationlog.sc, human_bowler.sc, REALM2)], 
            gzfile(here('output', 'turnover_w_covariates.csv.gz')), row.names = FALSE)
  write.csv(scalingall, here('output','turnover_w_covariates_scaling.csv'), row.names = FALSE)
}
```

# A bit more prep for visualizing covariate distributions

## Make a data.table summarized by rarefyID

``` r
# function for returning the slope of a regression and catching errors related to NAs
lmNAcoef <- function(y, x){
  if(sum(!is.na(y)) > 1){ # if enough data points to calculate a slope
    b <- coef(lm(y ~ x))[2]
    return(b)
  } else {
    return(NA_real_)
  }
}
btaveall <- btall[, .(Jtuslope = lmNAcoef(Jtu, duration), Jbetaslope = lmNAcoef(Jbeta, duration), Hornslope = lmNAcoef(Horn, duration), 
                    minyrBT = min(year1), maxyrBT = max(year2), nyrBT = length(Jtu) + 1,
                    REALM = unique(REALM), REALM2 = unique(REALM2), rarefyID_y = mean(rarefyID_y),
                    tempchange = mean(tempchange), tempave = mean(tempave),
                    tempave_metab = mean(tempave_metab), seas = mean(seas), microclim = mean(microclim),
                    mass_mean_weight = mean(mass_mean_weight), speed_mean_weight = mean(speed_mean_weight),
                    lifespan_mean_weight = mean(lifespan_mean_weight), consfrac = mean(consfrac),
                    endofrac = mean(endofrac), Nspp = mean(Nspp), thermal_bias = mean(thermal_bias),
                    npp = mean(npp), veg = mean(veg), human_bowler = mean(human_bowler),
                    human_venter = mean(human_venter), human_halpern = mean(human_halpern),
                    tempave.sc = mean(tempave.sc), tempave_metab.sc = mean(tempave_metab.sc), 
                    tempchange.sc = mean(tempchange.sc), tempchange_abs.sc = mean(tempchange_abs.sc),
                    seas.sc = mean(seas.sc), microclim.sc = mean(microclim.sc),
                    mass.sc = mean(mass.sc), speed.sc = mean(speed.sc), lifespan.sc = mean(lifespan.sc),
                    consumerfrac.sc = mean(consumerfrac.sc), endothermfrac.sc = mean(endothermfrac.sc),
                    nspp.sc = mean(nspp.sc), thermal_bias.sc = mean(thermal_bias.sc),
                    npp.sc = mean(npp.sc), veg.sc = mean(veg.sc), human_bowler.sc = mean(human_bowler.sc),
                    human_footprint.sc = mean(human_footprint.sc)),
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
invisible(btaveall[, hist(mass_mean_weight, main = 'Mass (g)', cex.main = cexmain)])
invisible(btaveall[, hist(speed_mean_weight, main = 'Speed (km/hr)', cex.main = cexmain)])
invisible(btaveall[, hist(lifespan_mean_weight, main = 'Lifespan (yr)', cex.main = cexmain)])
invisible(btaveall[, hist(tempave_metab, main = 'Metabolic temperature (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(consfrac, main = 'Consumers (fraction)', cex.main = cexmain)])
invisible(btaveall[, hist(endofrac, main = 'Endotherms (fraction)', cex.main = cexmain)])
invisible(btaveall[, hist(tempave, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(tempchange, main = 'Temperature trend (°C/yr)', cex.main = cexmain)]) # all the raw data
invisible(btaveall[, hist(seas, main = 'Seasonality (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(microclim, main = 'Microclimates (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(Nspp, main = 'Species richness', cex.main = cexmain)])
invisible(btaveall[, hist(thermal_bias, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(npp, main = 'Net primary productivity', cex.main = cexmain)])
invisible(btaveall[, hist(veg, main = 'Vegetation index', cex.main = cexmain)])
invisible(btaveall[, hist(human_bowler, main = 'Human impact score (Bowler)', cex.main = cexmain)])
invisible(btaveall[, hist(human_venter, main = 'Human impact score (Venter)', cex.main = cexmain)])
invisible(btaveall[, hist(human_halpern, main = 'Human impact score (Halpern)', cex.main = cexmain)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20unscaled-1.png)<!-- -->

## Scaled covariates

``` r
# histograms to examine
cexmain = 0.6
par(mfrow = c(5,4))
invisible(btaveall[, hist(tempave.sc, main = 'Environmental temperature (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(tempave_metab.sc, main = 'Metabolic temperature (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(seas.sc, main = 'Seasonality (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(microclim.sc, main = 'log Microclimates (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(tempchange.sc, main = 'Temperature trend (°C/yr)', cex.main = cexmain)])
invisible(btaveall[, hist(tempchange_abs.sc, main = 'abs(Temperature trend) (°C/yr)', cex.main = cexmain)])
invisible(btaveall[, hist(mass.sc, main = 'log Mass (g)', cex.main = cexmain)])
invisible(btaveall[, hist(speed.sc, main = 'log Speed (km/hr)', cex.main = cexmain)])
invisible(btaveall[, hist(lifespan.sc, main = 'log Lifespan (yr)', cex.main = cexmain)])
invisible(btaveall[, hist(consumerfrac.sc, main = 'Consumers (fraction)', cex.main = cexmain)])
invisible(btaveall[, hist(endothermfrac.sc, main = 'Endotherms (fraction)', cex.main = cexmain)])
invisible(btaveall[, hist(nspp.sc, main = 'log Species richness', cex.main = cexmain)])
invisible(btaveall[, hist(thermal_bias.sc, main = 'Thermal bias (°C)', cex.main = cexmain)])
invisible(btaveall[, hist(npp.sc, main = 'log Net primary productivity', cex.main = cexmain)])
invisible(btaveall[, hist(veg.sc, main = 'log Vegetation index', cex.main = cexmain)])
invisible(btaveall[, hist(human_bowler.sc, main = 'log Human impact score (Bowler)', cex.main = cexmain)])
invisible(btaveall[, hist(human_footprint.sc, main = 'log Human impact score (Venter & Halpern)', cex.main = cexmain)])
```

![](assemble_turnover_covariates_files/figure-gfm/histograms%20scaled-1.png)<!-- -->

# Check correlations among variables

Pearson’s r is in the lower triangle

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
pairs(formula = ~ tempave.sc + tempave_metab.sc + seas.sc + microclim.sc + tempchange.sc + tempchange_abs.sc + mass.sc + speed.sc + lifespan.sc + consumerfrac.sc + endothermfrac.sc + nspp.sc + thermal_bias.sc + npp.sc + veg.sc + human_bowler.sc + human_footprint.sc, data = btaveall, gap = 1/10, cex = 0.2, col = '#00000022', lower.panel = panel.cor)
```

![](assemble_turnover_covariates_files/figure-gfm/pairs-1.png)<!-- -->

  - Mass and lifespan look tightly correlated, but r only 0.56…?
  - Tempave\_metab and lifespan don’t look tightly correlated, but r=
    -0.81
  - Tempave\_metab and speed don’t look tightly correlated, but r= -0.83
  - Lifespan and speed don’t look tightly correlated, but r = 0.73

# Compare covariates across realms

``` r
par(mfrow=c(5,3))
beanplot(rarefyID_y ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Latitude (degN)', ll = 0.05)
beanplot(tempave ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature (degC)', ll = 0.05)
beanplot(tempave_metab ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Metabolic Temperature (degC)', ll = 0.05, bw = 'nrd0') # nrd0 bandwidth to calculation gap
beanplot(seas ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Seasonality (degC)', ll = 0.05)
beanplot(microclim ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Microclimates (degC)', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(tempchange ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Temperature change (degC)', ll = 0.05)
beanplot(mass_mean_weight ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Mass (g)', ll = 0.05, log = 'y')
beanplot(speed_mean_weight +1 ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Speed (km/hr)', ll = 0.05, log = 'y')
beanplot(lifespan_mean_weight ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Lifespan (yr)', ll = 0.05, log = 'y')
#beanplot(consfrac ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Consumers (fraction)', ll = 0.05, log = '') # too sparse
#beanplot(endofrac ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Endotherms (fraction)', ll = 0.05, log = '') # too sparse
beanplot(Nspp ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Number of species', ll = 0.05, log = 'y')
beanplot(thermal_bias ~ REALM, data = btaveall[!is.na(thermal_bias),], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'Thermal bias (degC)', ll = 0.05)
beanplot(npp ~ REALM, data = btaveall, what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'NPP', ll = 0.05)
```

    ## log="y" selected

``` r
beanplot(veg ~ REALM, data = btaveall[REALM !='Marine',], what = c(1,1,1,1), col = c("#CAB2D6", "#33A02C", "#B2DF8A"), border = "#CAB2D6", ylab = 'veg', ll = 0.05)
```

![](assemble_turnover_covariates_files/figure-gfm/compare%20across%20realms-1.png)<!-- -->

  - Marine are in generally warmer locations (seawater doesn’t freeze)
  - Marine have much lower seasonality.
  - Marine and freshwater have some very small masses (plankton), but
    much of dataset is similar to terrestrial.
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
  - Some indication of less turnover for larger organisms (mass)
  - Higher turnover on land with higher seasonality?
  - More turnover for shorter-lived organisms?
  - No really clear differences among realms.

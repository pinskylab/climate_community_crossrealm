Calculate biodiversity urnover
================

# Set up

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

# richness
rich <- fread('output/richness_by_rarefyID.csv.gz') # number of species
```

Remove studies with only 1 species and trim to count data

``` r
length(unique(bt$rarefyID))
```

    ## [1] 53467

``` r
bt <- bt[rarefyID %in% rich[Nspp > 1, rarefyID] & BROAD_TYPE == 'count', ]
length(unique(bt$rarefyID))
```

    ## [1] 51578

# Plot turnover example

Michaelis-Menton or Exponential fits (2- or 3-parameter).

``` r
test <- bt[rarefyID == '339_1085477',]
test[, YEAR0 := YEAR - min(YEAR)]

# fit a 2- or 3-parameter Michaelis-Menton curve (latter includes an intercept)
modmm2 <- nls(Jtu_base ~ Vm * YEAR0/(K + YEAR0), data = test[YEAR0 != 0, ], 
              start = list(K = mean(test$YEAR0), Vm = 1),
              algorithm = 'port', # to allow bounds to be set
              lower = c(0, 0),
              upper = c(max(test$YEAR0), 1))
modmm3 <- nls(Jtu_base ~ c + Vm * YEAR0/(K + YEAR0), 
              data = test[YEAR0 != 0, ], 
              start = list(c = 0, K = mean(test$YEAR0), Vm = 1),
              algorithm = 'port',
              lower = c(0, 0, 0),
              upper = c(1, max(test$YEAR0), 1))

# fit a 2- or 3-parameter exponential decay
modexp2 <- nls(Jtu_base ~ Vm - exp(a * YEAR0)*Vm, data = test[YEAR0 != 0, ], 
              start = list(a = -0.5, Vm = 1),
              algorithm = 'port',
              lower = c(-Inf, 0),
              upper = c(0, 1))
modexp3 <- nls(Jtu_base ~ Vm - exp(a * (YEAR0 + o))*Vm, data = test[YEAR0 != 0, ], 
              start = list(a = -0.5, Vm = 1, o = 0),
              algorithm = 'port',
              lower = c(-Inf, 0, -Inf),
              upper = c(Inf, 1, Inf))

# make predictions
mml <- data.frame(YEAR0 = seq(0, diff(range(test$YEAR)), length.out = 100))
mml$YEAR <- mml$YEAR0 + min(test$YEAR)
mml$Jtu_base.mm2 <- predict(modmm2, newdata = mml)
mml$Jtu_base.mm3 <- predict(modmm3, newdata = mml)
mml$Jtu_base.exp2 <- predict(modexp2, newdata = mml)
mml$Jtu_base.exp3 <- predict(modexp3, newdata = mml)

# calc half-saturations
cmm2 <- coef(modmm2)
cmm3 <- coef(modmm3)
cexp2 <- coef(modexp2)
cexp3 <- coef(modexp3)

hsmm2 <- cmm2[['K']] + min(test$YEAR)
hsmm3 <- with(as.data.frame(t(cmm3)), K*(Vm-c)/(Vm+c)) + min(test$YEAR)
hsexp2 <- with(as.data.frame(t(cexp2)), log(1/2)/a) + min(test$YEAR)
hsexp3 <- with(as.data.frame(t(cexp3)), log(exp(a*o)/2)/a - o) + min(test$YEAR)

# get colors to match ggplot
hues = seq(15, 375, length = 5 + 1)
cols <- hcl(h = hues, l = 65, c = 100)[1:5]

# make plot
ggplot(test, aes(YEAR, Jtu_base, color = 'Linear')) +
    geom_point() +
    geom_smooth(data = test[YEAR != min(YEAR), ], 
                method = 'lm') + # a straight line fit
    geom_line(data = mml, aes(YEAR, Jtu_base.mm2, color = 'Michaelis-Menton 2-param')) + # MM
    geom_line(data = mml, aes(YEAR, Jtu_base.mm3, color = 'Michaelis-Menton 3-param'), linetype = 'dashed', size = 1) +
    geom_line(data = mml, aes(YEAR, Jtu_base.exp2, color = 'Exponential 2-param')) +
    geom_line(data = mml, aes(YEAR, Jtu_base.exp3, color = 'Exponential 3-param'), linetype = 'dashed', size = 1) +
    geom_vline(xintercept = c(hsmm2, hsmm3, hsexp2, hsexp3), color = cols[c(4,5,1,2)], size = c(0.5, 1, 0.5, 1), 
               linetype = c('solid', 'dashed', 'solid', 'dashed')) +
    labs(x = 'Year', y = 'Jaccard dissimilarity', color = 'Functional form')
```

![](calc_turnover_files/figure-gfm/plot%20turnover-1.png)<!-- -->

Vertical lines are the half-saturations.

# Calculate temporal trends (temporal turnover)

This only recalculates from scratch if temp/trendstemp.rds is not
available. \#\# Set up functions

``` r
# function to calc linear trend, removing first year with 0
calctrendrem0 <- function(y, YEAR, nm = 'y'){
  if(length(YEAR)>2){
     # remove year 1 comparison
    o <- order(YEAR)
    YEAR2 <- YEAR[o][2:length(YEAR)]
    y2 <- y[o][2:length(y)]
    
    if(sum(!is.na(y2)) >= 2){ # make sure enough values to fit a line
      mod <- lm(y2 ~ YEAR2) # fit line
      out <- list(y = coef(mod)[2], # coef for the slope
                  y_se = sqrt(diag(vcov(mod)))[2]) # SE
      names(out) <- c(nm, paste0(nm, '_se'))
      return(out)
      
    } else {
      out <- list(y = NA_real_, y_se = NA_real_)
      names(out) <- c(nm, paste0(nm, '_se'))
      return(out)
    }
    
  } else {
    out <- list(y = NA_real_, y_se = NA_real_)
    names(out) <- c(nm, paste0(nm, '_se'))
    return(out)
  }
  

}

# function to get distance from last to first year
calcfirstlast <- function(y, YEAR, nm = 'y'){
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
  return(out)
}

# calc half-saturation constant from exponential fit
# if pred = FALSE, return half-saturation
# if pred = TRUE, return prediction data.frame
calcexp <- function(y, YEAR, nm = 'y', pred = FALSE){
    if(length(YEAR)>2){
        o <- order(YEAR)
        YEAR2 <- YEAR[o][2:length(YEAR)] - YEAR[1]
        y2 <- y[o][2:length(y)]
        
        if(sum(!is.na(y2)) >= 1){ # make sure enough values to fit
            
            out <- tryCatch({
                
                modexp <- nls(y2 ~ Vm - Vm*exp(a * YEAR2), 
                              start = list(a = -0.5, Vm = 1),
                              algorithm = 'port',
                              lower = c(-Inf, 0),
                              upper = c(0, 1))
                
                if(!pred){
                  # calc half-saturations
                  cexp <- coef(modexp)
                  hsexp <- with(as.data.frame(t(cexp)), log(1/2)/a)
                  
                  out <- list(y = hsexp) # half-saturation
                  names(out) <- nm
                }
                if(pred){
                  # make predictions
                  out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
                  out$YEAR <- out$YEAR2 + YEAR[1]
                  out$pred <- predict(modexp, newdata = out)
                }
                out
                
            }, warning = function(w) {
              #print('warning')
              if(!pred){
                out <- list(y = NA_real_)
                names(out) <- nm
                return(out)
                
              }
              if(pred){
                out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
                out$YEAR <- out$YEAR2 + YEAR[1]
                out$pred <- NA
                return(out)
              }
            }, error = function(e) {
              #print('error')
              if(!pred){
                out <- list(y = NA_real_)
                names(out) <- nm
                return(out)
                
              }
              if(pred){
                out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
                out$YEAR <- out$YEAR2 + YEAR[1]
                out$pred <- NA
                return(out)
              }
            })
            return(out)
            
        } else {
          if(!pred){
            out <- list(y = NA_real_)
            names(out) <- nm
            return(out)
            
          }
          if(pred){
            out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
            out$YEAR <- out$YEAR2 + YEAR[1]
            out$pred <- NA
            return(out)
          }
        }
        
    } else {
      if(!pred){
        out <- list(y = NA_real_)
        names(out) <- nm
        return(out)
        
      }
      if(pred){
        out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
        out$YEAR <- out$YEAR2 + YEAR[1]
        out$pred <- NA
        return(out)
      }
    }
}

# calc half-saturation constant from Michaelis-Menton fits
calcmm <- function(y, YEAR, nm = 'y', pred = FALSE){
    if(length(YEAR)>2){
        o <- order(YEAR)
        YEAR2 <- YEAR[o][2:length(YEAR)] - YEAR[1]
        y2 <- y[o][2:length(y)]
        
        if(sum(!is.na(y2)) >= 1){ # make sure enough values to fit
          #print('fit model')
          out <- tryCatch({
            
            modmm <- nls(y2 ~ Vm * YEAR2/(K + YEAR2), 
                         start = list(K = mean(YEAR2), Vm = 1),
                         algorithm = 'port', 
                         lower = c(0, 0),
                         upper = c(1000, 1))
            
            if(!pred){
              #print('calc half-saturations')
              cmm <- coef(modmm)
              hsmm <- cmm[['K']]
              
              out <- list(y = hsmm)
              names(out) <- nm
            }
            if(pred){
              #print('make predictions')
              out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
              out$YEAR <- out$YEAR2 + YEAR[1]
              out$pred <- predict(modmm, newdata = out)
            }
            out
          }, warning = function(w) {
            #print('warning')
            if(!pred){
              out <- list(y = NA_real_)
              names(out) <- nm
              return(out)
              
            }
            if(pred){
              out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
              out$YEAR <- out$YEAR2 + YEAR[1]
              out$pred <- NA
              return(out)
            }
            
          }, error = function(e) {
            #print('error')
            if(!pred){
              out <- list(y = NA_real_)
              names(out) <- nm
              return(out)
            }
            if(pred){
              out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
              out$YEAR <- out$YEAR2 + YEAR[1]
              out$pred <- NA
              return(out)
            }
            
          })
          return(out)
          
        } else {
          if(!pred){
            out <- list(y = NA_real_)
            names(out) <- nm
            return(out)
          }
          if(pred){
            out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
            out$YEAR <- out$YEAR2 + YEAR[1]
            out$pred <- NA
            return(out)
          }
        }
        
    } else {
      if(!pred){
        out <- list(y = NA_real_)
        names(out) <- nm
        return(out)
      }
      if(pred){
        out <- data.frame(YEAR2 = seq(0, max(YEAR2), length.out = 100))
        out$YEAR <- out$YEAR2 + YEAR[1]
        out$pred <- NA
        return(out)
      }
    }
}
```

## Do calculations

And write out trendstemp.rds

``` r
setkey(bt, STUDY_ID, rarefyID, YEAR)

if(file.exists('temp/trendstemp.rds')){
    print('File already exists. Will not do calculations')
    trends <- readRDS('temp/trendstemp.rds')
} else {
  print('Calculating from scratch')
  
  # remove the next file (with standardized slopes) as well if we are recalculting this from scratch
  if(file.exists('temp/trendstemp2.rds')){
    file.remove('temp/trendstemp2.rds')
  }
  
  trends1 <- bt[, calctrendrem0(Jtu_base, YEAR, 'Jtutrendrem0'), 
                by = .(REALM, Biome, taxa_mod, STUDY_ID, 
                       rarefyID, rarefyID_x, rarefyID_y)] # calculate trend in Jaccard turnover without first year
  trends2 <- bt[, calctrendrem0(Jbeta_base, YEAR, 'Jbetatrendrem0'),
                by = .(rarefyID)]
  trends3 <- bt[, calctrendrem0(1-Horn_base, YEAR, 'Horntrendrem0'),
                by = .(rarefyID)]
  
  trends4 <- bt[, calcfirstlast(Jtu_base, YEAR, 'Jtulast'),
                by = .(rarefyID)]
  trends5 <- bt[, calcfirstlast(Jbeta_base, YEAR, 'Jbetalast'),
                by = .(rarefyID)]
  trends6 <- bt[, calcfirstlast(1-Horn_base, YEAR, 'Hornlast'),
                by = .(rarefyID)]
  
  trends7 <- bt[, calcexp(Jtu_base, YEAR, 'Jtuexp'),
                by = .(rarefyID)]
  trends8 <- bt[, calcexp(Jbeta_base, YEAR, 'Jbetaexp'),
                by = .(rarefyID)]
  trends9 <- bt[, calcexp(1-Horn_base, YEAR, 'Hornexp'),
                by = .(rarefyID)]
  
  trends10 <- bt[, calcmm(Jtu_base, YEAR, 'Jtumm'),
                by = .(rarefyID)]
  trends11 <- bt[, calcmm(Jbeta_base, YEAR, 'Jbetamm'),
                by = .(rarefyID)]
  trends12 <- bt[, calcmm(1-Horn_base, YEAR, 'Hornmm'),
                by = .(rarefyID)]
  
  nyrBT <-  bt[, .(nyrBT = length(YEAR),
                   minyrBT = min(YEAR),
                   maxyrBT = max(YEAR),
                   medianyrBT = median(YEAR),
                   meanyrBT = mean(YEAR)),
               by = .(rarefyID)] # number of years in time-series
  
  trends <- merge(trends1, trends2)
  trends <- merge(trends, trends3)
  trends <- merge(trends, trends4)
  trends <- merge(trends, trends5)
  trends <- merge(trends, trends6)
  trends <- merge(trends, trends7)
  trends <- merge(trends, trends8)
  trends <- merge(trends, trends9)
  trends <- merge(trends, trends10)
  trends <- merge(trends, trends11)
  trends <- merge(trends, trends12)
  trends <- merge(trends, nyrBT)
  
  saveRDS(trends, file = 'temp/trendstemp.rds')
  
}
```

    ## [1] "File already exists. Will not do calculations"

``` r
nrow(trends)
```

    ## [1] 51578

Add species richness

``` r
trends <- merge(trends, rich, all.x = TRUE) # species richness
```

## Plot every Jtu timeseries (a lot\!)

Not run during knitting

``` r
rids <- trends[nyrBT > 2, sort(unique(rarefyID))]
setkey(bt, rarefyID, YEAR)
print(paste(length(rids), ' rarefyIDs'))
filenum <- 1
plotnum <- 1
for(i in 1:length(rids)){
#for(i in 1:400){ # for testing
  jtutrendrem0 <- trends[rarefyID == rids[i], Jtutrendrem0]
  jtuexp <- trends[rarefyID == rids[i], Jtuexp]
  jtumm <- trends[rarefyID == rids[i], Jtumm]
  nspp <- trends[rarefyID == rids[i], Nspp]
  x <- bt[rarefyID == rids[i], YEAR]
  y <- bt[rarefyID == rids[i], Jtu_base]
  
  if(length(x) > 2 & nspp > 1){
    if(plotnum %% 400 == 1){
      if(plotnum >1) dev.off()
      png(file = paste0('figures/jtu_plots/jtu_plots', formatC(filenum, width = 3, format = 'd', flag = '0'), '.png'), 
          width = 36, height = 36, units = 'in', res = 100)
      par(mfrow=c(20,20), mai = c(0.4, 0.5, 0.5, 0.1))
      filenum <- filenum + 1
    }
    plot(x, y, main = paste('Jtu rem0:', signif(jtutrendrem0, 3), 'Jtu exp:', signif(jtuexp, 3), '\nJtu mm:', signif(jtumm, 3), 'Nspp:', nspp, '\nrID:', rids[i]), xlab = '', ylab = 'Jtu',
         cex.main = 0.7)
    abline(lm(y ~ x))
    if(plotnum %% 400 == 1){
      legend('topleft', legend = c('linear', 'rem0', 'exp', 'mm'), lty = 1, 
             col = c('black', 'red', 'blue', 'green'), cex = 0.5)
    }
    
    x2 <- x[2:length(x)]
    y2 <- y[2:length(y)]
    abline(lm(y2 ~ x2), col = 'red')
    
    x3 <- calcexp(y, x, pred = TRUE)
    lines(x3$YEAR, x3$pred, col = 'blue')
    
    x4 <- calcmm(y, x, pred = TRUE)
    lines(x4$YEAR, x4$pred, col = 'green')
    
    plotnum <- plotnum + 1
  }
}
  
dev.off()
```

## Standardize slope values against null model simulations

Set up null model sims

``` r
# remove some duplicates
sim_summary <- as.data.table(sim_summary)
sim_summary <- sim_summary[!duplicated(sim_summary), ]

# hist
hist(sim_summary[metric == 'Jaccard (total)', sd], breaks = 200, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds-1.png)<!-- -->

``` r
hist(sim_summary[metric == 'Nestedness', sd], breaks = 200, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds-2.png)<!-- -->

``` r
hist(sim_summary[metric == 'Morisita-Horn', sd], breaks = 200, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds-3.png)<!-- -->

``` r
# plot
sim_summary[metric == 'Jaccard (total)', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](calc_turnover_files/figure-gfm/plot%20sds-4.png)<!-- -->

``` r
sim_summary[metric == 'Nestedness', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](calc_turnover_files/figure-gfm/plot%20sds-5.png)<!-- -->

``` r
sim_summary[metric == 'Morisita-Horn', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](calc_turnover_files/figure-gfm/plot%20sds-6.png)<!-- -->

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

![](calc_turnover_files/figure-gfm/plot%20sds-7.png)<!-- -->

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

![](calc_turnover_files/figure-gfm/plot%20sds-8.png)<!-- -->

Merge null model

``` r
# merge
trends <- merge(trends, sim_summary[metric == 'Jaccard (total)', .(Jbeta_exp = expected, Jbeta_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
trends <- merge(trends, sim_summary[metric == 'Nestedness', .(Jtu_exp = expected, Jtu_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
trends <- merge(trends, sim_summary[metric == 'Morisita-Horn', .(Horn_exp = expected, Horn_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
```

Set a threshold and only standardize those timeseries with SD \> 0.001
Write out/load trendstemp2.rds

``` r
if(file.exists('temp/trendstemp2.rds')){
    print('File already exists. Will not write out')
  trends <- readRDS('temp/trendstemp2.rds')
  trends <- merge(trends, rich, all.x = TRUE) # species richness
} else {
# standardize, only where sd >0
  trends[, Jtutrendz := NA_real_]
  trends[, Jbetatrendz := NA_real_]
  trends[, Horntrendz := NA_real_]
  
  trends[Jtu_sd > 0.001, Jtutrendz := (Jtutrendrem0 - Jtu_exp)/Jtu_sd]
  trends[Jbeta_sd > 0.001, Jbetatrendz := (Jbetatrendrem0 - Jbeta_exp)/Jbeta_sd]
  trends[Horn_sd > 0.001, Horntrendz := (Horntrendrem0 - Horn_exp)/Horn_sd]
  
  
  print('writing out')
  saveRDS(trends[, .(rarefyID, REALM, Biome, taxa_mod, STUDY_ID, rarefyID_x, rarefyID_y,
                     Jtutrendrem0, Jtutrendrem0_se, Jbetatrendrem0, Jbetatrendrem0_se,
                     Horntrendrem0, Horntrendrem0_se, 
                     Jtutrendz, Jbetatrendz, Horntrendz,
                     Jtulast, Jbetalast, Hornlast,
                     Jtuexp, Jbetaexp, Hornexp,
                     Jtumm, Jbetamm, Hornmm,
                     nyrBT, minyrBT, maxyrBT, medianyrBT, meanyrBT)], file = 'temp/trendstemp2.rds')
}
```

    ## [1] "File already exists. Will not write out"

``` r
nrow(trends)
```

    ## [1] 51578

# Examine the turnover calculations

## How many values?

``` r
apply(trends[, .(Jtutrendrem0, Jbetatrendrem0, Horntrendrem0,
                 Jtutrendz, Jbetatrendz, Horntrendz,
                Jtulast, Jbetalast, Hornlast, 
                Jtuexp, Jbetaexp, Hornexp, Jtumm, Jbetamm, Hornmm)], MARGIN = 2, 
      function(x) sum(!is.na(x)))
```

    ##   Jtutrendrem0 Jbetatrendrem0  Horntrendrem0      Jtutrendz    Jbetatrendz 
    ##          38188          38188          38188          33381          34684 
    ##     Horntrendz        Jtulast      Jbetalast       Hornlast         Jtuexp 
    ##          34584          51578          51578          51578          26550 
    ##       Jbetaexp        Hornexp          Jtumm        Jbetamm         Hornmm 
    ##          25924          24769          38063          38096          38053

## Do some basic checks of the turnover calculations

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
    ##         Jtutrendz Jbetatrendz Horntrendz   Jtulast Jbetalast  Hornlast
    ##     1:         NA          NA         NA 0.4324324 0.5000000 0.3642417
    ##     2: 0.01012774    1.543874         NA 0.0000000 0.5000000 0.1541960
    ##     3:         NA          NA         NA 1.0000000 1.0000000 1.0000000
    ##     4:         NA          NA         NA 0.8000000 0.8571429 0.9177154
    ##     5:         NA          NA         NA 0.5000000 0.7142857 0.7026207
    ##    ---                                                                
    ## 51574:         NA          NA         NA 0.0000000 0.9285714 0.9456966
    ## 51575:         NA          NA         NA 1.0000000 1.0000000 1.0000000
    ## 51576:         NA          NA         NA 0.6666667 0.7500000 0.3265231
    ## 51577:         NA          NA         NA 0.4000000 0.8571429 0.9516536
    ## 51578:         NA          NA         NA 0.0000000 0.8750000 0.8025743
    ##            Jtuexp  Jbetaexp Hornexp      Jtumm   Jbetamm   Hornmm nyrBT
    ##     1: 0.02709627 0.2524249 1.22328  0.1253377 0.1942427 1.043563    31
    ##     2: 1.38629436 0.2599677      NA 15.5000000 0.1336071 0.000000    31
    ##     3:         NA        NA      NA         NA        NA       NA     2
    ##     4:         NA        NA      NA         NA        NA       NA     2
    ##     5:         NA        NA      NA         NA        NA       NA     2
    ##    ---                                                                 
    ## 51574:         NA        NA      NA  0.0000000 0.0000000 0.000000     3
    ## 51575:         NA        NA      NA         NA        NA       NA     2
    ## 51576:         NA        NA      NA         NA        NA       NA     2
    ## 51577:         NA        NA      NA         NA        NA       NA     2
    ## 51578:         NA        NA      NA         NA        NA       NA     2
    ##        minyrBT maxyrBT medianyrBT meanyrBT Nspp
    ##     1:    1981    2011     1996.0 1996.000   83
    ##     2:    1981    2011     1996.0 1996.000   15
    ##     3:    1985    1999     1992.0 1992.000    4
    ##     4:    1985    1993     1989.0 1989.000    7
    ##     5:    1985    1993     1989.0 1989.000    7
    ##    ---                                         
    ## 51574:    1934    1952     1950.0 1945.333   42
    ## 51575:    1934    1955     1944.5 1944.500   37
    ## 51576:    1973    1975     1974.0 1974.000    9
    ## 51577:    1973    1975     1974.0 1974.000   23
    ## 51578:    1971    1973     1972.0 1972.000   17

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

``` r
summary(trends$Jtuexp)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   0.002   0.394   1.386   5.011   3.889 648.434   25028

``` r
summary(trends$Jbetaexp)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   0.002   0.269   0.674   1.989   1.935 363.489   25654

``` r
summary(trends$Hornexp)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    ##    0.002    0.355    1.346    6.712    4.626 1683.947    26809

``` r
summary(trends$Jtumm)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   0.000   0.000   0.323   5.195   4.000 924.180   13515

``` r
summary(trends$Jbetamm)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   0.000   0.000   0.113   1.547   0.996 511.455   13482

``` r
summary(trends$Hornmm)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    ##    0.000    0.000    0.106    9.960    2.839 1000.000    13525

## Histograms of temporal change

Standardized slopes have very large and small values

``` r
x <- trends[, hist(Jtutrendrem0)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-1.png)<!-- -->

``` r
x <- trends[, hist(Jbetatrendrem0)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-2.png)<!-- -->

``` r
x <- trends[, hist(Horntrendrem0)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-3.png)<!-- -->

``` r
x <- trends[, hist(Jtutrendz, breaks = 100)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-4.png)<!-- -->

``` r
x <- trends[, hist(Jbetatrendz, breaks = 100)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-5.png)<!-- -->

``` r
x <- trends[, hist(Horntrendz, breaks = 100)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-6.png)<!-- -->

``` r
x <- trends[, hist(Jtulast)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-7.png)<!-- -->

``` r
x <- trends[, hist(Jbetalast)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-8.png)<!-- -->

``` r
x <- trends[, hist(Hornlast)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-9.png)<!-- -->

``` r
x <- trends[, hist(log10(Jtuexp), breaks = 20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-10.png)<!-- -->

``` r
x <- trends[, hist(log10(Jbetaexp), breaks = 20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-11.png)<!-- -->

``` r
x <- trends[, hist(log10(Hornexp), breaks = 20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-12.png)<!-- -->

``` r
x <- trends[, hist(log10(Jtumm), breaks = 20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-13.png)<!-- -->

``` r
x <- trends[, hist(log10(Jbetamm), breaks = 20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-14.png)<!-- -->

``` r
x <- trends[, hist(log10(Hornmm), breaks = 20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-15.png)<!-- -->

## Turnover calculations are correlated, though less so for Horn

``` r
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = 'pairwise.complete.obs')
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = 0.5) #, cex = cex.cor * r)
}
pairs(formula = ~ Jtutrendrem0 + Jbetatrendrem0 + Horntrendrem0 + 
        Jtutrendz + Jbetatrendz + Horntrendz + 
          Jtulast + Jbetalast + Hornlast + 
        log10(Jtuexp) + log10(Jbetaexp) + log10(Hornexp) +
        log10(Jtumm) + log10(Jbetamm) + log10(Hornmm), 
      data = trends, gap = 1/10, cex = 0.2, col = '#00000022', 
      lower.panel = panel.cor,
      upper.panel = panel.smooth)
```

![](calc_turnover_files/figure-gfm/pairs-1.png)<!-- --> Standardization
has the effect of up-weighting some of the smaller trends, as expected
Strange that Jbetatrendz and Jtutrendz are not really correlated \*last
correlated to \*trendrem0. Also constrained to a parallelogram Jtulast
is lower diagonal of Jbetalast, which makes sense Hornlast correlated
but with lots of scatter to Jbetalast

# Change compared to number of years and number of species in time-series

## Number of years

``` r
# number of year
trends[, summary(nyrBT)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   2.000   2.000   4.000   5.538   7.000  97.000

``` r
x <- trends[, hist(nyrBT)]
```

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20years-1.png)<!-- -->

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

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20years-2.png)<!-- -->

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

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20years-3.png)<!-- -->

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

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20years-4.png)<!-- -->

``` r
trends[, plot(nyrBT, Jtuexp, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Jbetaexp, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Hornexp, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20years-5.png)<!-- -->

``` r
trends[, plot(nyrBT, Jtumm, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Jbetamm, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(nyrBT, Hornmm, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20years-6.png)<!-- -->

## Duration

``` r
trends[, duration := maxyrBT - minyrBT + 1]
trends[, summary(duration)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    3.00    6.00   11.00   13.52   17.00   98.00

``` r
x <- trends[, hist(duration)]
```

![](calc_turnover_files/figure-gfm/change%20vs.%20duration-1.png)<!-- -->

``` r
par(mfrow=c(1,3))
trends[, plot(duration, Jtutrendrem0, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Jbetatrendrem0, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Horntrendrem0, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20duration-2.png)<!-- -->

``` r
trends[, plot(duration, Jtutrendz, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Jbetatrendz, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Horntrendz, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20duration-3.png)<!-- -->

``` r
trends[, plot(duration, Jtulast, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Jbetalast, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Hornlast, log = 'x', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20duration-4.png)<!-- -->

``` r
trends[, plot(duration, Jtuexp, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Jbetaexp, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Hornexp, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20duration-5.png)<!-- -->

``` r
trends[, plot(duration, Jtumm+1, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Jbetamm+1, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

``` r
trends[, plot(duration, Hornmm+1, log = 'xy', col = '#00000033')]; abline(h = 0)
```

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20duration-6.png)<!-- -->

## Number of species

``` r
trends[, summary(nyrBT)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   2.000   2.000   4.000   5.538   7.000  97.000

``` r
x <- trends[, hist(nyrBT)]
```

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20spp-1.png)<!-- -->

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

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20spp-2.png)<!-- -->

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

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20spp-3.png)<!-- -->

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

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20spp-4.png)<!-- -->

    ## NULL

``` r
trends[, plot(Nspp, Jtuexp, log = 'xy', col = '#00000033')]
```

    ## NULL

``` r
trends[, plot(Nspp, Jbetaexp, log = 'xy', col = '#00000033')]
```

    ## NULL

``` r
trends[, plot(Nspp, Hornexp, log = 'xy', col = '#00000033')]
```

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20spp-5.png)<!-- -->

    ## NULL

``` r
trends[, plot(Nspp, Jtumm, log = 'xy', col = '#00000033')]
```

    ## NULL

``` r
trends[, plot(Nspp, Jbetamm, log = 'xy', col = '#00000033')]
```

    ## NULL

``` r
trends[, plot(Nspp, Hornmm, log = 'xy', col = '#00000033')]
```

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20spp-6.png)<!-- -->

    ## NULL

## Average change compared to \#years

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

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20years-1.png)<!-- -->

``` r
ggplot(trends, aes(nyrBT, Jtutrendz, color = 'Jtu z')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetatrendz, color = 'Jbeta z')) +
    geom_smooth(aes(y = Horntrendz, color = 'Horn z')) +
    scale_x_log10() +
    labs(y = 'Standardize slope') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20years-2.png)<!-- -->

``` r
ggplot(trends, aes(nyrBT, Jtulast, color = 'Jtu last')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetalast, color = 'Jbeta last')) +
    geom_smooth(aes(y = Hornlast, color = 'Horn last')) +
    scale_x_log10() +
    labs(y = 'Dissimilarity') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20years-3.png)<!-- -->

``` r
ggplot(trends, aes(nyrBT, Jtuexp, color = 'Jtu exp')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetaexp, color = 'Jbeta exp')) +
    geom_smooth(aes(y = Hornexp, color = 'Horn exp')) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = 'Half-saturation year') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20years-4.png)<!-- -->

``` r
ggplot(trends, aes(nyrBT, Jtumm, color = 'Jtu mm')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetamm, color = 'Jbeta mm')) +
    geom_smooth(aes(y = Hornmm, color = 'Horn mm')) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = 'Half-saturation year') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20years-5.png)<!-- -->

## Average change compared to duration

``` r
trends[, duration := maxyrBT - minyrBT + 1]

ggplot(trends, aes(duration, Jtutrendrem0, color = 'Jtu trend')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetatrendrem0, color = 'Jbeta trend')) +
    geom_smooth(aes(y = Horntrendrem0, color = 'Horn trend')) +
    scale_x_log10() +
    labs(y = 'Slope') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20duration-1.png)<!-- -->

``` r
ggplot(trends, aes(duration, Jtutrendz, color = 'Jtu z')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetatrendz, color = 'Jbeta z')) +
    geom_smooth(aes(y = Horntrendz, color = 'Horn z')) +
    scale_x_log10() +
    labs(y = 'Standardize slope') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20duration-2.png)<!-- -->

``` r
ggplot(trends, aes(duration, Jtulast, color = 'Jtu last')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetalast, color = 'Jbeta last')) +
    geom_smooth(aes(y = Hornlast, color = 'Horn last')) +
    scale_x_log10() +
    labs(y = 'Dissimilarity') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20duration-3.png)<!-- -->

``` r
ggplot(trends, aes(duration, Jtuexp, color = 'Jtu exp')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetaexp, color = 'Jbeta exp')) +
    geom_smooth(aes(y = Hornexp, color = 'Horn exp')) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = 'Half-saturation year') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20duration-4.png)<!-- -->

``` r
ggplot(trends, aes(duration, Jtumm, color = 'Jtu mm')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetamm, color = 'Jbeta mm')) +
    geom_smooth(aes(y = Hornmm, color = 'Horn mm')) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = 'Half-saturation year') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20duration-5.png)<!-- -->

## Average change compared to \#species

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

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20spp-1.png)<!-- -->

``` r
ggplot(trends, aes(Nspp, Jtutrendz, color = 'Jtu z')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetatrendz, color = 'Jbeta z')) +
    geom_smooth(aes(y = Horntrendz, color = 'Horn z')) +
    scale_x_log10() +
    labs(y = 'Slope') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20spp-2.png)<!-- -->

``` r
ggplot(trends, aes(Nspp, Jtulast, color = 'Jtu last')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetalast, color = 'Jbeta last')) +
    geom_smooth(aes(y = Hornlast, color = 'Horn last')) +
    scale_x_log10() +
    labs(y = 'Dissimilarity') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20spp-3.png)<!-- -->

``` r
ggplot(trends, aes(Nspp, Jtuexp, color = 'Jtu exp')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetaexp, color = 'Jbeta exp')) +
    geom_smooth(aes(y = Hornexp, color = 'Horn exp')) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = 'Half-saturation year') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20spp-4.png)<!-- -->

``` r
ggplot(trends, aes(Nspp, Jtumm, color = 'Jtu mm')) +
    geom_smooth() +
    geom_smooth(aes(y = Jbetamm, color = 'Jbeta mm')) +
    geom_smooth(aes(y = Hornmm, color = 'Horn mm')) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = 'Half-saturation year') +
    geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20spp-5.png)<!-- -->

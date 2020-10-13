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

## Remove studies that don’t meet quality thresholds

Keep observations with: - \>6 individuals Keep studies with: - \>2
species - count data - \>2 years of samples

``` r
length(unique(bt$rarefyID))
```

    ## [1] 53467

``` r
# samples with >6 indivs
bt <- bt[N > 6, ]

# number of years per study (in remaining samples)
nyrs <- bt[, .(nyrBT = length(YEAR)), by = rarefyID]

# >2 species, count data, >2 yrs
bt <- bt[rarefyID %in% rich[Nspp > 2, rarefyID] & 
           BROAD_TYPE == 'count' & 
           rarefyID %in% nyrs[nyrBT > 2, rarefyID], ]


length(unique(bt$rarefyID))
```

    ## [1] 32581

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
                  if(is.infinite(hsexp)) hsexp <- NA_real_
                  
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

    ## [1] 32581

## Add species richness to trends

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

### Plot null model sims

#### Histograms

``` r
# remove some duplicates
sim_summary <- as.data.table(sim_summary)
sim_summary <- sim_summary[!duplicated(sim_summary), ]

# hist
hist(sim_summary[metric == 'Jaccard (total)', sd], breaks = 200, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20hist-1.png)<!-- -->

``` r
hist(sim_summary[metric == 'Nestedness', sd], breaks = 200, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20hist-2.png)<!-- -->

``` r
hist(sim_summary[metric == 'Morisita-Horn', sd], breaks = 200, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20hist-3.png)<!-- -->

#### SD vs. expected

``` r
# plot
sim_summary[metric == 'Jaccard (total)', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20vs.%20exp-1.png)<!-- -->

``` r
sim_summary[metric == 'Nestedness', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20vs.%20exp-2.png)<!-- -->

``` r
sim_summary[metric == 'Morisita-Horn', plot(expected, log10(sd+1))]
```

    ## NULL

``` r
abline(h = 0, col = 'red')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20vs.%20exp-3.png)<!-- -->

#### SD vs Nyrs

``` r
# vs #yrs
sim_summary[metric == 'Nestedness', plot(Nyrs, sd+0.001, log = 'y', main = 'Null SD declines with #yrs (J tu)')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20nyrs-1.png)<!-- -->

``` r
sim_summary[metric == 'Jaccard (total)', plot(Nyrs, sd+0.001, log = 'y', main = 'Null SD declines with #yrs (J total)')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20nyrs-2.png)<!-- -->

``` r
sim_summary[metric == 'Morisita-Horn', plot(Nyrs, sd+0.001, log = 'y', main = 'Null SD declines with #yrs (Horn)')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20nyrs-3.png)<!-- -->

#### SD vs duration

``` r
sim_summary <- merge(sim_summary, trends[, .(rarefyID, duration = maxyrBT - minyrBT)], all.x = TRUE, by = 'rarefyID')
# vs #yrs
sim_summary[metric == 'Nestedness', plot(duration, sd+0.001, log = 'y', main = 'Null SD declines with duration (J tu)')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20duration-1.png)<!-- -->

``` r
sim_summary[metric == 'Jaccard (total)', plot(duration, sd+0.001, log = 'y', main = 'Null SD declines with duration (J total)')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20duration-2.png)<!-- -->

``` r
sim_summary[metric == 'Morisita-Horn', plot(duration, sd+0.001, log = 'y', main = 'Null SD declines with duration (Horn)')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20duration-3.png)<!-- -->

#### SD vs Nspp

``` r
# vs #spp
sim_summary[metric == 'Nestedness', plot(Nspp, sd+0.001, log = 'y')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20nspp-1.png)<!-- -->

``` r
sim_summary[metric == 'Jaccard (total)', plot(Nspp, sd+0.001, log = 'y')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20nspp-2.png)<!-- -->

``` r
sim_summary[metric == 'Morisita-Horn', plot(Nspp, sd+0.001, log = 'y')]
```

    ## NULL

``` r
abline(h = 0.001, col = 'red')
abline(h = 0.001+0.001, col = 'green')
abline(h = 0.01+0.001, col = 'grey')
```

![](calc_turnover_files/figure-gfm/plot%20sds%20nspp-3.png)<!-- -->

### Merge null model

``` r
# remove some duplicates
sim_summary <- as.data.table(sim_summary)
sim_summary <- sim_summary[!duplicated(sim_summary), ]

# merge
trends <- merge(trends, sim_summary[metric == 'Jaccard (total)', .(Jbeta_exp = expected, Jbeta_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
trends <- merge(trends, sim_summary[metric == 'Nestedness', .(Jtu_exp = expected, Jtu_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
trends <- merge(trends, sim_summary[metric == 'Morisita-Horn', .(Horn_exp = expected, Horn_sd = sd, rarefyID)], by = 'rarefyID', all.x = TRUE)
```

### Set a threshold and only standardize those timeseries with SD \> 0.001

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

    ## [1] 32581

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
    ##          32581          32581          32581          29030          29132 
    ##     Horntrendz        Jtulast      Jbetalast       Hornlast         Jtuexp 
    ##          29032          32581          32581          32581          22523 
    ##       Jbetaexp        Hornexp          Jtumm        Jbetamm         Hornmm 
    ##          22026          20894          32483          32506          32460

## Do some basic checks of the turnover calculations

``` r
# basic checks
trends
```

    ##           rarefyID  REALM                      Biome      taxa_mod
    ##     1:  100_606491 Marine     Northern_European_Seas          Fish
    ##     2:  101_606491 Marine     Northern_European_Seas Invertebrates
    ##     3: 108_3948499 Marine Continental_High_Antarctic         Birds
    ##     4: 108_3949201 Marine Continental_High_Antarctic         Birds
    ##     5: 108_3949203 Marine Continental_High_Antarctic         Birds
    ##    ---                                                            
    ## 32577:  91_1620531 Marine     Northern_European_Seas         Birds
    ## 32578:  91_1620532 Marine     Northern_European_Seas         Birds
    ## 32579:  91_1621261 Marine     Northern_European_Seas         Birds
    ## 32580:  91_1624191 Marine     Northern_European_Seas         Birds
    ## 32581:  91_1627107 Marine     Northern_European_Seas         Birds
    ##        STUDY_ID rarefyID_x rarefyID_y  Jtutrendrem0 Jtutrendrem0_se
    ##     1:      100   -3.08000   51.14000  2.976230e-03     0.001309059
    ##     2:      101   -3.08000   51.14000  0.000000e+00     0.000000000
    ##     3:      108   61.95750  -64.10500 -2.380952e-02             NaN
    ##     4:      108   61.97167  -66.84500 -4.000000e-01             NaN
    ##     5:      108   61.98500  -66.65667 -1.666667e-01             NaN
    ##    ---                                                             
    ## 32577:       91   20.94600   55.98300 -5.583900e-02     0.050492580
    ## 32578:       91   20.89200   56.03600  3.925231e-17             NaN
    ## 32579:       91   21.01500   56.04000  3.140185e-16             NaN
    ## 32580:       91   20.14625   57.30000 -6.666667e-01             NaN
    ## 32581:       91   20.76667   57.37167  5.000000e-01             NaN
    ##        Jbetatrendrem0 Jbetatrendrem0_se Horntrendrem0 Horntrendrem0_se
    ##     1:   2.587992e-03      0.0009975709  0.0023147105     0.0016045472
    ##     2:   2.208307e-03      0.0009103225  0.0001735803     0.0009087936
    ##     3:  -1.515152e-02               NaN -0.0002810624              NaN
    ##     4:  -1.111111e-01               NaN -0.9619098672              NaN
    ##     5:  -2.380952e-02               NaN -0.0540243109              NaN
    ##    ---                                                                
    ## 32577:  -2.419604e-02      0.0267409270 -0.0226161165     0.0337895033
    ## 32578:   1.587302e-02               NaN  0.0556650638              NaN
    ## 32579:   3.140185e-16               NaN -0.0799274013              NaN
    ## 32580:  -2.857143e-01               NaN  0.5887908368              NaN
    ## 32581:   3.500000e-01               NaN -0.0785534094              NaN
    ##          Jtutrendz Jbetatrendz   Horntrendz   Jtulast Jbetalast   Hornlast
    ##     1:          NA          NA           NA 0.4324324 0.5000000 0.36424173
    ##     2:  0.01012774  1.54387383           NA 0.0000000 0.5000000 0.15419602
    ##     3: -0.16171471 -0.09097099  0.002943179 0.8571429 0.9090909 0.99831363
    ##     4: -2.77813454 -0.63547169 -8.326522938 0.0000000 0.5555556 0.03120719
    ##     5: -1.90394821 -0.30780072 -0.800247051 0.0000000 0.8571429 0.67585413
    ##    ---                                                                    
    ## 32577:          NA          NA           NA 0.4444444 0.5454545 0.89295821
    ## 32578:          NA          NA           NA 0.5000000 0.7777778 0.98641246
    ## 32579:          NA          NA           NA 0.8333333 0.8333333 0.90547285
    ## 32580:          NA          NA           NA 0.0000000 0.4285714 0.73493976
    ## 32581:          NA          NA           NA 0.5000000 0.7500000 0.79578150
    ##            Jtuexp  Jbetaexp    Hornexp      Jtumm   Jbetamm    Hornmm
    ##     1: 0.02709627 0.2524249 1.22327984  0.1253377 0.1942427 1.0435632
    ##     2: 1.38629436 0.2599677         NA 15.5000000 0.1336071 0.0000000
    ##     3:         NA        NA         NA  0.0000000 0.0000000 0.0000000
    ##     4: 0.02943714        NA 0.03175053  0.0000000 0.0000000 0.0000000
    ##     5:         NA        NA         NA  0.0000000 0.0000000 0.0000000
    ##    ---                                                               
    ## 32577:         NA        NA         NA  0.0000000 0.0000000 0.0000000
    ## 32578: 0.06280353 0.5546708 0.95977976  0.0000000 0.2790698 0.5206255
    ## 32579:         NA        NA         NA  0.0000000 0.0000000 0.0000000
    ## 32580: 0.03140673        NA 1.60747190  0.0000000 0.0000000 1.7519826
    ## 32581: 3.00000002 1.1476482         NA  3.7699967 1.0557483 0.0000000
    ##        nyrBT minyrBT maxyrBT medianyrBT meanyrBT Nspp
    ##     1:    31    1981    2011       1996 1996.000   83
    ##     2:    31    1981    2011       1996 1996.000   15
    ##     3:     3    1981    1989       1983 1984.333   11
    ##     4:     3    1983    1985       1984 1984.000   11
    ##     5:     3    1986    1995       1989 1990.000    9
    ##    ---                                               
    ## 32577:     7    1992    1999       1995 1995.429   17
    ## 32578:     3    1993    1999       1995 1995.667    9
    ## 32579:     3    1993    1995       1994 1994.000   13
    ## 32580:     3    1992    1994       1993 1993.000    8
    ## 32581:     3    1992    1994       1993 1993.000    9

``` r
summary(trends$Jtutrendrem0)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -1.000000 -0.005726  0.001484  0.007008  0.020226  1.000000

``` r
summary(trends$Jbetatrendrem0)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -1.000000 -0.004611  0.003798  0.006793  0.017766  1.000000

``` r
summary(trends$Horntrendrem0)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -1.000000 -0.009260  0.003069  0.009457  0.024740  1.000000

``` r
summary(trends$Jtutrendz)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ## -58.283  -0.311   0.063   0.358   0.983  61.472    3551

``` r
summary(trends$Jbetatrendz)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ## -31.591  -0.177   0.137   0.303   0.655  84.690    3449

``` r
summary(trends$Horntrendz)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ## -74.446  -0.439   0.134   0.604   1.343 104.775    3549

``` r
summary(trends$Jtulast)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.2222  0.4706  0.4479  0.6667  1.0000

``` r
summary(trends$Jbetalast)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.5000  0.6400  0.6329  0.7778  1.0000

``` r
summary(trends$Hornlast)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.1854  0.5362  0.5326  0.9020  1.0000

``` r
summary(trends$Jtuexp)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   0.002   0.528   1.386   5.522   4.326 648.434   10058

``` r
summary(trends$Jbetaexp)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   0.002   0.330   0.768   2.174   2.070 363.489   10555

``` r
summary(trends$Hornexp)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    ##    0.002    0.520    1.645    7.530    5.499 1683.947    11687

``` r
summary(trends$Jtumm)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    ##   0.0000   0.0000   0.4911   5.7149   4.5000 924.1801       98

``` r
summary(trends$Jbetamm)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    ##   0.0000   0.0000   0.1852   1.7254   1.1333 511.4554       75

``` r
summary(trends$Hornmm)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
    ##    0.0000    0.0000    0.2619   10.8118    3.4915 1000.0000       121

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
    ##   3.000   4.000   5.000   6.913   9.000  97.000

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
    ##    3.00    8.00   13.00   15.88   21.00   98.00

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
    ##   3.000   4.000   5.000   6.913   9.000  97.000

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
trends[, duration := maxyrBT - minyrBT]

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

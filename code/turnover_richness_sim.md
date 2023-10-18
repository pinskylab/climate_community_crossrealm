Simulations to test effects of random sampling and species richness
interactions
================

- [Simple scenario of species
  turnover](#simple-scenario-of-species-turnover)
  - [100% of species sampled](#100-of-species-sampled)
  - [50% of species sampled](#50-of-species-sampled)
    - [Calculate turnover rates](#calculate-turnover-rates)
    - [Plot turnover rates as lines](#plot-turnover-rates-as-lines)
    - [Plot turnover rates as
      densities](#plot-turnover-rates-as-densities)

Do assemblages with higher species richness have higher turnover rates?
Let’s test this with some simple simulations. If true, it might explain
higher turnover rates at higher temperatures, since species richness is
often higher in the tropics.

``` r
# set parameters for simulations
require(vegan)
```

    ## Loading required package: vegan

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-4

``` r
require(reshape2)
```

    ## Loading required package: reshape2

``` r
require(data.table)
```

    ## Loading required package: data.table

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:reshape2':
    ## 
    ##     dcast, melt

``` r
lown <- 100 # number of species in low richness treatment
highn <- 1000 # number of species in high richness treatment
ngen <- 10 # number of sampling periods
frac <- 0.05 # fraction of species added and removed each sampling period
prob <- 0.5 # probability of seeing a species that is present
reps <- 1000 # number of times to repeat the sampling
set.seed(1001)
```

# Simple scenario of species turnover

Add 5% new species and remove 5% of existing species each time step. A
high and a low richness treatment.

## 100% of species sampled

When we sample all species present, the turnover rates don’t depend on
species richness. They are exactly the same in both species richness
treatments.

``` r
# make the communities
low <- matrix(c(rep(1, lown), rep(0, 2*lown*ngen - lown)), nrow=2*lown, ncol = ngen) # rows are species. in first time step, half are present
high <- matrix(c(rep(1, highn), rep(0, 2*highn*ngen - highn)), nrow = 2*highn, ncol = ngen)
for(i in 2:10){
    low[,i] <- low[,i-1] * c(rep(0, (i-1)*frac*lown), rep(1, 2*lown - (i-1)*frac*lown)) + c(rep(0, lown + (i-2)*frac*lown), rep(1, frac*lown), rep(0, lown - (i-1)*frac*lown)) # remove frac existing species and add frac new species
    high[,i] <- high[,i-1] * c(rep(0, (i-1)*frac*highn), rep(1, 2*highn - (i-1)*frac*highn)) + c(rep(0, highn + (i-2)*frac*highn), rep(1, frac*highn), rep(0, highn - (i-1)*frac*highn))
}

# calculate dissimilarities
lowdist <- vegdist(t(low), method='jaccard', binary=TRUE) # vegan expects spp in columns
highdist <- vegdist(t(high), method='jaccard', binary=TRUE)

# examine distances
print('low species richness')
```

    ## [1] "low species richness"

``` r
lowdist
```

    ##            1         2         3         4         5         6         7
    ## 2  0.0952381                                                            
    ## 3  0.1818182 0.0952381                                                  
    ## 4  0.2608696 0.1818182 0.0952381                                        
    ## 5  0.3333333 0.2608696 0.1818182 0.0952381                              
    ## 6  0.4000000 0.3333333 0.2608696 0.1818182 0.0952381                    
    ## 7  0.4615385 0.4000000 0.3333333 0.2608696 0.1818182 0.0952381          
    ## 8  0.5185185 0.4615385 0.4000000 0.3333333 0.2608696 0.1818182 0.0952381
    ## 9  0.5714286 0.5185185 0.4615385 0.4000000 0.3333333 0.2608696 0.1818182
    ## 10 0.6206897 0.5714286 0.5185185 0.4615385 0.4000000 0.3333333 0.2608696
    ##            8         9
    ## 2                     
    ## 3                     
    ## 4                     
    ## 5                     
    ## 6                     
    ## 7                     
    ## 8                     
    ## 9  0.0952381          
    ## 10 0.1818182 0.0952381

``` r
print('high species richness')
```

    ## [1] "high species richness"

``` r
highdist
```

    ##            1         2         3         4         5         6         7
    ## 2  0.0952381                                                            
    ## 3  0.1818182 0.0952381                                                  
    ## 4  0.2608696 0.1818182 0.0952381                                        
    ## 5  0.3333333 0.2608696 0.1818182 0.0952381                              
    ## 6  0.4000000 0.3333333 0.2608696 0.1818182 0.0952381                    
    ## 7  0.4615385 0.4000000 0.3333333 0.2608696 0.1818182 0.0952381          
    ## 8  0.5185185 0.4615385 0.4000000 0.3333333 0.2608696 0.1818182 0.0952381
    ## 9  0.5714286 0.5185185 0.4615385 0.4000000 0.3333333 0.2608696 0.1818182
    ## 10 0.6206897 0.5714286 0.5185185 0.4615385 0.4000000 0.3333333 0.2608696
    ##            8         9
    ## 2                     
    ## 3                     
    ## 4                     
    ## 5                     
    ## 6                     
    ## 7                     
    ## 8                     
    ## 9  0.0952381          
    ## 10 0.1818182 0.0952381

## 50% of species sampled

This simulation adds sampling noise such that a random 50% of species
are observed. This is repeated many times to understand central
tendencies.

``` r
for(i in 1:reps){
    if(i %% 50 == 0) cat(paste0('#', i, " "))
    # make the observed communities with only X% of species observed
    lowobs <- low * sample(c(rep(1, prob*2*lown*ngen), rep(0, (1-prob)*2*lown*ngen)), 2*lown*ngen, replace=FALSE) # turn X% of obs to absent, randomly
    highobs <- high * sample(c(rep(1, prob*2*highn*ngen), rep(0, (1-prob)*2*highn*ngen)), 2*highn*ngen, replace=FALSE)
    
    
    # calculate dissimilarities
    lowobsdist <- as.matrix(vegdist(t(lowobs), method='jaccard', binary=TRUE)) # vegan expects spp in columns
    highobsdist <- as.matrix(vegdist(t(highobs), method='jaccard', binary=TRUE))
    
    # examine distances (first example only)
    if(i == 1){
        print('low species richness with obs error')
        print(lowobsdist)
        print('high species richness')
        print(highobsdist)
    }
    
    # convert to long format (only one triangle)
    dimnames(lowobsdist) <- list(1:ngen, 1:ngen)
    lownms <- t(combn(colnames(lowobsdist), 2))
    thislowobslong <- data.frame(lownms, dist=lowobsdist[lownms])
    thislowobslong$rep <- i
    
    dimnames(highobsdist) <- list(1:ngen, 1:ngen)
    highnms <- t(combn(colnames(highobsdist), 2))
    thishighobslong <- data.frame(highnms, dist=highobsdist[highnms])
    thishighobslong$rep <- i
    
    if(i == 1){
        lowobslong <- thislowobslong
        highobslong <- thishighobslong
    } else {
        lowobslong <- rbind(lowobslong, thislowobslong)
        highobslong <- rbind(highobslong, thishighobslong)
    }
}
```

    ## [1] "low species richness with obs error"
    ##            1         2         3         4         5         6         7
    ## 1  0.0000000 0.7236842 0.6666667 0.7468354 0.7333333 0.7750000 0.7738095
    ## 2  0.7236842 0.0000000 0.7051282 0.6849315 0.8051948 0.7974684 0.8636364
    ## 3  0.6666667 0.7051282 0.0000000 0.6623377 0.6800000 0.7710843 0.7977528
    ## 4  0.7468354 0.6849315 0.6623377 0.0000000 0.7466667 0.7402597 0.6923077
    ## 5  0.7333333 0.8051948 0.6800000 0.7466667 0.0000000 0.7083333 0.7594937
    ## 6  0.7750000 0.7974684 0.7710843 0.7402597 0.7083333 0.0000000 0.6533333
    ## 7  0.7738095 0.8636364 0.7977528 0.6923077 0.7594937 0.6533333 0.0000000
    ## 8  0.7710843 0.8352941 0.7674419 0.7530864 0.8024691 0.7341772 0.6538462
    ## 9  0.7875000 0.8250000 0.7375000 0.7532468 0.7222222 0.7500000 0.5915493
    ## 10 0.8313253 0.9195402 0.8522727 0.8148148 0.8354430 0.7662338 0.6301370
    ##            8         9        10
    ## 1  0.7710843 0.7875000 0.8313253
    ## 2  0.8352941 0.8250000 0.9195402
    ## 3  0.7674419 0.7375000 0.8522727
    ## 4  0.7530864 0.7532468 0.8148148
    ## 5  0.8024691 0.7222222 0.8354430
    ## 6  0.7341772 0.7500000 0.7662338
    ## 7  0.6538462 0.5915493 0.6301370
    ## 8  0.0000000 0.7307692 0.6800000
    ## 9  0.7307692 0.0000000 0.7631579
    ## 10 0.6800000 0.7631579 0.0000000
    ## [1] "high species richness"
    ##            1         2         3         4         5         6         7
    ## 1  0.0000000 0.6747759 0.7182741 0.7424623 0.7183272 0.7500000 0.7882070
    ## 2  0.6747759 0.0000000 0.6907895 0.6818182 0.7146433 0.7327478 0.7887668
    ## 3  0.7182741 0.6907895 0.0000000 0.7271523 0.7316456 0.7124183 0.7759494
    ## 4  0.7424623 0.6818182 0.7271523 0.0000000 0.6760000 0.7157895 0.7590674
    ## 5  0.7183272 0.7146433 0.7316456 0.6760000 0.0000000 0.6877419 0.7042802
    ## 6  0.7500000 0.7327478 0.7124183 0.7157895 0.6877419 0.0000000 0.6976127
    ## 7  0.7882070 0.7887668 0.7759494 0.7590674 0.7042802 0.6976127 0.0000000
    ## 8  0.8179775 0.7885514 0.7643468 0.7678133 0.7283800 0.6867008 0.7030848
    ## 9  0.8257840 0.8014440 0.7997528 0.8019925 0.7586634 0.7534766 0.7110519
    ## 10 0.8422857 0.8231132 0.8105134 0.7843632 0.7531017 0.7741935 0.7336815
    ##            8         9        10
    ## 1  0.8179775 0.8257840 0.8422857
    ## 2  0.7885514 0.8014440 0.8231132
    ## 3  0.7643468 0.7997528 0.8105134
    ## 4  0.7678133 0.8019925 0.7843632
    ## 5  0.7283800 0.7586634 0.7531017
    ## 6  0.6867008 0.7534766 0.7741935
    ## 7  0.7030848 0.7110519 0.7336815
    ## 8  0.0000000 0.6912145 0.7053571
    ## 9  0.6912145 0.0000000 0.6873315
    ## 10 0.7053571 0.6873315 0.0000000
    ## #50 #100 #150 #200 #250 #300 #350 #400 #450 #500 #550 #600 #650 #700 #750 #800 #850 #900 #950 #1000

``` r
lowobslong$X1 <- as.numeric(lowobslong$X1)
lowobslong$X2 <- as.numeric(lowobslong$X2)

highobslong$X1 <- as.numeric(highobslong$X1)
highobslong$X2 <- as.numeric(highobslong$X2)

lowobslong <- data.table(lowobslong)
highobslong <- data.table(highobslong)
```

### Calculate turnover rates

Calculate turnover rates as the slope of dissimilarity vs. time. We find
the same central tendencies in both high and low species richness
scenarios.

``` r
lowslopes <- lowobslong[, as.list(coef(lm(dist ~ I(X2-X1)))), by=rep]
highslopes <- highobslong[, as.list(coef(lm(dist ~ I(X2-X1)))), by=rep]
setnames(lowslopes, c('rep', 'int', 'slope'))
setnames(highslopes, c('rep', 'int', 'slope'))

print('low richness slopes')
```

    ## [1] "low richness slopes"

``` r
lowslopes[, summary(slope)]
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## 0.007309 0.016866 0.019464 0.019356 0.021881 0.030277

``` r
print('high richness slopes')
```

    ## [1] "high richness slopes"

``` r
highslopes[, summary(slope)]
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 0.01604 0.01866 0.01939 0.01940 0.02012 0.02323

### Plot turnover rates as lines

We can plot all the turnover rates as sloping slopes. We find more
variance in turnover rates for the lower richness scenario, though note
this doesn’t plot correctly in the github document (only one line is
shown when there should be 1000).

``` r
par(mfrow=c(1,2))
plot(c(0,10), c(0,1), col='white', xlab = 'temporal distance', ylab='dissimilarity', main = 'low richness')
lowslopes[, abline(a = int, b = slope, col = '#00000033'), by = rep]
```

    ## Empty data.table (0 rows and 1 cols): rep

``` r
plot(c(0,10), c(0,1), col='white', xlab = 'temporal distance', ylab='dissimilarity', main = 'high richness')
highslopes[, abline(a = int, b = slope, col = '#00000033'), by = rep]
```

![](turnover_richness_sim_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

    ## Empty data.table (0 rows and 1 cols): rep

### Plot turnover rates as densities

Another way to plot the turnover rates for high and low richness
scenarios. In addition to the same central tendencies irrespective of
richness, we do see more variance in turnover rates for the lower
richness scenario (black) compared to higher richness (red).

``` r
lowslopes[, plot(density(slope), ylim =c(0,500), main = 'Turnover rates')]
```

    ## NULL

``` r
highslopes[, lines(density(slope), col='red')]
```

![](turnover_richness_sim_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

    ## NULL

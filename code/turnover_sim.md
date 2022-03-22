Turnover simulations to test effect of duration
================

  - [Basic settings](#basic-settings)
  - [Neutral simulation with an even
    metacommunity](#neutral-simulation-with-an-even-metacommunity)
  - [Neutral simulation with an uneven
    metacommunity](#neutral-simulation-with-an-uneven-metacommunity)
  - [Plot slope vs. duration](#plot-slope-vs.-duration)

# Basic settings

``` r
## Basic settings
nburn = 1000 # length of burnin
ngen = 500 # number of simulated generations
nrep = 5 # number of replicates
nindiv = 100 # number of local individuals
nsamp = 100 # number of local individuals sampled (if == nindiv, then sample all of them)
prob = 0.1 # probability of replacing an individual from the regional pool rather than locally
```

# Neutral simulation with an even metacommunity

Each replicate has a burnin and then a simulation period. Do many
replicates. Plot dynamics from the first replicate.

``` r
set.seed(6)
metaeven = rep(1:100, 100) # metacommunity of 100 evenly distributed species

simburnineven <- vector('list', nrep)
simeven <- vector('list', nrep)
names(simeven) <- 1:nrep
for(i in 1:nrep){
    simburnineven[[i]] <- untb(start = as.count(sample(metaeven, nindiv)), prob = prob, D = 1, gens = nburn, keep = TRUE, meta = as.count(metaeven)) # burnin
    simeven[[i]] <- untb(start = simburnineven[[i]][nburn,], prob = prob, D=1, gens = ngen, keep = TRUE, meta = as.count(metaeven)) # simulation
}

# plot an example from the first
plot(as.count(metaeven), main = 'Regional SAD')
```

![](turnover_sim_files/figure-gfm/sim%20even-1.png)<!-- -->

``` r
plot(species.count(simburnineven[[1]]), type = 'b', main = 'Burnin: num spp through time') # evaluate burnin
```

![](turnover_sim_files/figure-gfm/sim%20even-2.png)<!-- -->

``` r
matplot(species.table(simburnineven[[1]]), type='l', lty=1, main = 'Burnin: spp dynamics')
```

![](turnover_sim_files/figure-gfm/sim%20even-3.png)<!-- -->

``` r
plot(count(simburnineven[[1]][1,]), main = 'Burnin: Initial local SAD')
```

![](turnover_sim_files/figure-gfm/sim%20even-4.png)<!-- -->

``` r
plot(count(simburnineven[[1]][nburn,]), main = 'Burnin: Final local SAD')
```

![](turnover_sim_files/figure-gfm/sim%20even-5.png)<!-- -->

``` r
plot(species.count(simeven[[1]]), type = 'b', main = 'Sim: num spp through time')
```

![](turnover_sim_files/figure-gfm/sim%20even-6.png)<!-- -->

``` r
matplot(species.table(simeven[[1]]), type='l', lty=1, main = 'Sim: spp dynamics')
```

![](turnover_sim_files/figure-gfm/sim%20even-7.png)<!-- -->

``` r
plot(count(simeven[[1]][ngen,]), main = 'Sim: Final local SAD')
```

![](turnover_sim_files/figure-gfm/sim%20even-8.png)<!-- -->

Sample from each annual community and calculate slope of Jaccard
dissimilarity among pairs of years for timeseries of 3:20 years.

``` r
# sample from the community
sampeven = lapply(simeven, FUN = function(x) t(apply(x, MARGIN = 1, FUN = sample, size = nsamp)))

# calculate dissimilarity for all pairwise communities
disteven <- lapply(sampeven, FUN = function(x) as.matrix(vegdist(tocommat(x), method='jaccard', binary=TRUE)))

# convert to long format
for(i in 1:length(disteven)){
    dimnames(disteven[[i]]) <- list(1:ngen, 1:ngen)
    xy <- t(combn(colnames(disteven[[i]]), 2))
    temp <- data.table(year1 = as.numeric(xy[,1]), year2 = as.numeric(xy[,2]), dist = as.numeric(disteven[[i]][xy]), rarefyID = i)
    if(i == 1){
        distlongeven <- temp
    } else {
        distlongeven <- rbind(distlongeven, temp)
    }
}

# calculate slopes. slow.
yrslist <- 3:100
for(yr in yrslist){
  cat(yr)
  temp <- distlongeven[, calctrendnsampsall(dist, year1, year2, numyrs = yr, nsamps = 3, 
                                            measure = 'Jtu', duration_group = paste0(yr, 'min3')), by = rarefyID]
  if(yr == min(yrslist)) trendseven = temp[!is.na(disstrend), ] # make a new dataset if first iteration through
  if(yr > min(yrslist)) trendseven = rbind(trendseven,temp[!is.na(disstrend), ]) # otherwise append
  
}
```

    ## 3456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100

Plot all pairwise dissimilarities and dissimilarities vs. year 1 for an
example

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](turnover_sim_files/figure-gfm/plot%20even-1.png)<!-- -->

# Neutral simulation with an uneven metacommunity

Use BCI data as the metacommunity.
![](turnover_sim_files/figure-gfm/sim%20SAD-1.png)<!-- -->![](turnover_sim_files/figure-gfm/sim%20SAD-2.png)<!-- -->![](turnover_sim_files/figure-gfm/sim%20SAD-3.png)<!-- -->![](turnover_sim_files/figure-gfm/sim%20SAD-4.png)<!-- -->![](turnover_sim_files/figure-gfm/sim%20SAD-5.png)<!-- -->![](turnover_sim_files/figure-gfm/sim%20SAD-6.png)<!-- -->![](turnover_sim_files/figure-gfm/sim%20SAD-7.png)<!-- -->![](turnover_sim_files/figure-gfm/sim%20SAD-8.png)<!-- -->

Sample from them and calculate Jaccard dissimilarity among pairs of
years.

    ## 3456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100

Plot all dissimilarities for an example

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](turnover_sim_files/figure-gfm/plot%20SAD%20example-1.png)<!-- -->

# Plot slope vs. duration

Top row has an example from one replicate. Middle row shows slopes from
all replicates. Bottom row is the average across all replicates. Left
column is even metacommunity, right column is BCI metacommunity (a
skewed SAD).

    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
    ## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

![](turnover_sim_files/figure-gfm/plot%20slopes-1.png)<!-- -->

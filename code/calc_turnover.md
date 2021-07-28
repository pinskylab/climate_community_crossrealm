Calculate biodiversity turnover as slope of dissimilarity vs. time
================

# Set up

# Load data

``` r
# BioTime
# biotime community dissimilarity data
load(here::here('data', 'biotime_blowes', 'all_pairs_beta.Rdata')) # load rarefied_beta_medians, which has all pairwise dissimilarities
bt <- data.table(rarefied_beta_medians); rm(rarefied_beta_medians)
bt[, year1 := as.numeric(year1)] # not sure why it gets read in as character
bt[, year2 := as.numeric(year2)]
bt[, Horn := 1-Hornsim] # convert similarity to dissimilarity
bt[, Hornsim := NULL]

# BioTime data type
load('data/biotime_blowes/time_series_data_type.Rdata') # loads rarefyID_type
bt <- merge(bt, rarefyID_type, by = 'rarefyID', all.x = TRUE) # merge with biotime

# biotime taxa category and other info
load(here::here('data', 'biotime_blowes', 'bt_malin.Rdata')) # load bt_malin
btinfo <- data.table(bt_malin); rm(bt_malin)
btinfo2 <- btinfo[!duplicated(rarefyID), .(rarefyID, rarefyID_x, rarefyID_y, Biome, taxa_mod, REALM, STUDY_ID)]
bt <- merge(bt, btinfo2, by = 'rarefyID', all.x = TRUE)

# richness
rich <- fread('output/richness_by_rarefyID.csv.gz') # number of species
```

## Examine number of species and individuals per sample

``` r
hist(log10(btinfo$N))
```

![](calc_turnover_files/figure-gfm/numspp%20per%20samp-1.png)<!-- -->

``` r
hist(log10(btinfo$S)); abline(v = log10(2), col = 'red') # line at 2 species
```

![](calc_turnover_files/figure-gfm/numspp%20per%20samp-2.png)<!-- -->

``` r
nrow(btinfo)
```

    ## [1] 295095

``` r
btinfo[S == 1, .N]
```

    ## [1] 37212

## Remove studies that don’t meet quality thresholds

Keep studies with: - \>2 species

``` r
length(unique(bt$rarefyID))
```

    ## [1] 59036

``` r
# number of years per study (in remaining samples)
nyrs <- bt[, .(nyrBT = length(unique(c(year1, year2)))), by = rarefyID]

# >2 species, >2 yrs
bt <- bt[rarefyID %in% rich[Nspp > 2, rarefyID] & 
           rarefyID %in% nyrs[nyrBT > 2, rarefyID], ]

length(unique(bt$rarefyID))
```

    ## [1] 38487

# Plot turnover example

``` r
test <- bt[rarefyID == '339_1085477',]
test[, dY := year2 - year1]
ggplot(test, aes(dY, Jtu)) +
    geom_point(alpha = 0.5) +
    geom_smooth(data = test, 
                method = 'lm') + # a straight line fit
    labs(x = 'Year', y = 'Jaccard turnover')
```

    ## `geom_smooth()` using formula 'y ~ x'

![](calc_turnover_files/figure-gfm/plot%20turnover-1.png)<!-- -->

``` r
# a declining 3-year slope
test <- bt[rarefyID == '112_529867' & year1>= 2002 & year2 <= 2004, ]
test[, dY := year2 - year1]
ggplot(test, aes(dY, Jtu)) +
    geom_point(alpha = 0.5) +
    geom_smooth(data = test, 
                method = 'lm') + # a straight line fit
    labs(x = 'Difference in years', y = 'Jaccard turnover')
```

    ## `geom_smooth()` using formula 'y ~ x'

![](calc_turnover_files/figure-gfm/plot%20turnover-2.png)<!-- -->

``` r
# a declining 20-year slope
test <- bt[rarefyID == '213_435199' & year1>= 1987 & year2 <= 2006, ]
test[, dY := year2 - year1]
ggplot(test, aes(dY, Jtu)) +
    geom_point(alpha = 0.5) +
    geom_smooth(data = test, 
                method = 'lm') + # a straight line fit
    labs(x = 'Difference in years', y = 'Jaccard turnover')
```

    ## `geom_smooth()` using formula 'y ~ x'

![](calc_turnover_files/figure-gfm/plot%20turnover-3.png)<!-- -->

# Calculate temporal trends (temporal turnover)

This only recalculates from scratch if temp/trendstemp.rds is not
available. \#\# Set up functions

``` r
# function to calc linear trend from all years of the time-series
calctrend <- function(y, year1, year2, nm = 'y'){
  dy <- year2 - year1
  if(length(dy)>1){
    mod <- lm(y ~ dy) # fit line
    out <- list(y = coef(mod)[2], # coef for the slope
                y_se = sqrt(diag(vcov(mod)))[2],
                year1 = min(year1), 
                year2 = max(year2))
    names(out) <- c(nm, paste0(nm, '_se'), paste0(nm, '_y1'), paste0(nm, '_y2'))
    return(out)
  } else{
    out <- list(y = NA_real_, y_se = NA_real_, year1 = NA_real_, year2 = NA_real_)
    names(out) <- c(nm, paste0(nm, '_se'), paste0(nm, '_y1'), paste0(nm, '_y2'))
    return(out)
  }
}

calctrendlast <- function(y, year1, year2, numyrs, nm = 'y'){
  yrs <- sort(unique(c(year1, year2)))
  dy <- diff(yrs) # find intervals between years
  rl <- rle(dy) # run length encoding
  end = cumsum(rl$lengths) # find ends of runs
  start = c(1, head(end, -1) + 1) # find starts of runs
  rlkeep <- which(rl$lengths >= numyrs & rl$values == 1) # find runs with desired length and 1 year intervals
  if(length(rlkeep)>0){
    rlkeep <- max(rlkeep) # find the latest run that meets our criteria
    start <- start[rlkeep] # start of this run
    end <- end[rlkeep] # end of this run
    dykeep <- rep(FALSE, length(dy)) # year intervals to keep
    dykeep[start:end] <- TRUE
    yrs2 <- yrs[c(dykeep, FALSE) | c(FALSE, dykeep)] # keep years involved in at least one interval to keep
    maxyr <- max(yrs2)
    yrs3 <- yrs2[yrs2 > maxyr - numyrs] # only last numyrs years
    ykeep <- year1 %in% yrs3 & year2 %in% yrs3 # keep values in pairwise comparisons for the years we want
    y <- y[ykeep] # trim the timeseries
    year1 <- year1[ykeep]
    year2 <- year2[ykeep]
    run <- TRUE
  } else{
    run <- FALSE
  }
  
  if(run){
    dy <- year2 - year1
    
    mod <- lm(y ~ dy) # fit line
    out <- list(y = coef(mod)[2], # coef for the slope
                y_se = sqrt(diag(vcov(mod)))[2],
                year1 = min(yrs3), 
                year2 = max(yrs3)) # SE
    names(out) <- c(nm, paste0(nm, '_se'), paste0(nm, '_y1'), paste0(nm, '_y2'))
    return(out)

  } else {
    out <- list(y = NA_real_, y_se = NA_real_, year1 = NA_real_, year2 = NA_real_)
    names(out) <- c(nm, paste0(nm, '_se'), paste0(nm, '_y1'), paste0(nm, '_y2'))
    return(out)
  }
  

}
```

## Do calculations

And write out trendstemp.rds

``` r
setkey(bt, STUDY_ID, rarefyID, year1,  year2)

if(file.exists('temp/trendstemp.rds')){
    print('File already exists. Will not do calculations')
    trends2 <- readRDS('temp/trendstemp.rds')
} else {
  print('Calculating from scratch')
  
  # 3-year trends
  trends1 <- bt[, calctrendlast(Jtu, year1, year2, 3, 'Jtutrend3'), 
                by = .(REALM, Biome, taxa_mod, STUDY_ID, 
                       rarefyID, rarefyID_x, rarefyID_y)] # calculate trend in Jaccard turnover
  trends2 <- bt[, calctrendlast(Jbeta, year1, year2, 3, 'Jbetatrend3'),
                by = .(rarefyID)]
  trends3 <- bt[!is.na(Horn), calctrendlast(Horn, year1, year2, 3, 'Horntrend3'),
                by = .(rarefyID)]
  
  # 5 year
  trends4 <- bt[, calctrendlast(Jtu, year1, year2, 5, 'Jtutrend5'), 
                by = .(rarefyID)]
  trends5 <- bt[, calctrendlast(Jbeta, year1, year2, 5, 'Jbetatrend5'),
                by = .(rarefyID)]
  trends6 <- bt[!is.na(Horn), calctrendlast(Horn, year1, year2, 5, 'Horntrend5'),
                by = .(rarefyID)]


  # 10 year
  trends7 <- bt[, calctrendlast(Jtu, year1, year2, 10, 'Jtutrend10'), 
                by = .(rarefyID)]
  trends8 <- bt[, calctrendlast(Jbeta, year1, year2, 10, 'Jbetatrend10'),
                by = .(rarefyID)]
  trends9 <- bt[!is.na(Horn), calctrendlast(Horn, year1, year2, 10, 'Horntrend10'),
                by = .(rarefyID)]
  
  # 20 year
  trends10 <- bt[, calctrendlast(Jtu, year1, year2, 20, 'Jtutrend20'), 
                by = .(rarefyID)]
  trends11 <- bt[, calctrendlast(Jbeta, year1, year2, 20, 'Jbetatrend20'),
                by = .(rarefyID)]
  trends12 <- bt[!is.na(Horn), calctrendlast(Horn, year1, year2, 20, 'Horntrend20'),
                by = .(rarefyID)]
  
  # All years available
  trends13 <- bt[, calctrend(Jtu, year1, year2, 'JtutrendAll'), by = .(rarefyID)]
  trends14 <- bt[, calctrend(Jbeta, year1, year2, 'JbetatrendAll'), by = .(rarefyID)]
  trends15 <- bt[!is.na(Horn), calctrend(Horn, year1, year2, 'HorntrendAll'), by = .(rarefyID)]
  
  trends <- merge(trends1, trends2, all = TRUE)
  trends <- merge(trends, trends3, all = TRUE)
  trends <- merge(trends, trends4, all = TRUE)
  trends <- merge(trends, trends5, all = TRUE)
  trends <- merge(trends, trends6, all = TRUE)
  trends <- merge(trends, trends7, all = TRUE)
  trends <- merge(trends, trends8, all = TRUE)
  trends <- merge(trends, trends9, all = TRUE)
  trends <- merge(trends, trends10, all = TRUE)
  trends <- merge(trends, trends11, all = TRUE)
  trends <- merge(trends, trends12, all = TRUE)
  trends <- merge(trends, trends13, all = TRUE)
  trends <- merge(trends, trends14, all = TRUE)
  trends <- merge(trends, trends15, all = TRUE)
  
  trends2 <- trends[!is.na(Jtutrend3) | !is.na(Jbetatrend3) | !is.na(Horntrend3), ]

  saveRDS(trends2, file = 'temp/trendstemp.rds')
  
}
```

    ## [1] "Calculating from scratch"

``` r
nrow(trends2)
```

    ## [1] 10785

## Add species richness to trends

``` r
trends2 <- merge(trends2, rich, all.x = TRUE) # species richness
```

## Plot every Jtu timeseries (a lot\!)

Not run during knitting

``` r
rids <- trends2[nyrBT > 2, sort(unique(rarefyID))]
setkey(bt, rarefyID, YEAR)
print(paste(length(rids), ' rarefyIDs'))
filenum <- 1
plotnum <- 1
for(i in 1:length(rids)){
#for(i in 1:400){ # for testing
  jtutrendrem0 <- trends2[rarefyID == rids[i], Jtutrendrem0]
  jtuexp <- trends2[rarefyID == rids[i], Jtuexp]
  jtumm <- trends2[rarefyID == rids[i], Jtumm]
  nspp <- trends2[rarefyID == rids[i], Nspp]
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

# Examine the turnover calculations

## How many values?

``` r
apply(trends2[, .(Jtutrend3, Jbetatrend3, Horntrend3,
                 Jtutrend5, Jbetatrend5, Horntrend5,
                Jtutrend10, Jbetatrend10, Horntrend10,
                Jtutrend20, Jbetatrend20, Horntrend20,
                JtutrendAll, JbetatrendAll, HorntrendAll)], MARGIN = 2, 
      function(x) sum(!is.na(x)))
```

    ##     Jtutrend3   Jbetatrend3    Horntrend3     Jtutrend5   Jbetatrend5 
    ##         10785         10785         10492          5563          5563 
    ##    Horntrend5    Jtutrend10  Jbetatrend10   Horntrend10    Jtutrend20 
    ##          5385          2183          2183          2120           494 
    ##  Jbetatrend20   Horntrend20   JtutrendAll JbetatrendAll  HorntrendAll 
    ##           494           483         10785         10785         10492

## Do some basic checks of the turnover calculations

``` r
# basic checks
trends2
```

    ##           rarefyID  REALM                      Biome      taxa_mod STUDY_ID
    ##     1:  100_606491 Marine     Northern_European_Seas          Fish      100
    ##     2:  101_606491 Marine     Northern_European_Seas Invertebrates      101
    ##     3: 108_4114762 Marine Continental_High_Antarctic         Birds      108
    ##     4: 108_4183957 Marine Continental_High_Antarctic         Birds      108
    ##     5: 108_4595757 Marine Southeast_Australian_Shelf         Birds      108
    ##    ---                                                                     
    ## 10781:  91_1619798 Marine     Northern_European_Seas         Birds       91
    ## 10782:  91_1619799 Marine     Northern_European_Seas         Birds       91
    ## 10783:  91_1620530 Marine     Northern_European_Seas         Birds       91
    ## 10784:  91_1620531 Marine     Northern_European_Seas         Birds       91
    ## 10785:  91_1624191 Marine     Northern_European_Seas         Birds       91
    ##        rarefyID_x rarefyID_y   Jtutrend3 Jtutrend3_se Jtutrend3_y1 Jtutrend3_y2
    ##     1:   -3.08000   51.14000  0.06031746 1.044729e-01         2009         2011
    ##     2:   -3.08000   51.14000  0.00000000 0.000000e+00         2009         2011
    ##     3:  109.54444  -65.44889  0.41666667 4.330127e-01         1994         1996
    ##     4:  139.82600  -65.22800  0.50000000 8.660254e-01         1994         1996
    ##     5:  146.01167  -44.73833 -0.16666667 2.634496e-16         1995         1997
    ##    ---                                                                         
    ## 10781:   21.03750   55.69375 -0.25000000 4.330127e-01         1994         1996
    ## 10782:   20.97727   55.75727  0.12500000 2.165064e-01         1997         1999
    ## 10783:   20.97500   55.90875  0.17460317 3.024216e-01         1993         1995
    ## 10784:   20.94600   55.98300  0.36666667 5.196152e-01         1993         1995
    ## 10785:   20.14625   57.30000 -0.16666667 2.886751e-01         1993         1995
    ##        Jbetatrend3 Jbetatrend3_se Jbetatrend3_y1 Jbetatrend3_y2 Horntrend3
    ##     1:  0.03419901     0.13594270           2009           2011  0.1586654
    ##     2:  0.01818182     0.15745916           2009           2011  0.1519895
    ##     3:  0.05952381     0.04123930           1994           1996 -0.4629824
    ##     4:  0.16666667     0.28867513           1994           1996  0.4999542
    ##     5: -0.17500000     0.04330127           1995           1997  0.2140151
    ##    ---                                                                    
    ## 10781: -0.13690476     0.09278844           1994           1996 -0.1311599
    ## 10782:  0.13888889     0.24056261           1997           1999  0.2825524
    ## 10783:  0.14318182     0.33460072           1993           1995  0.3741819
    ## 10784:  0.20000000     0.28867513           1993           1995  0.3730799
    ## 10785:  0.12142857     0.30929479           1993           1995  0.4072077
    ##        Horntrend3_se Horntrend3_y1 Horntrend3_y2   Jtutrend5 Jtutrend5_se
    ##     1:    0.07938981          2009          2011  0.01417706   0.02229779
    ##     2:    0.03255360          2009          2011  0.00000000   0.00000000
    ##     3:    0.04546305          1994          1996          NA           NA
    ##     4:    0.86594610          1994          1996          NA           NA
    ##     5:    0.28035861          1995          1997          NA           NA
    ##    ---                                                                   
    ## 10781:    0.03025215          1994          1996          NA           NA
    ## 10782:    0.46657341          1997          1999 -0.06547619   0.04861131
    ## 10783:    0.77766614          1993          1995          NA           NA
    ## 10784:    0.74460278          1993          1995          NA           NA
    ## 10785:    0.48435286          1993          1995          NA           NA
    ##        Jtutrend5_y1 Jtutrend5_y2 Jbetatrend5 Jbetatrend5_se Jbetatrend5_y1
    ##     1:         2007         2011  0.01829960     0.02121080           2007
    ##     2:         2007         2011  0.06491841     0.02735165           2007
    ##     3:           NA           NA          NA             NA             NA
    ##     4:           NA           NA          NA             NA             NA
    ##     5:           NA           NA          NA             NA             NA
    ##    ---                                                                    
    ## 10781:           NA           NA          NA             NA             NA
    ## 10782:         1995         1999 -0.05415973     0.04646827           1995
    ## 10783:           NA           NA          NA             NA             NA
    ## 10784:           NA           NA          NA             NA             NA
    ## 10785:           NA           NA          NA             NA             NA
    ##        Jbetatrend5_y2   Horntrend5 Horntrend5_se Horntrend5_y1 Horntrend5_y2
    ##     1:           2011  0.097426333    0.03629820          2007          2011
    ##     2:           2011 -0.006999522    0.02101022          2007          2011
    ##     3:             NA           NA            NA            NA            NA
    ##     4:             NA           NA            NA            NA            NA
    ##     5:             NA           NA            NA            NA            NA
    ##    ---                                                                      
    ## 10781:             NA           NA            NA            NA            NA
    ## 10782:           1999 -0.074264193    0.11135998          1995          1999
    ## 10783:             NA           NA            NA            NA            NA
    ## 10784:             NA           NA            NA            NA            NA
    ## 10785:             NA           NA            NA            NA            NA
    ##          Jtutrend10 Jtutrend10_se Jtutrend10_y1 Jtutrend10_y2 Jbetatrend10
    ##     1: -0.007288403   0.005205582          2002          2011  0.008630235
    ##     2:  0.001078972   0.003874275          2002          2011  0.004222293
    ##     3:           NA            NA            NA            NA           NA
    ##     4:           NA            NA            NA            NA           NA
    ##     5:           NA            NA            NA            NA           NA
    ##    ---                                                                    
    ## 10781:           NA            NA            NA            NA           NA
    ## 10782:           NA            NA            NA            NA           NA
    ## 10783:           NA            NA            NA            NA           NA
    ## 10784:           NA            NA            NA            NA           NA
    ## 10785:           NA            NA            NA            NA           NA
    ##        Jbetatrend10_se Jbetatrend10_y1 Jbetatrend10_y2 Horntrend10
    ##     1:     0.004356945            2002            2011 0.025979735
    ##     2:     0.006004550            2002            2011 0.001354416
    ##     3:              NA              NA              NA          NA
    ##     4:              NA              NA              NA          NA
    ##     5:              NA              NA              NA          NA
    ##    ---                                                            
    ## 10781:              NA              NA              NA          NA
    ## 10782:              NA              NA              NA          NA
    ## 10783:              NA              NA              NA          NA
    ## 10784:              NA              NA              NA          NA
    ## 10785:              NA              NA              NA          NA
    ##        Horntrend10_se Horntrend10_y1 Horntrend10_y2  Jtutrend20 Jtutrend20_se
    ##     1:    0.009168350           2002           2011 0.001301872   0.001359989
    ##     2:    0.003620485           2002           2011 0.002267440   0.001508569
    ##     3:             NA             NA             NA          NA            NA
    ##     4:             NA             NA             NA          NA            NA
    ##     5:             NA             NA             NA          NA            NA
    ##    ---                                                                       
    ## 10781:             NA             NA             NA          NA            NA
    ## 10782:             NA             NA             NA          NA            NA
    ## 10783:             NA             NA             NA          NA            NA
    ## 10784:             NA             NA             NA          NA            NA
    ## 10785:             NA             NA             NA          NA            NA
    ##        Jtutrend20_y1 Jtutrend20_y2 Jbetatrend20 Jbetatrend20_se Jbetatrend20_y1
    ##     1:          1992          2011  0.004917939     0.001018181            1992
    ##     2:          1992          2011  0.005892922     0.001494195            1992
    ##     3:            NA            NA           NA              NA              NA
    ##     4:            NA            NA           NA              NA              NA
    ##     5:            NA            NA           NA              NA              NA
    ##    ---                                                                         
    ## 10781:            NA            NA           NA              NA              NA
    ## 10782:            NA            NA           NA              NA              NA
    ## 10783:            NA            NA           NA              NA              NA
    ## 10784:            NA            NA           NA              NA              NA
    ## 10785:            NA            NA           NA              NA              NA
    ##        Jbetatrend20_y2  Horntrend20 Horntrend20_se Horntrend20_y1
    ##     1:            2011 0.0026249022    0.001928495           1992
    ##     2:            2011 0.0002789841    0.001070418           1992
    ##     3:              NA           NA             NA             NA
    ##     4:              NA           NA             NA             NA
    ##     5:              NA           NA             NA             NA
    ##    ---                                                           
    ## 10781:              NA           NA             NA             NA
    ## 10782:              NA           NA             NA             NA
    ## 10783:              NA           NA             NA             NA
    ## 10784:              NA           NA             NA             NA
    ## 10785:              NA           NA             NA             NA
    ##        Horntrend20_y2  JtutrendAll JtutrendAll_se JtutrendAll_y1 JtutrendAll_y2
    ##     1:           2011  0.002519604   0.0005271386           1981           2011
    ##     2:           2011  0.001844628   0.0006968803           1981           2011
    ##     3:             NA  0.020507246   0.0237470099           1990           2002
    ##     4:             NA  0.380000000   0.2374868417           1993           1996
    ##     5:             NA -0.150000000   0.1870828693           1994           1997
    ##    ---                                                                         
    ## 10781:             NA  0.099209486   0.0699589739           1993           1999
    ## 10782:             NA -0.029162688   0.0241569824           1992           1999
    ## 10783:             NA  0.069399626   0.0426281017           1992           1999
    ## 10784:             NA  0.031313418   0.0359748410           1992           1999
    ## 10785:             NA  0.200000000   0.2380476143           1992           1995
    ##        JbetatrendAll JbetatrendAll_se JbetatrendAll_y1 JbetatrendAll_y2
    ##     1:   0.003285291      0.000449605             1981             2011
    ##     2:   0.006713102      0.000810327             1981             2011
    ##     3:  -0.004143526      0.008255619             1990             2002
    ##     4:   0.132380952      0.088805375             1993             1996
    ##     5:   0.005303030      0.069984379             1994             1997
    ##    ---                                                                 
    ## 10781:   0.009302654      0.024901900             1993             1999
    ## 10782:   0.009590836      0.013622507             1992             1999
    ## 10783:   0.052660237      0.021213411             1992             1999
    ## 10784:   0.009932716      0.019876409             1992             1999
    ## 10785:   0.108571429      0.121741194             1992             1995
    ##         HorntrendAll HorntrendAll_se HorntrendAll_y1 HorntrendAll_y2 Nspp
    ##     1:  0.0044639799    0.0008210369            1981            2011   83
    ##     2: -0.0004329909    0.0003722616            1981            2011   15
    ##     3: -0.0485170937    0.0217421003            1990            2002   15
    ##     4:  0.2889130783    0.2069851260            1993            1996    7
    ##     5:  0.2483119073    0.1713067801            1994            1997   12
    ##    ---                                                                   
    ## 10781:  0.0754701872    0.0444102580            1993            1999   16
    ## 10782: -0.0069010373    0.0348588399            1992            1999   25
    ## 10783:  0.0386443798    0.0510251318            1992            1999   18
    ## 10784:  0.0059745794    0.0356418670            1992            1999   17
    ## 10785:  0.3345072230    0.1439054702            1992            1995    8

``` r
summary(trends2$Jtutrend3)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -1.000000 -0.090909  0.000000  0.009606  0.125000  1.000000

``` r
summary(trends2$Jbetatrend3)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -1.00000 -0.05263  0.01695  0.01590  0.12500  0.65233

``` r
summary(trends2$Horntrend3)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    ## -1.00000 -0.05304  0.00802  0.01344  0.14626  0.99984      293

## Histograms of temporal change

Standardized slopes have very large and small values

``` r
x <- trends2[, hist(Jtutrend3)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-1.png)<!-- -->

``` r
x <- trends2[, hist(Jbetatrend3)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-2.png)<!-- -->

``` r
x <- trends2[, hist(Horntrend3)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-3.png)<!-- -->

``` r
x <- trends2[, hist(Jtutrend5)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-4.png)<!-- -->

``` r
x <- trends2[, hist(Jbetatrend5)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-5.png)<!-- -->

``` r
x <- trends2[, hist(Horntrend5)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-6.png)<!-- -->

``` r
x <- trends2[, hist(Jtutrend10)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-7.png)<!-- -->

``` r
x <- trends2[, hist(Jbetatrend10)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-8.png)<!-- -->

``` r
x <- trends2[, hist(Horntrend10)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-9.png)<!-- -->

``` r
x <- trends2[, hist(Jtutrend20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-10.png)<!-- -->

``` r
x <- trends2[, hist(Jbetatrend20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-11.png)<!-- -->

``` r
x <- trends2[, hist(Horntrend20)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-12.png)<!-- -->

``` r
x <- trends2[, hist(JtutrendAll)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-13.png)<!-- -->

``` r
x <- trends2[, hist(JbetatrendAll)]
```

![](calc_turnover_files/figure-gfm/histograms%20of%20change-14.png)<!-- -->

``` r
x <- trends2[, hist(HorntrendAll)]
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
pairs(formula = ~ Jtutrend3 + Jbetatrend3 + Horntrend3 +
        Jtutrend5 + Jbetatrend5 + Horntrend5 +
        Jtutrend10 + Jbetatrend10 + Horntrend10 +
        Jtutrend20 + Jbetatrend20 + Horntrend20 +
        JtutrendAll + JbetatrendAll + HorntrendAll, 
      data = trends2, gap = 1/10, cex = 0.2, col = '#00000022', 
      lower.panel = panel.cor,
      upper.panel = panel.smooth)
```

![](calc_turnover_files/figure-gfm/pairs-1.png)<!-- -->

# Change compared to time-series characteristics

## Number of species

    ## NULL

    ## NULL

    ## NULL

    ## NULL

    ## NULL

    ## NULL

    ## NULL

    ## NULL

    ## NULL

    ## NULL

    ## NULL

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20spp-1.png)<!-- -->

    ## NULL

    ## NULL

    ## NULL

![](calc_turnover_files/figure-gfm/change%20vs.%20num%20spp-2.png)<!-- -->

## Average change compared to \#species

``` r
# number of species
ggplot(trends2, aes(Nspp, Jtutrend3, color = 'Jtu trend3')) +
  geom_smooth() +
  geom_smooth(aes(y = Jbetatrend3, color = 'Jbeta trend3')) +
  geom_smooth(aes(y = Horntrend3, color = 'Horn trend3')) +
  geom_smooth(aes(y = Jtutrend5, color = 'Jtu trend5')) +
  geom_smooth(aes(y = Jbetatrend5, color = 'Jbeta trend5')) +
  geom_smooth(aes(y = Horntrend5, color = 'Horn trend5')) +
  geom_smooth(aes(y = Jtutrend10, color = 'Jtu trend10')) +
  geom_smooth(aes(y = Jbetatrend10, color = 'Jbeta trend10')) +
  geom_smooth(aes(y = Horntrend10, color = 'Horn trend10')) +
  geom_smooth(aes(y = Jtutrend20, color = 'Jtu trend20')) +
  geom_smooth(aes(y = Jbetatrend20, color = 'Jbeta trend20')) +
  geom_smooth(aes(y = Horntrend20, color = 'Horn trend20')) +
  geom_smooth(aes(y = JtutrendAll, color = 'Jtu trendAll')) +
  geom_smooth(aes(y = JbetatrendAll, color = 'Jbeta trendAll')) +
  geom_smooth(aes(y = HorntrendAll, color = 'Horn trendAll')) +
  scale_x_log10() +
  labs(y = 'Slope') +
  geom_abline(intercept = 0, slope = 0)
```

![](calc_turnover_files/figure-gfm/ave%20change%20vs.%20num%20spp-1.png)<!-- -->

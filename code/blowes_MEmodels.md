Replicate and extend analyses from Blowes et al. 2019
================

# Load data

    ## here() starts at /local/home/malinp/climate_community_crossrealm

# Fit biome-taxon model

## On slopes from year 1 for all years

Finds a significantly positive slope. Not as high as in 2019 paper,
perhaps because year 1 self-comparison not included.

``` r
# published model
modJtuAlly1 <- glmmTMB(disstrend ~ 1 + (1|Biome/taxa_mod/STUDY_ID),
                data = trends[measure == 'Jtu' & duration_group == 'Ally1', ])

summary(modJtuAlly1)
```

    ##  Family: gaussian  ( identity )
    ## Formula:          disstrend ~ 1 + (1 | Biome/taxa_mod/STUDY_ID)
    ## Data: trends[measure == "Jtu" & duration_group == "Ally1", ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -55925.3 -55882.5  27967.7 -55935.3    38727 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                    Name        Variance  Std.Dev. 
    ##  STUDY_ID:(taxa_mod:Biome) (Intercept) 9.198e-05 9.591e-03
    ##  taxa_mod:Biome            (Intercept) 8.289e-05 9.104e-03
    ##  Biome                     (Intercept) 1.545e-10 1.243e-05
    ##  Residual                              1.377e-02 1.174e-01
    ## Number of obs: 38732, groups:  
    ## STUDY_ID:(taxa_mod:Biome), 399; taxa_mod:Biome, 164; Biome, 59
    ## 
    ## Dispersion estimate for gaussian family (sigma^2): 0.0138 
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.008564   0.001957   4.375 1.22e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Slopes by realm/biome

Essentially no difference among the realms. If anything, terrestrial is
fastest.

``` r
# extract
res <- ranef(modJtuAlly1)$cond$Biome
colnames(res) <- "re"
res$Biome = rownames(res)
realms <- trends[, .(REALM = unique(REALM)), by = Biome] # realm for each biome
res <- as.data.table(merge(res, realms, by = 'Biome'))
res[, intercept := re + as.numeric(fixef(modJtuAlly1))[1]] # add the fixed intercept to the REs

# by realm
res[, .(mean = mean(intercept), se = sd(intercept)/sqrt(.N), meanRE = mean(re), seRE = sd(re)/sqrt(.N)), by= REALM]
```

    ##          REALM        mean           se        meanRE         seRE
    ## 1:      Marine 0.008563691 1.637036e-09 -1.333878e-10 1.637036e-09
    ## 2: Terrestrial 0.008563693 2.096735e-09  1.531757e-09 2.096735e-09
    ## 3:  Freshwater 0.008563689 2.318680e-09 -2.428471e-09 2.318680e-09

## On slopes from all pairs for all years

Finds a significantly positive slope. Not as high as in 2019 paper,
perhaps because year 1 self-comparison not included.

``` r
# published model
modJtuAll <- glmmTMB(disstrend ~ 1 + (1|Biome/taxa_mod/STUDY_ID),
                data = trends[measure == 'Jtu' & duration_group == 'All', ])

summary(modJtuAll)
```

    ##  Family: gaussian  ( identity )
    ## Formula:          disstrend ~ 1 + (1 | Biome/taxa_mod/STUDY_ID)
    ## Data: trends[measure == "Jtu" & duration_group == "All", ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -91132.0 -91089.1  45571.0 -91142.0    38727 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                    Name        Variance  Std.Dev. 
    ##  STUDY_ID:(taxa_mod:Biome) (Intercept) 9.091e-05 9.535e-03
    ##  taxa_mod:Biome            (Intercept) 1.200e-04 1.095e-02
    ##  Biome                     (Intercept) 6.409e-10 2.532e-05
    ##  Residual                              5.539e-03 7.443e-02
    ## Number of obs: 38732, groups:  
    ## STUDY_ID:(taxa_mod:Biome), 399; taxa_mod:Biome, 164; Biome, 59
    ## 
    ## Dispersion estimate for gaussian family (sigma^2): 0.00554 
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.011372   0.001786   6.367 1.93e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Slopes by realm/biome

Essentially no difference among the realms. If anything, terrestrial is
fastest.

``` r
# extract
res <- ranef(modJtuAll)$cond$Biome
colnames(res) <- "re"
res$Biome = rownames(res)
realms <- trends[, .(REALM = unique(REALM)), by = Biome] # realm for each biome
res <- as.data.table(merge(res, realms, by = 'Biome'))
res[, intercept := re + as.numeric(fixef(modJtuAll))[1]] # add the fixed intercept to the REs

# by realm
res[, .(mean = mean(intercept), se = sd(intercept)/sqrt(.N), meanRE = mean(re), seRE = sd(re)/sqrt(.N)), by= REALM]
```

    ##          REALM       mean           se        meanRE         seRE
    ## 1:      Marine 0.01137218 8.092646e-09 -1.828668e-09 8.092646e-09
    ## 2: Terrestrial 0.01137219 1.184445e-08  5.588601e-09 1.184445e-08
    ## 3:  Freshwater 0.01137218 4.339448e-09  8.495960e-11 4.339448e-09

## On slopes from 5 consecutive years

Finds a significantly positive slope. Not as high as in 2019 paper,
perhaps because year 1 self-comparison not included.

``` r
# published model
modJtu5 <- glmmTMB(disstrend ~ 1 + (1|Biome/taxa_mod/STUDY_ID),
                data = trends[measure == 'Jtu' & duration_group == '5', ])

summary(modJtu5)
```

    ##  Family: gaussian  ( identity )
    ## Formula:          disstrend ~ 1 + (1 | Biome/taxa_mod/STUDY_ID)
    ## Data: trends[measure == "Jtu" & duration_group == "5", ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -14038.5 -14005.3   7024.2 -14048.5     5625 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                    Name        Variance  Std.Dev. 
    ##  STUDY_ID:(taxa_mod:Biome) (Intercept) 1.190e-11 3.450e-06
    ##  taxa_mod:Biome            (Intercept) 3.426e-05 5.853e-03
    ##  Biome                     (Intercept) 2.144e-13 4.631e-07
    ##  Residual                              4.811e-03 6.936e-02
    ## Number of obs: 5630, groups:  
    ## STUDY_ID:(taxa_mod:Biome), 208; taxa_mod:Biome, 108; Biome, 43
    ## 
    ## Dispersion estimate for gaussian family (sigma^2): 0.00481 
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.009790   0.001587   6.167 6.97e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Slopes by realm/biome

Essentially no difference among the realms. If anything, terrestrial is
fastest.

``` r
# extract
res <- ranef(modJtu5)$cond$Biome
colnames(res) <- "re"
res$Biome = rownames(res)
realms <- trends[, .(REALM = unique(REALM)), by = Biome] # realm for each biome
res <- as.data.table(merge(res, realms, by = 'Biome'))
res[, intercept := re + as.numeric(fixef(modJtu5))[1]] # add the fixed intercept to the REs

# by realm
res[, .(mean = mean(intercept), se = sd(intercept)/sqrt(.N), meanRE = mean(re), seRE = sd(re)/sqrt(.N)), by= REALM]
```

    ##          REALM        mean           se        meanRE         seRE
    ## 1:      Marine 0.009789614 5.081636e-12  1.000758e-12 5.081636e-12
    ## 2: Terrestrial 0.009789614 3.022008e-12  1.072424e-12 3.022008e-12
    ## 3:  Freshwater 0.009789614 2.146757e-12 -6.314673e-12 2.146757e-12

## On slopes from 10 consecutive years

Finds a significantly positive slope. Not as high as in 2019 paper,
perhaps because year 1 self-comparison not included.

``` r
# published model
modJtu10 <- glmmTMB(disstrend ~ 1 + (1|Biome/taxa_mod/STUDY_ID),
                data = trends[measure == 'Jtu' & duration_group == '10', ])

summary(modJtu10)
```

    ##  Family: gaussian  ( identity )
    ## Formula:          disstrend ~ 1 + (1 | Biome/taxa_mod/STUDY_ID)
    ## Data: trends[measure == "Jtu" & duration_group == "10", ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -11203.7 -11175.3   5606.9 -11213.7     2184 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                    Name        Variance  Std.Dev.
    ##  STUDY_ID:(taxa_mod:Biome) (Intercept) 4.888e-06 0.002211
    ##  taxa_mod:Biome            (Intercept) 1.502e-05 0.003875
    ##  Biome                     (Intercept) 4.018e-06 0.002005
    ##  Residual                              3.413e-04 0.018474
    ## Number of obs: 2189, groups:  
    ## STUDY_ID:(taxa_mod:Biome), 123; taxa_mod:Biome, 66; Biome, 30
    ## 
    ## Dispersion estimate for gaussian family (sigma^2): 0.000341 
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.005521   0.001125    4.91 9.12e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Slopes by realm/biome

Essentially no difference among the realms. If anything, marine is
fastest.

``` r
# extract
res <- ranef(modJtu10)$cond$Biome
colnames(res) <- "re"
res$Biome = rownames(res)
realms <- trends[, .(REALM = unique(REALM)), by = Biome] # realm for each biome
res <- as.data.table(merge(res, realms, by = 'Biome'))
res[, intercept := re + as.numeric(fixef(modJtu10))[1]] # add the fixed intercept to the REs

# by realm
res[, .(mean = mean(intercept), se = sd(intercept)/sqrt(.N), meanRE = mean(re), seRE = sd(re)/sqrt(.N)), by= REALM]
```

    ##          REALM        mean           se        meanRE         seRE
    ## 1:      Marine 0.005530257 2.237110e-04  9.085209e-06 2.237110e-04
    ## 2: Terrestrial 0.005518891 1.495499e-04 -2.280771e-06 1.495499e-04
    ## 3:  Freshwater 0.005490536 5.365838e-05 -3.063605e-05 5.365838e-05

## On slopes from 20 consecutive years

Finds a significantly positive slope. Not as high as in 2019 paper,
perhaps because year 1 self-comparison not included.

``` r
# published model
modJtu20 <- glmmTMB(disstrend ~ 1 + (1|Biome/taxa_mod/STUDY_ID),
                data = trends[measure == 'Jtu' & duration_group == '20', ])

summary(modJtu20)
```

    ##  Family: gaussian  ( identity )
    ## Formula:          disstrend ~ 1 + (1 | Biome/taxa_mod/STUDY_ID)
    ## Data: trends[measure == "Jtu" & duration_group == "20", ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  -4261.4  -4240.4   2135.7  -4271.4      490 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                    Name        Variance  Std.Dev. 
    ##  STUDY_ID:(taxa_mod:Biome) (Intercept) 4.171e-05 6.458e-03
    ##  taxa_mod:Biome            (Intercept) 1.128e-12 1.062e-06
    ##  Biome                     (Intercept) 2.385e-13 4.884e-07
    ##  Residual                              8.204e-06 2.864e-03
    ## Number of obs: 495, groups:  
    ## STUDY_ID:(taxa_mod:Biome), 48; taxa_mod:Biome, 29; Biome, 18
    ## 
    ## Dispersion estimate for gaussian family (sigma^2): 8.2e-06 
    ## 
    ## Conditional model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) 0.0042356  0.0009956   4.254  2.1e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Slopes by realm/biome

Essentially no difference among the realms. If anything, freshwater is
fastest.

``` r
# extract
res <- ranef(modJtu20)$cond$Biome
colnames(res) <- "re"
res$Biome = rownames(res)
realms <- trends[, .(REALM = unique(REALM)), by = Biome] # realm for each biome
res <- as.data.table(merge(res, realms, by = 'Biome'))
res[, intercept := re + as.numeric(fixef(modJtu20))[1]] # add the fixed intercept to the REs

# by realm
res[, .(mean = mean(intercept), se = sd(intercept)/sqrt(.N), meanRE = mean(re), seRE = sd(re)/sqrt(.N)), by= REALM]
```

    ##          REALM        mean           se        meanRE         seRE
    ## 1: Terrestrial 0.004235562 1.492518e-11  3.795197e-12 1.492518e-11
    ## 2:      Marine 0.004235562 2.543591e-11 -1.949495e-11 2.543591e-11
    ## 3:  Freshwater 0.004235562 1.108973e-11  1.000695e-11 1.108973e-11

Replicate and extend analyses from Antao et al. 2020
================

# Load data

    ## here() starts at /local/home/malinp/climate_community_crossrealm

# Fit model

``` r
mod0 <- glmmTMB(slope ~ 0 + new_sTempYear:REALM + TempGAMCoef:REALM + TempGAMCoef:new_sTempYear:REALM
                +(0 + TempGAMCoef|REALM/taxa_mod1) +(1|REALM/taxa_mod1/STUDY_ID), 
                disp = ~std.error*REALM,
                data = dat[model_id == 'logS_lm',])


#mod0Mar <- glmmTMB(slope ~ 0 + TempGAMCoef+new_sTempYear +TempGAMCoef:new_sTempYear +(0 + TempGAMCoef|taxa_mod1) +(1|taxa_mod1/STUDY_ID), 
#                disp = ~std.error,
#                data = datRichMar)

#mod0Terr <- glmmTMB(slope ~ 0 + TempGAMCoef+new_sTempYear +TempGAMCoef:new_sTempYear +(0 + TempGAMCoef|taxa_mod1) +(1|taxa_mod1/STUDY_ID), 
 #               disp = ~std.error,
 #               data = datRichTerr)

summary(mod0)
```

    ##  Family: gaussian  ( identity )
    ## Formula:          
    ## slope ~ 0 + new_sTempYear:REALM + TempGAMCoef:REALM + TempGAMCoef:new_sTempYear:REALM +  
    ##     (0 + TempGAMCoef | REALM/taxa_mod1) + (1 | REALM/taxa_mod1/STUDY_ID)
    ## Dispersion:             ~std.error * REALM
    ## Data: dat[model_id == "logS_lm", ]
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## -80528.1 -80408.4  40279.0 -80558.1    21485 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                     Name        Variance  Std.Dev. 
    ##  taxa_mod1.REALM            TempGAMCoef 4.220e-03 0.0649587
    ##  REALM                      TempGAMCoef 4.028e-07 0.0006347
    ##  STUDY_ID..taxa_mod1.REALM. (Intercept) 3.035e-04 0.0174222
    ##  taxa_mod1.REALM.1          (Intercept) 2.012e-05 0.0044856
    ##  REALM.1                    (Intercept) 1.188e-05 0.0034464
    ##  Residual                                      NA        NA
    ## Number of obs: 21500, groups:  
    ## taxa_mod1:REALM, 12; REALM, 2; STUDY_ID:(taxa_mod1:REALM), 156
    ## 
    ## Conditional model:
    ##                                              Estimate Std. Error z value
    ## new_sTempYear:REALMMarine                   0.0029022  0.0005907   4.913
    ## new_sTempYear:REALMTerrestrial              0.0011181  0.0007450   1.501
    ## REALMMarine:TempGAMCoef                     0.0938505  0.0331252   2.833
    ## REALMTerrestrial:TempGAMCoef               -0.0091852  0.0398106  -0.231
    ## new_sTempYear:REALMMarine:TempGAMCoef       0.0848609  0.0117952   7.195
    ## new_sTempYear:REALMTerrestrial:TempGAMCoef -0.0140901  0.0132704  -1.062
    ##                                            Pr(>|z|)    
    ## new_sTempYear:REALMMarine                  8.96e-07 ***
    ## new_sTempYear:REALMTerrestrial              0.13342    
    ## REALMMarine:TempGAMCoef                     0.00461 ** 
    ## REALMTerrestrial:TempGAMCoef                0.81753    
    ## new_sTempYear:REALMMarine:TempGAMCoef      6.27e-13 ***
    ## new_sTempYear:REALMTerrestrial:TempGAMCoef  0.28834    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Dispersion model:
    ##                            Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                -7.87426    0.01876  -419.8   <2e-16 ***
    ## std.error                  39.41112    0.47960    82.2   <2e-16 ***
    ## REALMTerrestrial           -0.60158    0.07103    -8.5   <2e-16 ***
    ## std.error:REALMTerrestrial 27.14306    2.46121    11.0   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Compare datasets

``` r
# merge
comb <- merge(dat, trends[duration_group == 'Ally1', ], by = c('STUDY_ID', 'rarefyID', 'REALM', 'rarefyID_y'))

# compare
comb[, plot(startYear, year1)]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

    ## NULL

``` r
comb[, plot(endYear, year2)]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

    ## NULL

``` r
comb[, plot(TempGAMCoef, temptrend)]; abline(0,1)
```

    ## NULL

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
comb[, plot(new_sTempYear, tempave.sc, col = c('green', 'blue')[(REALM=='Marine')+1])]; abline(0,1)
```

    ## NULL

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

``` r
comb[model_id == 'logS_lm' & measure == 'Jtu', plot(slope, disstrend, xlab = 'logS slope', ylab = 'Jtu slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'Gains_lm' & measure == 'Jtu', plot(slope, disstrend, xlab = 'Gains slope', ylab = 'Jtu slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-6.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'Losses_lm' & measure == 'Jtu', plot(slope, disstrend, xlab = 'Losses slope', ylab = 'Jtu slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-7.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'logS_lm' & measure == 'Jbeta', plot(slope, disstrend, xlab = 'logS slope', ylab = 'Jbeta slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-8.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'Gains_lm' & measure == 'Jbeta', plot(slope, disstrend, xlab = 'Gains slope', ylab = 'Jbeta slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-9.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'Losses_lm' & measure == 'Jbeta', plot(slope, disstrend, xlab = 'Losses slope', ylab = 'Jbeta slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-10.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'logS_lm' & measure == 'Horn', plot(slope, disstrend, xlab = 'logS slope', ylab = 'Horn slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-11.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'logN_lm' & measure == 'Horn', plot(slope, disstrend, xlab = 'logN slope', ylab = 'Horn slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-12.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'Gains_lm' & measure == 'Horn', plot(slope, disstrend, xlab = 'Gains slope', ylab = 'Horn slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-13.png)<!-- -->

    ## NULL

``` r
comb[model_id == 'Losses_lm' & measure == 'Horn', plot(slope, disstrend, xlab = 'Losses slope', ylab = 'Horn slope')]
```

![](antao_MEmodels_files/figure-gfm/unnamed-chunk-2-14.png)<!-- -->

    ## NULL

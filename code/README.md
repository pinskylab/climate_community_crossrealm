# code

The scripts in this analysis were run in four computer environments:
- A MacBook Pro 16 GB RAM with R version 4.0.3 and RStudio 1.3.1093
- A scientific workstation ('Annotate') with 36 cores (Intel Xeon CPU E5-2697 v4 @ 2.30GHz) running CentOS 7.6.1810 with R version 3.5.2 and RStudio 1.3.1093
- A scientific workstation ('Annotate2') with 56 cores (Intel Xeon CPU w9-3495X @ 4.8GHz) running Red Hat 8 with R version 4.3.1 and RStudio 2023.06.2 (fitting ordered beta GLMMs)
- A compute cluster with 64-bit x86 processors running Rocky Linux 9.1 with SLURM job management and R version 4.2.2 (Step 2 below).

The libraries included:
- bbmle 1.0.23.1
- beanplot 1.2
- betapart 1.6
- data.table 1.11.8
- dggridR
- dplyr 1.1.0
- DHARMa 0.3.3.0
- ggplot2 3.3.5
- ggpubr 0.4.0
- glmmTMB 1.0.2.1
- grid 3.5.2
- gridExtra 2.3
- here 0.3.3
- iNext
- lattice 0.20-45
- lazyeval 0.2.2
- lme4 1.1-27.1
- maps 3.3.0
- mgcv 1.8.26
- ncdf4 1.17
- nlme 3.1-137
- performance 0.7.0
- permut 0.9-7
- purrr 1.0.1
- raster 2.7-15
- RColorBrewer 1.1-2
- rcompanion 2.3.27
- reshape2 1.4.4
- rstudioapi 0.13
- scales 1.0.0
- sdmpredictors 0.2.9
- taxize 0.9.99
- tibble 3.2.0
- tidyr 1.3.0
- tidyverse 1.2.1
- vegan 2.6-4

In addition, we use custom functions in `util.R`.

We expect install time will take a couple hours on a standard desktop computer, with the slowest part being package installation.

The code is organized in six main steps. Run time for each step should be less than an hour, except for steps 2 and 4. Each model in step 4 takes up to 48 hrs to fit on a single core of a scientific workstation.

## 1. Examine alternative statistical approaches
- `duration_sim.R`: simulate data and fit models. Writes `output/simulated_ts.csv.gz`
- `duration_sim.Rmd`: examine simulation results
- `turnover_richness_sim.Rmd`: examine the effects of species richness on turnover rates

## 2. Prep dissimilarity data
- `01_Study_to_Grid.R`: grids the observations to 96 km2 hexagons
- `02_all_pairs_cluster.R`: calculate rarefied dissimilarities among years within each time series. Runs in SLURM.
- `03_collate_resamps_cluster.R`: script to combine rarefied resamples and calculate median dissimilarity. Runs in SLURM. Writes `bt-rarefy-collate.Rdata`

## 3. Prep turnover data
- `assemble_microclimates.R`: extract microclimate variability for each time series. Writes `output/microclimates.csv.gz`
- `assemble_temp.Rmd`: calculate temperature averages and trends. Writes `output/temperature_byrarefyID.csv.gz`
- `extract_human_impact.R`: extract human impact values for each time series. Writes `output/humanimpact_by_rarefyID.csv.gz`
- `extract_richness.R`: extract species richness data for each time series. Writes `output/richness_by_rarefyID.csv.gz`
- `assemble_turnover_covariates.Rmd`: put dissimilarity data together with covariates. Writes `output/turnover_w_covariates_scaling.csv`
- `calc_turnover.R`: calculate temporal turnover as the slope of dissimilarity vs. year. Writes `output/slope.csv.gz`
- `sample_global_temp.Rmd`: make a representative sample of global temperature trends. Writes `output/temperature_trends_sampled.csv.gz`

## 4. Fit models
- `fit_turnover_GLMM.R`: fit one mixed effects model at a time (specified as an argument). Writes `temp/mod*.rds`
- `fit_turnover_GLMM.sh`: bash script to initiate multiple threads, one to fit each model specified as an argument. Calls `fit_turnover_GLMM.R`
- `fit_turnover_GLMMaltlink.R`: fit models with Gaussian errors and linear link. Writes `temp/mod*_lin.rds`
- `fit_turnover_GLMMaltlink.sh`: bash script to initiate multiple threads, one to fit each model specified as an argument. Calls `fit_turnover_GLMMaltlink.R`
- `fit_turnover_GLMMlong.R`: fit mixed effects models to time series >= 7 years long. Writes `temp/mod*LongJtu.rds`

## 5. Make predictions from the models
- `pred_GLMM.R`: make predictions from the models for plotting, including turnover rates and sensitivity of turnover to temperature change. Takes the model name as an argument. Writes `temp/preds_modOBsdTMERtsRealmtsigninitAllJtu.rds` and `preds_modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds` (dissimilarity), `temp/slopes_modOBsdTMERtsRealmtsigninitAllJtu.rds` and `temp/slopes_modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds` (turnover rate), and `temp/sensitivity_modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds` (sensitivity of turnover rates to temperature change)
- `pred_GLMM.sh`: Shell script to spawn multiple instances of `pred_GLMM.R`
- `pred_GLMMcov.R`: make predictions for the microcolimate and human impact covariate models. Writes `temp/preds_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu_marterr.rds` (dissimilarity), `temp/slopes_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu_marterr.rds` (turnover rates), and `temp/sensitivity_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu_marterr.rds` (sensitivity)
- `pred_GLMMaltlink.R`: for predictions from models with Gaussian errors
- `pred_GLMMaltlink.sh`: shell script to accompany R script
- `pred_GLMMaltlinkcov.R`: for predictions from models with covariates and Gaussian errors
- `pred_GLMMLong.R`: make predictions from the models fit to longer time series. Writes `temp/preds_*Long*.rds`, `temp/slopes_*Long*.rds` and `temp/sensitivity_*Long*.rds`
- `make_residuals.R`: make DHARMa standardized residuals from models for checking model assumptions

## 6. Reshuffling and downsampling sensitivity analyses
- `fit_pred_turnover_GLMM_downsamp.R`: fit models to downsampled data and predict from them. Set up to run with arguments to specify the model and range of random seeds. Writes `temp/'*_boot*.rds`, `preds_*_boot*.rds`, `slopes_*_boot*.rds`, and `sensitivity_*_boot*.rds`.
- `fit_pred_turnover_GLMM_downsamp.sh`: bash script to spawn multiple instances of `fit_pred_turnover_GLMM_downsamp.R`, each with a range of random seeds. Takes arguments.
- `read_downsamp.R`: Reads in and repackages the output from `fit_pred_turnover_GLMM_downsamp.R`. Writes `output/downsampTchange.csv.gz` and `output/downsampTave.csv.gz`.
- `fit_turnover_GLMM_reshuffle.R`: fit the Tchange x Year x Realm model with reshuffled Tchange. Set up to be run with arguments to specify the range of random seeds.
- `fit_turnover_GLMM_reshuffle.sh`: bash script to spawn multiple instances of `fit_turnover_GLMM_reshuffle.R`, each with a range of random seeds.
- `read_reshuffle.R`: Reads in and repackages the output from `fit_turnover_GLMM_reshuffle.R`. Writes `output/coeftable_reshuffle_modOBsdTMERtsRealmtsigninitAllJtu.csv`.

## 7. Make outputs
- `figures_for_paper.R`: figures and tables and stats for paper
- `figures_with_altlink.R`: make figures and tables with alternative link function
- `examine_altlink_mods.R`: make alternative figures and tables as a sensitivity analysis with models that have linear link and Gaussian errors

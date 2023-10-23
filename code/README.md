# code

The scripts in this analysis were run in three computer environments:
- A MacBook Pro 16 GB RAM with R version 4.0.3 and RStudio 1.3.1093
- A scientific workstation ('Annotate') with 72 CPUs (Intel Xeon CPU E5-2697 v4 @ 2.30GHz) running CentOS 7.6.1810 with R version 3.5.2 and RStudio 1.3.1093
- A high performance computing environment with SLURM job management (Step 2 below).

The libraries included:
- bbmle 1.0.23.1
- beanplot 1.2
- betapart
- data.table 1.11.8
- dggridR
- DHARMa 0.3.3.0
- ggplot2 3.3.5
- ggpubr 0.4.0
- glmmTMB 1.0.2.1
- grid 3.5.2
- gridExtra 2.3
- here 0.3.3
- iNext
- lazyeval
- lme4 1.1-27.1
- maps 3.3.0
- mgcv 1.8.26
- ncdf4 1.17
- nlme 3.1-137
- performance 0.7.0
- purrr
- raster 2.7-15
- RColorBrewer 1.1-2
- rcompanion 2.3.27
- reshape2 1.4.4
- rstudioapi 0.13
- scales 1.0.0
- sdmpredictors 0.2.9
- taxize 0.9.99
- tibble
- tidyr
- tidyverse 1.2.1
- vegan 2.6-4

In addition, we use custom functions here:

- `util.R`: some utility functions used in multiple scripts

We expect install time will take a couple hours on a standard desktop computer.

The code is organized in six main steps. Run time for each step should be less than an hour, except for steps 2 and 4. Each model in step 4 takes up to 48 hrs to fit,

## 1. Examine alternative statistical approaches
- `duration_sim.R`: simulate data and fit models. Writes `output/simulated_ts.csv.gz`
- `duration_sim.Rmd`: examine simulation results
- `turnover_richness_sim.Rmd`: examine the effects of species richness on turnover rates

## 2. Prep dissimilarity data
- `01_Study_to_Grid.R`: grids the observations to 96 km2 hexagons
- `02_all_pairs_cluster.R`: calculate rarefied dissimilarities among years within each time series. Runs in SLURM.
- `03_collate_resamps_cluster.R`: script to combine rarefied resamples and calculate median dissimilarity. Runs in SLURM.

## 3. Prep turnover data
- `assemble_microclimates.R`: extract microclimate variability for each time series. Writes `output/microclimates.csv.gz`
- `assemble_temp.Rmd`: calculate temperature averages and trends. Writes `output/temperature_byrarefyID.csv.gz`
- `extract_human_impacts.R`: extract human impact values for each time series. Writes `output/humanimpact_by_rarefyID.csv.gz`
- `extract_richness.R`: extract species richness data for each time series. Writes `output/richness_by_rarefyID.csv.gz`
- `assemble_turnover_covariates.Rmd`: put dissimilarity data together with covariates. Writes `output/turnover_w_covariates_scaling.csv`
- `calc_turnover.R`: calculate temporal turnover as the slope of dissimilarity vs. year. Writes `output/slope.csv.gz`
- `sample_global_temp.Rmd`: make a representative sample of global temperature trends. Writes `output/temperature_trends_sampled.csv.gz`

## 4. Fit models
- `turnover_GLMM_fit.R`: fit one mixed effects model (specified as an argument). Writes `temp/mod*.rds`
- `turnover_GLMM_fit.sh`: bash script to initiate multiple threads, one to fit each model specified as an argument. Calls `turnover_GLMM_fit.R`
- `turnover_GLMMgainloss_fit.R`: fit mixed effects models with the relative proportions of species gains and losses as a covariate (a sensitivity analysis). Writes `temp/mod*InitGainLossAll*.rds`
- `turnover_GLMMlong_fit.R`: fit mixed effects models to time series >= 7 years long. Writes `temp/mod*LongJtu.rds`

## 5. Make predictions from the models
- `pred_GLMMmodrawXAllJtu.R`: make predictions from the models for plotting, including turnover rates and sensitivity of turnover to temperature change. Writes `temp/preds_modsdTRealmtsigninitAllJtu.rds` and `preds_rawTsdTTRealmtsigninit.rds` (dissimilarity), `temp/slopes_modsdTRealmtsigninitAllJtu.rds` and `slopes_rawTsdTTRealmtsigninit.rds` (turnover rate), and `sensitivity_rawTsdTTRealmtsigninit.rds` (sensitivity of turnover rates to temperature change)
- `pred_modrawXAllJtu.sh`: Shell script to spawn multiple instances of `pred_GLMMmodrawXAllJtu.R`
- `pred_GLMMmodrawCovariate.R`: make predictions for the covariate models. Writes `temp/preds_rawTsdTTRealmtsignCovariateInit.rds` (dissimilarity), `temp/slopes_rawTsdTTRealmtsignCovariateInit.rds` (turnover rates), and `temp/sensitivity_rawTsdTTRealmtsignCovariateInit.rds` (sensitivity)
- `pred_GLMMmodrawXLongJtu.R`: make predictions from the models fit to longer time series. Writes `temp/preds_*Long*.rds`, `temp/slopes_*Long*.rds` and `temp/sensitivity_*Long*.rds`

## 6. Make outputs for communication
- `figures_for_paper.R`: figures and tables and stats for paper

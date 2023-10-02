# code

The scripts in this analysis were run in two computer environments:
- A MacBook Pro 16 GB RAM with R version 4.0.3 and RStudio 1.3.1093
- A scientific workstation with 72 CPUs (Intel Xeon CPU E5-2697 v4 @ 2.30GHz) running CentOS 7.6.1810 with R version 3.5.2 and RStudio 1.3.1093

The libraries included:
- tidyverse 1.2.1
- raster 2.7-15
- sdmpredictors 0.2.9
- ggplot2 3.3.5
- maps 3.3.0
- data.table 1.11.8
- ncdf4 1.17
- here 0.3.3
- mgcv 1.8.26
- beanplot 1.2
- gridExtra 2.3
- grid 3.5.2
- scales 1.0.0
- glmmTMB 1.0.2.1
- ggpubr 0.4.0
- nlme 3.1-137
- RColorBrewer 1.1-2
- rcompanion 2.3.27
- taxize 0.9.99
- rstudioapi 0.13
- lavaan 0.6-9
- lme4 1.1-27.1
- bbmle 1.0.23.1
- DHARMa 0.3.3.0
- performance 0.7.0

In addition, we use custom functions here:

- util.R: some utility functions used in multiple scripts


The code is organized in five main steps:
## 1. Examine alternative statistical approaches
- duration_sim.R: simulate data and fit models
- duration_sim.Rmd: examine results

## 2. Prep data
- assemble_microclimates.R: extract microclimate variability for each time series
- assemble_temp.Rmd: calculate temperature averages and trends
- extract_human_impacts.R: extract human impact values for each time series
- extract_richness.R: extract species richness data for each time series
- assemble_turnover_covariates.Rmd: put dissimilarity data together with covariates
- calc_turnover.R: calculate temporal turnover as the slope of dissimilarity vs. year
- sample_global_temp.Rmd: make a representative sample of global temperature trends

## 3. Fit models
- turnover_GLMM_fit.R: fit one mixed effects model (specified as an argument)
- turnover_GLMM_fit.sh: bash script to initiate multiple threads, one to fit each model specified as an argument. Calls turnover_GLMM_fit.R.
- turnover_GLMMgainloss_fit.R: fit mixed effects models with the relative proportions of species gains and losses as a covariate (a sensitivity analysis)

## 4. Make predictions from the models
- pred_GLMMmodrawXAllJtu.R: make predictions from the models for plotting, including turnover rates and sensitivity of turnover to temperature change
- pred_modrawXAllJtu.sh: Shell script to spawn multiple instances of pred_GLMMmodrawXAllJtu.R
- pred_GLMMmodrawCovariate.R: make predictions for the covariate models

## 5. Make outputs for communication
- figures_for_paper.R: figures and tables and stats for paper

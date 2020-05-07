
###Laura Antao
###06-05-2020

##first script to send to Malin using brms models to analyse turnover~Year trends from the BioTIME gridded data
##note the models can take several days to run

# DISCLAIMER: code is under development

##======================================================================
## The biodiversity data used originates from pre-processing steps from Blowes et al (https://github.com/sChange-workshop/BioGeo-BioDiv-Change)

## Code below follows from script in https://github.com/pinskylab/climate_community_crossrealm/blob/master/code/turnover_vs_temperature.Rmd


## Input is object with turnover change estimates, temperature change estimates, and metadata for the time series
## Will also add "microclimate" estimates from each grid cells used in the biodiversity data processing (accompanying script)


##======================================================================

library(tidyverse)
library(brms)  ##v 2.6.0
# library(ggpubr)
# library(ggthemes)
# library(viridis)
# library(ggExtra)
# library(cowplot)

# setwd("~/Project_collabs/Turnover_temperature_Pinsky/climate_community_crossrealm-master")

## Load data and run code as in "turnover_vs_temperature.Rmd"

##======================================================================

## then replicate the models at the end of the script (move from lmer to brms)

# fit models
# mods <- vector('list', 0)
# mods[[1]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr.sc + REALM + (1|STUDY_ID))]#, weights = 1/(Jtutrend_se_nz^2))]
# mods[[2]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc + REALM + (1|STUDY_ID))] 
# mods[[3]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*REALM + (1|STUDY_ID))]
# mods[[4]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*seas_comb.sc + (1|STUDY_ID))]
# mods[[5]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*tempave_comb.sc + (1|STUDY_ID))]
# mods[[6]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*npp.sc + (1|STUDY_ID))]
# mods[[7]] <- trends[i, lmer(Jtutrend ~ temptrend_comb_allyr_abs.sc*mass_geomean.sc + (1|STUDY_ID))]


###data filtering steps (suggestions!!)

##remove cells (=rarefyID) with <5 years sampled
##this should also allow to deal with the majority (all??) of the "perfect fits"
names(trends)


# trends <- trends %>%
#   filter(nyrBT >5)  ##would have 18114 records



##overall model strcuture examples to account for the data structure??
##              Jtutrend ~ temptrend_comb_allyr_abs.sc * REALM + (1 | STUDY_ID)
##              Jtutrend ~ temptrend_comb_allyr_abs.sc * REALM + (1 | Taxa/STUDY_ID)   ##because we have multiple per taxa group?
##              Jtutrend ~ temptrend_comb_allyr_abs.sc * REALM + (temptrend_comb_allyr_abs.sc | Taxa/STUDY_ID)  ##varying responses per taxa?

##once decide on random effects, models with different co-variates numbered ~ above
##              [1] Jtutrend ~ temptrend_comb_allyr_abs.sc * REALM + (1 | STUDY_ID)
##              [2] Jtutrend ~ temptrend_comb_allyr_abs.sc * seas_comb.sc + (1|STUDY_ID)
##              [3] Jtutrend ~ temptrend_comb_allyr_abs.sc * tempave_comb.sc + (1|STUDY_ID)
##              [4] Jtutrend ~ temptrend_comb_allyr_abs.sc * npp.sc + (1|STUDY_ID)
##              [5] Jtutrend ~ temptrend_comb_allyr_abs.sc * mass_geomean.sc + (1|STUDY_ID)



###create formula for brms models
formula1 <- bf(Jtutrend | se(Jtutrend_se, sigma = TRUE) ~ temptrend_comb_allyr_abs.sc * REALM + (1 | STUDY_ID))

##the syntax using se() allows to specify standard errors of the observations, thus allowing to perform meta-analysis
##https://rdrr.io/cran/brms/man/brmsformula.html
##https://vuorre.netlify.app/post/2016/09/29/meta-analysis-is-a-special-case-of-bayesian-multilevel-modeling/


formula2 <- bf(Jtutrend | se(Jtutrend_se, sigma = TRUE) ~ temptrend_comb_allyr_abs.sc * seas_comb.sc + (1 | STUDY_ID))

formula3 <- bf(Jtutrend | se(Jtutrend_se, sigma = TRUE) ~ temptrend_comb_allyr_abs.sc * tempave_comb.sc + (1 | STUDY_ID))

#...

############################################################################
############################################################################

##the code below is for the simplest model [1]

##run model
mod_brms_1 <- brm(formula= formula1,
                     data= trends,  ##na.omit(trends)
                     control = list(adapt_delta = 0.99),
                     iter = 4000)

##control is used to improve sampler behaviour
##the remaining arguments were used as default; model as such uses non-informative flat priors
##default for chains=4; iter= 2000; warmup= iter/2

##we used 8k iter for the temperature analysis and running models for each realm separately



summary(mod_brms_1)  ##rhat = 1 indicates the chains have converged
plot(mod_brms_1)

##save model output object
##save(mod_brms_1, file="mod_brms_1.Rdata")



###it is possible to define priors, as well as different settings for the runs

# e.g. using get_prior() to see all the params for which priors can be defined
get_prior(formula1, data = trends)


##example from Shane's code
##	set some weakly regularising priors...
hier_prior <- c(set_prior(prior = 'normal(0,1)', class='b', coef='cYEAR'), 	# global slope
                set_prior(prior = 'normal(0,2)', class='Intercept', coef=''), 		# global intercept
                set_prior(prior = 'cauchy(0,2)', class='sd'),							# group-level intercepts and slopes
                set_prior(prior = 'lkj(2)', class='cor'))

Jtu_norm_BTSRfyID <- brm(bf(Jtu_base ~ cYEAR + (cYEAR|Biome/taxa_mod/STUDY_ID/rarefyID), 
                            family = brmsfamily('gaussian')),	
                         data= rarefied_medians,
                         prior=hier_prior,
                         inits = '0',
                         iter = 2000,
                         cores = 4,
                         chains = 4)



###second model (start the runs in pararell)
mod_brms_2 <- brm(formula= formula2,
                  data= trends,  ##na.omit(trends)
                  control = list(adapt_delta = 0.99),
                  iter = 4000)

############################################################################
############################################################################
##for model selection
##using WAIC and/or LOOIC

#https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/overfitting-regularization-and-information-criteria.html
#https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html


library(loo)

waic(mod_brms_1, mod_brms_2)
loo(mod_brms_1, mod_brms_2)


##might be useful to add the IC values to the model outputs, as this takes a long time to calculate
mod_brms_1 <- add_criterion(mod_brms_1, "waic")
mod_brms_2 <- add_criterion(mod_brms_2, "waic")

#add_criterion() also calculates "loo" - replace argument

# and compare the WAIC estimates
w <- loo_compare(mod_brms_1, mod_brms_2, criterion = "waic")









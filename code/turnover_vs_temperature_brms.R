##======================================================================
# Laura Antao and Malin Pinsky

## note the models can take several days to run

## The biodiversity data used originates from pre-processing steps from Blowes et al (https://github.com/sChange-workshop/BioGeo-BioDiv-Change)

##======================================================================

library(brms)  ##v 2.12.0
library(data.table)

###############
# Load the data
###############

trends <- fread('output/turnover_w_covariates.csv.gz')

#################
# Prep the data
#################
# trim to only data with some temperature change
# important since sign of temperature change is a variable
trends[tempchange == 0, .N] # number to remove
trends <- trends[tempchange != 0, ] # also removes any NA values

# set realm order
trends[, REALM := factor(REALM, levels = c('Freshwater', 'Marine', 'Terrestrial'), ordered = FALSE)]

# set up sign of temperature change
trends[, tsign := factor(sign(tempchange))]

# realm that combines Terrestrial and Freshwater, for interacting with human impact
trends[, REALM2 := REALM]
levels(trends$REALM2) = list(TerrFresh = "Freshwater", TerrFresh = "Terrestrial", Marine = "Marine")

# group Marine invertebrates/plants in with All
trends[, taxa_mod2 := taxa_mod]
trends[taxa_mod == 'Marine invertebrates/plants', taxa_mod2 := 'All']

# calculate duration
trends[, duration := year2 - year1]

#add a comparison id
trends[, compID := paste0(rarefyID, '_', year1, '_', year2)]

## Log-transform some variables, then center and scale. 
trends[, tempave.sc := scale(tempave)]
trends[, tempave_metab.sc := scale(tempave_metab)]
trends[, seas.sc := scale(seas)]
trends[, microclim.sc := scale(log(microclim))]
trends[, tempchange.sc := scale(tempchange, center = FALSE)] # do not center
trends[, tempchange_abs.sc := scale(abs(tempchange), center = FALSE)] # do not center, so that 0 is still 0 temperature change
trends[, mass.sc := scale(log(mass_mean_weight))]
trends[, speed.sc := scale(log(speed_mean_weight+1))]
trends[, lifespan.sc := scale(log(lifespan_mean_weight))]
trends[, consumerfrac.sc := scale(consfrac)]
trends[, endothermfrac.sc := scale(endofrac)]
trends[, nspp.sc := scale(log(Nspp))]
trends[, thermal_bias.sc := scale(thermal_bias)]
trends[, npp.sc := scale(log(npp))]
trends[, veg.sc := scale(log(veg+1))]
trends[, duration.sc := scale(log(duration))]
trends[, human_bowler.sc := scale(log(human_bowler+1)), by = REALM2] # separate scaling by realm
trends[REALM2 == 'TerrFresh', human_footprint.sc := scale(log(human_venter+1))]
trends[REALM2 == 'Marine', human_footprint.sc := scale(log(human_halpern))]


#################
## Set up models
#################
###create formula for brms models
rf0 <- bf(Jbeta ~ tempchange_abs.sc)
rf1 <- bf(Jbeta ~ tempchange_abs.sc + (1 | taxa_mod2), family = zero_one_inflated_beta())
rf2 <- bf(Jbeta ~ tempchange_abs.sc + (1 | taxa_mod2) + (1 | STUDY_ID))
rf3 <- bf(Jbeta ~ tempchange_abs.sc + (1 | taxa_mod2) + (1 | STUDY_ID) + (1 | rarefyID))

#...

############################################################################
# Fit models
############################################################################

##the code below is for the simplest model [1]

##run model
##control is used to improve sampler behavior
##the remaining arguments were used as default; model as such uses non-informative flat priors
##default for chains=4; iter= 2000; warmup= iter/2
modrf0 <- brm(formula= rf0, data= trends[1:10000,])
modrf1 <- brm(formula= rf1, data= trends)

modrf1 <- brm(formula= rf1,
                     data= trends,
                     control = list(adapt_delta = 0.99),
                     iter = 4000)


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









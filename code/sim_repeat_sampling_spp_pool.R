# source code to do turnover calculations
rm(list=ls())

library(tidyverse)

source('~/Dropbox/1current/BioTime/BioGeo-BioDiv-Change/R/06_rarefysamplesturnoverbnh.R')

set.seed(403)

# simulation parameters
n_sim <- 500

# this code loads the data and determines the characteristics of the empirical 
# distribution of time series lengths 

# load('~/Dropbox/BiogeoBioTIME/rarefied_medians.Rdata')
# cell_count <- rarefied_medians %>%
#   group_by(Biome, taxa_mod) %>%
#   dplyr::summarise(n_cells = n_distinct(rarefyID)) %>%
#   ungroup() 
# 
# ##	rejoin
# rarefied_medians <- left_join(cell_count, rarefied_medians, by=c('Biome', 'taxa_mod'))
# 
# # get parameters of poisson lognormal for duration
# duration_poilog <- poilog::poilogMLE(rarefied_medians %>% 
#                     group_by(rarefyID) %>%
#                     # distinct(num_years) %>%
#                     distinct(duration) %>% .$duration)
# 
# ##	filter to count data and biome/taxa combinations with >3 cells
# # these are the data that the model was fit to
# rarefied_medians <- rarefied_medians %>%
#   filter(BROAD_TYPE=='count' & n_cells > 3)

# these are the parameters that describe the empirical distribution of time series lengths
time_series_length <- floor(rlnorm(n_sim, meanlog = 2.3, sdlog = 0.65))
time_series_length <- time_series_length[time_series_length>1]  
# make sure we have enough of these values for each simulation
if(length(time_series_length) < n_sim){
  n_sim = length(time_series_length)
}

# want to examine for species pool dependence
regional_pool <- c(20, 50, 100, 200, 400, 800)

# number of samples for each simulation (max could be higher, but is set to 20
# to work for the smallest species pool)
sample_effort <- floor(runif(n_sim, min = 2, max = 20))

# initialise df to store slope estimates
jtu_slopes <- tibble()
jne_slopes <- tibble()

for(pool in 1:length(regional_pool)){
  jtu_slopes_temp <- c()
  jne_slopes_temp <- c()
  for(sim in 1:n_sim){
    print(paste('simulation ', sim, 'in ', n_sim, 'for species pool', pool, 'of', length(regional_pool)))
    # initialise empty vectors for this simulation  
    Year = c()
    SampleID = c()
    Species = c()
    Abundance = c()
  
    for(y in 1:time_series_length[sim]){
      # random sampling effort
      S = sample_effort[sim]
      # create year covariate for modelling 
      Year = c(Year,rep(y,S))
      # assign sampleID
      SampleID = c(SampleID,rep(1,S))
      # sample from regional pool
      Species = c(Species,sample(regional_pool[pool],S))
      # presence only data: abundance==1
      Abundance = c(Abundance,rep(1,S))
    }
    
    # call function to calculate temporal turnover metrics 
    # (function also does sample-based rarefaction to standardise sampling effort, not required here)
    a <- rarefysamplesturnoverbnh(Year, SampleID, Species, Abundance, resamps = 1)
    
    # fit models to turnover and nestedness; extract the slope coefficient
    jtu_slopes_temp = c(jtu_slopes_temp, with(a[[1]], coef(lm(Jtu ~ Year))['Year']))
    jne_slopes_temp = c(jne_slopes_temp, with(a[[1]], coef(lm(Jne ~ Year))['Year']))
  }
  
  jtu_slopes = bind_rows(jtu_slopes,
                         tibble(metric = 'Turnover',
                                spp_pool = regional_pool[pool],
                                slope = jtu_slopes_temp
                                )
  )
  jne_slopes = bind_rows(jne_slopes,
                         tibble(metric = 'Nestedness',
                                spp_pool = regional_pool[pool],
                                slope = jne_slopes_temp
                                )
  )
}


slopes = bind_rows(jtu_slopes,
                  jne_slopes) %>% 
  mutate(spp_pool = paste0('Species pool = ', spp_pool))

# order for plot
slopes$spp_pool <- factor(slopes$spp_pool,
                          levels = c("Species pool = 20",
                                     "Species pool = 50",
                                     "Species pool = 100",
                                     "Species pool = 200",
                                     "Species pool = 400",
                                     "Species pool = 800")
)
# calculate quantiles of simulated slopes
sim_quantiles <- slopes %>% 
  group_by(spp_pool, metric) %>% 
  summarise(median = median(slope, na.rm = T),
            lower = quantile(slope, probs = 0.025, na.rm = T),
            upper = quantile(slope, probs = 0.975, na.rm = T))

# ditch nestedness, there is nothing going on
ggplot() +
  facet_grid(spp_pool~metric) +
  geom_histogram(data = slopes %>% filter(metric=='Turnover'),
                 aes(x = slope), binwidth = 0.005) +
  geom_vline(data = sim_quantiles %>% filter(metric=='Turnover'),
             aes(xintercept = median),
             lty = 1) +
  geom_vline(data = sim_quantiles %>% filter(metric=='Turnover'),
             aes(xintercept = upper),
             lty = 2) +
  geom_vline(data = sim_quantiles %>% filter(metric=='Turnover'),
             aes(xintercept = lower),
             lty = 2) +
  geom_text(data = slopes %>% 
              group_by(spp_pool, metric) %>% 
              summarise(lower = quantile(slope, probs = 0.025, na.rm = T),
                        upper = quantile(slope, probs = 0.975, na.rm = T)) %>% 
              filter(metric=='Turnover'),
            aes(x = -0.2, y = 200,
                label = paste0('2.5% quantile = ', round(lower,3))),
            size = 3
  ) +
  geom_text(data = slopes %>% 
                group_by(spp_pool, metric) %>% 
                summarise(lower = quantile(slope, probs = 0.025, na.rm = T),
                          upper = quantile(slope, probs = 0.975, na.rm = T)) %>% 
                filter(metric=='Turnover'),
              aes(x = 0.2, y = 200,
                  label = paste0('97.5% quantile = ', round(upper,3))),
              size = 3) +  
  labs(x = 'Slope estimate',
       y = 'Number of simulations') +
  theme_bw()

ggsave('~/Dropbox/1current/conceptual/beta_spp_pool_sims/sim_results.png',
       width = 150, 
       height = 200, 
       units = 'mm')



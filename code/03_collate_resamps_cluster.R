##============================================================
##      script to combine rarefied resamples and calculate medians
##      for analysis
library(tidyverse)
##============================================================
rm(list=ls())
setwd('/data/idiv_chase/sablowes/biotime/resamps/with-gains-losses/')
##      set the pattern to load the files to be compiled
filelist = dir(pattern="bt-rarefy.+Rdata")
nSimul = length(filelist)

##      initialise data.frames to store all rarefied resamples
rarefy_abund_alpha <- data_frame()
rarefy_abund_beta <- data_frame()
rarefy_biomass_alpha <- data_frame()
rarefy_biomass_beta <- data_frame()
rarefy_pres_alpha <- data_frame()
rarefy_pres_beta <- data_frame()

for (iSimul in 1:nSimul) {
  ##      load a rarefied sample
  load(filelist[iSimul])
  print(iSimul)
  ##      add a column describing measurement type
  alpha_abundance <- rarefy_abund[[1]] %>% mutate(measurement_type = 'abundance')
  beta_abundance <- rarefy_abund[[2]] %>% mutate(measurement_type = 'abundance')
  
  alpha_biomass <- rarefy_biomass[[1]] %>% mutate(measurement_type = 'biomass')
  beta_biomass <- rarefy_biomass[[2]] %>% mutate(measurement_type = 'biomass')
  
  alpha_pres <- rarefy_pres[[1]] %>% mutate(measurement_type = 'presence')    
  beta_pres <- rarefy_pres[[2]] %>% mutate(measurement_type = 'presence')    
  
  ##      collate rarefied samples
  rarefy_abund_alpha <- bind_rows(rarefy_abund_alpha, alpha_abundance)
  rarefy_abund_beta <- bind_rows(rarefy_abund_beta, beta_abundance)
  
  rarefy_biomass_alpha <- bind_rows(rarefy_biomass_alpha, alpha_biomass)
  rarefy_biomass_beta <- bind_rows(rarefy_biomass_beta, beta_biomass)
  
  rarefy_pres_alpha <- bind_rows(rarefy_pres_alpha, alpha_pres)
  rarefy_pres_beta <- bind_rows(rarefy_pres_beta, beta_pres)
}

##      put them all together
alpha_metrics <- bind_rows(rarefy_abund_alpha, rarefy_biomass_alpha, rarefy_pres_alpha)
beta_metrics <- bind_rows(rarefy_abund_beta, rarefy_biomass_beta, rarefy_pres_beta)

##      pull out new metadata
new_meta <- alpha_metrics %>%
  distinct(rarefyID, SamplePool, SampleN, num_years, duration, startYear, endYear)


##      calculate the medians for all the metrics
rarefied_alpha_medians <- ungroup(alpha_metrics) %>%
  group_by(rarefyID, YEAR) %>%
  dplyr::summarise(
    N = median(N),
    N_int = round(median(N)),
    S = median(S),    # as above?
    S_int = round(median(S)),
    PIE = median(PIE),
    ENSPIE = median(ENSPIE),
    ENSPIE_int = round(median(ENSPIE)),
    pielou = median(pielou),
    Hill1 = median(Hill1)) %>%
  ungroup()

rarefied_beta_medians <- ungroup(beta_metrics) %>%
  # combine pair of years into single variable for grouping
  unite(yr_pair, c(YEAR1, YEAR2)) %>% 
  group_by(rarefyID, yr_pair) %>%
  dplyr::summarise(
    Jbeta = median(Jbeta),
    Jtu = median(Jtu),
    Jne = median(Jne),
    Hornsim = median(Hornsim),
    gains = median(gains),
    losses = median(losses)) %>%
  ungroup() %>%
  separate(yr_pair, c('YEAR1', 'YEAR2'))

##      recombine with new metadata
rarefied_alpha_medians <- inner_join(new_meta, rarefied_alpha_medians, by='rarefyID')
rarefied_beta_medians <- inner_join(new_meta, rarefied_beta_medians, by='rarefyID')

##      save
save(rarefied_alpha_medians,
     rarefied_beta_medians, 
     file=Sys.getenv('OFILE'))

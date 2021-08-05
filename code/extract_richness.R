# extract richness for each rarefyID

require(data.table)


# BioTime studies
load(here::here('data', 'biotime_blowes', 'all_pairs_beta.Rdata')) # load rarefied_beta_medians, which has all pairwise dissimilarities
bt <- data.table(rarefied_beta_medians); rm(rarefied_beta_medians)

# BioTime species
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has species list per study, including some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp

length(setdiff(bt$rarefyID, btspp$rarefyID)) # not missing any

btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID)], by = 'rarefyID', all.y = TRUE) # trim to focal studies

# BioTime abundance
load('data/biotime_blowes/bt_grid_spp_list_abund.Rdata') # loads abundance data
btabund <- as.data.table(bt_grid_spp_list_abund); rm(bt_grid_spp_list_abund, bt_grid_spp_list_abund_year) # rename and delete unneeded df

length(setdiff(bt$rarefyID, btabund$rarefyID)) # not missing any

btspp <- merge(btspp, btabund[, .(rarefyID, Species, N_bar)], by = c('rarefyID', 'Species'), all.x = TRUE)

btspp[, length(unique(Species))] # 27876 species
nrow(btspp) #1142602

# calculate richness in 3 ways: total species, total species with abund > 0 or NA, or total species with abund > 0 and not NA
rich <- btspp[, .(Nspp = .N, 
                  NsppNAorNonZero = sum(N_bar > 0 | is.na(N_bar)),
                  NsppNonZero = sum(N_bar > 0)), by = rarefyID]

rich

rich[Nspp != NsppNAorNonZero, ] # where trimming to abund >0 makes a difference. 630
summary(rich[Nspp != NsppNAorNonZero, ])

rich[Nspp != NsppNonZero, ] # where trimming to abund >0 and !NA makes a difference. 629
summary(rich[Nspp != NsppNonZero, ])

rich[, hist(Nspp)] # highly right-skewed
rich[, summary(Nspp)] # 1 to 1427
rich[, sum(Nspp == 1)] # 596
rich[, sum(Nspp <= 2)] # 3255

# how many missing timeseries? These had no species.
length(setdiff(bt$rarefyID, rich$rarefyID)) # 0

# write out
write.csv(rich, file = gzfile('output/richness_by_rarefyID.csv.gz'), row.names = FALSE)

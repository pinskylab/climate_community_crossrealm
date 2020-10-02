# extract richness for each rarefyID

require(data.table)


# BioTime studies
load('data/biotime_blowes/bt_malin.Rdata') # load bt_malin, has our focal studies
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt

# BioTime species
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has species list per study, including some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID)], by = 'rarefyID', all.y = TRUE) # trim to focal studies

# BioTime abundance
load('data/biotime_blowes/bt_grid_spp_list_abund.Rdata') # loads abundance data
btabund <- as.data.table(bt_grid_spp_list_abund); rm(bt_grid_spp_list_abund, bt_grid_spp_list_abund_year) # rename and delete unneeded df
btspp <- merge(btspp, btabund[, .(rarefyID, Species, N_bar)], by = c('rarefyID', 'Species'), all.x = TRUE)

btspp[, length(unique(Species))] # 26289 species
nrow(btspp) #1034700

# remove species with 0 abundance
btspp <- btspp[N_bar > 0, ]

btspp[, length(unique(Species))] # 23805 species
nrow(btspp) #1010338

# calculate richness
rich <- btspp[, .(Nspp = .N), by = rarefyID]

rich

rich[, hist(Nspp)] # highly right-skewed
rich[, summary(Nspp)] # 1 to 1427
rich[, sum(Nspp == 1)] # 457
rich[, sum(Nspp <= 2)] # 2697

# write out
write.csv(rich, file = gzfile('output/richness_by_rarefyID.csv.gz'), row.names = FALSE)

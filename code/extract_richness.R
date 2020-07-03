# extract richness for each rarefyID

require(data.table)


# BioTime species
load('data/biotime_blowes/bt_malin.Rdata') # load bt_malin, has our focal studies
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID)], by = 'rarefyID') # trim to focal studies

btspp[, length(unique(Species))] # 26289 species
nrow(btspp) #1034700

# calculate richness
rich <- btspp[, .(Nspp = .N), by = rarefyID]

rich

rich[, hist(Nspp)] # highly right-skewed
rich[, summary(Nspp)] # 1 to 1427
rich[, sum(Nspp == 1)] # 454
rich[, sum(Nspp <= 2)] # 2677

# write out
write.csv(rich, file = gzfile('output/richness_by_rarefyID.csv.gz'), row.names = FALSE)

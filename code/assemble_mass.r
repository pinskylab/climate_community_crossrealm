# assemble the mass data extracted from various databases

##############
## functions
##############
require(data.table)

geomean = function(x){ # geometric mean. remove NAs and x<0
    exp(sum(log(x[x > 0 & !is.na(x)])) / length(x[x > 0 & !is.na(x)]))
}

geosd <- function(x){ # geometric standard deviation
    exp(sd(log(x[x > 0 & !is.na(x)])))
}

##############
# Load data
##############

## read in trait data
vectraits <- fread('output/mass_BioTrait.csv.gz', drop = 1) # mass in g
try <- fread('output/mass_tryplants.csv.gz', drop = 1) # mass in g
gateway <- fread('output/mass_gateway.csv.gz', drop = 1) # appears to be mass in g
fishbase <- fread('output/mass_fishbase.csv.gz', drop = 1) # mass in g
sealifebase <- fread('output/mass_sealifebase.csv.gz', drop = 1) # mass in g
eltonbirds <- fread('output/mass_eltonbirds.csv.gz', drop = 1) # mass in g
eltonmammals <- fread('output/mass_eltonmammals.csv.gz', drop = 1) # mass in g

# BioTime change (for taxa_mod) and species
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID') # add taxa_mod to spp list

btspp[, length(unique(Species))] # 26289 species

#################
# Assemble data
#################

## add traits to BioTime
btmass <- merge(btspp, vectraits[, .(Species = gsub(' ', '_', Species), rarefyID, mass_vectraits = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, try[, .(Species, rarefyID, mass_try = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, gateway[, .(Species = gsub(' ', '_', Species), rarefyID, mass_gateway = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, fishbase[, .(Species, rarefyID, mass_fishbase = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, sealifebase[, .(Species, rarefyID, mass_sealifebase = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, eltonbirds[, .(Species, rarefyID, mass_eltonbirds = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, eltonmammals[, .(Species, rarefyID, mass_eltonmammals = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)

#######################################################
# pairwise comparisons among values from each source
#######################################################
btmass[!is.na(mass_vectraits) & !is.na(mass_try), .(.N, length(unique(Species)))] # 74, 14
    btmass[!is.na(mass_vectraits) & !is.na(mass_try), ][!duplicated(Species), plot(mass_vectraits, mass_try, log='xy')]; abline(0,1) # anti-correlated!
    btmass[!is.na(mass_vectraits) & !is.na(mass_try), table(taxa_mod)] # all Plants
    btmass[!is.na(mass_vectraits) & !is.na(mass_try), ][!duplicated(Species), ]
btmass[!is.na(mass_vectraits) & !is.na(mass_gateway), .(.N, length(unique(Species)))] # 102872 600
    btmass[!is.na(mass_vectraits) & !is.na(mass_gateway), ][!duplicated(Species), plot(mass_vectraits, mass_gateway, log='xy')]; abline(0,1) # correlated
    btmass[!is.na(mass_vectraits) & !is.na(mass_gateway), table(taxa_mod)] # a wide mix
btmass[!is.na(mass_vectraits) & !is.na(mass_fishbase), .(.N, length(unique(Species)))] # 77130 107
    btmass[!is.na(mass_vectraits) & !is.na(mass_fishbase), ][!duplicated(Species), plot(mass_vectraits, mass_fishbase, log='xy')]; abline(0,1) # uncorrelated
    btmass[!is.na(mass_vectraits) & !is.na(mass_fishbase), table(taxa_mod)] # a mix
btmass[!is.na(mass_vectraits) & !is.na(mass_sealifebase), .(.N, length(unique(Species)))] # 13041 19
    btmass[!is.na(mass_vectraits) & !is.na(mass_sealifebase), ][!duplicated(Species), plot(mass_vectraits, mass_sealifebase, log='xy')]; abline(0,1) # correlated
    btmass[!is.na(mass_vectraits) & !is.na(mass_sealifebase), table(taxa_mod)] # a mix
    btmass[!is.na(mass_vectraits) & !is.na(mass_sealifebase), ][!duplicated(Species), ]
btmass[!is.na(mass_vectraits) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 10515, 58
    btmass[!is.na(mass_vectraits) & !is.na(mass_eltonbirds), ][!duplicated(Species), plot(mass_vectraits, mass_eltonbirds, log='xy')]; abline(0,1) # correlated
    btmass[!is.na(mass_vectraits) & !is.na(mass_eltonbirds), table(taxa_mod)] # All, Birds
    btmass[!is.na(mass_vectraits) & !is.na(mass_eltonbirds), ][!duplicated(Species), .(Species, taxa_mod, mass_vectraits, mass_eltonbirds)]
btmass[!is.na(mass_vectraits) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 5, 2
    btmass[!is.na(mass_vectraits) & !is.na(mass_eltonmammals), ] # Mus_musculus, Mus_sp


btmass[!is.na(mass_try) & !is.na(mass_gateway), .(.N, length(unique(Species)))] # 1, 1
    btmass[!is.na(mass_try) & !is.na(mass_gateway), ] # Bougainvillea_glabra 
btmass[!is.na(mass_try) & !is.na(mass_fishbase), .(.N, length(unique(Species)))] # 0
btmass[!is.na(mass_try) & !is.na(mass_sealifebase), .(.N, length(unique(Species)))] # 0
btmass[!is.na(mass_try) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 1, 1
    btmass[!is.na(mass_try) & !is.na(mass_eltonbirds), ] # Magnolia_Warbler. Elton is correct, TRY is wrong
btmass[!is.na(mass_try) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 0

btmass[!is.na(mass_gateway) & !is.na(mass_fishbase), .(.N, length(unique(Species)))] # 266428 943
    btmass[!is.na(mass_gateway) & !is.na(mass_fishbase), ][!duplicated(Species), plot(mass_gateway, mass_fishbase, log='xy')]; abline(0,1) # correlated
    btmass[!is.na(mass_gateway) & !is.na(mass_fishbase), table(taxa_mod)] # a mix
btmass[!is.na(mass_gateway) & !is.na(mass_sealifebase), .(.N, length(unique(Species)))] # 50745 172
    btmass[!is.na(mass_gateway) & !is.na(mass_sealifebase), ][!duplicated(Species), plot(mass_gateway, mass_sealifebase, log='xy')]; abline(0,1) # correlated
    btmass[!is.na(mass_gateway) & !is.na(mass_sealifebase), table(taxa_mod)] # a mix
btmass[!is.na(mass_gateway) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 66048, 345
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonbirds), ][!duplicated(Species), plot(mass_gateway, mass_eltonbirds, log='xy')]; abline(0,1) # correlated
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonbirds), table(taxa_mod)] # All, Benthos, Birds, Plants
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonbirds) & taxa_mod != "Birds", ][!duplicated(Species), .(Species, taxa_mod, mass_gateway, mass_eltonbirds)]
btmass[!is.na(mass_gateway) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 5880 64
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonmammals), ][!duplicated(Species), plot(mass_gateway, mass_eltonmammals, log='xy')]; abline(0,1) # correlated
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonmammals), table(taxa_mod)] # All, Birds, Mammals
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonmammals), ][!duplicated(Species), .(Species, taxa_mod, mass_gateway, mass_eltonmammals)]

btmass[!is.na(mass_fishbase) & !is.na(mass_sealifebase), .(.N, length(unique(Species)))] # 0
btmass[!is.na(mass_fishbase) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 0
btmass[!is.na(mass_fishbase) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 0

btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 66502 214
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonbirds), ][!duplicated(Species), plot(mass_sealifebase, mass_eltonbirds, log='xy')] # correlated
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonbirds), table(taxa_mod)] # All, Benthos
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonbirds), ][!duplicated(Species), .(Species, taxa_mod, mass_sealifebase, mass_eltonbirds)]
btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 10932 60
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonmammals), ][!duplicated(Species), plot(mass_sealifebase, mass_eltonmammals, log='xy')] # correlated
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonmammals), table(taxa_mod)] # All
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonmammals), ][!duplicated(Species), .(Species, taxa_mod, mass_sealifebase, mass_eltonmammals)]

btmass[!is.na(mass_eltonbirds) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 0
    
    
# merge mass values based on a decision tree
btmass[, mass := NA_real_]
btmass[, mass_source := NA_character_]
btmass[!is.na(mass_eltonbirds), ':='(mass = mass_eltonbirds, mass_source = 'Elton birds')]
btmass[!is.na(mass_eltonmammals), ':='(mass = mass_eltonmammals, mass_source = 'Elton mammals')]
btmass[!is.na(mass_fishbase), ':='(mass = mass_fishbase, mass_source = 'Fishbase')]
btmass[!is.na(mass_sealifebase) & is.na(mass), ':='(mass = mass_sealifebase, mass_source = 'Sealifebase')]
btmass[!is.na(mass_try) & is.na(mass) & taxa_mod == "Plants", ':='(mass = mass_try, mass_source = 'TRY')]
btmass[!is.na(mass_gateway) & is.na(mass), ':='(mass = mass_gateway, mass_source = 'GATEway')]
btmass[!is.na(mass_vectraits) & is.na(mass), ':='(mass = mass_vectraits, mass_source = 'VecTraits')]
btmass[!is.na(mass_try) & is.na(mass), ':='(mass = mass_try, mass_source = 'TRY')]

# check
btmass[(!is.na(mass_eltonbirds) | !is.na(mass_eltonmammals) | !is.na(mass_fishbase) | !is.na(mass_sealifebase) |
           !is.na(mass_try) | !is.na(mass_gateway) | !is.na(mass_vectraits)) & is.na(mass), .N] # 0: good
btmass[, length(unique(Species)), by = !is.na(mass)] # 9426 species with data, 16863 without

btmass[, .(n = length(unique(Species)), val = sum(!is.na(mass))), by = rarefyID][, hist(val/n)] # most rarefyIDs have >50% of species represented!
btmass[, .(n = length(unique(Species)), val = sum(!is.na(mass))), by = rarefyID][(val/n) < 0.5, ] # 1726 rarefyID with <50% of 53467 (3.2%)

btmass[, .(n = length(unique(Species)), val = sum(!is.na(mass))), by = 
           .(rarefyID, STUDY_ID)][, .(n = sum(n), val = sum(val)), by = STUDY_ID][, hist(val/n)] # a bit more than half the studies have >50% of species represented!
btmass[, .(n = length(unique(Species)), val = sum(!is.na(mass))), by = 
           .(rarefyID, STUDY_ID, taxa_mod)][, .(n = sum(n), val = sum(val)), by = .(STUDY_ID, taxa_mod)][(val/n) < 0.5, ] #127 studies of 332 have <50% of species represented
btmass[, .(n = length(unique(Species)), val = sum(!is.na(mass))), by = 
           .(rarefyID, STUDY_ID, taxa_mod)][, .(n = sum(n), val = sum(val)), by = .(STUDY_ID, taxa_mod)][(val/n) < 0.5, table(taxa_mod)] # 66 of studies <50% are plants, 34 inverts

# make a list of priority studies for further mass finding
# less thatn 50% species with data, at least 5 species in the study
priorities <- btmass[rarefyID %in% bt[, rarefyID], .(n = length(unique(Species)), val = sum(!is.na(mass))), by = 
           .(rarefyID, STUDY_ID, taxa_mod)][, .(n = sum(n), val = sum(val)), by = .(STUDY_ID, taxa_mod)][(val/n) < 0.5 & n > 4, STUDY_ID]
length(priorities) # 124 studies



#########
# Output
#########

# output mass values
btmass.out <- btmass[!is.na(mass), .(Species, rarefyID, REALM, STUDY_ID, taxa_mod, mass, mass_source)]
nrow(btmass)
nrow(btmass.out) # 814040

write.csv(btmass.out, gzfile('output/mass_byspecies.csv.gz'), row.names = FALSE)

# output average mass by rarefyID
#see how we're doing (per rarefyID: # species with data, mean mass, sd mass, geometric mean mass, geometric standard deviation mass)
btmass.sum <- btmass[, .(mass_mean = mean(mass, na.rm = TRUE), mass_sd = sd(mass, na.rm = TRUE), mass_geomean = geomean(mass),
                         mass_geosd = geosd(mass), nspp = length(unique(Species)), nspp_wdata = sum(!is.na(mass))), 
                     by = .(STUDY_ID, rarefyID, REALM, taxa_mod)]
btmass.sum[is.nan(mass_mean), mass_mean := NA_real_]
btmass.sum[is.nan(mass_geomean), mass_geomean := NA_real_]
nrow(btmass.sum) # 53467
setkey(btmass.sum, STUDY_ID, rarefyID)
btmass.sum

write.csv(btmass.sum, gzfile('output/mass_byrarefyID.csv.gz'))

# output list of priority species
priorities.out <- btmass[STUDY_ID %in% priorities & rarefyID %in% bt[, rarefyID] & is.na(mass), .(Species, REALM, taxa_mod)][!duplicated(cbind(Species, REALM, taxa_mod)), .(REALM, taxa_mod, Species)]
setkey(priorities.out, REALM, taxa_mod, Species)
nrow(priorities.out) # 12907
priorities.out

write.csv(priorities.out, gzfile('output/mass_prioritymissing.csv.gz'))

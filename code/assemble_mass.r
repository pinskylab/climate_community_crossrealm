# assemble the mass data extracted from various databases

require(data.table)

## read in trait data
vectraits <- fread('output/mass_BioTrait.csv.gz', drop = 1)
try <- fread('output/mass_tryplants.csv.gz', drop = 1)
gateway <- fread('output/mass_gateway.csv.gz', drop = 1)
fishbase <- fread('output/mass_fishbase.csv.gz', drop = 1)
sealifebase <- fread('output/mass_sealifebase.csv.gz', drop = 1)
eltonbirds <- fread('output/mass_eltonbirds.csv.gz', drop = 1)
eltonmammals <- fread('output/mass_eltonmammals.csv.gz', drop = 1)

# BioTime change (for taxa_mod) and species
load('data/biotime_blowes/bt.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID') # add taxa_mod to spp list

## add traits to BioTime
btmass <- merge(btspp, vectraits[, .(Species = gsub(' ', '_', Species), rarefyID, mass_vectraits = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, try[, .(Species, rarefyID, mass_try = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, gateway[, .(Species = gsub(' ', '_', Species), rarefyID, mass_gateway = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, fishbase[, .(Species, rarefyID, mass_fishbase = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, sealifebase[, .(Species, rarefyID, mass_sealifebase = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, eltonbirds[, .(Species, rarefyID, mass_eltonbirds = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
btmass <- merge(btmass, eltonmammals[, .(Species, rarefyID, mass_eltonmammals = mass)], by = c('Species', 'rarefyID'), all.x = TRUE)

# compare among values
btmass[!is.na(mass_vectraits) & !is.na(mass_try), .(.N, length(unique(Species)))] # 74, 14
    btmass[!is.na(mass_vectraits) & !is.na(mass_try), ][!duplicated(Species), plot(mass_vectraits, mass_try, log='xy')] # anti-correlated!
    btmass[!is.na(mass_vectraits) & !is.na(mass_try), table(taxa_mod)] # all Plants
    btmass[!is.na(mass_vectraits) & !is.na(mass_try), ][!duplicated(Species), ]
btmass[!is.na(mass_vectraits) & !is.na(mass_gateway), .(.N, length(unique(Species)))] # 136731, 676
    btmass[!is.na(mass_vectraits) & !is.na(mass_gateway), ][!duplicated(Species), plot(mass_vectraits, mass_gateway, log='xy')] # correlated
    btmass[!is.na(mass_vectraits) & !is.na(mass_gateway), table(taxa_mod)] # a wide mix
btmass[!is.na(mass_vectraits) & !is.na(mass_fishbase), .(.N, length(unique(Species)))] # 57053, 110
    btmass[!is.na(mass_vectraits) & !is.na(mass_fishbase), ][!duplicated(Species), plot(mass_vectraits, mass_fishbase, log='xy')] # uncorrelated
    btmass[!is.na(mass_vectraits) & !is.na(mass_fishbase), table(taxa_mod)] # all Fish
btmass[!is.na(mass_vectraits) & !is.na(mass_sealifebase), .(.N, length(unique(Species)))] # 5237, 9
    btmass[!is.na(mass_vectraits) & !is.na(mass_sealifebase), ][!duplicated(Species), plot(mass_vectraits, mass_sealifebase, log='xy')] # correlated
    btmass[!is.na(mass_vectraits) & !is.na(mass_sealifebase), table(taxa_mod)] # All, Benthos, Inverts
    btmass[!is.na(mass_vectraits) & !is.na(mass_sealifebase), ][!duplicated(Species), ]
btmass[!is.na(mass_vectraits) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 10515, 58
    btmass[!is.na(mass_vectraits) & !is.na(mass_eltonbirds), ][!duplicated(Species), plot(mass_vectraits, mass_eltonbirds, log='xy')] # mildly correlated
    btmass[!is.na(mass_vectraits) & !is.na(mass_eltonbirds), table(taxa_mod)] # All, Birds
    btmass[!is.na(mass_vectraits) & !is.na(mass_eltonbirds), ][!duplicated(Species), .(Species, taxa_mod, mass_vectraits, mass_eltonbirds)]
btmass[!is.na(mass_vectraits) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 5, 2
    btmass[!is.na(mass_vectraits) & !is.na(mass_eltonmammals), ] # Mus_musculus, Mus_sp


btmass[!is.na(mass_try) & !is.na(mass_gateway), .(.N, length(unique(Species)))] # 1, 1
    btmass[!is.na(mass_try) & !is.na(mass_gateway), ] # Bougainvillea_glabra 
btmass[!is.na(mass_try) & !is.na(mass_fishbase), .(.N, length(unique(Species)))] # 0
btmass[!is.na(mass_try) & !is.na(mass_sealifebase), .(.N, length(unique(Species)))] # 0
btmass[!is.na(mass_try) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 1, 1
    btmass[!is.na(mass_try) & !is.na(mass_eltonbirds), ] # Magnolia_Warbler
btmass[!is.na(mass_try) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 0

btmass[!is.na(mass_gateway) & !is.na(mass_fishbase), .(.N, length(unique(Species)))] # 207991, 1123
    btmass[!is.na(mass_gateway) & !is.na(mass_fishbase), ][!duplicated(Species), plot(mass_gateway, mass_fishbase, log='xy')] # correlated
    btmass[!is.na(mass_gateway) & !is.na(mass_fishbase), table(taxa_mod)] # Fish
btmass[!is.na(mass_gateway) & !is.na(mass_sealifebase), .(.N, length(unique(Species)))] # 11067, 88
    btmass[!is.na(mass_gateway) & !is.na(mass_sealifebase), ][!duplicated(Species), plot(mass_gateway, mass_sealifebase, log='xy')] # correlated
    btmass[!is.na(mass_gateway) & !is.na(mass_sealifebase), table(taxa_mod)] # All, Benthos, Inverts
btmass[!is.na(mass_gateway) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 66048, 345
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonbirds), ][!duplicated(Species), plot(mass_gateway, mass_eltonbirds, log='xy')] # correlated
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonbirds), table(taxa_mod)] # All, Benthos, Birds, Plants
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonbirds) & taxa_mod != "Birds", ][!duplicated(Species), .(Species, taxa_mod, mass_gateway, mass_eltonbirds)]
btmass[!is.na(mass_gateway) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 5630, 62
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonmammals), ][!duplicated(Species), plot(mass_gateway, mass_eltonmammals, log='xy')] # correlated
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonmammals), table(taxa_mod)] # All, Birds, Mammals
    btmass[!is.na(mass_gateway) & !is.na(mass_eltonmammals), ][!duplicated(Species), .(Species, taxa_mod, mass_gateway, mass_eltonmammals)]

btmass[!is.na(mass_fishbase) & !is.na(mass_sealifebase), .(.N, length(unique(Species)))] # 0
btmass[!is.na(mass_fishbase) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 0
btmass[!is.na(mass_fishbase) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 0

btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonbirds), .(.N, length(unique(Species)))] # 4102, 26
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonbirds), ][!duplicated(Species), plot(mass_sealifebase, mass_eltonbirds, log='xy')] # correlated
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonbirds), table(taxa_mod)] # All, Benthos
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonbirds), ][!duplicated(Species), .(Species, taxa_mod, mass_sealifebase, mass_eltonbirds)]
btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 4077, 21
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonmammals), ][!duplicated(Species), plot(mass_sealifebase, mass_eltonmammals, log='xy')] # correlated
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonmammals), table(taxa_mod)] # All
    btmass[!is.na(mass_sealifebase) & !is.na(mass_eltonmammals), ][!duplicated(Species), .(Species, taxa_mod, mass_sealifebase, mass_eltonmammals)]

btmass[!is.na(mass_eltonbirds) & !is.na(mass_eltonmammals), .(.N, length(unique(Species)))] # 0
    
    
# merge mass values based on a decision tree
btmass[, mass := mass_eltonbirds]
btmass[!is.na(mass_eltonmammals), mass := mass_eltonmammals]
btmass[!is.na(mass_fishbase), mass := mass_fishbase]
btmass[!is.na(mass_sealifebase) & is.na(mass), mass := mass_sealifebase]
btmass[!is.na(mass_try) & is.na(mass) & taxa_mod == "Plants", mass := mass_try]
btmass[!is.na(mass_gateway) & is.na(mass), mass := mass_gateway]
btmass[!is.na(mass_vectraits) & is.na(mass), mass := mass_vectraits]
btmass[!is.na(mass_try) & is.na(mass), mass := mass_try]

# check
btmass[(!is.na(mass_eltonbirds) | !is.na(mass_eltonmammals) | !is.na(mass_fishbase) | !is.na(mass_sealifebase) |
           !is.na(mass_try) | !is.na(mass_gateway) | !is.na(mass_vectraits)) & is.na(mass), .N] # 0: good
btmass[, length(unique(Species)), by = !is.na(mass)] # 9455 species out of 17022 (about half)

btmass[, .(n = length(unique(Species)), val = sum(!is.na(mass))), by = rarefyID][, hist(val/n)] # most have >50% of species represented!

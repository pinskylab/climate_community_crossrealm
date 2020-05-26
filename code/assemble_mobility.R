# calculate mobility index for each species, average to community time-series
# uses maximum speed from Hirt et al. 2017 Nature E&E
# needs output from assemble_mass.r

require(data.table)

##############
# Load data
##############
mass <- fread('output/mass_byspecies.csv.gz') # assembled data on mass, in g
phylo <- fread('output/taxonomy.csv.gz') # phylogeny info

# BioTime change (for taxa_mod) and species
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID') # add taxa_mod to spp list and trims to spp in bt

#################
# process data
#################

# merge taxonomy and mass
phylo[, original_name := gsub(' ', '_', original_name)] # add _ back to name
phylo[, ':='(gateway.masses = NULL, V1 = NULL, X = NULL, NA. = NULL)] # remove extraneous columns
mob <- merge(mass[, .(original_name = Species, rarefyID, REALM, STUDY_ID, taxa_mod, mass, mass_source)], phylo, by = 'original_name', all.x = TRUE)
mob <- merge(btspp[, .(rarefyID, original_name = Species, REALM, STUDY_ID, taxa_mod)], 
             mob[, .(original_name, rarefyID, mass, mass_source, kingdom, phylum, class, order, family, genus, species)], 
             all.x = TRUE, by = c("rarefyID", "original_name"))

# make a pseudo-genus from original_name
mob[, original_genus := vapply(strsplit(original_name, "_"), `[`, 1, FUN.VALUE=character(1))]


# classify by movement mode
mob[, table(class, mass_source)] # to examine

mob[, mode := NA_character_]
mob[mass_source == 'Fishbase', mode := 'swimming']
mob[class %in% c('Actinopterygii', 'Cephalaspidomorphi', 'Elasmobranchii', 'Holocephali', 'Cephalopoda'), mode := 'swimming'] # fishes
mob[class == 'Aves' & order == 'Sphenisciformes', mode := 'swimming'] # penguins
mob[order %in% c('Cetacea', 'Sirenia') | family %in% c('Odobenidae', 'Otariidae', 'Phocidae') | genus %in% c('Enhydra'), mode := 'swimming'] # whales, manatees, seals, sea otter
mob[phylum %in% c('Rotifera', 'Chaetognatha', 'Ctenophora'), mode := 'swimming'] # rotifers, arrow worms, comb jellies
mob[class %in% c('Branchiopoda', 'Ostracoda', 'Scyphozoa', 'Appendicularia', 'Thaliacea'), mode := 'swimming'] # fairy shrimp, seed shrimp, true jellies, tunicates
mob[order %in% c('Calanoida', 'Cyclopoida', 'Harpacticoida', 'Siphonostomatoida'), mode := 'swimming'] # copepods
mob[order %in% c('Amphipoda', 'Euphausiacea', 'Mysida', 'Siphonophorae'), mode := 'swimming'] # amphipods, krill, mysid shrimp, siphonophores
mob[original_genus == 'Acartia', mode := 'swimming'] # certain copepods

mob[mass_source == 'Elton birds' & (order != 'Sphenisciformes' | is.na(order)), mode := 'flying']
mob[class == 'Aves' & order != 'Sphenisciformes', mode := 'flying']
mob[class == 'Insecta' & !(order %in% c('Plecoptera', 'Trichoptera')), mode := 'flying'] # except stoneflies and caddisflies
mob[original_genus %in% c('Cecidomyiidae'), mode := 'flying'] # gall midges
mob[REALM == 'Terrestrial' & class == 'Insecta' & order %in% c('Plecoptera', 'Trichoptera'), mode := 'crawling'] # stoneflies and caddisflies

mob[mass_source == 'Elton mammals' & !(order %in% c('Cetacea', 'Sirenia')) & !(family %in% c('Odobenidae', 'Otariidae', 'Phocidae')), mode := 'running']
mob[class %in% c('Reptilia', 'Amphibia'), mode := 'running']
mob[class %in% c('Arachnida', 'Entognatha'), mode := 'running'] # spiders, Entognaths
mob[REALM == 'Terrestrial' & order %in% c('Isopoda'), mode := 'running'] # terrestrial isopods
mob[REALM == 'Terrestrial' & original_genus %in% c('Isopoda'), mode := 'running'] # terrestrial isopods
mob[original_genus %in% c('Acari', 'Cryptorhynchinae'), mode := 'running'] # mites & ticks, weevils
mob[REALM == 'Freshwater' & class == 'Insecta' & order %in% c('Plecoptera', 'Trichoptera'), mode := 'running'] # stoneflies and caddisflies
mob[class %in% c('Merostomata', 'Pycnogonida'), mode := 'running'] # horseshoe crabs, sea spiders
mob[order %in% c('Cumacea', 'Decapoda', 'Stomatopoda'), mode := 'running'] # arthropods in mud/sand, crabs/lobsters, mantis shrimp
mob[REALM == 'Marine' & order %in% c('Isopoda'), mode := 'running'] # marine isopods
mob[REALM == 'Marine' & original_genus %in% c('Isopoda'), mode := 'running'] # marine isopods

mob[phylum %in% c( 'Annelida', 'Echinodermata', 'Sipuncula', 'Cephalorhyncha', 'Nemertea'), mode := 'crawling'] # segmented worms, echinoderms, sipunculid worms, ribbon worms
mob[class %in% c('Gastropoda'), mode := 'crawling'] # snails and slugs
mob[class %in% c('Polyplacophora', 'Priapulida'), mode := 'crawling'] # chitons, priapulid worms
mob[original_genus %in% c('Turbellaria', 'Valvata'), mode := 'crawling'] # flatworms (some can swim...), snails

mob[mass_source == 'TRY', mode := 'stationary'] # plants
mob[kingdom %in% c('Bacteria', 'Chromista', 'Fungi', 'Plantae'), mode := 'stationary']
mob[original_genus %in% c('Chaetoceros', 'Chlamydomonas'), mode := 'stationary'] # diatoms, green algae
mob[phylum %in% c('Bryozoa', 'Porifera'), mode := 'stationary'] # bryozoans, sponges
mob[class %in% c('Bivalvia', 'Scaphopoda', 'Ascidiacea', 'Anthozoa'), mode := 'stationary'] # clams/oysters, tusk shells, ascidians, anemones/corals
mob[order %in% c('Sessilia', 'Anthoathecata', 'Leptothecata'), mode := 'stationary'] # barnacles, hydrozoans, 


# examine what is missing
mob[is.na(mode) & !duplicated(species), table(phylum)] # nothing missing that has a scientific name match
mob[is.na(mode) & !duplicated(original_name), table(mass_source)]


# calculate max speed (km/hr) from allometry
# Table S4 in Hirt et al. 2017 Nature E&E
mob[, mass_kg := mass * 1000] # convert to kg from g
mob[, speed := NA_real_]
mob[mode == 'flying', speed := 142.8 * mass_kg^0.24 * (1 - exp(-2.4 * mass_kg^(-0.72)))]
mob[mode == 'running', speed := 25.5 * mass_kg^0.26 * (1 - exp(-22 * mass_kg^(-0.6)))]
mob[mode == 'crawling', speed := 0.1] # set crawlers to something below the slowest runners
mob[mode == 'swimming', speed := 11.2 * mass_kg^0.36 * (1 - exp(-19.5 * mass_kg^(-0.56)))]
mob[mode == 'stationary', speed := 0]


# check coverage
mob[, .(n = length(unique(original_name)), val = sum(!is.na(speed))), by = rarefyID][, hist(val/n)] # most have >50% of species represented!
mob[, .(n = length(unique(original_name)), val = sum(!is.na(speed))), by = rarefyID][(val/n) < 0.5, ] # 4159 rarefyID with <50%



# make a list of priority studies for further mobility finding
# have mass, but less than 50% species with mobility, at least 5 species in the study with mass
priorities <- mob[!is.na(mass), .(n = length(unique(original_name)), val = sum(!is.na(speed))), by = 
                         .(rarefyID, STUDY_ID, taxa_mod)][, .(n = sum(n), val = sum(val)), by = .(STUDY_ID, taxa_mod)][(val/n) < 0.5 & n > 4, STUDY_ID]
length(priorities) # 13 studies
priorities.out <- mob[STUDY_ID %in% priorities & is.na(speed), .(original_name, original_genus, REALM, taxa_mod)][!duplicated(cbind(original_name, REALM, taxa_mod)), .(REALM, taxa_mod, original_name, original_genus)] # 552 species
nrow(priorities.out) # 1492 taxa

priorities.out[, .(n = length(original_name)), by = original_genus][order(n, decreasing = TRUE)]





############
# output
############

# output speed values
mob.out <- mob[!is.na(speed), .(original_name, rarefyID, REALM, STUDY_ID, taxa_mod, speed, mode)]
nrow(mob)
nrow(mob.out) # 765848

write.csv(mob.out, gzfile('output/speed_byspecies.csv.gz'), row.names = FALSE)


# output average speed by rarefyID
#see how we're doing (per rarefyID: # species with data, mean mass, sd mass, geometric mean mass, geometric standard deviation mass)
mob.sum <- mob[, .(speed_mean = mean(speed, na.rm = TRUE), speed_sd = sd(speed, na.rm = TRUE), speed_geomean = geomean(speed),
                   speed_geosd = geosd(speed), nspp = length(unique(original_name)), nspp_wdata = sum(!is.na(speed))), 
                     by = .(STUDY_ID, rarefyID, REALM, taxa_mod)]
mob.sum[is.nan(speed_mean), speed_mean := NA_real_]
mob.sum[is.nan(speed_geomean), speed_geomean := NA_real_]
nrow(mob.sum) # 53467
setkey(mob.sum, STUDY_ID, rarefyID)
mob.sum

write.csv(mob.sum, gzfile('output/speed_byrarefyID.csv.gz'))


# output list of priority species
setkey(priorities.out, REALM, taxa_mod, original_genus, original_name)
nrow(priorities.out) # 1492
priorities.out

write.csv(priorities.out, gzfile('output/speed_prioritymissing.csv.gz'))



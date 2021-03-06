# calculate mobility index for each species, average to community time-series
# uses maximum speed from Hirt et al. 2017 Nature E&E
# needs output from assemble_mass.r

require(data.table)
source('code/util.R')

##############
# Load data
##############
mass <- fread('output/mass_byspecies.csv.gz') # assembled data on mass, in g
phylo <- fread('output/taxonomy.csv.gz') # phylogeny info
crawl <- fread('data/crawling/Crawling_speeds.csv') # crawling speeds from the literature

# BioTime change (for taxa_mod) and species
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID') # add taxa_mod to spp list and trims to spp in bt

# Add abundance information
load('data/biotime_blowes/bt_grid_spp_list_abund.Rdata') # loads abundance data
btabund <- as.data.table(bt_grid_spp_list_abund); rm(bt_grid_spp_list_abund, bt_grid_spp_list_abund_year) # rename and delete unneeded df
btspp <- merge(btspp, btabund[, .(rarefyID, Species, N_bar)], by = c('rarefyID', 'Species'), all.x = TRUE)

btspp[, length(unique(Species))] # 26289 species
nrow(btspp) #1034700

#################
# process data
#################

# merge taxonomy and mass
phylo[, original_name := gsub(' ', '_', original_name)] # add _ back to name
phylo[, ':='(gateway.masses = NULL, V1 = NULL, X = NULL, NA. = NULL)] # remove extraneous columns
mob <- merge(mass[, .(original_name = Species, rarefyID, REALM, STUDY_ID, taxa_mod, mass, mass_source)], 
             phylo, by = 'original_name', all.x = TRUE)
mob <- merge(btspp[, .(rarefyID, original_name = Species, REALM, STUDY_ID, taxa_mod, N_bar)], 
             mob[, .(original_name, rarefyID, mass, mass_source, kingdom, phylum, class, order, family, genus, species)], 
             all.x = TRUE, by = c("rarefyID", "original_name"))

# make a pseudo-genus from original_name
mob[, original_genus := vapply(strsplit(original_name, "_"), `[`, 1, FUN.VALUE=character(1))]


# classify by movement mode
mob[, table(class, mass_source)] # to examine

mob[, mode := NA_character_]

# swimming
mob[mass_source == 'Fishbase', mode := 'swimming']
mob[class %in% c('Actinopterygii', 'Cephalaspidomorphi', 'Elasmobranchii', 'Holocephali', 'Cephalopoda'), mode := 'swimming'] # fishes
mob[class == 'Aves' & order == 'Sphenisciformes', mode := 'swimming'] # penguins
mob[order %in% c('Cetacea', 'Sirenia') | family %in% c('Odobenidae', 'Otariidae', 'Phocidae') | genus %in% c('Enhydra'), mode := 'swimming'] # whales, manatees, seals, sea otter
mob[phylum %in% c('Rotifera', 'Chaetognatha', 'Ctenophora'), mode := 'swimming'] # rotifers, arrow worms, comb jellies
mob[class %in% c('Branchiopoda', 'Ostracoda', 'Scyphozoa', 'Appendicularia', 'Thaliacea'), mode := 'swimming'] # fairy shrimp, seed shrimp, true jellies, tunicates
mob[order %in% c('Calanoida', 'Cyclopoida', 'Harpacticoida', 'Siphonostomatoida'), mode := 'swimming'] # copepods
mob[order %in% c('Amphipoda', 'Euphausiacea', 'Mysida', 'Siphonophorae'), mode := 'swimming'] # amphipods, krill, mysid shrimp, siphonophores
mob[original_genus == 'Acartia', mode := 'swimming'] # certain copepods
mob[original_genus == 'Delphinidae', mode := 'swimming'] # dolphins
mob[original_genus == 'Sphyrna', mode := 'swimming'] # hammerhead sharks

# flying
mob[mass_source == 'Elton birds' & (order != 'Sphenisciformes' | is.na(order)), mode := 'flying']
mob[class == 'Aves' & order != 'Sphenisciformes', mode := 'flying']
mob[class == 'Insecta' & !(order %in% c('Plecoptera', 'Trichoptera')), mode := 'flying'] # except stoneflies and caddisflies
mob[original_genus %in% c('Cecidomyiidae'), mode := 'flying'] # gall midges
mob[REALM == 'Terrestrial' & class == 'Insecta' & order %in% c('Plecoptera', 'Trichoptera'), mode := 'crawling'] # stoneflies and caddisflies
mob[original_name == 'Unknown_Yellow-rumped-Warbler', mode := 'flying'] # bird
mob[is.na(mode) & taxa_mod == 'Birds', mode := 'flying'] # probably birds that can fly, too

# running
mob[mass_source == 'Elton mammals' & !(order %in% c('Cetacea', 'Sirenia')) & !(family %in% c('Odobenidae', 'Otariidae', 'Phocidae')), mode := 'running']
mob[class %in% c('Reptilia', 'Amphibia') & !(family %in% c('Colubridae', 'Elapidae', 'Natricidae', 'Pygopodidae', 'Typhlopidae')), mode := 'running'] # reptiles and amphibians run, except snakes/legless lizards
mob[class %in% c('Arachnida', 'Entognatha'), mode := 'running'] # spiders, Entognaths
mob[REALM == 'Terrestrial' & order %in% c('Isopoda'), mode := 'running'] # terrestrial isopods
mob[REALM == 'Terrestrial' & original_genus %in% c('Isopoda'), mode := 'running'] # terrestrial isopods
mob[original_genus %in% c('Acari', 'Cryptorhynchinae'), mode := 'running'] # mites & ticks, weevils
mob[REALM == 'Freshwater' & class == 'Insecta' & order %in% c('Plecoptera', 'Trichoptera'), mode := 'running'] # stoneflies and caddisflies
mob[class %in% c('Merostomata', 'Pycnogonida'), mode := 'running'] # horseshoe crabs, sea spiders
mob[order %in% c('Cumacea', 'Decapoda', 'Stomatopoda'), mode := 'running'] # arthropods in mud/sand, crabs/lobsters, mantis shrimp
mob[REALM == 'Marine' & order %in% c('Isopoda'), mode := 'running'] # marine isopods
mob[REALM == 'Marine' & original_genus %in% c('Isopoda'), mode := 'running'] # marine isopods
mob[original_genus %in% c('Scoloplos'), mode := 'crawling'] # polychaetes
mob[is.na(mode) & taxa_mod == 'Amphibians', mode := 'running'] # probably amphibians that can run, too
mob[is.na(mode) & taxa_mod == 'Reptiles', mode := 'running'] # probably reptiles that can run, too

# crawling
mob[phylum %in% c( 'Annelida', 'Echinodermata', 'Sipuncula', 'Cephalorhyncha', 'Nemertea'), mode := 'crawling'] # segmented worms, echinoderms, sipunculid worms, ribbon worms
mob[class %in% c('Gastropoda'), mode := 'crawling'] # snails and slugs
mob[order == 'Squamata' & family %in% c('Colubridae', 'Elapidae', 'Natricidae', 'Pygopodidae', 'Typhlopidae'), mode := 'crawling'] # snakes and legless lizards
mob[class %in% c('Polyplacophora', 'Priapulida'), mode := 'crawling'] # chitons, priapulid worms
mob[original_genus %in% c('Turbellaria', 'Valvata', 'Turbonilla', 'Odostomia', 'Epitonium'), mode := 'crawling'] # flatworms (some can swim...), snails
mob[original_genus %in% c('Turbonilla', 'Odostomia', 'Epitonium', 'Notomastus'), mode := 'crawling'] # sea snails, bristleworms

# stationary
mob[mass_source == 'TRY', mode := 'stationary'] # plants
mob[kingdom %in% c('Bacteria', 'Chromista', 'Fungi', 'Plantae'), mode := 'stationary']
mob[original_genus %in% c('Chaetoceros', 'Chlamydomonas', 'Diatom', 'Thalassiosira'), mode := 'stationary'] # diatoms, green algae
mob[phylum %in% c('Bryozoa', 'Porifera'), mode := 'stationary'] # bryozoans, sponges
mob[class %in% c('Bivalvia', 'Scaphopoda', 'Ascidiacea', 'Anthozoa'), mode := 'stationary'] # clams/oysters, tusk shells, ascidians, anemones/corals
mob[order %in% c('Sessilia', 'Anthoathecata', 'Leptothecata'), mode := 'stationary'] # barnacles, hydrozoans, 
mob[original_name %in% c('Pandora_gouldiana', 'Pandora_inflata', 'Pandora_wardiana', 'Pandora_inornata', 'Pandora_NA', 
                         'Pandora_trilineata', 'Pandora_bushiana', 'Pandora_filosa'), mode := 'flying'] # bivalves, 
mob[is.na(mode) & taxa_mod == 'Plant', mode := 'stationary'] # probably plants


# examine what is missing
mob[is.na(mode) & !duplicated(species), table(phylum)] # nothing missing that has a scientific name match
mob[is.na(mode) & !duplicated(original_name), table(mass_source)] # many GATEway species remain unclassified
tab <- mob[is.na(mode) & !duplicated(original_name) & !is.na(mass_source), table(original_genus, mass_source)]
tail(tab[order(rowSums(tab[,1:2]), decreasing = FALSE),]) # the most common genera missing

# calculate allometry for crawling
cmod <- crawl[, lm(log10(`max. speed [km/h]`) ~ log10(`body mass [kg]`))]
crawl_coef <- coef(cmod)
crawl_coef

# calculate max speed (km/hr) from allometry
# Table S4 in Hirt et al. 2017 Nature E&E
mob[, mass_kg := mass / 1000] # convert to kg from g
mob[, speed := NA_real_]
mob[mode == 'flying', speed := 142.8 * mass_kg^0.24 * (1 - exp(-2.4 * mass_kg^(-0.72)))]
mob[mode == 'running', speed := 25.5 * mass_kg^0.26 * (1 - exp(-22 * mass_kg^(-0.6)))]
mob[mode == 'crawling', speed := 10^(crawl_coef[1] + crawl_coef[2] * log10(mass_kg))]
mob[mode == 'swimming', speed := 11.2 * mass_kg^0.36 * (1 - exp(-19.5 * mass_kg^(-0.56)))]
mob[mode == 'stationary', speed := 0]

# examine the calculations
mob[!is.na(mode), .(min = min(speed, na.rm = TRUE), max = max(speed, na.rm = TRUE)), by = mode]

# check coverage
mob[, .(n = length(unique(original_name)), val = sum(!is.na(speed))), by = rarefyID][, hist(val/n)] # most have >50% of species represented!
mob[, .(n = length(unique(original_name)), val = sum(!is.na(speed))), by = rarefyID][(val/n) < 0.5, ] # 3857 rarefyID with <50%



# make a list of priority studies for further mobility finding
# have mass, but less than 50% species with mobility, at least 5 species in the study with mass
priorities <- mob[!is.na(mass), .(n = length(unique(original_name)), val = sum(!is.na(speed))), by = 
                         .(rarefyID, STUDY_ID, taxa_mod)][, .(n = sum(n), val = sum(val)), by = .(STUDY_ID, taxa_mod)][(val/n) < 0.5 & n > 4, STUDY_ID]
length(priorities) # 8 studies
priorities.out <- mob[STUDY_ID %in% priorities & is.na(speed), .(original_name, original_genus, REALM, taxa_mod)][!duplicated(cbind(original_name, REALM, taxa_mod)), .(REALM, taxa_mod, original_name, original_genus)] # 552 species
nrow(priorities.out) # 767 taxa

priorities.out[, .(n = length(original_name)), by = original_genus][order(n, decreasing = TRUE)]





############
# output
############

# output speed values
mob.out <- mob[!is.na(speed), .(original_name, rarefyID, REALM, STUDY_ID, taxa_mod, speed, mode)]
nrow(mob)
nrow(mob.out) # 782382

write.csv(mob.out, gzfile('output/speed_byspecies.csv.gz'), row.names = FALSE)


# output average speed by rarefyID
#see how we're doing (per rarefyID: # species with data, mean mass, sd mass, geometric mean mass, geometric standard deviation mass)
mob.sum <- mob[, .(speed_mean = mean(speed, na.rm = TRUE), speed_sd = sd(speed, na.rm = TRUE), 
                   speed_geomean = geomean(speed), speed_geosd = geosd(speed), 
                   speed_mean_weight = meanwt(speed, N_bar), speed_sd_weight = sdwt(speed, N_bar),
                   nspp = length(unique(original_name)), nspp_wdata = sum(!is.na(speed))), 
                     by = .(STUDY_ID, rarefyID, REALM, taxa_mod)]
mob.sum[is.nan(speed_mean), speed_mean := NA_real_]
mob.sum[is.nan(speed_geomean), speed_geomean := NA_real_]
mob.sum[is.nan(speed_mean_weight), speed_mean_weight := NA_real_]
summary(mob.sum) # summary
nrow(mob.sum) # 53467
setkey(mob.sum, STUDY_ID, rarefyID)
mob.sum

write.csv(mob.sum, gzfile('output/speed_byrarefyID.csv.gz'))


# output list of priority species
setkey(priorities.out, REALM, taxa_mod, original_genus, original_name)
nrow(priorities.out) # 767
priorities.out

write.csv(priorities.out, gzfile('output/speed_prioritymissing.csv.gz'))



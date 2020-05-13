################################################################
#######  extract body mass from the Elton Traits database  #####
######   birds and mammals                                 #####
################################################################

###########
# Packages
###########
require(data.table)

############
# Load and merge data
############
# BioTime change, useful for taxa_mod
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt

# BioTime species
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp

# add taxa_mod to spp list
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID')
nrow(btspp) # 1034700

# load Elton traits birds. mass in g
birds <- fread('dataDL/eltontraits/BirdFuncDat.txt')
setnames(birds, "BodyMass-Value", "mass")

# load Elton traits mammals. mass in g
mams <- fread('dataDL/eltontraits/MamFuncDat.txt')
setnames(mams, "BodyMass-Value", "mass")

# load gbif species names
phylo = fread('output/taxonomy.csv.gz', drop = 1)

# merge spp with gbif names
btspp2 <- merge(btspp, phylo[, .(class_gbif = class, genus_gbif = genus, Species_gbif = species, Species = gsub(' ', '_', original_name))], all.x = TRUE, by = 'Species')

#########################
# Find BM for birds
#########################

# merge with bird mass on biotime name
btbirds <- merge(btspp2, birds[, .(Species = gsub(' ', '_', Scientific), spp_elton = Scientific, mass_species_bt = mass)], by = 'Species', all.x = TRUE)

# merge with bird mass on gbif name
btbirds <- merge(btbirds, birds[, .(Species_gbif = Scientific, mass_species_gbif = mass)], by = 'Species_gbif', all.x = TRUE)
btbirds[!is.na(mass_species_gbif) & is.na(mass_species_bt), .(.N, length(unique(Species)))] # adds 4716 rows of data across 63 species

# match on common name in Elton
btbirds <- merge(btbirds, birds[, .(Species = gsub(' ', '_', English), spp_elton_comname = Scientific, mass_comname = mass)], by = 'Species', all.x = TRUE)

btbirds[!is.na(mass_comname), ] # look at the matches
btbirds[!is.na(mass_comname), length(unique(Species))] # 121

# look at the matches
btbirds[!is.na(mass_species_bt) | !is.na(mass_species_gbif) | !is.na(mass_comname), ]
btbirds[!is.na(mass_species_bt) | !is.na(mass_species_gbif) | !is.na(mass_comname), length(unique(Species))] # 1166
btbirds[!is.na(mass_species_bt) | !is.na(mass_species_gbif) | !is.na(mass_comname), sort(unique(taxa_mod))] # some are not labeled as birds by BioTime
btbirds[(!is.na(mass_species_bt) | !is.na(mass_species_gbif) | !is.na(mass_comname)) & taxa_mod != 'Birds', .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # in Benthos, Fish, All, Mammals. But all seem correct.

# compare matches on biotime and gbif names. no overlaps with common name to plot
btbirds[, plot(mass_species_bt, mass_species_gbif, log='xy')]; abline(0,1) # on 1:1 line: good!

# combine bt, gbif, and common name masses
btbirds[, mass_species := mass_species_bt]
btbirds[is.na(mass_species_bt) & !is.na(mass_species_gbif), mass_species := mass_species_gbif]
btbirds[is.na(mass_species) & !is.na(mass_comname), mass_species := mass_comname]
btbirds[!is.na(mass_species), length(unique(Species))] # 1166. good that it matches from before.


# check for missing
btbirds[taxa_mod == 'Birds' & is.na(mass_species), length(unique(Species))] # 833 species
btbirds[taxa_mod == 'Birds' & is.na(mass_species), sort(unique(Species))] # mostly common names and genus names, some 4-letter codes. Also a seal.



# get genus-level values
btbirds[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))] # not necessarily the genus, if "species name" is not a binomial
birds[, genus := vapply(strsplit(Scientific, " "), `[`, 1, FUN.VALUE=character(1))]
btbirdsgenu <- btbirds[(taxa_mod == 'Birds' | class_gbif == 'Aves'), .(genus = sort(unique(c(genus, genus_gbif))))] # list of potential bird genera (ignore potential birds in other categories)

btbirdsgenu[, mass_genus := NA_real_] 
for(i in 1:nrow(btbirdsgenu)){ # a brute-force for loop approach, but fast enough
    k <- which(btbirdsgenu$genus[i] == birds$genus)
    if(length(k) > 0){
        btbirdsgenu$mass_genus[i] <- mean(birds$mass[k])
    }
}
btbirdsgenu[!is.na(mass_genus),] # found 382 values
btbirdsgenu[!is.na(mass_genus), sort(unique(genus))] # check by eye that these look like genera

# merge genus-values with species list
btbirds <- merge(btbirds, btbirdsgenu[, .(genus, mass_genus_bt = mass_genus)], by = 'genus', all.x = TRUE) # match on bt genus
btbirds <- merge(btbirds, btbirdsgenu[, .(genus_gbif = genus, mass_genus_gbif = mass_genus)], by = 'genus_gbif', all.x = TRUE) # match on gbif genus

# examine cases where bt and gbif genus values do not match
btbirds[mass_genus_bt != mass_genus_gbif, ][!duplicated(Species), .(Species, Species_gbif, genus, genus_gbif, mass_genus_bt, mass_genus_gbif)]

# choose from species-, common- or genus-level mass values
btbirds[, mass := mass_species]
btbirds[is.na(mass_species), mass := mass_comname] # use common where species-level value is NA
btbirds[is.na(mass_species) & is.na(mass_comname), mass := mass_genus_gbif] # use genus where species- and common-level values are NA. prefer gbif genus
btbirds[is.na(mass_species) & is.na(mass_comname) & is.na(mass_genus_gbif), mass := mass_genus_bt] # bt genus as last choice
btbirds[is.nan(mass), mass := NA] # NaN to NA

btbirds[!is.na(mass), length(unique(Species))] # 1310 species now with data

btbirds[mass_genus_gbif != mass_species & !is.na(mass_genus_gbif), plot(mass_genus_gbif, mass_species)] # compare species and genus. pretty well correlated
btbirds[mass_genus_bt != mass_species & !is.na(mass_genus_bt), plot(mass_genus_bt, mass_species)] # compare species and genus. pretty well correlated

# merge Elton scientific names (from scientific vs. common-name merging)
btbirds[is.na(spp_elton) & !is.na(spp_elton_comname), spp_elton := spp_elton_comname]
btbirds[, spp_elton_comname := NULL]

# trim to species with data
nrow(btbirds) # 1,034,700
btbirds.out <- btbirds[!is.na(mass), ]
nrow(btbirds.out) # 211,269

### write out
write.csv(btbirds.out, file = gzfile('output/mass_eltonbirds.csv.gz'))


#########################
# Find BM for mammals
#########################
# merge with mamals mass on biotime name
btmams <- merge(btspp2, mams[, .(Species = gsub(' ', '_', Scientific), spp_elton = Scientific, mass_species_bt = mass)], by = 'Species', all.x = TRUE)

# merge with mass on gbif name
btmams <- merge(btmams, mams[, .(Species_gbif = Scientific, mass_species_gbif = mass)], by = 'Species_gbif', all.x = TRUE)
btmams[!is.na(mass_species_gbif) & is.na(mass_species_bt), .(.N, length(unique(Species)))] # adds 25 rows of data across 20 species

# don't match on common name, since don't see any in the list

# look at the matches
btmams[!is.na(mass_species_bt)| !is.na(mass_species_gbif), ]
btmams[!is.na(mass_species_bt) | !is.na(mass_species_gbif), length(unique(Species))] # 238
btmams[!is.na(mass_species_bt) | !is.na(mass_species_gbif), sort(unique(taxa_mod))] # some are not labeled as mammals by BioTime
btmams[(!is.na(mass_species_bt) | !is.na(mass_species_gbif)) & taxa_mod != 'Mammals', .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # in Birds, Fish, All. But all seem correct.

# compare matches on biotime and gbif names
btmams[, plot(mass_species_bt, mass_species_gbif, log='xy')]; abline(0,1) # on 1:1 line: good!

# combine bt and gbif masses
btmams[, mass_species := mass_species_bt]
btmams[is.na(mass_species_bt) & !is.na(mass_species_gbif), mass_species := mass_species_gbif]
btmams[!is.na(mass_species), length(unique(Species))] # 238. good that it matches from before.

# check for missing
btmams[taxa_mod == 'Mammals' & is.na(mass_species), length(unique(Species))] # 35
btmams[taxa_mod == 'Mammals' & is.na(mass_species), sort(unique(Species))] # mostly species names and genus names.

# get genus-level values
btmams[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))] # not necessarily the genus, if "species name" is not a binomial
mams[, genus := vapply(strsplit(Scientific, " "), `[`, 1, FUN.VALUE=character(1))]
btmamsgenu <- btmams[taxa_mod == 'Mammals' | class_gbif == 'Mammalia', .(genus = sort(unique(c(genus, genus_gbif))))] # list of potential mammal genera (ignore potential mammals in other categories)

btmamsgenu[, mass_genus := NA_real_] 
for(i in 1:nrow(btmamsgenu)){ # a brute-force for loop approach, but fast
    k <- which(btmamsgenu$genus[i] == mams$genus)
    if(length(k) > 0){
        btmamsgenu$mass_genus[i] <- mean(mams$mass[k])
    }
}
btmamsgenu[!is.na(mass_genus),] # found 144 values
btmamsgenu[!is.na(mass_genus), sort(unique(genus))] # check by eye that these look like genera

# merge genus-values with species list
btmams <- merge(btmams, btmamsgenu[, .(genus, mass_genus_bt = mass_genus)], by = 'genus', all.x = TRUE) # match on bt genus
btmams <- merge(btmams, btmamsgenu[, .(genus_gbif = genus, mass_genus_gbif = mass_genus)], by = 'genus_gbif', all.x = TRUE) # match on gbif genus

# examine cases where bt and gbif genus values do not match
btmams[mass_genus_bt != mass_genus_gbif, ][!duplicated(Species), .(Species, Species_gbif, genus, genus_gbif, mass_genus_bt, mass_genus_gbif)] # nothing

# choose from species- or genus-level values
btmams[, mass := mass_species]
btmams[is.na(mass_species), mass := mass_genus_gbif] # use genus where species-level values are NA. prefer gbif
btmams[is.na(mass_species) & is.na(mass_genus_gbif), mass := mass_genus_bt] # bt genus as last choice
btmams[is.nan(mass), mass := NA] # NaN to NA

btmams[!is.na(mass), length(unique(Species))] # 257 species now with data

btmams[mass_genus_gbif != mass_species & !is.na(mass_genus_gbif), plot(mass_genus_gbif, mass_species, log = 'xy')] # compare species and genus. pretty well correlated

# check for missing
btmams[taxa_mod == 'Mammals' & is.na(mass), length(unique(Species))] # 21


# trim to species with data
nrow(btmams) # 1,034,700
btmams.out <- btmams[!is.na(mass), ]
nrow(btmams.out) # 12097

### write out. mass in g
write.csv(btmams.out, file = gzfile('output/mass_eltonmammals.csv.gz'))

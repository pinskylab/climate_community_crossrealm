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
# BioTime change, useful to taxa_mod
load('data/biotime_blowes/bt.Rdata')
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

#########################
# Find BM for birds
#########################
# merge with bird traits
btbirds <- merge(btspp, birds[, .(Species = gsub(' ', '_', Scientific), spp_elton = Scientific, mass_species = mass)], by = 'Species', all.x = TRUE)

# look at the matches
btbirds[!is.na(mass_species), ]
btbirds[!is.na(mass_species), length(unique(Species))] # 982
btbirds[!is.na(mass_species), sort(unique(taxa_mod))] # some are not labeled as birds by BioTime
btbirds[!is.na(mass_species) & taxa_mod != 'Birds', .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # in Benthos, Fish, All, Mammals. But all seem correct.

# check for missing
btbirds[taxa_mod == 'Birds' & is.na(mass_species), length(unique(Species))] # 1013
btbirds[taxa_mod == 'Birds' & is.na(mass_species), sort(unique(Species))] # mostly common names and genus names, some 4-letter codes. Also a seal.

# match on common name
btbirds <- merge(btbirds, birds[, .(Species = gsub(' ', '_', English), spp_elton_comname = Scientific, mass_comname = mass)], by = 'Species', all.x = TRUE)

# look at the matches
btbirds[!is.na(mass_comname), ]
btbirds[!is.na(mass_comname), length(unique(Species))] # 121


# get genus-level values
btbirds[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))] # not necessarily the genus, if "species name" is not a binomial
birds[, genus := vapply(strsplit(Scientific, " "), `[`, 1, FUN.VALUE=character(1))]
btbirdsgenu <- btbirds[!duplicated(genus) & taxa_mod == 'Birds', .(genus)] # list of potential bird genera (ignore potential birds in other categories)
btbirdsgenu[, mass_genus := NA_real_] 

for(i in 1:nrow(btbirdsgenu)){ # a brute-force for loop approach, but fast
    k <- which(btbirdsgenu$genus[i] == birds$genus)
    if(length(k) > 0){
        btbirdsgenu$mass_genus[i] <- mean(birds$mass[k])
    }
}
btbirdsgenu[!is.na(mass_genus),] # found 369 values
btbirdsgenu[!is.na(mass_genus), sort(unique(genus))] # check by eye that these look like genera

# merge genus-values with species list
btbirds <- merge(btbirds, btbirdsgenu, by = 'genus', all.x = TRUE)

# choose from species-, common- or genus-level mass values
btbirds[, mass := mass_species]
btbirds[is.na(mass_species), mass := mass_comname] # use common where species-level value is NA
btbirds[is.na(mass_species) & is.na(mass_comname), mass := mass_genus] # use genus where species- and common-level values are NA
btbirds[is.nan(mass), mass := NA] # NaN to NA

btbirds[!is.na(mass), length(unique(Species))] # 1277 species now with data

btbirds[mass_genus != mass_species & !is.na(mass_genus), plot(mass_genus, mass_species)] # compare species and genus. pretty well correlated

# merge Elton scientific names (from scientific vs. common-name merging)
btbirds[is.na(spp_elton) & !is.na(spp_elton_comname), spp_elton := spp_elton_comname]
btbirds[, spp_elton_comname := NULL]

# trim to species with data
nrow(btbirds) # 1,034,700
btbirds.out <- btbirds[!is.na(mass), ]
nrow(btbirds.out) # 209,663

### write out
write.csv(btbirds.out, file = gzfile('output/mass_eltonbirds.csv.gz'))

#########################
# Find BM for mammals
#########################
# merge with mamals traits
btmams <- merge(btspp, mams[, .(Species = gsub(' ', '_', Scientific), spp_elton = Scientific, mass_species = mass)], by = 'Species', all.x = TRUE)

# look at the matches
btmams[!is.na(mass_species), ]
btmams[!is.na(mass_species), length(unique(Species))] # 218
btmams[!is.na(mass_species), sort(unique(taxa_mod))] # some are not labeled as mammals by BioTime
btmams[!is.na(mass_species) & taxa_mod != 'Mammals', .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # in Birds, Fish, All. But all seem correct.

# check for missing
btmams[taxa_mod == 'Mammals' & is.na(mass_species), length(unique(Species))] # 55
btmams[taxa_mod == 'Mammals' & is.na(mass_species), sort(unique(Species))] # mostly species names and genus names.

# don't match on common name, since don't see any in the list

# get genus-level values
btmams[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))] # not necessarily the genus, if "species name" is not a binomial
mams[, genus := vapply(strsplit(Scientific, " "), `[`, 1, FUN.VALUE=character(1))]
btmamsgenu <- btmams[!duplicated(genus) & taxa_mod == 'Mammals', .(genus)] # list of potential mammal genera (ignore potential mammals in other categories)
btmamsgenu[, mass_genus := NA_real_] 

for(i in 1:nrow(btmamsgenu)){ # a brute-force for loop approach, but fast
    k <- which(btmamsgenu$genus[i] == mams$genus)
    if(length(k) > 0){
        btmamsgenu$mass_genus[i] <- mean(mams$mass[k])
    }
}
btmamsgenu[!is.na(mass_genus),] # found 89 values
btmamsgenu[!is.na(mass_genus), sort(unique(genus))] # check by eye that these look like genera

# merge genus-values with species list
btmams <- merge(btmams, btmamsgenu, by = 'genus', all.x = TRUE)

# choose from species- or genus-level values
btmams[, mass := mass_species]
btmams[is.na(mass_species), mass := mass_genus] # use genus where species-level values are NA
btmams[is.nan(mass), mass := NA] # NaN to NA

btmams[!is.na(mass), length(unique(Species))] # 240 species now with data

btmams[mass_genus != mass_species & !is.na(mass_genus), plot(mass_genus, mass_species)] # compare species and genus. pretty well correlated


# trim to species with data
nrow(btmams) # 1,034,700
btmams.out <- btmams[!is.na(mass), ]
nrow(btmams.out) # 11,604

### write out. mass in g
write.csv(btmams.out, file = gzfile('output/mass_eltonmammals.csv.gz'))

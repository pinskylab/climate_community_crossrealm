################################################################
#######  extract body mass from the Fishbase database  #########
################################################################

###########
# Packages
###########
require(data.table)
require(rfishbase)

############
# Load and merge data
############
# BioTime change, useful to taxa_mod
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt

# BioTime species
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp

# load gbif species names
phylo = fread('output/taxonomy.csv.gz', drop = 1)

# merge spp with gbif names
btspp2 <- merge(btspp, phylo[, .(class_gbif = class, genus_gbif = genus, Species_gbif = species, Species = gsub(' ', '_', original_name))], all.x = TRUE, by = 'Species')


# add taxa_mod to spp list and trim to ones we will search
btspp2 <- merge(btspp2, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], all.x = TRUE, by = 'rarefyID')
nrow(btspp2)

finfish <- c('Actinopterygii', 'Myxini', 'Cephalaspidomorphi', 'Chondrichthyes', 'Elasmobranchii', 'Holocephali') # classes of finfish
btsppfish <- btspp2[class_gbif %in% finfish | is.na(class_gbif)] # from fishbase
nrow(btsppfish) # 964040 rows
btsppfishu <- btsppfish[!duplicated(Species), .(Species, REALM, taxa_mod, class_gbif, genus_gbif, Species_gbif)]
nrow(btsppfishu) # 12800 species

btsppsl <- btspp2[(REALM == 'Marine' & taxa_mod %in% c('All', 'Benthos', 'Invertebrates', 'Marine invertebrates/plants')) |
                     !(class_gbif %in% finfish),] # from sealifebase
nrow(btsppsl) # 1046793 rows
btsppslu <- btsppsl[!duplicated(Species), .(Species, REALM, taxa_mod, class_gbif, genus_gbif, Species_gbif)]
nrow(btsppslu) # 26319 species


#########################
# Find BM from Fishbase
#########################
# validate species names from bt
btsppfishu[, spp_fb := paste0(validate_names(gsub('_', ' ', Species)), collapse = ";"), by = seq_len(NROW(btsppfishu))] # validate names. Separate multiple names with ;
btsppfishu[, sum(spp_fb == '')] # 7938 number of blanks
btsppfishu[, sum(spp_fb != '')] # 4862 number of validated names
btsppfishu[gsub('_', ' ', Species) != spp_fb & spp_fb != '', .(Species, spp_fb)] # check where names don't match exactly. 509 non-blanks
btsppfishu[grep(';', spp_fb), .(Species, spp_fb)] # multiple names returned. 6.

# manually fix the multiple names (only 5)
btsppfishu[, spp_fb2 := spp_fb] # new column with manually fixed FB names (for 5)
btsppfishu[spp_fb2 == '', spp_fb2 := NA]
    bt[rarefyID %in% btspp[Species == 'Epigonus_macrops', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[spp_fb == 'Epigonus robustus;Epigonus macrops', spp_fb2 := 'Epigonus macrops'] # since in NW Atlantic
btsppfishu[spp_fb == 'Hemitripterus americanus;Hemitripterus americanus', spp_fb2 := 'Hemitripterus americanus'] # identical
    bt[rarefyID %in% btspp[Species == 'Centriscus_cristatus', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[spp_fb == 'Centriscus cristatus;Notopogon lilliei', spp_fb2 := 'Centriscus cristatus'] # since in NW Australian shelf
btsppfishu[spp_fb == 'Cociella crocodilus;Cociella crocodilus', spp_fb2 := 'Cociella crocodilus'] # identical
btsppfishu[spp_fb == 'Parajulis poecilepterus;Parajulis poecilepterus', spp_fb2 := 'Parajulis poecilepterus'] # identical
btsppfishu[spp_fb == 'Chaetodon melannotus;Chaetodon melannotus', spp_fb2 := 'Chaetodon melannotus'] # identical

# validate species name from gbif
btsppfishu[, spp_fb_gbif := paste0(validate_names(Species_gbif), collapse = ";"), by = seq_len(NROW(btsppfishu))]
btsppfishu[, sum(spp_fb_gbif == '')] # 7880 number of blanks
btsppfishu[, sum(spp_fb_gbif != '')] # 4920 number of validated names
btsppfishu[, sum(spp_fb_gbif != '' & spp_fb == '')] # 118 new names found by using gbif
btsppfishu[Species_gbif != spp_fb & spp_fb_gbif != '', .(Species, spp_fb_gbif)] # check where names don't match exactly. 164 non-blanks
btsppfishu[grep(';', spp_fb_gbif), .(Species, spp_fb_gbif)] # multiple names returned. 4

# manually fix the multiple names from gbif search
btsppfishu[, spp_fb2_gbif := spp_fb_gbif] # new column with manually fixed FB names (for 4)
btsppfishu[spp_fb2_gbif == '', spp_fb2_gbif := NA]
    bt[rarefyID %in% btspp[Species == 'Epigonus_macrops', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[spp_fb_gbif == 'Epigonus robustus;Epigonus macrops', spp_fb2_gbif := 'Epigonus macrops'] # since in NW Atlantic
btsppfishu[spp_fb_gbif == 'Hemitripterus americanus;Hemitripterus americanus', spp_fb2_gbif := 'Hemitripterus americanus'] # identical
    bt[rarefyID %in% btspp[Species == 'Centriscus_cristatus', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[spp_fb_gbif == 'Centriscus cristatus;Notopogon lilliei', spp_fb2_gbif := 'Centriscus cristatus'] # since in NW Australian shelf
btsppfishu[spp_fb_gbif == 'Cociella crocodilus;Cociella crocodilus', spp_fb2_gbif := 'Cociella crocodilus'] # identical
btsppfishu[grep(';', spp_fb2_gbif), .(Species, spp_fb2_gbif)] # multiple names fixed


# combine bt and gbif names
btsppfishu[!is.na(spp_fb2) & !is.na(spp_fb2_gbif) & spp_fb2 != spp_fb2_gbif, .(Species, genus_gbif, Species_gbif, spp_fb2, spp_fb2_gbif)] # only 7 where bt and gbif do not match after running through fishbase validation

btsppfishu[, spp_fb3 := spp_fb2_gbif] # new column with manually fixed FB names (for 7)
btsppfishu[is.na(spp_fb3) & !is.na(spp_fb2), spp_fb3 := spp_fb2]
btsppfishu[spp_fb3 == '', spp_fb3 := NA]

bt[rarefyID %in% btspp2[Species == 'Inopsetta_ischyra', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[Species == 'Inopsetta_ischyra', spp_fb3 := spp_fb2_gbif] # since in Arctic as well

bt[rarefyID %in% btspp2[Species == 'Epinephelus_striatus', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[Species == 'Epinephelus_striatus', spp_fb3 := spp_fb2] # since in NW Atlantic

bt[rarefyID %in% btspp2[Species == 'Stomias_ferox', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[Species == 'Stomias_ferox', spp_fb3 := spp_fb2] # prefer FB over GBIF

bt[rarefyID %in% btspp2[Species == 'Apogon_striatus', rarefyID], ][!duplicated(rarefyID),] # check location. none returned from gridded BT dataset
btsppfishu[Species == 'Apogon_striatus', spp_fb3 := spp_fb2] # prefer FB over GBIF

bt[rarefyID %in% btspp2[Species == 'Platybelone_platyura', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[Species == 'Platybelone_platyura', spp_fb3 := spp_fb2] # since NW pacific

bt[rarefyID %in% btspp2[Species == 'Tylosurus_melanotus', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[Species == 'Tylosurus_melanotus', spp_fb3 := spp_fb2] # since NW pacific

bt[rarefyID %in% btspp2[Species == 'Hemigymnus_fasciatus', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[Species == 'Hemigymnus_fasciatus', spp_fb3 := spp_fb2] # prefer FB

btsppfishu[is.na(spp_fb3) & (!is.na(spp_fb2) | !is.na(spp_fb2_gbif)), ] # fb3 not missing any: good

# get species-level masses
btsppfishu[!is.na(spp_fb3), mass_species := as.numeric(species(spp_fb3, fields=c("Weight"))), by = seq_len(NROW(btsppfishu[!is.na(spp_fb3)]))] # in grams



# get genus-level values
btsppfishu[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))] # from the biotime name
btsppfishu[, genus_fb3 := vapply(strsplit(spp_fb3, " "), `[`, 1, FUN.VALUE=character(1))] # from the fishbase-processed name
btsppfishgenu <- btsppfishu[!is.na(mass_species), ][!duplicated(cbind(genus, genus_fb3)), .(genus, genus_gbif, genus_fb3)] # list of genera where we don't have a species value
    btsppfishgenu[!is.na(genus) & !is.na(genus_fb3) & genus != genus_fb3,]
nrow(btsppfishgenu) # 507

getgenmass <- function(x){ # function to return average weight of species in a genus
    spps <- species_list(Genus = x)
    mean(species(spps, fields = 'Weight')$Weight, na.rm=TRUE)
}
btsppfishgenu[, mass_genus_bt := getgenmass(genus), by = seq_len(NROW(btsppfishgenu))] # find based on Biotime genus
btsppfishgenu[!is.na(mass_genus_bt),] # found 481 values

btsppfishgenu[!is.na(genus_fb3), mass_genus_fb3 := getgenmass(genus_fb3), by = seq_len(NROW(btsppfishgenu[!is.na(genus_fb3),]))] # based on fishbase-checked name
btsppfishgenu[!is.na(mass_genus_fb3),] # found 507 values
btsppfishgenu[!is.na(mass_genus_fb3) & is.na(mass_genus_bt),] # found 26 new values
btsppfishgenu[is.na(mass_genus_fb3) & !is.na(mass_genus_bt),] # found all biotime genus values as well

# merge genus-level values
btsppfishgenu[mass_genus_bt != mass_genus_fb3,] # 26 do not match
btsppfishgenu[, mass_genus := mass_genus_fb3] # prefer the fishbase-checked genera
btsppfishgenu[is.na(mass_genus) & !is.na(mass_genus_bt), mass_genus := mass_genus_bt] # prefer the fishbase-checked genera


# merge genus-values with species list
btsppfishu2 <- merge(btsppfishu, btsppfishgenu[, .(genus, genus_fb3, mass_genus)], by = c('genus', 'genus_fb3'), all.x = TRUE)

# choose from species- or genus-level values
btsppfishu2[, mass := mass_species]
btsppfishu2[is.na(mass_species), mass := mass_genus] # use genus where species-level value is NA
btsppfishu2[is.nan(mass), mass := NA] # NaN to NA

btsppfishu2[is.na(mass), length(unique(Species))] # 10876 species now with data

btsppfishu2[mass_genus != mass_species, plot(mass_genus, mass_species, log='xy')] # compare species and genus. correlated.


# merge back with full list of species and studies
# check this!
btsppfish <- merge(btsppfish, btsppfishu2[, .(Species, spp_fb3, mass_genus, mass_species, mass)], by = 'Species', all.x = TRUE)

# trim to species with data
nrow(btsppfish)
btsppfish.out <- btsppfish[!is.na(mass), ]
nrow(btsppfish.out) # 515055 rows

### write out
write.csv(btsppfish.out, file = gzfile('output/mass_fishbase.csv.gz'))
          



#########################
# Find BM from Sealifebase
#########################
# validate species names
btsppslu[, spp_fb := paste0(unique(validate_names(gsub('_', ' ', Species), server = 'sealifebase')), 
                            collapse = ";"), by = seq_len(NROW(btsppslu))] # validate names. Separate multiple names with ;. Slow.
btsppslu[, sum(spp_fb == '')] # 18906 number of blanks
btsppslu[, sum(spp_fb != '')] # 7413 number of validated names
btsppslu[gsub('_', ' ', Species) != spp_fb & spp_fb != '' & spp_fb != 'NA', .(Species, spp_fb)] # check where names don't match exactly. 518 non-blanks


# manually fix the multiple names
# skip this and read in file of choices if already done
# for each displayed name, go to https://www.sealifebase.ca/summary/Genus-species.html for the options and make a decision based on location
# if no obvious choice, pick name closest to the BioTime name
# tofix <- btsppslu[grep(';', spp_fb), .(Species, spp_fb)]; tofix # multiple names returned. 21
# tofix[, spp_fb2 := NA_character_]
# 
# for(i in 1:nrow(tofix)){
#     print(paste0("Studies with ", tofix$Species[i], " were conducted in these locations:"))
#     print(bt[rarefyID %in% btspp[Species == tofix$Species[i], rarefyID], ][!duplicated(rarefyID), .(Biome, STUDY_ID, rarefyID_x, rarefyID_y)]) # check location
#     sppchoices <- unlist(strsplit(tofix$spp_fb[i], split = ';'))
#     print(paste0('For ', tofix$Species[i], ", choose ", paste(1:length(sppchoices), sppchoices, collapse = ',')))
#     n <- NA
#     while(is.na(n)){
#         n <- as.integer(readline(prompt = "Pick one: "))
#         if(!is.na(n)) tofix$spp_fb2[i] <- sppchoices[n]
#     }
# }
# tofix # examine choices made

# write out choices made for multiple names
# write.csv(tofix, file = 'data/biotime_blowes/name_validation_sealifebase.csv')

# validate species names from gbif
btsppslu[, spp_fb_gbif := paste0(unique(validate_names(Species_gbif, server = 'sealifebase')), 
                                 collapse = ";"), by = seq_len(NROW(btsppslu))] # validate names. Separate multiple names with ;. Slow.
btsppslu[, sum(spp_fb_gbif == '')] # 19592 number of blanks
btsppslu[, sum(spp_fb_gbif != '')] # 6727 number of validated names
btsppslu[, sum(spp_fb_gbif != '' & spp_fb == '')] # 302 new names found by gbif
btsppslu[Species_gbif != spp_fb_gbif & spp_fb_gbif != '' & spp_fb_gbif != 'NA', .(Species, spp_fb)] # check where names don't match exactly. 143 non-blanks


# manually fix the multiple names returned for gbif names
# skip this and read in file of choices if already done
# for each displayed name, go to https://www.sealifebase.ca/summary/Genus-species.html for the options and make a decision based on location
# if no obvious choice, pick name closest to the BioTime name
# tofix_gbif <- btsppslu[grep(';', spp_fb_gbif), .(Species, spp_fb_gbif)]; tofix_gbif # multiple names returned. 10
# tofix_gbif[, spp_fb2_gbif := NA_character_]
# 
# btsppslu[, spp_fb2_gbif := spp_fb_gbif] # new column with manually fixed FB names
# btsppslu[spp_fb2_gbif == '', spp_fb2_gbif := NA]
# 
# for(i in 1:nrow(tofix_gbif)){
#     print(paste0("Studies with ", tofix_gbif$Species[i], " were conducted in these locations:"))
#     print(bt[rarefyID %in% btspp[Species == tofix_gbif$Species[i], rarefyID], ][!duplicated(rarefyID), .(Biome, STUDY_ID, rarefyID_x, rarefyID_y)]) # check location
#     sppchoices <- unlist(strsplit(tofix_gbif$spp_fb_gbif[i], split = ';'))
#     print(paste0('For ', tofix_gbif$Species[i], ", choose ", paste(1:length(sppchoices), sppchoices, collapse = ',')))
#     n <- NA
#     while(is.na(n)){
#         n <- as.integer(readline(prompt = "Pick one: "))
#         if(!is.na(n)) tofix_gbif$spp_fb2_gbif[i] <- sppchoices[n]
#     }
# }
# tofix_gbif # examine choices made


# write out choices made for multiple names
# write.csv(tofix_gbif, file = 'data/biotime_blowes/name_gbif_validation_sealifebase.csv')


# read in name choices made to resolve multiple valid FB names
tofix <- fread('data/biotime_blowes/name_validation_sealifebase.csv', header = TRUE)
tofix_gbif <- fread('data/biotime_blowes/name_gbif_validation_sealifebase.csv', header = TRUE)

# resolve multiple names from FB
btsppslu[, spp_fb2 := spp_fb] # new column with manually fixed FB names
btsppslu[spp_fb2 == '', spp_fb2 := NA]
for(i in 1:nrow(tofix)){
    btsppslu[spp_fb == tofix$spp_fb[i], spp_fb2 := tofix$spp_fb2[i]]
}

btsppslu[, spp_fb2_gbif := spp_fb_gbif] # new column with manually fixed FB names
btsppslu[spp_fb2_gbif == '', spp_fb2_gbif := NA]
for(i in 1:nrow(tofix_gbif)){
    btsppslu[spp_fb_gbif == tofix_gbif$spp_fb_gbif[i], spp_fb2_gbif := tofix_gbif$spp_fb2_gbif[i]]
}

btsppslu[grep(';', spp_fb2), .(Species, spp_fb2)] # 0: no more multiple names: good
btsppslu[grep(';', spp_fb2_gbif), .(Species, spp_fb2_gbif)] # also 0


# combine bt and gbif names
btsppslu[!is.na(spp_fb2) & !is.na(spp_fb2_gbif) & spp_fb2 != spp_fb2_gbif, .(Species, genus_gbif, Species_gbif, spp_fb2, spp_fb2_gbif)] # 447 where bt and gbif do not match after running through fishbase validation

btsppslu[, spp_fb3 := spp_fb2] # new column with combined name. Prefer the direct BT -> sealifebase name, without GBIF in the middle
btsppslu[is.na(spp_fb2) & !is.na(spp_fb2_gbif), .N] # gbif adds 302 names
btsppslu[is.na(spp_fb2) & !is.na(spp_fb2_gbif), spp_fb3 := spp_fb2_gbif]
btsppslu[spp_fb3 == '', spp_fb3 := NA]

btsppfishu[is.na(spp_fb3) & (!is.na(spp_fb2) | !is.na(spp_fb2_gbif)), ] # fb3 not missing any: good


# get species-level masses
btsppslu[!is.na(spp_fb3), mass_species := as.numeric(species(spp_fb3, fields=c("Weight"), server = 'sealifebase')), by = seq_len(NROW(btsppslu[!is.na(spp_fb3)]))] # in grams
btsppslu[, sum(!is.na(mass_species))] # 319

# get genus-level values
btsppslu[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))] # from the biotime name
btsppslu[, genus_fb3 := vapply(strsplit(spp_fb3, " "), `[`, 1, FUN.VALUE=character(1))] # from the sealifebase-processed name
btsppslgenu <- btsppslu[!is.na(mass_species), ][!duplicated(cbind(genus, genus_fb3)), .(genus, genus_gbif, genus_fb3)] # list of genera where we don't have a species value
    btsppslgenu[!is.na(genus) & !is.na(genus_fb3) & genus != genus_fb3,] # where biotime and fishbase genera differ
nrow(btsppslgenu) # 147


getgenmassSLB <- function(x){ # function to return average weight of species in a genus
    spps <- species_list(Genus = x, server = 'sealifebase')
    mean(species(spps, fields = 'Weight', server = 'sealifebase')$Weight, na.rm=TRUE)
}
btsppslgenu[, mass_genus_bt := getgenmassSLB(genus), by = seq_len(NROW(btsppslgenu))]
btsppslgenu[!is.na(mass_genus_bt),] # 143 values

btsppslgenu[!is.na(genus_fb3), mass_genus_fb3 := getgenmassSLB(genus_fb3), by = seq_len(NROW(btsppslgenu[!is.na(genus_fb3),]))] # based on fishbase-checked name
btsppslgenu[!is.na(mass_genus_fb3),] # found 147 values
btsppslgenu[!is.na(mass_genus_fb3) & is.na(mass_genus_bt),] # found 4 new values
btsppslgenu[is.na(mass_genus_fb3) & !is.na(mass_genus_bt),] # found all biotime genus values as well

# merge genus-level values (from BT genus or from FB-checked genus)
btsppslgenu[mass_genus_bt != mass_genus_fb3,] # 13 do not match
btsppslgenu[, mass_genus := mass_genus_fb3] # prefer the fishbase-checked genera
btsppslgenu[is.na(mass_genus) & !is.na(mass_genus_bt), mass_genus := mass_genus_bt] # then try the BT genera


# merge genus-values with species list
btsppslu2 <- merge(btsppslu, btsppslgenu[, .(genus, genus_fb3, mass_genus)], by = c('genus', 'genus_fb3'), all.x = TRUE)

# choose from species- or genus-level values
btsppslu2[, mass := mass_species]
btsppslu2[is.na(mass_species) | mass_species == 0, mass := mass_genus] # use genus where species-level value is NA or 0
btsppslu2[is.nan(mass), mass := NA] # NaN to NA

btsppslu2[is.na(mass), length(unique(Species))] # 25954 species now with data

btsppslu2[ mass <= 0,] # 0, good

btsppslu2[mass_genus != mass_species, plot(mass_genus, mass_species, log='xy')] # compare species and genus


# merge back with full list of species and studies
btsppsl <- merge(btsppsl, btsppslu2[, .(Species, spp_fb3, mass_genus, mass_species, mass)], by = 'Species', all.x = TRUE)

# trim to species with data
nrow(btsppsl) # 1046793
btsppsl.out <- btsppsl[!is.na(mass), ]
nrow(btsppsl.out) # 174532 rows

### write out. mass in g
write.csv(btsppsl.out, file = gzfile('output/mass_sealifebase.csv.gz'))

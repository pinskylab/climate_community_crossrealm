# Decide whether species are ectotherms, birds, or mammals (for body temp calculations)
# Also classify as producer vs. consumer
# Then calculation fraction of each rarefyID, weighting by abundance

require(data.table)
#################
# Read in data
#################
# biotime community similarity data
load('data/biotime_blowes/bt_malin.Rdata') # load bt_malin
bt <- data.table(bt_malin); rm(bt_malin)

# biotime species lists
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID') # add taxa_mod to spp list

btspp[, length(unique(Species))] # 26289 species

# Add abundance information
load('data/biotime_blowes/bt_grid_spp_list_abund.Rdata') # loads abundance data (bt_grid_spp_list_abund, bt_grid_spp_list_abund_year)
btabund <- as.data.table(bt_grid_spp_list_abund) # rename 
rm(bt_grid_spp_list_abund, bt_grid_spp_list_abund_year) # delete unneeded df

# phylogeny info
phylo <- fread('output/taxonomy.csv.gz')
phylo[, ':='(gateway.masses = NULL, V1 = NULL, X = NULL, NA. = NULL)] # remove extraneous columns


###########
# Prep data
###########
# merge abundance and bt
btspp <- merge(btspp, btabund[, .(rarefyID, Species, N_bar)], by = c('rarefyID', 'Species'), all.x = TRUE)

# merge taxonomy and bt
phylo[, Species := gsub(' ', '_', original_name)] # add _ back to name
endo <- merge(btspp[, .(rarefyID, Species, REALM, STUDY_ID, taxa_mod, N_bar)], phylo, 
              all.x = TRUE, by = "Species")

##############
# Classify
##############
# classify by bird vs. mammal vs. ectotherm. Assume ecto if not otherwise known.
endo[, metab := NA_character_]
endo[class == 'Aves', metab := 'bird']
endo[class == 'Mammalia', metab := 'mamm']
endo[is.na(metab) & is.na(class) & taxa_mod == 'Birds', metab := 'bird'] # use taxa_mod if phylum not available
endo[is.na(metab) & is.na(class) & taxa_mod == 'Mammals', metab := 'mamm']
endo[is.na(metab), metab := 'ecto']

# examine metabolic classification
endo[, table(taxa_mod, metab)]

endo[taxa_mod == 'Mammals' & metab != 'mamm', ] # seem correct
endo[taxa_mod == 'Birds' & metab != 'bird', ][!duplicated(Species)] # seem correct
endo[taxa_mod == 'Fish' & metab != 'ecto', ][!duplicated(Species)] # seem correct
endo[taxa_mod == 'Plant' & metab != 'ecto', ][!duplicated(Species)] # wrong

# fix known errors
endo[Species == 'Phyllanthus_emblica', metab := 'ecto']


# classify as consumer vs. producer
endo[, feeding := NA_character_]
endo[kingdom %in% c('Animalia', 'Bacteria', 'Fungi', 'Protozoa'), feeding := 'consumer']
endo[kingdom == 'Bacteria' & phylum == 'Cyanobacteria', feeding := 'producer'] # photosynthetic bacteria
endo[kingdom %in% c('Plantae', 'Chromista'), feeding := 'producer']
endo[is.na(feeding) & taxa_mod %in% c('Amphibians', 'Birds', 'Fish', 'Invertebrates', 'Mammals', 'Reptiles'), feeding := 'consumer']
endo[is.na(feeding) & taxa_mod %in% c('Plant'), feeding := 'producer']

# examine feeding classification
endo[, table(taxa_mod, feeding)]

endo[taxa_mod == 'Birds' & feeding != 'consumer', ][!duplicated(Species)] # seem correct
endo[taxa_mod == 'Invertebrates' & feeding != 'consumer', ][!duplicated(Species)] # seem correct
endo[taxa_mod == 'Plant' & feeding != 'producer', ][!duplicated(Species)] # seem correct


# write out
write.csv(endo[, .(Species, rarefyID, metab, feeding)], file = gzfile('output/metab_feed_byspecies.csv.gz'), row.names = FALSE)


#############################################
# Average endo/ecto within rarefyIDs
#############################################

# fraction endotherms by rarefyID
endofrac <- endo[, .(N_sum = sum(N_bar, na.rm = TRUE),
                     nspp = length(unique(Species)), nspp_wdata = sum(!is.na(metab))), 
                 by = .(STUDY_ID, rarefyID, REALM, taxa_mod)]
endofrac <- merge(endofrac, endo[metab %in% c('bird', 'mamm'), 
                                 .(N_sum_endo = sum(N_bar, na.rm = TRUE),
                                   nspp_endo = length(unique(Species))), 
                                 by = .(rarefyID)], by = 'rarefyID', all.x = TRUE)
endofrac[is.na(N_sum_endo), N_sum_endo := 0] # if no endotherms, set sum of endotherm abundance to 0
endofrac[is.na(nspp_endo), nspp_endo := 0] # if no endotherms, set # of endotherm species to 0
endofrac[, endofrac := N_sum_endo/N_sum]
endofrac[is.na(endofrac), endofrac := nspp_endo/nspp] # if no abundance data, assume all species have the same abundance
nrow(endofrac) # 53467
setkey(endofrac, STUDY_ID, rarefyID)
endofrac
summary(endofrac)

# write out
write.csv(endofrac, gzfile('output/endofrac_byrarefyID.csv.gz'), row.names = FALSE)


#############################################
# Average consumer/producer within rarefyIDs
#############################################

# fraction consumers by rarefyID
consfrac <- endo[, .(N_sum = sum(N_bar, na.rm = TRUE),
                   nspp = length(unique(Species)), nspp_wdata = sum(!is.na(feeding))), 
               by = .(STUDY_ID, rarefyID, REALM, taxa_mod)]
consfrac <- merge(consfrac, endo[feeding == 'consumer', 
                                 .(N_sum_cons = sum(N_bar, na.rm = TRUE),
                                   nspp_cons = length(unique(Species))), 
                                 by = .(rarefyID)], by = 'rarefyID', all.x = TRUE)
consfrac[is.na(N_sum_cons), N_sum_cons := 0] # if no consumers, set sum of consumer abundance to 0
consfrac[is.na(nspp_cons), nspp_cons := 0] # if no consumers, set # of consumer species to 0
consfrac[, consfrac := N_sum_cons/N_sum]
consfrac[is.na(consfrac), consfrac := nspp_cons/nspp] # if no abundance data, assume all species have the same abundance
nrow(consfrac) # 53467
setkey(consfrac, STUDY_ID, rarefyID)
consfrac


# write out
write.csv(consfrac, gzfile('output/consfrac_byrarefyID.csv.gz'), row.names = FALSE)

# Decide whether species are ectotherms, birds, or mammals (for body temp calculations)
# Also classify as producer vs. consumer

require(data.table)

# biotime community similarity data
load('data/biotime_blowes/bt_malin.Rdata') # load bt_malin
bt <- data.table(bt_malin); rm(bt_malin)

# biotime species lists
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID') # add taxa_mod to spp list

btspp[, length(unique(Species))] # 26289 species

# phylogeny info
phylo <- fread('output/taxonomy.csv.gz')
phylo[, ':='(gateway.masses = NULL, V1 = NULL, X = NULL, NA. = NULL)] # remove extraneous columns



# merge taxonomy and bt
phylo[, Species := gsub(' ', '_', original_name)] # add _ back to name
endo <- merge(btspp[, .(rarefyID, Species, REALM, STUDY_ID, taxa_mod)], phylo, 
              all.x = TRUE, by = "Species")

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

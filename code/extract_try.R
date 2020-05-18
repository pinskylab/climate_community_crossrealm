################################################################
#######  extract body mass from the TRY Traits database    #####
######   plants                                            #####
################################################################

###########
# Packages
###########
require(data.table)
require(stringi)

############
# Load and merge data
############
# BioTime change, useful to taxa_mod
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt

# BioTime species
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp

# add taxa_mod to spp list
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID')
nrow(btspp) # 1034700
btspp[, genus := tolower(vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1)))] # genus name in lower case
btspp[, speciesnm := tolower(vapply(strsplit(Species, "_"), `[`, 2, FUN.VALUE=character(1)))] # genus name in lower case

# load gbif species names
phylo = fread('output/taxonomy.csv.gz', drop = 1)

# merge spp with gbif names
btspp2 <- merge(btspp, phylo[, .(class_gbif = class, genus_gbif = tolower(genus), Species_gbif = tolower(species), Species = gsub(' ', '_', original_name))], all.x = TRUE, by = 'Species')
btspp2[, speciesnm_gbif := tolower(vapply(strsplit(Species_gbif, " "), `[`, 2, FUN.VALUE=character(1)))]

# load TRY traits plants
trypl <- fread('dataDL/try/8916.txt', encoding = 'Latin-1')


# clean species names from TRY traits
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), AccSpeciesName] # examine the non-ASCII ones
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), stri_enc_mark(AccSpeciesName)] # "latin1" encoding
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), stri_trans_general(AccSpeciesName, "Latin-ASCII")]

trypl[, AccSpeciesName_clean := AccSpeciesName] # create a clean species name
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), AccSpeciesName_clean := gsub("[^[:alnum:][:blank:]?&/\\-]", "x", AccSpeciesName)] # remove non-Latin characters
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), AccSpeciesName_clean := gsub("\\s+", " ", AccSpeciesName_clean)] # remove double spaces
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), AccSpeciesName_clean := gsub(" x ", " x", AccSpeciesName_clean)] # turn x to an x before species name
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), AccSpeciesName_clean] # examine the non-ASCII ones
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), AccSpeciesName_clean := stri_encode(AccSpeciesName_clean, to = "ASCII")] # encode as ASCII
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), AccSpeciesName_clean] # examine the non-ASCII ones
trypl[which(stri_enc_mark(AccSpeciesName)!='ASCII'), stri_enc_mark(AccSpeciesName_clean)] # encoding is now ASCII


###############################
# Find mass for plants in TRY
################################
# examine (optional)
trypl[, sort(unique(TraitName))] # "Dispersal distance", "Dispersal syndrome", "Plant biomass and allometry: Plant dry mass", "Plant biomass and allometry: Plant mass", "Plant fresh mass", "Plant lifespan (longevity)"         
trypl[, unique(TraitName), by = TraitID]

trypl[TraitName == "Plant biomass and allometry: Plant mass", sort(unique(OrigValueStr))][1:100] # numbers as strings
trypl[TraitID == 3451, sum(is.na(as.numeric(unique(OrigValueStr))))] # 0 convert to NA
trypl[TraitID == 3451, sort(unique(OrigUnitStr))] # all kg
trypl[TraitID == 3451, sort(unique(ValueKindName))] # all Single
trypl[TraitID == 3451, sort(unique(StdValue))] # numbers! These have been standardized by TRY. Use this
trypl[TraitID == 3451, sort(unique(UnitName))] # all g
trypl[TraitID == 3451, .(AccSpeciesName_clean, StdValue)][!duplicated(AccSpeciesName_clean), sum(!is.na(StdValue))] # 277 species
trypl[TraitID == 3451, .(AccSpeciesName_clean, OrigValueStr)][!duplicated(AccSpeciesName_clean), sum(!is.na(OrigValueStr))] # 277 species

trypl[TraitName == "Plant biomass and allometry: Plant dry mass", sort(unique(OrigValueStr))][1:100] # numbers as strings
trypl[TraitID == 700, sum(is.na(as.numeric(unique(OrigValueStr))))] # 2 convert to NA
trypl[TraitID == 700, .(val = unique(OrigValueStr))][is.na(as.numeric(val)), ] # 2 convert to NA: <NA> and .
trypl[TraitID == 700 & OrigValueStr %in% c('<NA>', '.'), ] # look at them
trypl[TraitID == 700, sort(unique(OrigUnitStr))] # g or mg
trypl[TraitID == 700, sort(unique(ValueKindName))] # Mean or Single
trypl[TraitID == 700, sort(unique(StdValue))] # numbers! These have been standardized by TRY. Use this
trypl[TraitID == 700, sort(unique(UnitName))] # all g
trypl[TraitID == 700, .(AccSpeciesName_clean, StdValue)][!duplicated(AccSpeciesName_clean), sum(!is.na(StdValue))] # 163 species. Not 238 as reported by TRY
trypl[TraitID == 700, .(AccSpeciesName_clean, OrigValueStr)][!duplicated(AccSpeciesName_clean), sum(!is.na(OrigValueStr))] # 163 species

trypl[TraitName == "Plant fresh mass", sort(unique(OrigValueStr))][1:100] # numbers as strings
trypl[TraitID == 3407, sum(is.na(as.numeric(unique(OrigValueStr))))] # 0 convert to NA
trypl[TraitID == 3407, sort(unique(OrigUnitStr))] # g
trypl[TraitID == 3407, sort(unique(ValueKindName))] # Single
trypl[TraitID == 3407, sort(unique(StdValue))] # numbers! These have been standardized by TRY. Use this
trypl[TraitID == 3407, sort(unique(UnitName))] # all g
trypl[TraitID == 3407, .(AccSpeciesName_clean, StdValue)][!duplicated(AccSpeciesName_clean), sum(!is.na(StdValue))] # 20 species.
trypl[TraitID == 3407, .(AccSpeciesName_clean, OrigValueStr)][!duplicated(AccSpeciesName_clean), sum(!is.na(OrigValueStr))] # 20 species

# summarize to one trait per species
trypl[, genus := tolower(vapply(strsplit(AccSpeciesName_clean, " "), `[`, 1, FUN.VALUE=character(1)))] # genus name in lower case
trypl[, speciesnm := tolower(vapply(strsplit(AccSpeciesName_clean, " "), `[`, 2, FUN.VALUE=character(1)))] # species

trymass <- trypl[TraitID %in% c(700, 3407, 3451), .(mean = mean(StdValue, na.rm = TRUE), n = sum(!is.na(StdValue))), 
                 by = .(TraitID, TraitName, genus, speciesnm)]
trymass[is.nan(mean), mean := NA] # convert nan to NA
trymass[!is.na(mean) & TraitID == 700, ] # examine. 163 species for dry mass
trymass[!is.na(mean) & TraitID == 3407, ] # examine. 20 species for fresh mass
trymass[!is.na(mean) & TraitID == 3451, ] # examine. 277 species for mass

# calculate genera averages
trymassgenu <- trymass[, .(mean = mean(mean), n = sum(n)), by = .(TraitID, TraitName, genus)] # genus-level averages

#####################################
# merge with plant traits by species
#####################################
# merge based on biotime name
btplants <- merge(btspp2, trymass[TraitID == 700, .(genus, speciesnm, mass_dry_species_bt = mean)], by = c("genus", "speciesnm"), all.x = TRUE)
btplants <- merge(btplants, trymass[TraitID == 3407, .(genus, speciesnm, mass_fresh_species_bt = mean)], by = c("genus", "speciesnm"), all.x = TRUE)
btplants <- merge(btplants, trymass[TraitID == 3451, .(genus, speciesnm, mass_species_bt = mean)], by = c("genus", "speciesnm"), all.x = TRUE)

# merge based on gbif name
btplants <- merge(btplants, trymass[TraitID == 700, .(genus_gbif = genus, speciesnm_gbif = speciesnm, mass_dry_species_gbif = mean)], by = c("genus_gbif", "speciesnm_gbif"), all.x = TRUE)
btplants <- merge(btplants, trymass[TraitID == 3407, .(genus_gbif = genus, speciesnm_gbif = speciesnm, mass_fresh_species_gbif = mean)], by = c("genus_gbif", "speciesnm_gbif"), all.x = TRUE)
btplants <- merge(btplants, trymass[TraitID == 3451, .(genus_gbif = genus, speciesnm_gbif = speciesnm, mass_species_gbif = mean)], by = c("genus_gbif", "speciesnm_gbif"), all.x = TRUE)


# look at the matches
btplants[!is.na(mass_dry_species_bt), ]
btplants[!is.na(mass_fresh_species_bt), ]
btplants[!is.na(mass_species_bt), ]
btplants[!is.na(mass_dry_species_gbif), ]
btplants[!is.na(mass_fresh_species_gbif), ]
btplants[!is.na(mass_species_gbif), ]

btplants[!is.na(mass_dry_species_bt), length(unique(Species))] # 71 species
btplants[!is.na(mass_fresh_species_bt), length(unique(Species))] # 16
btplants[!is.na(mass_species_bt), length(unique(Species))] # 67
btplants[!is.na(mass_dry_species_gbif), length(unique(Species))] # 72 species
btplants[!is.na(mass_fresh_species_gbif), length(unique(Species))] # 14
btplants[!is.na(mass_species_gbif), length(unique(Species))] # 73

btplants[!is.na(mass_dry_species_gbif) & is.na(mass_dry_species_bt), length(unique(Species))] # 4 new species from gbif
btplants[!is.na(mass_fresh_species_gbif) & is.na(mass_fresh_species_bt), length(unique(Species))] # 0 new
btplants[!is.na(mass_species_gbif) & is.na(mass_species_bt), length(unique(Species))] # 8 new

btplants[mass_dry_species_bt != mass_dry_species_gbif, ] # 0 cases where bt and gbif do not match
btplants[mass_fresh_species_bt != mass_fresh_species_gbif, ] # 0
btplants[mass_species_bt != mass_species_gbif, ] # 0


btplants[!is.na(mass_dry_species_bt) | !is.na(mass_fresh_species_bt) | !is.na(mass_species_bt) |
           !is.na(mass_dry_species_gbif) | !is.na(mass_fresh_species_gbif) | !is.na(mass_species_gbif), sort(unique(taxa_mod))] # some are not labeled as Plant by BioTime
btplants[(!is.na(mass_dry_species_bt) | !is.na(mass_fresh_species_bt) | !is.na(mass_species_bt) |
            !is.na(mass_dry_species_gbif) | !is.na(mass_fresh_species_gbif) | !is.na(mass_species_gbif)) & taxa_mod != 'Plant', 
         .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # one in Benthos. Looks correct.

# check for missing
btplants[taxa_mod == 'Plant' & (is.na(mass_dry_species_bt) & is.na(mass_fresh_species_bt) & is.na(mass_species_bt) &
                                  is.na(mass_dry_species_gbif) & is.na(mass_fresh_species_gbif) & is.na(mass_species_gbif)), length(unique(Species))] # 7280
btplants[taxa_mod == 'Plant' & (is.na(mass_dry_species_bt) & is.na(mass_fresh_species_bt) & is.na(mass_species_bt) &
                                  is.na(mass_dry_species_gbif) & is.na(mass_fresh_species_gbif)), 
         sort(unique(Species))] # mostly species names and genus names. Some common names.


# merge genus-values with species list based on biotime genus
btplants <- merge(btplants, trymassgenu[TraitID == 700, .(genus, mass_dry_genus_bt = mean)], by = 'genus', all.x = TRUE)
btplants <- merge(btplants, trymassgenu[TraitID == 3407, .(genus, mass_fresh_genus_bt = mean)], by = 'genus', all.x = TRUE)
btplants <- merge(btplants, trymassgenu[TraitID == 3451, .(genus, mass_genus_bt = mean)], by = 'genus', all.x = TRUE)

# merge genus-values with species list based on gbif genus
btplants <- merge(btplants, trymassgenu[TraitID == 700, .(genus_gbif = genus, mass_dry_genus_gbif = mean)], by = 'genus_gbif', all.x = TRUE)
btplants <- merge(btplants, trymassgenu[TraitID == 3407, .(genus_gbif = genus, mass_fresh_genus_gbif = mean)], by = 'genus_gbif', all.x = TRUE)
btplants <- merge(btplants, trymassgenu[TraitID == 3451, .(genus_gbif = genus, mass_genus_gbif = mean)], by = 'genus_gbif', all.x = TRUE)


# examine genus values
btplants[!is.na(mass_dry_genus_bt), .(Species, genus, mass_dry_species_bt, mass_dry_genus_bt)]
btplants[!is.na(mass_fresh_genus_bt), .(Species, genus, mass_fresh_species_bt, mass_fresh_genus_bt)]
btplants[!is.na(mass_genus_bt), .(Species, genus, mass_species_bt, mass_genus_bt)]

btplants[!is.na(mass_dry_genus_gbif), .(Species, genus, mass_dry_species_gbif, mass_dry_genus_gbif)]
btplants[!is.na(mass_fresh_genus_gbif), .(Species, genus, mass_fresh_species_gbif, mass_fresh_genus_gbif)]
btplants[!is.na(mass_genus_gbif), .(Species, genus, mass_species_gbif, mass_genus_gbif)]

btplants[!is.na(mass_dry_genus_bt) & is.na(mass_dry_species_bt), .N] # genus adds 937 values
btplants[!is.na(mass_fresh_genus_bt) & is.na(mass_fresh_species_bt), .N] # 145
btplants[!is.na(mass_genus_bt) & is.na(mass_species_bt), .N] # 1726

btplants[!is.na(mass_dry_genus_gbif) & is.na(mass_dry_genus_bt), .N] # gbif adds 40 genera
btplants[!is.na(mass_fresh_genus_gbif) & is.na(mass_fresh_genus_bt), .N] # 5
btplants[!is.na(mass_genus_gbif) & is.na(mass_genus_bt), .N] # 76

btplants[mass_dry_genus_bt != mass_dry_genus_gbif, ] # 2 cases (2 species in 1 genus) where BT and GBIF genus masses do not match
btplants[mass_fresh_genus_bt != mass_fresh_genus_gbif, ] # 0
btplants[mass_genus_bt != mass_genus_gbif, ] # 5

# choose from species- or genus-level mass values and from BioTime vs. GBIF names. 
# Prefer GBIF species, then BT species, then GBIF genus
btplants[, mass_dry := mass_dry_species_gbif] # prefer GBIF species
btplants[is.na(mass_dry), mass_dry := mass_dry_species_bt] # then BT species
btplants[is.na(mass_dry), mass_dry := mass_dry_genus_gbif] # then GBIF genus
btplants[is.na(mass_dry), mass_dry := mass_dry_genus_bt] # then BT genus

btplants[, mass_fresh := mass_fresh_species_gbif]
btplants[is.na(mass_fresh), mass_fresh := mass_fresh_species_bt]
btplants[is.na(mass_fresh), mass_fresh := mass_fresh_genus_gbif]
btplants[is.na(mass_fresh), mass_fresh := mass_fresh_genus_bt]

btplants[, mass := mass_species_gbif]
btplants[is.na(mass), mass := mass_species_bt]
btplants[is.na(mass), mass := mass_genus_gbif]
btplants[is.na(mass), mass := mass_genus_bt]

btplants[!is.na(mass_dry), length(unique(Species))] # 557 species now with data
btplants[!is.na(mass_fresh), length(unique(Species))] # 108 species now with data
btplants[!is.na(mass), length(unique(Species))] # 868 species now with data

btplants[mass_dry_genus_bt != mass_dry_species_bt & !is.na(mass_dry_genus_bt), plot(mass_dry_genus_bt, mass_dry_species_bt, log='xy')] # compare species and genus. pretty well correlated
btplants[mass_dry_genus_gbif != mass_dry_species_gbif & !is.na(mass_dry_genus_gbif), plot(mass_dry_genus_gbif, mass_dry_species_gbif, log='xy')] # compare species and genus. pretty well correlated

btplants[mass_fresh_genus_bt != mass_fresh_species_bt & !is.na(mass_fresh_genus_bt), plot(mass_fresh_genus_bt, mass_fresh_species_bt, log='xy')] # compare species and genus. maybe correlated
btplants[mass_fresh_genus_gbif != mass_fresh_species_gbif & !is.na(mass_fresh_genus_gbif), plot(mass_fresh_genus_gbif, mass_fresh_species_gbif, log='xy')] # compare species and genus. maybe correlated

btplants[mass_genus_bt != mass_species_bt & !is.na(mass_genus_bt), plot(mass_genus_bt, mass_species_bt, log='xy')] # compare species and genus. mabye correlated
btplants[mass_genus_gbif != mass_species_gbif & !is.na(mass_genus_gbif), plot(mass_genus_gbif, mass_species_gbif, log='xy')] # compare species and genus. mabye correlated

# compare among the different mass values
btplants[!is.na(mass_dry) & !is.na(mass_fresh), plot(mass_dry, mass_fresh, log='xy')] # not well correlated
btplants[!is.na(mass_dry) & !is.na(mass), plot(mass_dry, mass, log='xy')] # not well correlated
btplants[!is.na(mass_fresh) & !is.na(mass), plot(mass_fresh, mass, log='xy')] # 0


# trim to species with data
nrow(btplants) # 1,034,700
btplants.out <- btplants[!is.na(mass_dry) | !is.na(mass_fresh) | !is.na(mass), ]
nrow(btplants.out) # 2855

### write out. mass in g
write.csv(btplants.out, file = gzfile('output/mass_tryplants.csv.gz'))



#########################
# Find longevity for plants
#########################
# examine
trypl[, sort(unique(TraitName))] # "Dispersal distance", "Dispersal syndrome", "Plant biomass and allometry: Plant dry mass", "Plant biomass and allometry: Plant mass", "Plant fresh mass", "Plant lifespan (longevity)"         
trypl[, unique(TraitName), by = TraitID]
trypl[TraitID == 59, sort(unique(OrigValueStr))][1:100] # numbers as strings, but some have > or , or a range
trypl[TraitID == 59, sort(unique(OrigUnitStr))] # blank, or various strings implying year
trypl[TraitID == 59, sort(unique(ValueKindName))] # High, Low, Maximum, Mean, or Single
trypl[TraitID == 59, sort(unique(StdValue))] # numbers! These have been standardized by TRY. Use this where possible.
trypl[TraitID == 59, sort(unique(UnitName))] # blank or years
trypl[TraitID == 59, .(AccSpeciesName_clean, StdValue)][!duplicated(AccSpeciesName_clean), sum(!is.na(StdValue))] # 508 species.
trypl[TraitID == 59, .(AccSpeciesName_clean, OrigValueStr)][!duplicated(AccSpeciesName_clean), sum(!is.na(OrigValueStr))] # 24315 species. Closer to TRY report of 24712

trypl[TraitID == 59 & is.na(StdValue), sort(unique(OrigValueStr))] # 456 options. some are numbers. some are ranges with - or ,

# parse longevity strings
getlong <- function(x){
    out <- suppressWarnings(as.numeric(x))
    if(!is.na(out)){
        return(out)
    }

    out <- suppressWarnings(as.numeric(gsub('>', '', x)))
    if(!is.na(out)){
        return(out)
    }
    
    if(grepl(',', x)){
        x2 <- suppressWarnings(as.numeric(unlist(strsplit(x, split = ','))))
        out <- mean(x2)
        if(!is.na(out)){
            return(out)
        }
    }

    if(grepl('-', x)){
        x2 <- suppressWarnings(as.numeric(unlist(strsplit(x, split = '-'))))
        out <- mean(x2)
        if(!is.na(out)){
            return(out)
        }
    }
    
    if(x %in% c("annual-winter annual", "ephemeral", "summer annuals", "winter annual", "winter annuals")){
        return(0.5)
    }

    if(x %in% c("always annual", "ann", "annual", "Annual", "annual/ephemeral", "annual (occ. Perennial)", "annuals",
                "annual/woody", "annual-winter annual, biennial", "winter annual-biennial")){
        return(1)
    }

    if(x %in% c("annual/bieenial", "annual-biennial", "annual/biennial", "annual, Biennial", "Annual, Biennial", "annual/biennual",
                "annual/bisannual", "annual/shortperennial", "annual/ short perennial", "sometimes annual, always biennial")){
        return(1.5)
    }
    
    if(x %in% c("always biennial", "always biennial, always pluriennial-hapaxanthic", "always biennial, always pluriennial-hapaxanthic, always pluriennial-pollakanthic",
                "always biennial, always pluriennial-hapaxanthic, sometimes pluriennial-pollakanthic",  "always biennial, always pluriennial-pollakanthic",
                "Bennial", "biannual", "biasannual", "biennial", "Biennial")){
        return(2)
    }
    
    
    if(x %in% c("poly-annuals < 5 years (short-lived perennials)")){
        return(3)
    }
    
    if(x %in% c("poly-annuals 5-50 years (medium-lived perennials)")){
        return(27)
    }
    
    return(NA)
}

trypl[TraitID == 59 & is.na(StdValue), OrigValueStr_parse := vapply(OrigValueStr, FUN = getlong, FUN.VALUE = numeric(1))]
trypl[TraitID == 59 & !is.na(OrigValueStr_parse), 
      .(StdValue, OrigValueStr, OrigValueStr_parse)][!duplicated(cbind(StdValue, OrigValueStr, OrigValueStr_parse)),] # examine results

# combine StdValue and OrigValueStr_parse
trypl[TraitID == 59, longevity := StdValue]
trypl[TraitID == 59 & is.na(StdValue), longevity := OrigValueStr_parse]

# summarize to one trait per species
trylong <- trypl[TraitID == 59, .(mean = mean(longevity, na.rm = TRUE), n = sum(!is.na(longevity))), 
                 by = .(TraitID, TraitName, AccSpeciesName_clean)]
trylong[is.nan(mean), mean := NA] # convert nan to NA
trylong[!is.na(mean), ] # examine. 7413 species

trylong[, genus := tolower(vapply(strsplit(AccSpeciesName_clean, " "), `[`, 1, FUN.VALUE=character(1)))] # genus name in lower case
trylong[, speciesnm := tolower(vapply(strsplit(AccSpeciesName_clean, " "), `[`, 2, FUN.VALUE=character(1)))] # species

trylonggenu <- trylong[, .(mean = mean(mean), n = sum(n)), by = .(TraitID, TraitName, genus)] # genus-level averages

# merge with plant traits by species
btlong <- merge(btspp, trylong[TraitID == 59, .(genus, speciesnm, spp_try = AccSpeciesName_clean, longevity_species = mean)], by = c("genus", "speciesnm"), all.x = TRUE)

# look at the matches
btlong[!is.na(longevity_species), ]

btlong[!is.na(longevity_species), length(unique(Species))] # 660 species

btlong[!is.na(longevity_species), sort(unique(taxa_mod))] # some are not labeled as All, Benthos, or Invertebrates by BioTime
btlong[(!is.na(longevity_species)) & taxa_mod != 'Plant', 
         .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # Look correct.

# check for missing
btlong[taxa_mod == 'Plant' & is.na(longevity_species), length(unique(Species))] # 6840
btlong[taxa_mod == 'Plant' & is.na(longevity_species), 
         sort(unique(Species))] # mostly spcies names and genus names


# merge genus-values with species list
btlong <- merge(btlong, trylonggenu[TraitID == 59, .(genus, longevity_genus = mean)], by = 'genus', all.x = TRUE)

# examine genus values
btlong[!is.na(longevity_genus), .(Species, genus, longevity_species, longevity_genus)]

# choose from species- or genus-level mass values
btlong[, longevity := longevity_species]
btlong[is.na(longevity_species), longevity := longevity_genus] # use common where genus-level value is NA

btlong[!is.na(longevity), length(unique(Species))] # 831 species now with data

btlong[longevity_genus != longevity_species & !is.na(longevity_genus), plot(longevity_genus, longevity_species, log='xy')] # compare species and genus. pretty well correlated

# trim to species with data
nrow(btlong) # 1,035186
btlong.out <- btlong[!is.na(longevity), ]
nrow(btlong.out) # 6306

### write out
write.csv(btlong.out, file = gzfile('output/longevity_tryplants.csv.gz'))


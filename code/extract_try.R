################################################################
#######  extract body mass from the TRY Traits database    #####
######   plants                                            #####
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
btspp[, genus := tolower(vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1)))] # genus name in lower case
btspp[, speciesnm := tolower(vapply(strsplit(Species, "_"), `[`, 2, FUN.VALUE=character(1)))] # genus name in lower case

# load TRY traits plants
trypl <- fread('dataDL/try/8916.txt')


#########################
# Find mass for plants
#########################
# examine
trypl[, sort(unique(TraitName))] # "Dispersal distance", "Dispersal syndrome", "Plant biomass and allometry: Plant dry mass", "Plant biomass and allometry: Plant mass", "Plant fresh mass", "Plant lifespan (longevity)"         
trypl[, unique(TraitName), by = TraitID]

trypl[TraitName == "Plant biomass and allometry: Plant mass", sort(unique(OrigValueStr))][1:100] # numbers as strings
trypl[TraitID == 3451, sum(is.na(as.numeric(unique(OrigValueStr))))] # 0 convert to NA
trypl[TraitID == 3451, sort(unique(OrigUnitStr))] # all kg
trypl[TraitID == 3451, sort(unique(ValueKindName))] # all Single
trypl[TraitID == 3451, sort(unique(StdValue))] # numbers! These have been standardized by TRY. Use this
trypl[TraitID == 3451, sort(unique(UnitName))] # all g
trypl[TraitID == 3451, .(AccSpeciesName, StdValue)][!duplicated(AccSpeciesName), sum(!is.na(StdValue))] # 277 species
trypl[TraitID == 3451, .(AccSpeciesName, OrigValueStr)][!duplicated(AccSpeciesName), sum(!is.na(OrigValueStr))] # 277 species

trypl[TraitName == "Plant biomass and allometry: Plant dry mass", sort(unique(OrigValueStr))][1:100] # numbers as strings
trypl[TraitID == 700, sum(is.na(as.numeric(unique(OrigValueStr))))] # 2 convert to NA
trypl[TraitID == 700, .(val = unique(OrigValueStr))][is.na(as.numeric(val)), ] # 2 convert to NA: <NA> and .
trypl[TraitID == 700 & OrigValueStr %in% c('<NA>', '.'), ] # look at them
trypl[TraitID == 700, sort(unique(OrigUnitStr))] # g or mg
trypl[TraitID == 700, sort(unique(ValueKindName))] # Mean or Single
trypl[TraitID == 700, sort(unique(StdValue))] # numbers! These have been standardized by TRY. Use this
trypl[TraitID == 700, sort(unique(UnitName))] # all g
trypl[TraitID == 700, .(AccSpeciesName, StdValue)][!duplicated(AccSpeciesName), sum(!is.na(StdValue))] # 163 species. Not 238 as reported by TRY
trypl[TraitID == 700, .(AccSpeciesName, OrigValueStr)][!duplicated(AccSpeciesName), sum(!is.na(OrigValueStr))] # 163 species

trypl[TraitName == "Plant fresh mass", sort(unique(OrigValueStr))][1:100] # numbers as strings
trypl[TraitID == 3407, sum(is.na(as.numeric(unique(OrigValueStr))))] # 0 convert to NA
trypl[TraitID == 3407, sort(unique(OrigUnitStr))] # g
trypl[TraitID == 3407, sort(unique(ValueKindName))] # Single
trypl[TraitID == 3407, sort(unique(StdValue))] # numbers! These have been standardized by TRY. Use this
trypl[TraitID == 3407, sort(unique(UnitName))] # all g
trypl[TraitID == 3407, .(AccSpeciesName, StdValue)][!duplicated(AccSpeciesName), sum(!is.na(StdValue))] # 20 species.
trypl[TraitID == 3407, .(AccSpeciesName, OrigValueStr)][!duplicated(AccSpeciesName), sum(!is.na(OrigValueStr))] # 20 species

# summarize to one trait per species
trymass <- trypl[TraitID %in% c(700, 3407, 3451), .(mean = mean(StdValue, na.rm = TRUE), n = sum(!is.na(StdValue))), 
                 by = .(TraitID, TraitName, AccSpeciesName)]
trymass[is.nan(mean), mean := NA] # convert nan to NA
trymass[!is.na(mean) & TraitID == 700, ] # examine. 163 species for dry mass
trymass[!is.na(mean) & TraitID == 3407, ] # examine. 20 species for fresh mass
trymass[!is.na(mean) & TraitID == 3451, ] # examine. 277 species for mass

trymass[, genus := tolower(vapply(strsplit(AccSpeciesName, " "), `[`, 1, FUN.VALUE=character(1)))] # genus name in lower case
trymass[, speciesnm := tolower(vapply(strsplit(AccSpeciesName, " "), `[`, 2, FUN.VALUE=character(1)))] # species

trymassgenu <- trymass[, .(mean = mean(mean), n = sum(n)), by = .(TraitID, TraitName, genus)] # genus-level averages

# merge with plant traits by species
btplants <- merge(btspp, trymass[TraitID == 700, .(genus, speciesnm, spp_try_mass_dry = AccSpeciesName, mass_dry_species = mean)], by = c("genus", "speciesnm"), all.x = TRUE)
btplants <- merge(btplants, trymass[TraitID == 3407, .(genus, speciesnm, spp_try_mass_fresh = AccSpeciesName, mass_fresh_species = mean)], by = c("genus", "speciesnm"), all.x = TRUE)
btplants <- merge(btplants, trymass[TraitID == 3451, .(genus, speciesnm, spp_try_mass = AccSpeciesName, mass_species = mean)], by = c("genus", "speciesnm"), all.x = TRUE)

# look at the matches
btplants[!is.na(mass_dry_species), ]
btplants[!is.na(mass_fresh_species), ]
btplants[!is.na(mass_species), ]

btplants[!is.na(mass_dry_species), length(unique(Species))] # 71 species
btplants[!is.na(mass_fresh_species), length(unique(Species))] # 16
btplants[!is.na(mass_species), length(unique(Species))] # 67

btplants[!is.na(mass_dry_species) | !is.na(mass_fresh_species) | !is.na(mass_species), sort(unique(taxa_mod))] # some are not labeled as Benthos by BioTime
btplants[(!is.na(mass_dry_species) | !is.na(mass_fresh_species) | !is.na(mass_species)) & taxa_mod != 'Plant', 
         .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # one in Benthos. Looks correct.

# check if we can merge the TRY species name columns from the 3 data types
btplants[spp_try_mass_dry != spp_try_mass_fresh, .(spp_try_mass_dry, spp_try_mass_fresh)] # 0 don't match
btplants[spp_try_mass_dry != spp_try_mass, .(spp_try_mass_dry, spp_try_mass)] # 1 spp doesn't match, has variant included. shouldn't merge
btplants[spp_try_mass_fresh != spp_try_mass, .(spp_try_mass_fresh, spp_try_mass)] # 0 don't match. overall

# check for missing
btplants[taxa_mod == 'Plant' & (is.na(mass_dry_species) & is.na(mass_fresh_species) & is.na(mass_species)), length(unique(Species))] # 7292
btplants[taxa_mod == 'Plant' & (is.na(mass_dry_species) & is.na(mass_fresh_species) & is.na(mass_species)), 
         sort(unique(Species))] # mostly spcies names and genus names. Some common names.


# merge genus-values with species list
btplants <- merge(btplants, trymassgenu[TraitID == 700, .(genus, mass_dry_genus = mean)], by = 'genus', all.x = TRUE)
btplants <- merge(btplants, trymassgenu[TraitID == 3407, .(genus, mass_fresh_genus = mean)], by = 'genus', all.x = TRUE)
btplants <- merge(btplants, trymassgenu[TraitID == 3451, .(genus, mass_genus = mean)], by = 'genus', all.x = TRUE)

# examine genus values
btplants[!is.na(mass_dry_genus), .(Species, genus, mass_dry_species, mass_dry_genus)]
btplants[!is.na(mass_fresh_genus), .(Species, genus, mass_fresh_species, mass_fresh_genus)]
btplants[!is.na(mass_genus), .(Species, genus, mass_species, mass_genus)]

# choose from species- or genus-level mass values
btplants[, mass_dry := mass_dry_species]
btplants[is.na(mass_dry_species), mass_dry := mass_dry_genus] # use common where genus-level value is NA

btplants[, mass_fresh := mass_fresh_species]
btplants[is.na(mass_fresh_species), mass_fresh := mass_fresh_genus] # use common where genus-level value is NA

btplants[, mass := mass_species]
btplants[is.na(mass_species), mass := mass_genus] # use common where genus-level value is NA

btplants[!is.na(mass_dry), length(unique(Species))] # 534 species now with data
btplants[!is.na(mass_fresh), length(unique(Species))] # 104 species now with data
btplants[!is.na(mass), length(unique(Species))] # 813 species now with data

btplants[mass_dry_genus != mass_dry_species & !is.na(mass_dry_genus), plot(mass_dry_genus, mass_dry_species, log='xy')] # compare species and genus. pretty well correlated
btplants[mass_fresh_genus != mass_fresh_species & !is.na(mass_fresh_genus), plot(mass_fresh_genus, mass_fresh_species, log='xy')] # compare species and genus. maybe correlated
btplants[mass_genus != mass_species & !is.na(mass_genus), plot(mass_genus, mass_species, log='xy')] # compare species and genus. mabye correlated

# compare among the different mass values
btplants[!is.na(mass_dry) & !is.na(mass_fresh), plot(mass_dry, mass_fresh, log='xy')] # not well correlated
btplants[!is.na(mass_dry) & !is.na(mass), plot(mass_dry, mass, log='xy')] # not well correlated
btplants[!is.na(mass_fresh) & !is.na(mass), plot(mass_fresh, mass, log='xy')] # 0


# trim to species with data
nrow(btplants) # 1,034,700
btplants.out <- btplants[!is.na(mass_dry) | !is.na(mass_fresh) | !is.na(mass), ]
nrow(btplants.out) # 2746

### write out
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
trypl[TraitID == 59, .(AccSpeciesName, StdValue)][!duplicated(AccSpeciesName), sum(!is.na(StdValue))] # 508 species.
trypl[TraitID == 59, .(AccSpeciesName, OrigValueStr)][!duplicated(AccSpeciesName), sum(!is.na(OrigValueStr))] # 24315 species. Closer to TRY report of 24712

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
trypl[TraitID == 59 & !is.na(OrigValueStr_parse), .(StdValue, OrigValueStr, OrigValueStr_parse)][!duplicated(cbind(StdValue, OrigValueStr, OrigValueStr_parse)),]

# combine StdValue and OrigValueStr_parse

# summarize to one trait per species
trylong <- trypl[TraitID %in% c(59), .(mean = mean(StdValue, na.rm = TRUE), n = sum(!is.na(StdValue))), 
                 by = .(TraitID, TraitName, AccSpeciesName)]
trylong[is.nan(mean), mean := NA] # convert nan to NA
trylong[!is.na(mean), ] # examine. 582 species

trylong[, genus := tolower(vapply(strsplit(SpeciesName, " "), `[`, 1, FUN.VALUE=character(1)))] # genus name in lower case
trylong[, speciesnm := tolower(vapply(strsplit(SpeciesName, " "), `[`, 2, FUN.VALUE=character(1)))] # species

trylonggenu <- trylong[, .(mean = mean(mean), n = sum(n)), by = .(TraitID, TraitName, genus)] # genus-level averages

# merge with plant traits by species
btplants <- merge(btspp, trylong[TraitID == 700, .(genus, speciesnm, spp_try_mass_dry = SpeciesName, mass_dry_species = mean)], by = c("genus", "speciesnm"), all.x = TRUE)
btplants <- merge(btplants, trylong[TraitID == 3407, .(genus, speciesnm, spp_mass_fresh = SpeciesName, mass_fresh_species = mean)], by = c("genus", "speciesnm"), all.x = TRUE)
btplants <- merge(btplants, trylong[TraitID == 3451, .(genus, speciesnm, spp_try_mass = SpeciesName, mass_species = mean)], by = c("genus", "speciesnm"), all.x = TRUE)

# look at the matches
btplants[!is.na(mass_dry_species), ]
btplants[!is.na(mass_fresh_species), ]
btplants[!is.na(mass_species), ]

btplants[!is.na(mass_dry_species), length(unique(Species))] # 72 species
btplants[!is.na(mass_fresh_species), length(unique(Species))] # 16
btplants[!is.na(mass_species), length(unique(Species))] # 68

btplants[!is.na(mass_dry_species) | !is.na(mass_fresh_species) | !is.na(mass_species), sort(unique(taxa_mod))] # some are not labeled as Benthos by BioTime
btplants[(!is.na(mass_dry_species) | !is.na(mass_fresh_species) | !is.na(mass_species)) & taxa_mod != 'Plant', 
         .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # one in Benthos. Looks correct.

# check for missing
btplants[taxa_mod == 'Plant' & (is.na(mass_dry_species) & is.na(mass_fresh_species) & is.na(mass_species)), length(unique(Species))] # 7290
btplants[taxa_mod == 'Plant' & (is.na(mass_dry_species) & is.na(mass_fresh_species) & is.na(mass_species)), 
         sort(unique(Species))] # mostly spcies names and genus names. Some common names.


# merge genus-values with species list
btplants <- merge(btplants, trylonggenu[TraitID == 700, .(genus, mass_dry_genus = mean)], by = 'genus', all.x = TRUE)
btplants <- merge(btplants, trylonggenu[TraitID == 3407, .(genus, mass_fresh_genus = mean)], by = 'genus', all.x = TRUE)
btplants <- merge(btplants, trylonggenu[TraitID == 3451, .(genus, mass_genus = mean)], by = 'genus', all.x = TRUE)

# examine genus values
btplants[!is.na(mass_dry_genus), .(Species, genus, mass_dry_species, mass_dry_genus)]
btplants[!is.na(mass_fresh_genus), .(Species, genus, mass_fresh_species, mass_fresh_genus)]
btplants[!is.na(mass_genus), .(Species, genus, mass_species, mass_genus)]

# choose from species- or genus-level mass values
btplants[, mass_dry := mass_dry_species]
btplants[is.na(mass_dry_species), mass_dry := mass_dry_genus] # use common where genus-level value is NA

btplants[, mass_fresh := mass_fresh_species]
btplants[is.na(mass_fresh_species), mass_fresh := mass_fresh_genus] # use common where genus-level value is NA

btplants[, mass := mass_species]
btplants[is.na(mass_species), mass := mass_genus] # use common where genus-level value is NA

btplants[!is.na(mass_dry), length(unique(Species))] # 530 species now with data
btplants[!is.na(mass_fresh), length(unique(Species))] # 104 species now with data
btplants[!is.na(mass), length(unique(Species))] # 818 species now with data

btplants[mass_dry_genus != mass_dry_species & !is.na(mass_dry_genus), plot(mass_dry_genus, mass_dry_species, log='xy')] # compare species and genus. pretty well correlated
btplants[mass_fresh_genus != mass_fresh_species & !is.na(mass_fresh_genus), plot(mass_fresh_genus, mass_fresh_species, log='xy')] # compare species and genus. maybe correlated
btplants[mass_genus != mass_species & !is.na(mass_genus), plot(mass_genus, mass_species, log='xy')] # compare species and genus. mabye correlated

# compare among the different mass values
btplants[!is.na(mass_dry) & !is.na(mass_fresh), plot(mass_dry, mass_fresh, log='xy')] # not well correlated
btplants[!is.na(mass_dry) & !is.na(mass), plot(mass_dry, mass, log='xy')] # not well correlated
btplants[!is.na(mass_fresh) & !is.na(mass), plot(mass_fresh, mass, log='xy')] # 0


# trim to species with data
nrow(btplants) # 1,034,700
btplants.out <- btplants[!is.na(mass_dry) | !is.na(mass_fresh) | !is.na(mass), ]
nrow(btplants.out) # 2780

### write out
write.csv(btplants.out, file = gzfile('output/mass_tryplants.csv.gz'))


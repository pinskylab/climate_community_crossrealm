########################################################
#######  extract body mass from the LEDA database  #####
######   Plants                                    #####
########################################################

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

# load LEDA traits
# downloaded from https://uol.de/en/landeco/research/leda/data-files/ on 2020-03-13
# KLEYER, M., BEKKER, R.M., KNEVEL, I.C., BAKKER, J.P, THOMPSON, K., SONNENSCHEIN, M., POSCHLOD, P., VAN GROENENDAEL, 
# J.M., KLIMES, L., KLIMESOVÁ, J., KLOTZ, S., RUSCH, G.M., HERMY, M., ADRIAENS, D., BOEDELTJE, G., BOSSUYT, B., DANNEMANN, A., 
# ENDELS, P., GÖTZENBERGER, L., HODGSON, J.G., JACKEL, A-K., KÜHN, I., KUNZMANN, D., OZINGA, W.A., RÖMERMANN, C., STADLER, M., 
# SCHLEGELMILCH, J., STEENDAM, H.J., TACKENBERG, O., WILMANN, B., CORNELISSEN, J.H.C., ERIKSSON, O., GARNIER, E., PECO, B. (2008):
# The LEDA Traitbase: A database of life-history traits of Northwest European flora. Journal of Ecology 96: 1266-1274. 
lifespan <- fread('dataDL/ledatraits/plant_life_span.txt')
dispersal <- fread('dataDL/ledatraits/dispersal_type.txt')

# check for duplicate rows
lifespan[, .(uniq = length(unique(`SBS name`)), total = .N)] # many duplicate names
dispersal[, .(uniq = length(unique(`SBS name`)), total = .N)] # same

lifespan[`SBS name` %in% lifespan[duplicated(`SBS name`), `SBS name`], .(`SBS name`, `plant lifespan`)]
numlifespans <- lifespan[`SBS name` %in% lifespan[duplicated(`SBS name`), `SBS name`], .(nvals = length(unique(`plant lifespan`))), by = `SBS name`][nvals > 1,] # some are >1
numdisps <- dispersal[`SBS name` %in% dispersal[duplicated(`SBS name`), `SBS name`], .(nvals = length(unique(`gen. dispersal type`))), by = `SBS name`][nvals > 1,] # some are >1

# remove any speces with >1 value
nrow(lifespan) # 4187
lifespan <- lifespan[!(`SBS name` %in% numlifespans$`SBS name`), .(`SBS name`, `plant lifespan`)]
nrow(lifespan) # 1966

nrow(dispersal) # 16480
dispersal <- dispersal[!(`SBS name` %in% numdisps$`SBS name`), .(`SBS name`, `gen. dispersal type`)]
nrow(dispersal) # 6490

# collapse to one row per species
lifespan <- lifespan[!duplicated(`SBS name`), ]
nrow(lifespan) # 1561

dispersal <- dispersal[!duplicated(`SBS name`), ]
nrow(dispersal) # 6229


#########################
# Merge plant trait data
#########################
# merge traits
btplants <- merge(btspp, lifespan[, .(Species = gsub(' ', '_', `SBS name`), spp_leda_lifespan = `SBS name`, lifespan_cat_spp = `plant lifespan`)], by = 'Species', all.x = TRUE)
btplants <- merge(btplants, dispersal[, .(Species = gsub(' ', '_', `SBS name`), spp_leda_dispersal = `SBS name`, 
                                          dispersal_type_spp = `gen. dispersal type`)], by = 'Species', all.x = TRUE)

# look at the matches
btplants[!is.na(lifespan_cat_spp), ]
btplants[!is.na(lifespan_cat_spp), length(unique(Species))] # 204
btplants[!is.na(lifespan_cat_spp), sort(unique(taxa_mod))] # some are labeled as Benthos by BioTime
btplants[!is.na(lifespan_cat_spp) & taxa_mod != 'Plant', .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # in Benthos. But all seem correct.

btplants[!is.na(dispersal_type_spp), ]
btplants[!is.na(dispersal_type_spp), length(unique(Species))] # 289
btplants[!is.na(dispersal_type_spp), sort(unique(taxa_mod))] # some are labeled as All, Benthos, or Invertebrates by BioTime
btplants[!is.na(dispersal_type_spp) & taxa_mod != 'Plant', .(Species, taxa_mod)][!duplicated(cbind(Species, taxa_mod)),] # All seem correct.


# check for missing
btplants[taxa_mod == 'Plant' & is.na(lifespan_cat_spp), length(unique(Species))] # 7226
btplants[taxa_mod == 'Plant' & is.na(lifespan_cat_spp), sort(unique(Species))] # mostly species names and genus names.

# don't match on common name, since don't see any in the list

# get genus-level values (most common category)
btplants[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))] # not necessarily the genus, if "species name" is not a binomial
lifespan[, genus := vapply(strsplit(`SBS name`, " "), `[`, 1, FUN.VALUE=character(1))]
dispersal[, genus := vapply(strsplit(`SBS name`, " "), `[`, 1, FUN.VALUE=character(1))]
btplantsgenu <- btplants[!duplicated(genus) & taxa_mod == 'Plant', .(genus)] # list of potential mammal genera (ignore potential mammals in other categories)
btplantsgenu[, lifespan_cat_genus := NA_character_] 
btplantsgenu[, dispersal_type_genus := NA_character_] 

for(i in 1:nrow(btplantsgenu)){ # a brute-force for loop approach, but fast
    k1 <- which(btplantsgenu$genus[i] == lifespan$genus)
    if(length(k1) > 0){
        btplantsgenu$lifespan_cat_genus[i] <- names(sort(table(lifespan$`plant lifespan`[k1]),decreasing=TRUE)[1]) # will pick one if there is a tie
    }
    k2 <- which(btplantsgenu$genus[i] == dispersal$genus)
    if(length(k2) > 0){
        btplantsgenu$dispersal_type_genus[i] <- names(sort(table(dispersal$`gen. dispersal type`[k2]),decreasing=TRUE)[1])
    }
}
btplantsgenu[!is.na(lifespan_cat_genus),] # found 298 values
btplantsgenu[!is.na(lifespan_cat_genus), sort(unique(genus))] # check by eye that these look like genera

btplantsgenu[!is.na(dispersal_type_genus),] # found 408 values
btplantsgenu[!is.na(dispersal_type_genus), sort(unique(genus))] # check by eye that these look like genera

# merge genus-values with species list
btplants <- merge(btplants, btplantsgenu, by = 'genus', all.x = TRUE)

# choose from species- or genus-level values
btplants[, lifespan_cat := lifespan_cat_spp]
btplants[is.na(lifespan_cat_spp), lifespan_cat := lifespan_cat_genus] # use genus where species-level values are NA
btplants[is.nan(lifespan_cat), lifespan_cat := NA] # NaN to NA

btplants[, dispersal_type := dispersal_type_spp]
btplants[is.na(dispersal_type_spp), dispersal_type := dispersal_type_genus] # use genus where species-level values are NA
btplants[is.nan(dispersal_type), dispersal_type := NA] # NaN to NA

btplants[!is.na(lifespan_cat), length(unique(Species))] # 1305 species now with data
btplants[!is.na(dispersal_type), length(unique(Species))] # 1644 species now with data

btplants[!is.na(lifespan_cat_genus), table(lifespan_cat_spp, lifespan_cat_genus)] # compare species and genus. pretty well correlated (most values on the diagonal)
btplants[!is.na(dispersal_type_genus), table(dispersal_type_spp, dispersal_type_genus)] # pretty well correlated

# merge LEDA scientific names (from lifespan vs. dispersal merging)
btplants[, spp_leda := spp_leda_lifespan]
btplants[is.na(spp_leda_lifespan), spp_leda := spp_leda_dispersal]
btplants[, spp_leda_lifespan := NULL]
btplants[, spp_leda_dispersal := NULL]


# trim to species with data
nrow(btplants) # 1,034,700
btplants.out <- btplants[!is.na(lifespan_cat) | !is.na(dispersal_type), ]
nrow(btplants.out) # 5553


### write out
write.csv(btplants.out, file = gzfile('output/traits_ledaplants.csv.gz'))

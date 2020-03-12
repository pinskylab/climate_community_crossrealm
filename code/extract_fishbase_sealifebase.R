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
load('data/biotime_blowes/bt.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt

# BioTime species
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp

# add taxa_mod to spp list and trim to ones we will search
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], all.x = TRUE, by = 'rarefyID')
nrow(btspp)

btsppfish <- btspp[REALM %in% c('Marine', 'Freshwater') & taxa_mod == 'Fish',] # from fishbase
nrow(btsppfish) # 454431 rows
btsppfishu <- btsppfish[!duplicated(Species), .(Species, REALM, taxa_mod)]
nrow(btsppfishu) # 4700 species

btsppsl <- btspp[REALM == 'Marine' & taxa_mod %in% c('All', 'Benthos', 'Invertebrates', 'Marine invertebrates/plants'),] # from sealifebase
nrow(btsppsl) # 328412 rows
btsppslu <- btsppsl[!duplicated(Species), .(Species, REALM, taxa_mod)]
nrow(btsppslu) # 9511 species


#########################
# Find BM from Fishbase
#########################
# validate species names
btsppfishu[, spp_fb := paste0(validate_names(gsub('_', ' ', Species)), collapse = ";"), by = seq_len(NROW(btsppfishu))] # validate names. Separate multiple names with ;
btsppfishu[, sum(spp_fb == '')] # 996. number of blanks
btsppfishu[, sum(spp_fb != '')] # 3704. number of validated names
btsppfishu[gsub('_', ' ', Species) != spp_fb & spp_fb != '', .(Species, spp_fb)] # check where names don't match exactly. 372 non-blanks
btsppfishu[grep(';', spp_fb), .(Species, spp_fb)] # multiple names returned. 5.

# manually fix the multiple names (only 5)
btsppfishu[, spp_fb2 := spp_fb] # new column with manually fixed FB names (for 5)
btsppfishu[spp_fb2 == '', spp_fb2 := NA]
    bt[rarefyID %in% btspp[Species == 'Epigonus_macrops', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[spp_fb == 'Epigonus robustus;Epigonus macrops', spp_fb2 := 'Epigonus macrops'] # since in NW Atlantic
btsppfishu[spp_fb == 'Hemitripterus americanus;Hemitripterus americanus', spp_fb2 := 'Hemitripterus americanus'] # identical
    bt[rarefyID %in% btspp[Species == 'Centriscus_cristatus', rarefyID], ][!duplicated(rarefyID),] # check location
btsppfishu[spp_fb == 'Centriscus cristatus;Notopogon lilliei', spp_fb2 := 'Centriscus cristatus'] # since in NW Australian shelf
btsppfishu[spp_fb == 'Parajulis poecilepterus;Parajulis poecilepterus', spp_fb2 := 'Parajulis poecilepterus'] # identical
btsppfishu[spp_fb == 'Chaetodon melannotus;Chaetodon melannotus', spp_fb2 := 'Chaetodon melannotus'] # identical

# get species-level masses
btsppfishu[!is.na(spp_fb2), mass_species := as.numeric(species(spp_fb2, fields=c("Weight"))), by = seq_len(NROW(btsppfishu[!is.na(spp_fb2)]))] # in grams

# get genus-level values
btsppfishu[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))]
btsppfishgenu <- btsppfishu[!duplicated(genus), .(genus)] # list of genera

getgenmass <- function(x){ # function to return average weight of species in a genus
    spps <- species_list(Genus = x)
    mean(species(spps, fields = 'Weight')$Weight, na.rm=TRUE)
}
btsppfishgenu[, mass_genus := getgenmass(genus), by = seq_len(NROW(btsppfishgenu))]
btsppfishgenu[!is.na(mass_genus),] # found 490 values

# merge genus-values with species list
btsppfishu <- merge(btsppfishu, btsppfishgenu, by = 'genus', all.x = TRUE)

# choose from species- or genus-level values
btsppfishu[, mass := mass_species]
btsppfishu[is.na(mass_species), mass := mass_genus] # use genus where species-level value is NA
btsppfishu[is.nan(mass), mass := NA] # NaN to NA

btsppfishu[is.na(mass), length(unique(Species))] # 2720 species now with data

btsppfishu[mass_genus != mass_species, plot(mass_genus, mass_species)] # compare species and genus

# merge back with full list of species and studies
# check this!
btsppfish <- merge(btsppfish, btsppfishu[, .(Species, spp_fb2, mass_genus, mass_species, mass)], by = 'Species', all.x = TRUE)

# trim to species with data
nrow(btsppfish)
btsppfish.out <- btsppfish[!is.na(mass), ]
nrow(btsppfish.out)

### write out
write.csv(btsppfish.out, file = gzfile('output/mass_fishbase.csv.gz'))
          
#########################
# Find BM from Sealifebase
#########################
# validate species names
btsppslu[, spp_fb := paste0(validate_names(gsub('_', ' ', Species), server = 'sealifebase'), collapse = ";"), by = seq_len(NROW(btsppslu))] # validate names. Separate multiple names with ;
btsppslu[, sum(spp_fb == '')] # 4614 number of blanks
btsppslu[, sum(spp_fb != '')] # 4897 number of validated names
btsppslu[gsub('_', ' ', Species) != spp_fb & spp_fb != '' & spp_fb != 'NA', .(Species, spp_fb)] # check where names don't match exactly. 372 non-blanks
btsppslu[grep(';', spp_fb), .(Species, spp_fb)] # multiple names returned. 18.

# manually fix the multiple names (only 5)
btsppslu[, spp_fb2 := spp_fb] # new column with manually fixed FB names (for 5)
btsppslu[spp_fb2 == '', spp_fb2 := NA]
bt[rarefyID %in% btspp[Species == 'Esperiopsis_fucorum', rarefyID], ][!duplicated(rarefyID),] # check location
btsppslu[spp_fb == 'Esperiopsis fucorum;Amphilectus fucorum', spp_fb2 := 'Esperiopsis fucorum'] # unclear, so picked one

bt[rarefyID %in% btspp[Species == 'Leuconia_nivea', rarefyID], ][!duplicated(rarefyID),] # check location
btsppslu[spp_fb == 'Leuconia johnstoni;Leuconia nivea', spp_fb2 := 'Leuconia nivea'] # unclear, so picked on

bt[rarefyID %in% btspp[Species == 'Cancer_irroratus', rarefyID], ][!duplicated(rarefyID),] # check location
btsppslu[spp_fb == 'Cancer irroratus;Cancer plebejus', spp_fb2 := 'Cancer irroratus'] # picked because in NW Atlantic

bt[rarefyID %in% btspp[Species == 'Typosyllis_aciculata', rarefyID], ][!duplicated(rarefyID),] # check location
btsppslu[spp_fb == 'Syllis aciculata;Typosyllis aciculata', spp_fb2 := 'Syllis aciculata'] # picked because in NW Atlantic

bt[rarefyID %in% btspp[Species == 'Nematonereis_unicornis', rarefyID], ][!duplicated(rarefyID),] # check location
btsppslu[spp_fb == 'Nematonereis hebes;Lysidice unicornis', spp_fb2 := 'Nematonereis hebes'] # picked because in NW Atlantic

bt[rarefyID %in% btspp[Species == 'Haploscoloplos_fragilis', rarefyID], ][!duplicated(rarefyID),] # check location
btsppslu[spp_fb == 'Haploscoloplos fragilis;Leitoscoloplos fragilis', spp_fb2 := 'Leitoscoloplos fragilis'] # picked because in NW Atlantic

bt[rarefyID %in% btspp[Species == 'Tellina_listeri', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Tellina listeri;Tellinella listeri', spp_fb2 := 'Tellina listeri'] # not clear, so picked on

btsppslu[spp_fb == 'Chiton squamosus;Chiton squamosus', spp_fb2 := 'Chiton squamosus'] # simple duplicate

bt[rarefyID %in% btspp[Species == 'Pitar_cordatus', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Pitarenus cordatus;Pitar cordatus', spp_fb2 := 'Pitar cordatus'] # both in NW atlantic, so picked one

bt[rarefyID %in% btspp[Species == 'Murex_fulvescens', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Hexaplex fulvescens;Murex fulvescens', spp_fb2 := 'Murex fulvescens'] # both in NW atlantic, so picked one

bt[rarefyID %in% btspp[Species == 'Beroe_cucumis', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Hormiphora cucumis;Beroe cucumis', spp_fb2 := 'Beroe cucumis'] # picked because in Atlantic

bt[rarefyID %in% btspp[Species == 'Siphonodentalium_verrilli', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Pulsellum verrilli;Siphonodentalium verrilli', spp_fb2 := 'Pulsellum verrilli'] # picked because in Atlantic. Other species not in SLB

btsppslu[spp_fb == 'Inodrillia carpenteri;Inodrillia carpenteri', spp_fb2 := 'Inodrillia carpenteri'] # simple duplicate

bt[rarefyID %in% btspp[Species == 'Dodecaceria_fewkesi', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Dodecaceria fewkesi;Dodecaceria fistulicola', spp_fb2 := 'Dodecaceria fewkesi'] # in NE Pacific

bt[rarefyID %in% btspp[Species == 'Cercodemas_anceps', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Pentacta anceps;Cercodemas anceps', spp_fb2 := 'Cercodemas anceps'] # in SW Australia

bt[rarefyID %in% btspp[Species == 'Ergalatax_contracta', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Cronia contracta;Ergalatax contracta', spp_fb2 := 'Cronia contracta'] # more obviously in Australia

bt[rarefyID %in% btspp[Species == 'Chirona_tenuis', rarefyID], ][!duplicated(rarefyID),]
btsppslu[spp_fb == 'Striatobalanus tenuis;Chirona tenuis', spp_fb2 := 'Chirona tenuis'] # both in W Pac, so picked original name

btsppslu[spp_fb == 'Haliplanella lineata;Haliplanella lineata;Haliplanella lineata;Haliplanella lineata', spp_fb2 := 'Haliplanella lineata'] # duplicate


# get species-level masses
btsppslu[!is.na(spp_fb2), mass_species := as.numeric(species(spp_fb2, fields=c("Weight"), server = 'sealifebase')), by = seq_len(NROW(btsppslu[!is.na(spp_fb2)]))] # in grams
btsppslu[, sum(!is.na(mass_species))] # 62

# get genus-level values
btsppslu[, genus := vapply(strsplit(Species, "_"), `[`, 1, FUN.VALUE=character(1))]
btsppslgenu <- btsppslu[!duplicated(genus), .(genus)] # list of genera

getgenmassSLB <- function(x){ # function to return average weight of species in a genus
    spps <- species_list(Genus = x, server = 'sealifebase')
    mean(species(spps, fields = 'Weight', server = 'sealifebase')$Weight, na.rm=TRUE)
}
btsppslgenu[, mass_genus := getgenmassSLB(genus), by = seq_len(NROW(btsppslgenu))]
btsppslgenu[!is.na(mass_genus),] # 58 values

# merge genus-values with species list
btsppslu <- merge(btsppslu, btsppslgenu, by = 'genus', all.x = TRUE)

# choose from species- or genus-level values
btsppslu[, mass := mass_species]
btsppslu[is.na(mass_species), mass := mass_genus] # use genus where species-level value is NA
btsppslu[is.nan(mass), mass := NA] # NaN to NA

btsppslu[is.na(mass), length(unique(Species))] # 2720 species now with data

btsppslu[mass_genus != mass_species, plot(mass_genus, mass_species)] # compare species and genus


# merge back with full list of species and studies
btsppsl <- merge(btsppsl, btsppslu[, .(Species, spp_fb2, mass_genus, mass_species, mass)], by = 'Species', all.x = TRUE)

# trim to species with data
nrow(btsppsl)
btsppsl.out <- btsppsl[!is.na(mass), ]
nrow(btsppsl.out)

### write out
write.csv(btsppsl.out, file = gzfile('output/mass_sealifebase.csv.gz'))

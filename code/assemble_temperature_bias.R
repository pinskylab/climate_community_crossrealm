# assemble the community temperature bias (following Stuart-Smith et al. 2015 Nature)

##############
## functions
##############
require(data.table)
source('code/util.R') # for picklongest()

##############
# Load data
##############

# environmental temperature data
temperature <- fread('output/temperature_byrarefyID.csv.gz')

# species temperature index data
sti <- fread('data/taxa_burrows/stibysprealm.csv', drop = 1)
sti <- sti[!is.na(X.50.), ] # trim out species without data

# BioTime change (for taxa_mod) and species
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID') # add taxa_mod to spp list, trim to rarefyIDs we use
rm(bt)

btspp[, length(unique(Species))] # 26289 species

# load gbif species names
phylo = fread('output/taxonomy.csv.gz', drop = 1) # from match_names_with_gbif.R

# merge spp with gbif names
btspp <- merge(btspp, phylo[, .(class_gbif = class, genus_gbif = genus, 
                                Species_gbif = species, Species = gsub(' ', '_', original_name))], all.x = TRUE, by = 'Species')

#################
# Assemble data
#################
# merge on biotime name
btsti <- merge(btspp, sti[, .(REALM, Species = gsub(' ', '_', speciesname), sti_bt = X.50., stimax_bt = X.90.)], by = c('Species', 'REALM'), all.x = TRUE)

# merge on gbif name
btsti <- merge(btsti, sti[, .(REALM, Species_gbif = speciesname, sti_gbif = X.50., stimax_gbif = X.90.)], by = c('Species_gbif', 'REALM'), all.x = TRUE)
btsti[!is.na(sti_gbif) & is.na(sti_bt), .(.N, length(unique(Species)))] # adds 3101 rows of data across 304 species


# look at the matches
btsti[!is.na(sti_bt) | !is.na(sti_gbif), ]
btsti[!is.na(sti_bt) | !is.na(sti_gbif), length(unique(Species))] # 9686 (of 26289)

# compare matches on biotime and gbif names
btsti[!is.na(sti_bt) & !is.na(sti_gbif) & !duplicated(Species), plot(sti_bt, sti_gbif)]; abline(0,1); abline(-10,1, lty=2); abline(10,1, lty=2) # mostly on 1:1 line, but some are far off
btsti[abs(sti_bt - sti_gbif) > 10,][!duplicated(cbind(Species, sti_bt, sti_gbif)), .(Species, Species_gbif, REALM, taxa_mod, sti_bt, sti_gbif)] # 13 cases were >10degC off

# combine bt and gbif stis
btsti[, sti := sti_bt]
btsti[is.na(sti_bt) & !is.na(sti_gbif), sti := sti_gbif] # fill in secondarily with gbif matches
btsti[, stimax := stimax_bt]
btsti[is.na(stimax_bt) & !is.na(stimax_gbif), stimax := stimax_gbif] # fill in secondarily with gbif matches

btsti[!is.na(sti), length(unique(Species))] # 9686 good that it matches from before.
btsti[!is.na(stimax), length(unique(Species))] # 9686 good that it matches from before.

# calculate average STI (=CTI) by rarefyID
#see how we're doing (per rarefyID: # species with data, mean sti, sd sti, geometric mean sti, geometric standard deviation sti)
btcti <- btsti[, .(cti = mean(sti, na.rm = TRUE), sti_sd = sd(sti, na.rm = TRUE),
                   ctimax = mean(stimax, na.rm = TRUE), stimax_sd = sd(stimax, na.rm = TRUE),
                   nspp = length(unique(Species)), nspp_wdata = sum(!is.na(sti))), 
                   by = .(STUDY_ID, rarefyID, REALM, taxa_mod)]
btcti[is.nan(cti), cti := NA_real_]
btcti[is.nan(ctimax), ctimax := NA_real_]
nrow(btcti) # 53467
setkey(btcti, STUDY_ID, rarefyID)
btcti

btcti[, hist(nspp_wdata/nspp)] # fraction of species in each study with data. most >50%

# add environmental temperature to CTI, calculate bias
tempave <- temperature[, .(tempave = picklongest(year1, year2, tempave)), by = rarefyID] # pick the tempave for the full time series duration
btcti <- merge(btcti, tempave[, .(rarefyID, tempave)], all.x = TRUE, by = 'rarefyID')
btcti[, thermal_bias := cti - tempave]
btcti[, thermal_biasmax := ctimax - tempave]


# plot cti vs. temperature
# increases, but not as clear as in Stuart-Smith et al. 2015 Nature
btcti[, plot(tempave, cti)]; abline(0,1, col = 'red')
btcti[, plot(tempave, ctimax)]; abline(0,1, col = 'red')

# make a list of priority species for further sti finding
# studies with less than 50% species with data, at least 5 species in the study
priority_studies <- btsti[, .(n = length(unique(Species)), val = sum(!is.na(sti))), 
                    by = .(rarefyID, STUDY_ID, taxa_mod)][, .(n = sum(n), val = sum(val)), 
                                                          by = .(STUDY_ID, taxa_mod)][(val/n) < 0.5 & n > 4, STUDY_ID]
length(priority_studies) # 126 studies

priority_spp <- btsti[STUDY_ID %in% priority_studies & is.na(sti), 
                         .(Species, REALM, taxa_mod)][!duplicated(cbind(Species, REALM, taxa_mod)), .(REALM, taxa_mod, Species)]
setkey(priority_spp, REALM, taxa_mod, Species)
nrow(priority_spp) # 10447
priority_spp



#########
# Output
#########

# output CTI and bias by rarefyID
btcti[, ':='(tempave = NULL, STUDY_ID = NULL, REALM = NULL, taxa_mod = NULL)] # remove columns
write.csv(btcti, gzfile('output/cti_byrarefyID.csv.gz'), row.names = FALSE)

# output list of priority species
write.csv(priority_spp, gzfile('output/sti_prioritymissing.csv.gz'), row.names = FALSE)

# calculate lifespan index for each species, average to community
# uses natural mortality calculation from McCoy & Gillooly 2008 Ecology Letters
# needs output from assemble_temperature.Rmd, assemble_mass.R, and match_names_with_gbif.R

require(data.table)

geomean = function(x){ # geometric mean. remove NAs and x<0
    exp(sum(log(x[x > 0 & !is.na(x)])) / length(x[x > 0 & !is.na(x)]))
}

geosd <- function(x){ # geometric standard deviation
    exp(sd(log(x[x > 0 & !is.na(x)])))
}

##############
# Load data
##############
mass <- fread('output/mass_byspecies.csv.gz') # assembled data on mass, in g
temperature <- fread('output/temperature_byrarefyID.csv.gz') # data on environmental temperature
endo <- fread('output/metab_feed_byspecies.csv.gz') # endotherm vs. ectotherm classifications


# BioTime change (for taxa_mod) and species
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin); rm(bt_malin) # rename to bt
load('data/biotime_blowes/bt_grid_spp_list.Rdata') # loads bt_grid_spp_list. this has some studies not in bt.
btspp <- data.table(bt_grid_spp_list); rm(bt_grid_spp_list) # rename to btspp
btspp <- merge(btspp, bt[!duplicated(rarefyID), .(rarefyID, taxa_mod)], by = 'rarefyID') # add taxa_mod to spp list and trims to spp in bt

#################
# assemble data
#################

# merge bt with mass, endoecto, temperature
lsp <- merge(btspp[, .(rarefyID, Species, REALM, STUDY_ID, taxa_mod)], 
             mass[, .(rarefyID, Species, mass)], by = c('Species', 'rarefyID'), all.x = TRUE)
lsp <- merge(lsp, endo, by = c('Species', 'rarefyID'), all.x = TRUE)
lsp <- merge(lsp, temperature[, .(rarefyID, tempave)], by = 'rarefyID', all.x = TRUE)
             
# decide which temperature to use
# use bird and mammal body temps following McCoy & Gillooly
# deg Kelvin
lsp[, temp_spp := tempave + 273.15]
lsp[metab == 'bird', temp_spp := 40 + 273.15] 
lsp[metab == 'mamm', temp_spp := 38 + 273.15]

# decide which activation energy to use
# follow McCoy & Gillooly
lsp[, E := NA_real_]
lsp[feeding == 'consumer', E := 0.65]
lsp[feeding == 'producer', E := 0.32]

###############################################################
# calculate lifespan from allometry and temperature
# Equation in Fig. 3 of McCoy & Gillooly 2008
###############################################################

k = 8.62 * 10^-5 # Boltzmann's constant
T20 = 293 # standardization temperature
lsp[, Ztrans := -0.22 * log(1/4 * mass) - 1.3] # y-axis in Fig. 3 of M&G2008 (in corrigendum). Use dry mass = wet mass/4 per Supplement Table 1.
lsp[, Z := exp(Ztrans)/exp(E/k*(1/temp_spp - 1/T20))] # natural mortality rate (yr-1) from Fig. 3 transformation
lsp[, lifespan := 1/Z] # lifespan in years



# check coverage
lsp[, .(n = length(unique(Species)), val = sum(!is.na(lifespan))), by = rarefyID][, hist(val/n)] # most have >50% of species represented!
lsp[, .(n = length(unique(Species)), val = sum(!is.na(lifespan))), by = rarefyID][(val/n) < 0.5, ] # 4139 rarefyID with <50%

# check data
lsp[, hist(lifespan)]

# average lifespan by rarefyID
#also see how we're doing (per rarefyID): # species with data, mean mass, sd mass, geometric mean mass, geometric standard deviation mass)
lsp.sum <- lsp[, .(lifespan_mean = mean(lifespan, na.rm = TRUE), lifespan_sd = sd(lifespan, na.rm = TRUE), lifespan_geomean = geomean(lifespan),
                   lifespan_geosd = geosd(lifespan), nspp = length(unique(Species)), nspp_wdata = sum(!is.na(lifespan))), 
               by = .(STUDY_ID, rarefyID, REALM, taxa_mod)]
lsp.sum[is.nan(lifespan_mean), lifespan_mean := NA_real_]
lsp.sum[is.nan(lifespan_geomean), lifespan_geomean := NA_real_]
nrow(lsp.sum) # 53467
setkey(lsp.sum, STUDY_ID, rarefyID)
lsp.sum


############
# output
############

write.csv(lsp.sum, gzfile('output/lifespan_byrarefyID.csv.gz'), row.names = FALSE)



# Assemble temperature dataset

require(data.table)

### Read in biotime data
load('data/biotime_blowes/bt.Rdata')
bt <- data.table(bt_malin); rm(bt_malin)

# trim to locations
btterra <- bt[!duplicated(rarefyID), ][REALM %in% c('Terrestrial', 'Freshwater'), .(rarefyID, lon = rarefyID_x, lat = rarefyID_y)]
btaqua <- bt[!duplicated(rarefyID), ][REALM == 'Marine', .(rarefyID, lon = rarefyID_x, lat = rarefyID_y)]

### Add in MODIS data

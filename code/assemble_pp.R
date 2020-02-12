# Assemble temperature dataset

require(data.table)

### Read in biotime data
load('data/biotime_blowes/bt.Rdata')
bt <- data.table(bt_malin); rm(bt_malin)

# trim to locations
btterra <- bt[!duplicated(STUDY_ID), ][REALM %in% c('Terrestrial', 'Freshwater'), .(STUDY_ID, lon = rarefyID_x, lat = rarefyID_y)]
btaqua <- bt[!duplicated(STUDY_ID), ][REALM == 'Marine', .(STUDY_ID, lon = rarefyID_x, lat = rarefyID_y)]

### Add in MODIS data

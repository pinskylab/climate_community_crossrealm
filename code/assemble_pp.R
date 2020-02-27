# Assemble temperature dataset

require(data.table)
require(MODISTools)

### Read in biotime data
load('data/biotime_blowes/bt.Rdata')
bt <- data.table(bt_malin); rm(bt_malin)

# trim to locations
#btterra <- bt[!duplicated(rarefyID), ][REALM %in% c('Terrestrial', 'Freshwater'), .(site_name = rarefyID, lon = rarefyID_x, lat = rarefyID_y)]
btaqua <- bt[!duplicated(rarefyID), ][REALM == 'Marine', .(site_name = rarefyID, lon = rarefyID_x, lat = rarefyID_y)]

### Add in MODIS data on land
# tempterra <- mt_batch_subset(df = btterra,
#                               product = "MOD17A3H",
#                               band = "Npp_500m",
#                               internal = TRUE,
#                               start = "2004-01-01",
#                               end = "2018-12-31",
#                               out_dir = "temp")
# missing some processing steps from Helmut here
# write.csv(modis.terra, file = 'data/modis_npp/NPPterr.csv')

### Add in MODIS data in the ocean
for(i in 642:nrow(btaqua)){
    print(paste0(i, ' of ', nrow(btaqua)))
    temp <- tryCatch({
        mt_subset(lat = btaqua[i, lat],
                  lon = btaqua[i, lon],
                  site_name = btaqua[i, site_name],
                  product = "MOD17A3H",
                  band = "Npp_500m",
                  internal = TRUE,
                  start = "2004-01-01",
                  end = "2018-12-31",
                  out_dir = "temp",
                  progress = FALSE)
    }, warning = function(w) {
        print(paste0('i=', i, ': ', w))
        return(data.frame(value = numeric(0)))
    }, error = function(e) {
        print(paste0('i=', i, ': ', e))
        return(data.frame(value = numeric(0)))
    }, finally = {
        
    })
    
    temp$value[temp$value > 32700] <- NA # valid range is up to 32700
    tempsum <- btaqua[i, .(latitude = lat, longitude = lon, site = site_name)]
    tempsum[, ':='(N = nrow(temp), N.na = sum(is.na(temp$value)), npp.mean = mean(temp$value, na.rm = TRUE), npp.sd = sd(temp$value, na.rm = TRUE))]
    tempsum[, npp.cv := npp.sd/npp.mean]
    
    if(i == 1) modis.aqua <- tempsum
    if(i > 1) modis.aqua <- rbind(modis.aqua, tempsum)
}

write.csv(modis.aqua, file = 'data/modis_npp/NPPaqua.csv')


######################
# basic plots

library(ggplot2)
library(RColorBrewer)
modis.terra <- fread('data/modis_npp/NPPTerr.csv', drop = 1)

ggplot(modis.terra, aes(longitude, latitude, color = npp.mean)) +
    borders("world", size = 0.1) +
    geom_point(size = 0.5) +
    scale_color_gradientn(colours = brewer.pal(9, 'YlOrRd'))

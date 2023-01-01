# extract human impact data

require(data.table)
require(raster) # for geotiff
require(ggplot2) # for plotting

# Load Bowler anthropogenic threats
atc <- stack("dataDL/bowler_atcs/Figure_6.grd") # raster stack. Cumulative is all of them together
names(atc) # Climate_change, Human_use, Human_population, Pollution, Alien_potential, Cumulative
atc_cum <- stackApply(atc, indices = c(1,2,2,2,2,1), fun = sum, na.rm = TRUE)[[2]] # sum up non-climate and non-cumulative drivers in layer 2 and return only layer 2
#plot(atc_cum)

# BioTime locations
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin)[!duplicated(rarefyID), .(rarefyID, REALM, rarefyID_x, rarefyID_y)]; rm(bt_malin) # rename to bt

btcoords <- bt[, .(rarefyID_x, rarefyID_y)]
coordinates(btcoords) = ~ rarefyID_x + rarefyID_y # add coordinates
proj4string(btcoords) <- CRS('+init=epsg:4326') # set coordinate system to WGS84 latlong

# Extract from Bowler ATCs
bt3 <- spTransform(btcoords, proj4string(atc_cum)) # convert points to ATC coord system
atcval <- extract(atc_cum, bt3)
bt[, atc := atcval]


# plot maps
ggplot(bt, aes(rarefyID_x, rarefyID_y, color = atc)) +
    geom_point(size = 0.3)

# write out
write.csv(bt[, .(rarefyID, atc)], file = gzfile('output/humanimpact_by_rarefyID.csv.gz'), row.names = FALSE)

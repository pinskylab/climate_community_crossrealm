# extract human impact data

require(data.table)
require(raster) # for geotiff
require(ggplot2) # for plotting

# Load human footprint (Venter and Halpern)
hfp <- raster("dataDL/terrestrial_footprint/HFP1993.tif") # terrestrial human footprint (Venter)
himp <- raster('dataDL/human_impact_marine/model.tif') # marine human impact (Halpern)
proj4string(himp) <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs') # since WGS 1984 Mollweide
extent(himp)
#plot(himp)

# Load Bowler anthropogenic threats
atc <- stack("dataDL/bowler_atcs/Figure_6.grd") # raster stack. Cumulative is all of them together
names(atc)
atc_cum <- stackApply(atc, indices = c(1,2,2,2,2,1), fun = sum, na.rm = TRUE)[[2]] # sum up non-climate drivers
#plot(atc_cum)

# BioTime locations
load('data/biotime_blowes/bt_malin.Rdata')
bt <- data.table(bt_malin)[!duplicated(rarefyID), .(rarefyID, REALM, rarefyID_x, rarefyID_y)]; rm(bt_malin) # rename to bt

btcoords <- bt[, .(rarefyID_x, rarefyID_y)]
coordinates(btcoords) = ~ rarefyID_x + rarefyID_y # add coordinates
proj4string(btcoords) <- CRS('+init=epsg:4326') # set coordinate system to WGS84 latlong


# Extract from HFP
bt1 <- spTransform(btcoords, proj4string(hfp)) # convert points to HFP coord system
hfpval <- extract(hfp, bt1)
bt[, hfp := hfpval]


# Extract from HIMP
# Not working on western and southern hemisphere values
# bt2 <- spTransform(btcoords, proj4string(himp)) # convert to HIMP coord system
# extent(bt2) # doesn't match himp
# extent(himp) 
# 
# himpval <- extract(himp, bt2)
# bt[, himp := himpval]
# 
# bt[, summary(himp)]


# Extract from Bowler ATCs
bt3 <- spTransform(btcoords, proj4string(atc_cum)) # convert points to HFP coord system
atcval <- extract(atc_cum, bt3)
bt[, atc := atcval]


# plot maps
ggplot(bt, aes(rarefyID_x, rarefyID_y, color = hfp)) +
    geom_point(size = 0.3)

ggplot(bt, aes(rarefyID_x, rarefyID_y, color = himp)) +
    geom_point(size = 0.3)

ggplot(bt, aes(rarefyID_x, rarefyID_y, color = atc)) +
    geom_point(size = 0.3)

# plot pairwise
ggplot(bt, aes(hfp, atc)) +
    geom_point(size = 0.3, alpha = 0.5)

ggplot(bt, aes(himp, atc)) +
    geom_point(size = 0.3, alpha = 0.5)

# stats
bt[, cor.test(hfp, atc, use = 'pairwise.complete.obs')]

# write out Bowler
write.csv(bt[, .(rarefyID, atc)], file = gzfile('output/humanimpact_by_rarefyID.csv.gz'))

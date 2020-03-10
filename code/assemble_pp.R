# Assemble temperature dataset

require(data.table)
library(ggplot2)
library(RColorBrewer)

###########################
# Read in and average NPP
###########################

### Read in biotime data for location
load('data/biotime_blowes/bt.Rdata'); bt <- data.table(bt_malin); rm(bt_malin)

# trim to locations
bt <- bt[!duplicated(rarefyID), .(rarefyID, lon = rarefyID_x, lat = rarefyID_y)]

### read in NPP data from land/ocean merge
# units are mgC m-2 day-1
files <- list.files('dataDL/ocean_productivity/', pattern = '.csv')
for(i in 1:length(files)){
    print(i)
    temp <- fread(paste0('dataDL/ocean_productivity/', files[i]))
    if(i == 1){
        npp <- array(data = NA, dim = c(nrow(temp), ncol(temp), length(files)))
        npp[,,i] <- as.matrix(temp)
    }
    if(i > 1) npp[,,i] <- as.matrix(temp)
                 
}

# turn missing to NA
npp[npp == -9999] <- NA

# average by location
nppave <- apply(npp, MARGIN = c(1,2), FUN = mean, na.rm = TRUE)
dim(nppave)

# label rows and columns by lat/lon
nppave <- data.table(nppave)
setnames(nppave, 1:ncol(nppave), as.character(seq(-180, 180 - 360/ncol(nppave), by = 360/ncol(nppave)) + 360/ncol(nppave)/2))
nppave[, lat := seq(90, -90 + 180/nrow(nppave), by = -180/nrow(nppave)) - 180/nrow(nppave)/2]

# long format
nppavel <- melt(nppave, id.vars = 'lat', measure.vars = 1:(ncol(nppave)-1), variable.name = 'lon', value.name = 'npp')
nppavel[, lon := as.numeric(as.character(lon))]

# plot: slow to plot because so many points
ggplot(nppavel, aes(lon, lat, color = npp)) +
    borders("world", size = 0.1) +
    geom_point(size = 0.2) +
    scale_color_gradientn(colours = brewer.pal(9, 'YlOrRd'))

ggsave('figures/npp.png')

# add npp to biotime locations
# npp data are on 1/12 deg grids, upper-left grid corner is at -180°W, 90°N
# round to 5 decimal points to ensure accurate merging
bt[, latgrid := round(floor(lat*12)/12 + 1/24, 3)]
bt[, longrid := round(floor(lon*12)/12 + 1/24, 3)]
nppavel[, latgrid := round(lat, 3)]
nppavel[, longrid := round(lon, 3)]

bt <- merge(bt, nppavel[, .(latgrid, longrid, npp)], by = c('latgrid', 'longrid'), all.x = TRUE)

bt[is.na(npp), .N] # 2180
bt[is.na(npp), hist(latgrid)]
ggplot(bt[is.na(npp),], aes(longrid, latgrid)) +
    borders("world", size = 0.1) +
    geom_point(size = 0.5)

# for the missing data points, fill in as average from surrounding 25 cells
missing <- bt[, which(is.na(npp))]
length(missing)
for(i in 1:length(missing)){
    if(i %% 100 == 0) print(i)
    lats <- bt[missing[i], round(latgrid + c(2/12, 1/12, 0, -1/12, -2/12), 3)]
    lons <- bt[missing[i], round(longrid + c(2/12, 1/12, 0, -1/12, -2/12), 3)]
    bt[missing[i], npp := nppavel[latgrid %in% lats & longrid %in% lons, mean(npp, na.rm = TRUE)]]
}

bt[is.na(npp), .N] # 153 (761 if 9x9)
bt[is.na(npp), hist(latgrid)]
ggplot(bt[is.na(npp),], aes(longrid, latgrid)) +
    borders("world", size = 0.1) +
    geom_point(size = 0.5)

# turn NaNs to NAs
bt[is.nan(npp), npp := NA]

# drop lat/lon columns
bt[, ':='(lat = NULL, lon = NULL, latgrid = NULL, longrid = NULL)]

# write out
write.csv(bt, file = gzfile('output/npplandocean.csv.gz'), row.names = FALSE)




######################
# basic plots


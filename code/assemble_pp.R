# Assemble temperature dataset

require(data.table)

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

# remove NAs
nppavel <- nppavel[!is.na(npp),]

# write out
write.csv(nppavel, file = gzfile('data/npp/landoceannpp.csv.gz'), row.names = FALSE)



######################
# basic plots

library(data.table)
library(ggplot2)
library(RColorBrewer)
nppavel <- fread('data/npp/landoceannpp.csv.gz')

# slow to plot because so many points
ggplot(nppavel, aes(lon, lat, color = npp)) +
    borders("world", size = 0.1) +
    geom_point(size = 0.2) +
    scale_color_gradientn(colours = brewer.pal(9, 'YlOrRd'))

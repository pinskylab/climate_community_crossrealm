###Laura Antao
###06-05-2020; updated 13-05-2020

##script to send to Malin extracting temperature values from Wordclim to use in turnover~Year analysis from the BioTIME gridded data


# DISCLAIMER: code is under development
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##



##18-05-2020
##new approach - get values within radius around each central lat-long

library(tidyverse)
library(raster)
library(sdmpredictors)
library(ggplot2) # for a plot at the end
library(maps) # for a plot at the end

##load raster files as before
tempbio1 <- raster("wc2.0_bio_30s_01.tif")   ##downloaded from site
temp_sstmean <- load_layers("BO_sstmean")

# load BioTime change
load('data/biotime_blowes/bt_malin.Rdata')
trends <- bt_malin[!duplicated(bt_malin$rarefyID), c('REALM', 'STUDY_ID', 'rarefyID', 'rarefyID_x', 'rarefyID_y')]; rm(bt_malin) # rename to bt

##for terrestrial & freshwater
##get coords from rarefyIDs ####
CoordsTer <- trends %>%
  filter(REALM=="Terrestrial"| REALM=="Freshwater") %>%
  dplyr::select(rarefyID, rarefyID_x, rarefyID_y) %>%
  unique()


ptsTer <- SpatialPoints(CoordsTer[,2:3])

##extract values ####
##setting a buffer with 6Km radius (~113km2, so a bit larger than grid cell of ~96km2)
start_time <- Sys.time()

tempvaluesTer <- extract(tempbio1, ptsTer, buffer=6000, fun=sd)

end_time <- Sys.time()

time<- end_time - start_time  ##~25 minutes to extract

Temp_sd_Ter <- cbind.data.frame(coordinates(ptsTer), tempvaluesTer) %>%
  mutate(rarefyID = CoordsTer$rarefyID)

#write.csv(Temp_sd_Ter, "Temp_sd_Ter.csv", row.names = F)

nrow(Temp_sd_Ter[is.na(Temp_sd_Ter$tempvaluesTer),])  ##4 records without sd values



##for marine ####
##get coords from rarefyIDs ####
CoordsMar <- trends %>%
  filter(REALM=="Marine") %>%
  dplyr::select(rarefyID, rarefyID_x, rarefyID_y) %>%
  unique()

ptsMar <- SpatialPoints(CoordsMar[,2:3])


##extract values ####
##setting a buffer with 6km radius
start_time <- Sys.time()

tempvaluesMar <- extract(temp_sstmean, ptsMar, buffer=6000, fun=sd)

end_time <- Sys.time()

time<- end_time - start_time  ##~4h to extract!! with error first time... all data lost...
# Error in dimnames(x) <- dn : 
#   length of 'dimnames' [2] not equal to array extent
# In addition: Warning message:
#   In matrix(unlist(cv, use.names = FALSE), nrow = np, byrow = TRUE) :
#   data length [94701] is not a sub-multiple or multiple of the number of rows [49036]



###--> turning into a loop so not all the values are lost!! and adding try() to skip the errors
tempvaluesMar<- matrix(NA,length(ptsMar),2)

start_time <- Sys.time()

u<-1

for (u in 1:length(ptsMar)) {
  
  try({
    
    aux<- extract(temp_sstmean, ptsMar[u], buffer=6000)  ##this is a list
    
    tempvaluesMar[u,]<-c(sd(aux[[1]]), length(aux[[1]]))
    ##how many points were used to calculate the sd value also being saved
    
    #cat(u)
    u<- u+1
    
  })
  
}

end_time <- Sys.time()

#write.csv(tempvaluesMar, "tempvaluesMar.csv", row.names = F)


Temp_sd_Mar <- cbind.data.frame(coordinates(ptsMar), tempvaluesMar) %>%
  mutate(rarefyID = CoordsMar$rarefyID) %>%
  rename(Temp_sd ="1", numberPts_sd = "2")

#write.csv(Temp_sd_Mar, "Temp_sd_Mar.csv", row.names = F)


nrow(Temp_sd_Mar[is.na(Temp_sd_Mar$Temp_sd),])  ##12590 records without sd values (~25%...)




# save(trends, CoordsTer, CoordsMar, ptsTer, ptsMar,
#      tempvaluesTer, tempvaluesMar,
#      Temp_sd_Ter, Temp_sd_Mar,
#      file = "temp_extract_sd.Rdata")




##combine both dataframes to then join to the main data
Temp_sd6_ALL<- rbind(Temp_sd_Ter%>%
                      rename(Temp_sd6km ="tempvaluesTer"),
                    Temp_sd_Mar[,-4]%>%
                      rename(Temp_sd6km ="Temp_sd"))








##added 19-05-2020 ####
##to use a larger buffer hoping to overcome having so many NAs

##for terrestrial
##setting a buffer with 20Km radius (~1256.64 km2)
start_time <- Sys.time()

tempvaluesTer20 <- extract(tempbio1, ptsTer, buffer=20000, fun=sd, na.rm=T)

end_time <- Sys.time()

time<- end_time - start_time  ##~20 minutes to extract

Temp_sd_Ter20 <- cbind.data.frame(coordinates(ptsTer), tempvaluesTer20) %>%
  mutate(rarefyID = CoordsTer$rarefyID)

#write.csv(Temp_sd_Ter20, "Temp_sd_Ter20.csv", row.names = F)

nrow(Temp_sd_Ter20[is.na(Temp_sd_Ter20$tempvaluesTer),])  ##2 records without sd values




##for marine
tempvaluesMar20<- matrix(NA,length(ptsMar),2)

start_time <- Sys.time()

u<-1

for (u in 1:length(ptsMar)) {
  
  try({
    
    aux<- raster::extract(temp_sstmean, ptsMar[u], buffer=20000)  ##this is a list
    
    tempvaluesMar20[u,]<-c(sd(aux[[1]], na.rm = T), length(aux[[1]]))
    ##how many points were used to calculate the sd value also being saved
    
    cat(u)
    u<- u+1
    
  })
  
}

end_time <- Sys.time()
time<- end_time - start_time    ##~3h to extract

#write.csv(tempvaluesMar20, "tempvaluesMar20.csv", row.names = F)


Temp_sd_Mar20 <- cbind.data.frame(coordinates(ptsMar), tempvaluesMar20) %>%
  mutate(rarefyID = CoordsMar$rarefyID) %>%
  rename(Temp_sd ="1", numberPts_sd = "2")

#write.csv(Temp_sd_Mar20, "Temp_sd_Mar20.csv", row.names = F)


nrow(Temp_sd_Mar20[is.na(Temp_sd_Mar20$Temp_sd),])  ##1203 records without sd values (~2%...)






##combine both dataframes to then join to the main data
##object not saved in case still want to modify individual data
Temp_sd20_ALL<- rbind(Temp_sd_Ter20%>%
                        rename(Temp_sd20km ="tempvaluesTer20"),
                      Temp_sd_Mar20[,-4]%>%
                        rename(Temp_sd20km ="Temp_sd"))

#save(trends, CoordsTer, CoordsMar, ptsTer, ptsMar,
     
     #buffer=6
#     tempvaluesTer, tempvaluesMar, Temp_sd_Ter, Temp_sd_Mar,
     
     #buffer=20
#     tempvaluesTer20, Temp_sd_Ter20, tempvaluesMar20, Temp_sd_Mar20,
#     file = "temp_extract_sd.Rdata")


# combine 6km and 20km buffers
Temp_sd <- merge(Temp_sd6_ALL[, c('rarefyID', 'rarefyID_x', 'rarefyID_y', 'Temp_sd6km')], Temp_sd20_ALL[, c('rarefyID', 'Temp_sd20km')], by = 'rarefyID', all = TRUE)

# write out
write.csv(Temp_sd, gzfile('output/microclimates.csv.gz'))

# read in, if needed
#Temp_sd <- read.csv(gzfile('output/microclimates.csv.gz'))

# make a plot of 20 km SD
world <- map_data('world')
ggplot(world, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = 'lightgray', color = 'white', size = 0.1) +
  geom_point(data = Temp_sd, aes(rarefyID_x, rarefyID_y, color = Temp_sd20km, group = NA), size = 0.2, alpha = 0.4) +
  scale_color_gradientn(colours = c('blue', 'white', 'red'))

# make a plot of missing data
world <- map_data('world')
ggplot(world, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = 'lightgray', color = 'white', size = 0.1) +
  geom_point(data = Temp_sd, aes(rarefyID_x, rarefyID_y, color = !is.na(Temp_sd20km), group = NA), size = 0.2, alpha = 0.4)

# look at missing data by realm
Temp_sd2 <- merge(Temp_sd, trends[, c('rarefyID', 'REALM')], by= 'rarefyID')

ggplot(world, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = 'lightgray', color = 'white', size = 0.1) +
  geom_point(data = subset(Temp_sd2, REALM=='Marine'), aes(rarefyID_x, rarefyID_y, color = !is.na(Temp_sd20km), group = NA), size = 0.2, alpha = 0.4) +
  coord_cartesian(xlim = c(-130, -100), ylim = c(30,45)) # what are they doing on land?



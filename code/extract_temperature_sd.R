
###Laura Antao
###06-05-2020

##script to send to Malin extracting temperature values from Wordclim to use in turnover~Year analysis from the BioTIME gridded data


# DISCLAIMER: code is under development

##======================================================================

library(tidyverse)
library(raster)
library(dggridR)
#library(scales)

tempbio1 <- raster("wc2.0_bio_30s_01.tif")   ##downloaded from site

bio1Extent <- extent(tempbio1)

##extract each pixel value for tempbio1 into a dataframe
bio1_extract <- raster::extract(tempbio1, bio1Extent, df=TRUE, cellnumbers=TRUE)

##create a data frame with the coordinates of each cell and then bind
tempbio1_coords <- as.data.frame(xyFromCell(tempbio1, bio1_extract[,2]))
tempbio1_ALL <- cbind(tempbio1_coords, bio1_extract[,3])



## steps to link to the grid
load("~BioTIME_grid_filtered_011017.Rdata") 
##this is the Rdata object that Shane sent us, keeping only the dgg object to use the same specifications
##dggs_type=="ISEA3H"; dggs_res_spec==12; precision==7

rm(bt_grid_filtered)


##replicating the same processing steps as for the biodiversity data (https://github.com/sChange-workshop/BioGeo-BioDiv-Change/blob/master/R/01_Study_to_Grid.R)
##i.e. assign the corresponding grid cells for all temperature observations
tempbio1_ALL <- tempbio1_ALL %>%
  mutate(cell = dgtransform(dgg, lat=lat_to_grid, lon=lon_to_grid))  ###   !!!!! have to change the names of the coord variables here once we get them!!!!! ###


##recover the cell ID from the original rarefyIDs, filter only the relevant temperature values and merge
trends1<- trends %>%
  separate(., rarefyID, into = c("STUDY_ID1", "cell"), remove=F) %>%
  dplyr::select(-STUDY_ID1)


tempbio1_data <- tempbio1_ALL %>%
  filter(cell %in% trends1$cell) %>%
  group_by(cell) %>%
  summarise(sd_temp= sd(xxxxxxx))   ###xxxxxxx== name of the variable


###join data to the original trends object
trends<- left_join(trends,
                   select(tempbio1_data, cell, sd_temp),  ##to avoid duplicating the coordinates
                   by= "cell")


##======================================================================

##then repeat for marine data from Bio-Oracle
library(sdmpredictors)










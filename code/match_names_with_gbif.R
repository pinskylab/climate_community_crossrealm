rm(list = ls())
library(taxize)
library(tidyverse)
library(rstudioapi)
library(data.table)
setwd(dirname(getActiveDocumentContext()$path))


# select the databse to use
db = 'gbif'

trace(get_gbifid, edit=TRUE)
# replace line 91 by: 
# take = which(df$status == 'ACCEPTED' & df$matchtype == 'EXACT')

# it's possible to include correspondances at higher ranks, but the classification function behavve weirdly in such cases
# if the species is not found, it seems that it return kingdom names to put so;ething at for the species name
# there is also a possibility that several species are found, leading to issues. 
# so adding the next lines is not recomended for now
# if (length(take) == 0) {
#   take = which(df$status == "ACCEPTED" &
#                  df$matchtype == "HIGHERRANK")
# }
# I also commented lines printing messages for performance
#  !!! these changes are not permanant, so it needs to be redone every run.... !!!!

species.list = c('Vulpes vulpes')

get.taxize.classification = function(species.list, db){
  # get the txonomic levels returned by db (sshould exists a better apporach...)
  ranks = classification('Vulpes vulpes', db = db)[[1]][,2]
  gbif.list = classification(species.list, db = db)
  gbif2<-gbif.list %>% 
    map(1)
  
  foo = function(x, ranks){
    if (is.na(x) | length(x) == 1){
      x = rep(NA, 7)
    }
    return(x)
  }
  gbif3 = lapply(gbif2, foo, ranks)
  gbif.df = as.data.frame(do.call(rbind.data.frame, gbif3))
  names(gbif.df) = ranks
  gbif.df$original_name = species.list
  
  return(gbif.df)
}

load('climate_community_crossrealm/data/biotime_blowes/bt_grid_spp_list.Rdata')
species.level = which(tstrsplit(species.list, ' ')[[2]] != 'sp')
phylo2 = get.taxize.classification(species.list[species.level[1:34]], db)

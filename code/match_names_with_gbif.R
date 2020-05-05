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
# take = which(df$status == "ACCEPTED" & df$matchtype == 
#                "EXACT")
# if (length(take) > 1) {
#   print(df)
#   take = take[1]
#   sink("multiple_accepted_names.txt", append = TRUE)
#   cat('asked name: ', sciname[i], '\n')
#   print(df)
#   cat('-------------------------------------------------------------\n')
#   sink()
# }

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

# species.list = c('Vulpes vulpes')

get.taxize.classification = function(species.list, db){
  # get the txonomic levels returned by db (should exists a better apporach...)
  ranks = classification('Vulpes vulpes', db = db)[[1]][,2]
  sink("multiple_accepted_names.txt")
  sink()
  gbif.list = classification(species.list, db = db)
    gbif2<-gbif.list %>% 
    map(1)
  
  # gbif.names = names(lapply(gbif.list, names))
  
  foo = function(x, ranks){
    # if species was not found or if the taxonomic classification is not complete 
    # (sometimes classification only returns kingdom and genus)
    if (is.na(x) | length(x) < 7){
      x = rep(NA, 7)
    }
    return(x)
  }
  gbif3 = lapply(gbif2, foo, ranks)
  gbif.df = as.data.frame(do.call(rbind.data.frame, gbif3))
  names(gbif.df) = ranks
  # gbif.df$original_name = gbif.names
  gbif.df$original_name = species.list
  return(gbif.df)
}

# classification(species.list[species.level[5]], db = db)

load('../data/biotime_blowes/bt_grid_spp_list.Rdata')
species.names = gsub('_', ' ', unique(bt_grid_spp_list$Species))

species.level = which(tstrsplit(species.names, ' ')[[2]] != 'sp')

phylo2 = get.taxize.classification(species.names[species.level], db)

tot.differences = sum(phylo2$species != phylo2$original_name  & !is.na(phylo2$species))

write.csv(phylo2, file=gzfile('../output/taxonomy.csv.gz'))



  # --------------------------Extract BM from gateway for new names  --------------------------------------------------   
# phylo2 = read.csv('../output/taxonomy.csv.gz')
db= read.csv('../dataDL/GATEWay_brose/283_2_FoodWebDataBase_2018_12_10.csv', header = TRUE, sep = ",")

# get species list form both dataset
species.gateway = unique(c(as.character(db$con.taxonomy), as.character(db$res.taxonomy)))
res.mass = db[,c("res.taxonomy", "res.mass.mean.g.")]; names(res.mass) = c('species.names', 'mass')
con.mass = db[,c("con.taxonomy", "con.mass.mean.g.")]; names(con.mass) = c('species.names', 'mass')
db.masses = unique(rbind.data.frame(res.mass, con.mass))

# some species bodymasses have been coded as -999 (I guess in absence of information). That might change species' average bodymasses...
db.masses = subset(db.masses, mass>0)

average_species = function(species, gateway){
    species.bms = gateway$mass[which(gateway$species.names %in% species)]
    if (length(species.bms) == 0) {
        return(NA)
    }else{
        return(mean(species.bms, na.rm = TRUE))
    }
}

bms.species = sapply(phylo2$species, average_species, db.masses, USE.NAMES = FALSE)

phylo2$gateway.masses = bms.species
sum(!is.na(bms.species))
sum(!is.na(bms.species) & !(phylo2$species %in% phylo2$original_name))
sum(!(phylo2$species %in% phylo2$original_name))
write.csv(phylo2, file=gzfile('../output/taxonomy.csv.gz'))


    # --------------------------adding information in mass.gaeway.csv  --------------------------------------------------

bm.gateway = read.csv("../output/mass_gateway.csv.gz")
bm.gateway$name_used = "Blowes"

gate.in.phylo =  match(bm.gateway$Species, phylo2$original_name)
# checks
# cbind(bm.gateway$mass_species, phylo2$gateway.masses[gate.in.phylo])

# get the BM from phylo2 in the order corresponding to the bm.gateway
species.m.from.phylo2 = phylo2$gateway.masses[gate.in.phylo]
bm.gateway$gbif_names = phylo2$species[gate.in.phylo]

# reaplce only values for which information is missing using non corrected names
to.replace = is.na(bm.gateway$mass_species) & !is.na(species.m.from.phylo2)
# sum(to.replace)
# sum(!(bm.gateway$Species %in% bm.gateway$gbif_names)  & !is.na(species.m.from.phylo2))

# replacing values
bm.gateway$mass_species[to.replace] = species.m.from.phylo2[to.replace]

# using genus estimated masses only for the remaining missing species masses
bm.gateway$mass = bm.gateway$mass_species
bm.gateway$mass[is.na(bm.gateway$mass)] = bm.gateway$mass_genus[is.na(bm.gateway$mass)]

# information on the name used
bm.gateway$name_used[to.replace] = 'gbif'

write.csv(bm.gateway, file=gzfile('../output/mass_gateway.csv.gz'), row.names = FALSE)



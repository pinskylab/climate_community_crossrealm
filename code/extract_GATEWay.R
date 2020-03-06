################################################################
#######  extract BM from Uli's gateway database  ###############
################################################################
rm(list = ls())
# make relative paths
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

db= read.csv('../dataDL/GATEWay_brose/283_2_FoodWebDataBase_2018_12_10.csv', header = TRUE, sep = ",")
load('../data/biotime_blowes/bt_grid_spp_list.Rdata')

# get species list form both dataset
species.gateway = unique(c(as.character(db$con.taxonomy), as.character(db$res.taxonomy)))
genus.gateway = strsplit(species.gateway, ' ')
genus.gateway = sapply(genus.gateway, `[[`, 1)

species.blowes = unique(bt_grid_spp_list$Species)
species.blowes.space = gsub('_', ' ', species.blowes)
genus.blowes = strsplit(species.blowes.space, ' ')
genus.blowes = sapply(genus.blowes, `[[`, 1)


# check for dbs 'proximity'
sum(species.blowes.space %in% species.gateway)
sum(species.gateway %in% species.blowes.space)

sum(genus.blowes %in% genus.gateway)



res.mass = db[,c("res.taxonomy", "res.mass.mean.g.")]; names(res.mass) = c('species.names', 'mass')
con.mass = db[,c("con.taxonomy", "con.mass.mean.g.")]; names(con.mass) = c('species.names', 'mass')

db.masses = unique(rbind.data.frame(res.mass, con.mass))

genus.gateway = strsplit(as.character(db.masses$species.names), ' ')
genus.gateway = sapply(genus.gateway, `[[`, 1)

db.masses$genus.names = genus.gateway
db.masses = subset(db.masses, mass<=0)

#####################################################
############## incorporate BMs ######################
#####################################################

blowes.bms = rep(NA, length(species.blowes))

# 1) make attribution based on exact species match
average_species = function(species, gateway){
    species.bms = gateway$mass[which(gateway$species.names == species)]
    if (length(species.bms) == 0) {
        return(NA)
    }else{
        return(mean(species.bms, na.rm = TRUE))
    }
}


blowes.bms.species = sapply(species.blowes.space, average_species, db.masses, USE.NAMES = FALSE)

# 2) make attribution based on exact genus match (if no info from species level)

average_genus = function(genus, gateway){
    genus.bms = gateway$mass[which(gateway$genus.names == genus)]
    if (length(genus.bms) == 0) {
        return(NA)
    }else{
        return(mean(genus.bms, na.rm = TRUE))
    }
}

blowes.bms.genus = sapply(genus.blowes, average_genus, db.masses, USE.NAMES = FALSE)

# -------------------------------------------------------------------------------------------
# -------------------------------------- DONE -----------------------------------------------
# -------------------------------------------------------------------------------------------



# association with priority (first species level, then genus level)

bms.final = blowes.bms.species
# if species level biomass is na, take the genus level information (even if NA)
bms.final = bms.final[!is.na(bms.final)] = blowes.bms.genus[!is.na(bms.final)]



####################### add biomass information to the data: ######################################
bt_grid_spp_list$mass = bms.final




# -------------------------------------------------------------------------------------------
# -------------------------------------- DONE -----------------------------------------------
# -------------------------------------------------------------------------------------------


#####################################################
############## Fuzzy matching  ######################
#####################################################


# 4) check for potential typos // renaming in the databases
# the agrep function seems bugged. use adist to manually make distance comparisons afterward.
is.it.almost.in.genus = function(g.name, gateway, max.distance){
    dists = adist(g.name, gateway)
    # i might need to add a parameter to use the length of the string (to be more permissive with large strings)
    a = which(dists < max.distance)
    if (length(a)>0){
        return(cbind(rep(g.name), gateway[a]))
    }
    else{
        return(c(NULL))
    }
}


# use parrallel computing as the procedure can be long
library("future.apply")
nb.threads = 6

# don't know if there is a default behaviour here, but a moderate value for chunks might increase speed
chunk.size = 50 

# maximum distance for the fuzzy matching
max.dists = 3

plan(multiprocess, workers = nb.threads)
almost.in = future_sapply(genus.blowes, FUN = is.it.almost.in.genus, db.masses$genus.names, max.distance = max.dists,  future.chunk.size = chunk.size)
# should be possible to directly obtain an array from future_sapply, but here it returns a list.... (maybe beacuse of return(NULL) in the function)
# so the list needs to be changed in an array
almost.in = do.call(rbind.data.frame, almost.in)
row.names(almost.in) = NULL
names(almost.in) = c("blowes", "gateway")
almost.in = subset(almost.in, !(as.character(blowes) == gateway))
almost.in = unique(almost.in)

library(rbenchmark)
benchmark(is.it.almost.in.genus(genus.blowes[1], db.masses$genus.names), replications = 1000)

yy = adist(genus.blowes[1], db.masses$genus.names)
xx = agrep(genus.blowes[1], db.masses$genus.names, max.distance = 3)
tt = aregexec(genus.blowes[1], db.masses$genus.names, max.distance = 1)
length(xx)


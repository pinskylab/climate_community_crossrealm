################################################################
#######  extract BM from Uli's gateway database  ###############
################################################################
rm(list = ls())
# make relative paths
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

db= read.csv('dataDL/GATEWay_brose/283_2_FoodWebDataBase_2018_12_10.csv', header = TRUE, sep = ",")
load('data/biotime_blowes/bt_grid_spp_list.Rdata')

# get species list form both dataset
species.gateway = unique(c(as.character(db$con.taxonomy), as.character(db$res.taxonomy)))
genus.gateway = strsplit(species.gateway, ' ')
genus.gateway = sapply(genus.gateway, `[[`, 1)

species.blowes = bt_grid_spp_list$Species
species.blowes.space = gsub('_', ' ', species.blowes)
genus.blowes = strsplit(species.blowes.space, ' ')
genus.blowes = sapply(genus.blowes, `[[`, 1)





res.mass = db[,c("res.taxonomy", "res.mass.mean.g.")]; names(res.mass) = c('species.names', 'mass')
con.mass = db[,c("con.taxonomy", "con.mass.mean.g.")]; names(con.mass) = c('species.names', 'mass')

db.masses = unique(rbind.data.frame(res.mass, con.mass))

genus.gateway = strsplit(as.character(db.masses$species.names), ' ')
genus.gateway = sapply(genus.gateway, `[[`, 1)

db.masses$genus.names = genus.gateway

# some species bodymasses have been coded as -999 (I guess in absence of information). That might change species' average bodymasses...
db.masses = subset(db.masses, mass>0)


# check for dbs 'proximity'
sum(species.blowes.space %in% db.masses$species.names)
sum(species.blowes.space %in% species.gateway)
db.masses[db.masses$species.names %in% species.blowes.space,]

sum(genus.blowes %in% genus.gateway)


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

# a bit brutal to make that on duplicated species, but still ok
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


# association with priority (first species level, then genus level)

bms.final = blowes.bms.species
# if species level biomass is na, take the genus level information (even if NA)
bms.final[is.na(bms.final)] = blowes.bms.genus[is.na(bms.final)]

# number of values from a match at the species level
sum(!is.na(blowes.bms.species))
# number of values from a match at the species level OR at genus level
sum(!is.na(bms.final))


####################### add biomass information to the data: ######################################
bt_grid_spp_list$mass_gateway = bms.final

# -------------------------------------------------------------------------------------------
# -------------------------------------- DONE -----------------------------------------------
# -------------------------------------------------------------------------------------------

final = data.frame(Species = species.blowes.space,
                  rarefyID = bt_grid_spp_list$rarefyID,
                  REALM = bt_grid_spp_list$REALM,
                  STUDY_ID = bt_grid_spp_list$STUDY_ID,
                  spp_fb2 = species.blowes.space,
                  mass_genus = blowes.bms.genus,
                  mass_species = blowes.bms.species,
                  mass = bms.final)
final = subset(final, !is.na(mass))

write.csv(final, file = gzfile("output/mass_gateway.csv.gz"))

#####################################################
############## Fuzzy matching  ######################
#####################################################


# 4) check for potential typos // renaming in the databases
# the agrep function seems bugged. use adist to manually make distance comparisons afterward.
is.it.almost.in.genus = function(g.name, gateway, max.distance){
    dists = adist(g.name, gateway)
    # i might need to add a parameter to use the length of the string (to be more permissive with large strings)
    a = which(dists < max.distance & dists > 0)
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

# don't know if there is a default behaviour here, but a moderate value for chunks should increase speed
chunk.size = 50 
# (it seems to reduce cpu (each core is not working at 100%) usage though, so maybe not the optimal approach?)

# maximum distance for the fuzzy matching
max.dists = 3

plan(multiprocess, workers = nb.threads)
almost.in = future_sapply(unique(genus.blowes), FUN = is.it.almost.in.genus, unique(db.masses$genus.names), max.distance = max.dists,  future.chunk.size = chunk.size)
# should be possible to directly obtain an array from future_sapply, but here it returns a list.... (maybe beacuse of return(NULL) in the function)
# so the list needs to be changed in an array
almost.in = do.call(rbind.data.frame, almost.in)
row.names(almost.in) = NULL
names(almost.in) = c("blowes", "gateway")
# simplify results (are several species of same genus in gateway)
almost.in = unique(almost.in)

library(rbenchmark)
benchmark(is.it.almost.in.genus(genus.blowes[1], db.masses$genus.names), replications = 1000)


#########################################################################################################################
# --------------------------Extract BM from gateway for new names  --------------------------------------------------   
# uses new scientific names found by match_names_with_gbif.R
#########################################################################################################################
phylo2 = read.csv('output/taxonomy.csv.gz')
db= read.csv('dataDL/GATEWay_brose/283_2_FoodWebDataBase_2018_12_10.csv', header = TRUE, sep = ",")

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


# --------------------------adding information in mass.gaeway.csv  --------------------------------------------------

bm.gateway = read.csv("output/mass_gateway.csv.gz")
bm.gateway$name_used = "Blowes"

gate.in.phylo =  match(bm.gateway$Species, phylo2$original_name)
# checks
# cbind(bm.gateway$mass_species, phylo2$gateway.masses[gate.in.phylo])

# get the BM from phylo2 in the order corresponding to the bm.gateway
species.m.from.phylo2 = phylo2$gateway.masses[gate.in.phylo]
bm.gateway$gbif_names = phylo2$species[gate.in.phylo]

# repalce only values for which information is missing using non corrected names
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

write.csv(bm.gateway, file=gzfile('output/mass_gateway.csv.gz'), row.names = FALSE)




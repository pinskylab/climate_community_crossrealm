################################################################
#######  extract BM for Uli's gateway database  3###############
################################################################

# make relative paths
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

db= read.csv('../dataDL/GATEWay_brose/283_2_FoodWebDataBase_2018_12_10.csv', header = TRUE, sep = ",")
load('../data/biotime_blowes/bt_grid_spp_list.Rdata')

# get species list form both dataset
species.gateway = unique(c(as.character(db$con.taxonomy), as.character(db$res.taxonomy)))
species.blowes = unique(bt_grid_spp_list$Species)
species.blowes.space = gsub('_', ' ', species.blowes)
# check for dbs 'proximity'
sum(species.blowes.space %in% species.gateway)
sum(species.gateway %in% species.blowes.space)



#####################################################
####### check for spelling diffeences ###############
#####################################################
is.it.almost.in = function(sp.name, gateway){
    a = agrep(sp.name, species.gateway)
    if (length(a)>0){
      return(cbind(rep(sp.name), a, species.gateway[a]))
    }
    else{
        return(c(NULL))
    }
}


# can be long to run (parrallel version below)
# almost.in = lapply(species.blowes.space, FUN = is.it.almost.in, species.gateway)
# almost.in = do.call(rbind.data.frame, almost.in)
# names(almost.in) = c("blowes", "index.gateway", "gateway")
# almost.in = subset(almost.in, !(as.character(blowes) == gateway))


# same but using multiple cores
library("future.apply")
plan(multiprocess, workers = 6)
almost.in = future_sapply(species.blowes.space, FUN = is.it.almost.in, species.gateway)
almost.in = do.call(rbind.data.frame, almost.in)
row.names(almost.in) = NULL
names(almost.in) = c("blowes", "index.gateway", "gateway")
almost.in = subset(almost.in, !(as.character(blowes) == gateway))

################################################################
#######  extract BM from the VecTraits database  ###############
################################################################
rm(list = ls())
# make relative paths
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

db = read.csv('../dataDL/VecTraits/vectraits.csv', header = TRUE)
dim(db)
load('../data/biotime_blowes/bt_grid_spp_list.Rdata')

species.blowes = bt_grid_spp_list$Species
species.blowes.space = gsub('_', ' ', species.blowes)
genus.blowes = strsplit(species.blowes.space, ' ')
genus.blowes = sapply(genus.blowes, `[[`, 1) # don't make a unique on that, important to keep correspondences with species list


# unique(db$published_data.interactor1sizeunitsi)

grep('1', names(db), value = TRUE)
grep('mass', names(db), value = TRUE)

# create a databse with mass (or with size in kg) info for interactor 1 and 2
# then rbind them in db2
sp1 = db[, c('published_data.interactor1massunitsi', 
             'published_data.interactor1sizeunit',
             'published_data.interactor1sizeunitsi',
             "published_data.interactor1size",
             "published_data.interactor1sizesi",
             "published_data.interactor1massvaluesi",
             'published_data.interactor1genus',
             "published_data.interactor1species"
                             )]
names(sp1) = c('mass_units_i', 'size_units', 'size_units_i', 'size', 'sizei', 'massi', 'genus', 'species')


sp2 = db[, c('published_data.interactor2massunitsi', 
            'published_data.interactor2sizeunit',
            'published_data.interactor2sizeunitsi',
            "published_data.interactor2size",
            "published_data.interactor2sizesi",
            "published_data.interactor2massvaluesi",
            'published_data.interactor2genus',
            "published_data.interactor2species"
                              )]
names(sp2) = c('mass_units_i', 'size_units', 'size_units_i', 'size', 'sizei', 'massi', 'genus', 'species')

db2 = rbind.data.frame(sp1, sp2)

# look at all possible untis
c(unique(as.character(db2$mass_units_i)), '///', unique(as.character(db2$size_units)), '///', unique(as.character(db2$size_units_i)))

# keep lines with information on species bodymasses
to.keep = db2$mass_units_i %in% 'kilogram (wet body mass)' |
  db2$size_units %in% c('kilogram', 'kg') |
  db2$size_units_i %in% c('kilogram', 'kilograms')

db2 = subset(db2, to.keep)
db2 = unique(db2)
dim(db2)

species = paste(db2$genus, db2$species, sep = ' ')

db2[species %in% c("Thalassiosira weissflogii", "Pinus radiata", "Ditylum brightwellii"), ]
subset(db2, species %in% "Pinus radiata")
sort(unique(species))

c(unique(as.character(db2$mass_units_i)), unique(as.character(db2$size_units_i)))
c(unique(as.character(db2$size_units)))

get.species.info = function(blowes, vectraits, sp.list){
  acceptable.indices = c("kilogram (wet body mass)", "kilogram (wet tissue mass)", "kilogram", "kilograms")
  index = which(sp.list %in% blowes)
  # if no info found return null
  if (length(index) == 0){return(NA)}
  
  # If I have some measurements at the individual level, take them
  b = FALSE
  if (sum(!is.na(vectraits$mass_units_i[index])) > 0){
    i.mass = vectraits$massi[index]
    b = TRUE
  }
  if (sum(!is.na(vectraits$size_units_i[index])) > 0){
    i.mass2 =vectraits$sizei[index]
    b = TRUE
  }
  
  if (b) {return(mean(c(i.mass, i.mass2)))}
  
  # if this line is run, means that some info were found but not at individual level 
  return(mean(vectraits$size[index]))
}

################# species level extraction #############################
bms.species = vapply(species.blowes.space, FUN = get.species.info, FUN.VALUE = double(1), db2, species)

# % checks
sum(!is.na(bms.species))
bms.species[!is.na(bms.species)]
db2[db2$species %in% 'meyeri',]

bms.species[species.blowes.space %in% 'Protonemura  meyeri']
################# genus level extraction #############################
# 
db.genus = unique(db2[, - which(names(db2) == 'species')])
genus.vectraits = db.genus$genus
bms.genus = vapply(genus.blowes, FUN = get.species.info, FUN.VALUE = double(1), db.genus, genus.vectraits)

    # checks
which(!is.na(bms.genus))
sum(!is.na(bms.genus))
bms.genus[!is.na(bms.genus)]
db2[db2$genus %in% 'Pinus',]
blowes = genus.blowes[12]



# association with priority (first species level, then genus level)

bms.species = bms.species * 1000
bms.genus = bms.genus * 1000

bms.final = bms.species
# if species level biomass is na, take the genus level information (even if NA)
bms.final[is.na(bms.final)] = bms.genus[is.na(bms.final)]

# number of values from a match at the species level
sum(!is.na(bms.species))
# number of values from a match at the species level OR at genus level
sum(!is.na(bms.final))


# use of gbif corrected names

gbif = read.csv(gzfile("../output/taxonomy.csv.gz"))

gbif.names = paste(gbif$genus, gbif$species, sep = ' ')

xx = match(gbif$original_name, species.blowes.space)
# sum(xx != 1:length(species.blowes.space)) 
# sum(species.blowes.space[xx] == gbif$original_name)
# length(xx)

bms.gbif = vapply(gbif.names, FUN = get.species.info, FUN.VALUE = double(1), db2, species)
# length(bms.gbif) 
# length(gbif.names)

# is there something new?
# check if I have some info from gbif where I had no info with uncorrected names
bms.final[xx][is.na(bms.final[xx]) & !is.na(bms.gbif)]
length(is.na(bms.final[xx]))
# Nothing more, then stick to the old results

final = data.frame(Species = species.blowes.space,
                   rarefyID = bt_grid_spp_list$rarefyID,
                   REALM = bt_grid_spp_list$REALM,
                   STUDY_ID = bt_grid_spp_list$STUDY_ID,
                   spp_fb2 = species.blowes.space,
                   mass_genus = bms.genus,
                   mass_species = bms.species,
                   mass = bms.final)

final = subset(final, !is.na(mass))

write.csv(final, file = gzfile("../output/mass_BioTrait.csv.gz"))


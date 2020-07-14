#----------------------------------------------------------------------------------------------------------#
# 1. access required data
#----------------------------------------------------------------------------------------------------------#

load('data/biotime_blowes/bt_malin.Rdata') # load bt_malin

# find VCT layers (NOTE: dependent on MAS lab data structure)
f1 = list.files('/data/idiv_meyer/00_data/processed/VCF5KYR/', 'TreeCover.*.tif$', full.names=T)
f2 = list.files('/data/idiv_meyer/00_data/processed/VCF5KYR/', 'NonTreeVegetation.*.tif$', full.names=T)
f3 = list.files('/data/idiv_meyer/00_data/processed/VCF5KYR/', 'NonVegetated.*.tif$', full.names=T)

# read raster data
r1 = terra::rast(f1) # % tree cover
r2 = terra::rast(f2) # % non-tree veg. cover
r3 = terra::rast(f3) # % non-vegetated cover

#----------------------------------------------------------------------------------------------------------#
# 2. extract reference variables
#----------------------------------------------------------------------------------------------------------#

# estract coordinates (addresses format conflicts during VCT extraction)
xy = data.frame(x=bt_malin$rarefyID_x, y=bt_malin$rarefyID_y)

# unique years with biodiversity data
uy = unique(bt_malin$YEAR)

# unique years with VCT data 
years = unname(sapply(f1, function(i) as.numeric(substr(strsplit(basename(i), '[-|_]')[[1]][3], 1, 4))))

#----------------------------------------------------------------------------------------------------------#
# 3. sample each VCT variable for each sample, conditioned by the observation year
#----------------------------------------------------------------------------------------------------------#

v1 = v2 = v3 = vector('numeric', nrow(xy))

for (y in uy) {
    
    si = which(bt_malin$YEAR == y)
    ri = which(years == y)
    
    if (length(ri) > 0) {
        v1[si] = terra::extract(r1[[ri]], xy[si,])
        v2[si] = terra::extract(r2[[ri]], xy[si,])
        v3[si] = terra::extract(r3[[ri]], xy[si,])
    } else {
        v1[si] = NA
        v2[si] = NA
        v3[si] = NA
    }
    
}

bt_malin$'tree cover %' = v1
bt_malin$'non-tree veg. %' = v2
bt_malin$'non-veg. %' = v3

saveRDS(bt_malin, '/data/idiv_meyer/01_projects/Ruben/Malin/bt_malin_vct.rds')

#----------------------------------------------------------------------------------------------------------#
# 4. derive mean/sd of VCT variables for each unique sampling site
#----------------------------------------------------------------------------------------------------------#

uid = unique(bt_malin$rarefyID)

stats = do.call(rbind, lapply(uid, function(i) {
    
    ind = which(bt_malin$rarefyID == i)
    s1 = mean(bt_malin$'tree cover %'[ind], na.rm=T)
    s2 = sd(bt_malin$'tree cover %'[ind], na.rm=T)
    s3 = mean(bt_malin$'non-tree veg. %'[ind], na.rm=T)
    s4 = sd(bt_malin$'non-tree veg. %'[ind], na.rm=T)
    s5 = mean(bt_malin$'non-veg. %'[ind], na.rm=T)
    s6 = sd(bt_malin$'non-veg. %'[ind], na.rm=T)
    
    return(data.frame(id=i, s1=s1, s2=s2, s3=s3, s4=s4, s5=s5, s6=s6))
    
}))

colnames(stats) = c('rarefyID', 'tree cover % (mean)', 'tree cover % (sd)', 
                    'non-tree veg. % (mean)', 'non-tree veg. % (sd)', 
                    'non-veg. % (mean)', 'non-veg. % (sd)')

v0 = apply(stats[,2:7], 1, max)

stats[which(v0 == 0),2:7] = NA

saveRDS(stats, '/data/idiv_meyer/01_projects/Ruben/Malin/bt_malin_vct_uniqueID_stats.rds')

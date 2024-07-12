# Read in and summarize the downsampled model fits

### Functions -----------
library(data.table) # for handling large datasets
library(here)
source(here('code', 'util.R'))


# Read in Tchange slope files ------------
slopefiles <- list.files(path = 'temp', pattern = glob2rx('slopes_modTchangexYearxRealmJtu_boot*.rds'), full.names=TRUE) # from pred_GLMM_downsamp.R or fit_pred_turnover_GLMM_downsamp.R
n <- length(slopefiles)
n # number of downsample slopes made
n <- min(1000, n) # take first 1000
for(i in 1:n){
    if(i %% 20 == 0) cat(paste0(i,','))
    temp <- readRDS(slopefiles[i])
    if(i==1){ 
        slopesdownsamp <- cbind(data.table(boot =i), readRDS(slopefiles[i])) 
    } 
    else{
        slopesdownsamp <- rbind(slopesdownsamp, cbind(data.table(boot =i), 
                                                    readRDS(slopefiles[i])))
    }
}
slopesdownsamp <- slopesdownsamp[tempave ==15,] # no Tave effect, so trim out the various levels
slopesdownsamp[, tempave := NULL] # remove the column

# Write out
write.csv(slopesdownsamp, file = gzfile('output/downsampTchange.csv.gz'), row.names = FALSE)


# Read in Tave sensitivity files ------------
sensfiles <- list.files(path = 'temp', pattern = glob2rx('sensitivity_modTchangexTavexYearxRealmJtu_boot*.rds'), full.names=TRUE) # from fit_pred_turnover_GLMM_downsamp.R
n <- length(sensfiles)
n # number of downsample sensitivities made
n <- min(1000, n) # take first 1000
for(i in 1:n){
    if(i %% 20 == 0) cat(paste0(i,','))
    temp <- readRDS(sensfiles[i])
    if(i==1){ 
        temp <- readRDS(sensfiles[i])
        temp <- temp[tempave %in% c(0,25), ]
        senspred <- cbind(data.table(boot =i, type = 'downsamp'), temp) 
    } 
    else{
        temp <- readRDS(sensfiles[i])
        temp <- temp[tempave %in% c(0,25), ]
        senspred <- rbind(senspred, cbind(data.table(boot =i, type = 'downsamp'), temp))
    }
}

# Write out
write.csv(senspred, file = gzfile('output/downsampTave.csv.gz'), row.names = FALSE)

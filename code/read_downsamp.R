# Read in and summarize the downsampled model fits

### Functions -----------
library(data.table) # for handling large datasets
library(here)
source(here('code', 'util.R'))


# Read in files ------------
slopefiles <- list.files(path = 'temp', pattern = glob2rx('slopes_modOBsdTMERtsRealmtsigninitAllJtu_boot*.rds'), full.names=TRUE) # from pred_GLMM_downsamp.R or fit_pred_turnover_GLMM_downsamp.R
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
slopesdownsamp[, tempave := NULL]

# Write out ------------
write.csv(slopesdownsamp, file = gzfile('output/downsampTchange.csv.gz'), row.names = FALSE)

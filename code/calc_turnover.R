# Calculate biodiversity turnover as slope of dissimilarity vs. time

## Set up --------------------
if(Sys.info()['nodename'] == 'annotate.sebs.rutgers.edu'){
  library('mgcv', lib.loc = '/usr/lib64/R/library') # when running on Annotate. Need to load 1.8-26, not 1.8-33.
} else {
  library('mgcv')
}
library(data.table)
library(ggplot2)
library(gridExtra) # for grid.arrange
library(here)

# function to calc linear trend from all year pairs of the time-series
calctrend <- function(y, year1, year2, measure = 'y', duration_group = NA_character_){
  yrs <- sort(unique(c(year1, year2)))
  dy <- year2 - year1
  if(length(dy)>1){
    mod <- lm(y ~ dy) # fit line
    se <- suppressWarnings(sqrt(diag(vcov(mod)))[2])
    out <- list(disstrend = coef(mod)[2], # coef for the slope
                trendse = se,
                year1 = min(year1), 
                year2 = max(year2),
                measure = measure,
                duration_group = duration_group,
                nsamps = length(yrs))
    return(out)
  }
  if(length(dy)==1){
    out <- list(disstrend = y/dy, # coef for the slope
                trendse = NA_real_,
                year1 = min(year1), 
                year2 = max(year2),
                measure = measure,
                duration_group = duration_group,
                nsamps = length(yrs))
    return(out)
  }
  if(length(dy)<1){
    out <- list(disstrend = NA_real_, trendse = NA_real_, year1 = NA_real_, year2 = NA_real_,
                measure = measure, duration_group = duration_group, nsamps = NA_integer_)
    return(out)
  }
}



# Load data
# biotime community dissimilarity data
load(here::here('data', 'biotime_blowes', 'all_pairs_beta.Rdata')) # load rarefied_beta_medians, which has all pairwise dissimilarities
bt <- data.table(rarefied_beta_medians); rm(rarefied_beta_medians)
bt[, year1 := as.numeric(year1)] # not sure why it gets read in as character
bt[, year2 := as.numeric(year2)]
bt[, Horn := 1-Hornsim] # convert similarity to dissimilarity
bt[, Hornsim := NULL]

# BioTime data type
load(here('data', 'biotime_blowes', 'time_series_data_type.Rdata')) # loads rarefyID_type
bt <- merge(bt, rarefyID_type, by = 'rarefyID', all.x = TRUE) # merge with biotime

# biotime taxa category and other info
load(here::here('data', 'biotime_blowes', 'bt_malin.Rdata')) # load bt_malin
btinfo <- data.table(bt_malin); rm(bt_malin)
btinfo[, Nave := mean(N), by = rarefyID] # average number of individuals in the timeseries
btinfo2 <- btinfo[!duplicated(rarefyID), .(rarefyID, rarefyID_x, rarefyID_y, Nave, Biome, taxa_mod, REALM, STUDY_ID)]
bt <- merge(bt, btinfo2, by = 'rarefyID') # trims out rarefyIDs not in bt_malin

# richness
rich <- fread(here('output','richness_by_rarefyID.csv.gz')) # number of species



## QA/QC ---------------------------
# Keep studies with:
# - >=5 species
# - >=2 years
# - >=10 individuals on average.
length(unique(bt$rarefyID))
length(setdiff(bt$rarefyID, rich$rarefyID))


# number of years per study (in remaining samples)
nyrs <- bt[, .(nyrBT = length(unique(c(year1, year2)))), by = rarefyID]
length(setdiff(bt$rarefyID, nyrs$rarefyID))
nyrs[nyrBT < 2, .N]

# >=5 species, >=2 yrs, >=10 individuals
bt <- bt[rarefyID %in% rich[Nspp >= 5, rarefyID], ]
length(unique(bt$rarefyID))

bt <- bt[rarefyID %in% nyrs[nyrBT >= 2, rarefyID], ]
length(unique(bt$rarefyID))

bt <- bt[Nave >= 10 | is.na(Nave), ]
length(unique(bt$rarefyID))


setkey(bt, STUDY_ID, rarefyID, year1,  year2)



## Calcs with all pairs --------------
print('Calculating for output/slope.csv.gz')

# All years available using all year pairs
trends <- bt[, calctrend(Jtu, year1, year2, measure = 'Jtu', duration_group = 'All'), by = .(rarefyID)][!is.na(disstrend), ]
temp <- bt[, calctrend(Jbeta, year1, year2, measure = 'Jbeta', duration_group = 'All'), by = .(rarefyID)]
trends = rbind(trends, temp[!is.na(disstrend), ])
temp <- bt[!is.na(Horn), calctrend(Horn, year1, year2, measure = 'Horn', duration_group = 'All'), by = .(rarefyID)]
trends = rbind(trends, temp[!is.na(disstrend), ])

# Write out all pairs -----------------
write.csv(trends, file = gzfile(here('output', 'slope.csv.gz')), row.names = FALSE)


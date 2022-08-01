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


# function to calc linear trend from all years of the time-series
# only comparing back to the first year
calctrendy1 <- function(y, year1, year2, measure = 'y', duration_group = NA_character_){
  yrs <- sort(unique(c(year1, year2)))
  startyear = min(year1)
  dy <- year2 - year1
  keep = year1 == startyear
  if(sum(keep)>1){
    y2 <- y[keep]
    dy2 <- dy[keep]
    mod <- lm(y2 ~ dy2) # fit line
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
  if(sum(keep)==1){
    y2 <- y[keep]
    dy2 <- dy[keep]
    out <- list(disstrend = y2/dy2, # coef for the slope
                trendse = NA_real_,
                year1 = min(year1), 
                year2 = max(year2),
                measure = measure, 
                duration_group = duration_group,
                nsamps = length(yrs))
    return(out)
  } 
  if(sum(keep)<1){
    out <- list(disstrend = NA_real_, trendse = NA_real_, year1 = NA_real_, year2 = NA_real_,
                measure = measure, duration_group = duration_group, nsamps = NA_integer_)
    return(out)
  }
}

# function to calc linear trend from all pairs
# only using a sequence of annual samples that is numyrs long
# use the most recent sequence that fits this criterion
calctrendlast <- function(y, year1, year2, numyrs, measure = 'y', duration_group = NA_character_){
  # try to identify a suitable sequence of samples
  yrs <- sort(unique(c(year1, year2)))
  dy <- diff(yrs) # find intervals between years
  rl <- rle(dy) # run length encoding
  end = cumsum(rl$lengths) # find ends of runs
  start = c(1, head(end, -1) + 1) # find starts of runs
  rlkeep <- which(rl$lengths >= numyrs & rl$values == 1) # find runs with desired length and 1 year intervals
  
  # find the latest run that meets our criteria (if any exist)
  if(length(rlkeep)>0){
    rlkeep <- max(rlkeep) 
    start <- start[rlkeep] # start of this run
    end <- end[rlkeep] # end of this run
    dykeep <- rep(FALSE, length(dy)) # year intervals to keep
    dykeep[start:end] <- TRUE
    yrs2 <- yrs[c(dykeep, FALSE) | c(FALSE, dykeep)] # keep years involved in at least one interval to keep
    maxyr <- max(yrs2)
    yrs3 <- yrs2[yrs2 > maxyr - numyrs] # only last numyrs years
    ykeep <- year1 %in% yrs3 & year2 %in% yrs3 # keep values in pairwise comparisons for the years we want
    y <- y[ykeep] # trim the timeseries
    year1 <- year1[ykeep]
    year2 <- year2[ykeep]
    run <- TRUE
  } else{
    run <- FALSE
  }
  
  # calculate slope (if a suitable sequence of samples was found)
  if(run){
    dy <- year2 - year1
    yrs <- sort(unique(c(year1, year2)))
    
    mod <- invisible(lm(y ~ dy)) # fit line
    se <- suppressWarnings(sqrt(diag(vcov(mod)))[2])
    out <- list(disstrend = coef(mod)[2], # coef for the slope
                trendse = se,
                year1 = min(yrs3), 
                year2 = max(yrs3),
                measure = measure, 
                duration_group = duration_group,
                nsamps = length(yrs))
    return(out)
    
  } else {
    out <- list(disstrend = NA_real_, trendse = NA_real_, year1 = NA_real_, year2 = NA_real_,
                measure = measure, duration_group = duration_group, nsamps = NA_integer_)
    return(out)
  }
}

# function to calc linear trend from all pairs
# only using a sequence of samples that is numyrs samples long and has nsamps samples
# not necessarily annual samples
# use the most recent sequence that fits this criterion
calctrendnsamps <- function(y, year1, year2, numyrs, nsamps, 
                            measure = 'y', duration_group = NA_character_){
  if(nsamps > numyrs) stop('nsamps must be <= numyrs')
  if(length(y) != length(year1) | length(y) != length(year2)) stop('y, year1, and year2 must be the same length')
  yrs <- sort(unique(c(year1, year2)))
  
  # brute force search for a sequence of samples that match criteria
  i = length(yrs) # index for finding a suitable sequence of samples. start at the end.
  run <- FALSE # flag for whether we have found a suitable sequence of samples
  while(i > 0 & !run){ 
    dys <- yrs[i] - yrs
    j <- which(dys == numyrs-1)
    if(length(j) >0){
      proposedset <- yrs[j:i]
      if(length(proposedset) >= nsamps){
        ykeep <- year1 %in% proposedset & year2 %in% proposedset # keep values in pairwise comparisons for the years we want
        y <- y[ykeep] # trim the timeseries
        year1 <- year1[ykeep] # trim the initial years
        year2 <- year2[ykeep] # trim the final years
        run <- TRUE # mark that we should run the calcs
      }
    }
    if(run == FALSE){
      i <- i - 1 # try the next earliest sample
    }
  }
  
  if(run){
    dy <- year2 - year1
    yrs <- sort(unique(c(year1, year2)))
    
    mod <- lm(y ~ dy) # fit line
    se <- suppressWarnings(sqrt(diag(vcov(mod)))[2]) # standard error
    out <- list(disstrend = coef(mod)[2], # coef for the slope
                trendse = se, 
                year1 = min(proposedset), 
                year2 = max(proposedset),
                measure = measure, 
                duration_group = duration_group,
                nsamps = length(yrs)) # SE
    return(out)
    
  } else {
    out <- list(disstrend = NA_real_, trendse = NA_real_, year1 = NA_real_, year2 = NA_real_,
                measure = measure, duration_group = duration_group, nsamps = NA_integer_)
    return(out)
  }
}


# function to calc linear trend from all pairs
# only using a sequence of samples that is numyrs samples long and has nsamps samples
# not necessarily annual samples
# use all sequences that fits this criterion
calctrendnsampsall <- function(y, year1, year2, numyrs, nsamps, 
                               measure = 'y', duration_group = NA_character_){
  if(exists('out')) rm(out) # remove the output object, if it exists for some reason
  
  if(nsamps > numyrs) stop('nsamps must be <= numyrs')
  if(length(y) != length(year1) | length(y) != length(year2)) stop('y, year1, and year2 must be the same length')
  
  yrs <- sort(unique(c(year1, year2))) # list of years in the dataset
  
  # search for sequences of samples that match criteria
  run <- FALSE # flag for whether we have found a suitable sequence of samples
  for(i in length(yrs):1){ 
    dys <- yrs[i] - yrs # find difference from this year to all other years in the set
    j <- which(dys == numyrs-1) # indiex for the other year that matches our desired length
    if(length(j) >0){
      proposedset <- yrs[j:i] # set of years that match our duration criterion
      if(length(proposedset) >= nsamps){ # check if number of smaples is sufficient
        ykeep <- year1 %in% proposedset & year2 %in% proposedset # keep values in pairwise comparisons for the years we want
        thisy <- y[ykeep] # trim the timeseries for this calc
        thisyear1 <- year1[ykeep] # trim the initial years for this calc
        thisyear2 <- year2[ykeep] # trim the final years for this calc
        thisdy <- thisyear2 - thisyear1 # calculate temporal difference among the year pairs for this calc
        run <- TRUE # mark that we ran the calcs at least once and should return an answer
        
        mod <- lm(thisy ~ thisdy) # fit line
        se <- suppressWarnings(sqrt(diag(vcov(mod)))[2]) # standard error of the slope
        thisout <- data.table(disstrend = coef(mod)[2], # coef for the slope
                              trendse = se, 
                              year1 = min(proposedset), 
                              year2 = max(proposedset),
                              measure = measure, 
                              duration_group = duration_group,
                              nsamps = length(proposedset)) # SE
        if(exists('out')) out <- rbind(out, thisout) # append if the output object exists
        if(!exists('out')) out <- thisout # create the output object if it doesn't exist
      }
    }
  }
  
  if(run){
    return(out)
    
  } else {
    out <- data.table(disstrend = NA_real_, trendse = NA_real_, year1 = NA_real_, year2 = NA_real_,
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



## Calcs with first year or last subset --------------
# Trends of standardized length and annual sampling frequency
# Only use the last sequence of samples that meet criteria
# Add this to slopes using all pairs
print('Calculating for temp/trendstemp.rds')
yrslist <- 3:20
for(yr in yrslist){
  print(yr)
  temp <- bt[, calctrendlast(Jtu, year1, year2, numyrs = yr, measure = 'Jtu', 
                             duration_group = paste0(yr, 'annual')), by = .(rarefyID)]
  trends = rbind(trends,temp[!is.na(disstrend), ]) # append
  
  temp <- bt[, calctrendlast(Jbeta, year1, year2, numyrs = yr, measure = 'Jbeta',
                             duration_group = paste0(yr, 'annual')), by = .(rarefyID)]
  trends = rbind(trends, temp[!is.na(disstrend), ])
  temp <- bt[!is.na(Horn), calctrendlast(Horn, year1, year2, numyrs = yr, measure = 'Horn',
                                         duration_group = paste0(yr, 'annual')), by = .(rarefyID)]
  trends = rbind(trends, temp[!is.na(disstrend), ])
}

# Trends of standardized length and at least 3 years sampled
# Use last possible sequences of samples
yrslist <- 3:20
for(yr in yrslist){
  print(yr)
  temp <- bt[, calctrendnsamps(Jtu, year1, year2, numyrs = yr, nsamps = 3, measure = 'Jtu', 
                               duration_group = paste0(yr, 'min3')), by = .(rarefyID)]
  trends = rbind(trends,temp[!is.na(disstrend), ])
  
  temp <- bt[, calctrendnsamps(Jbeta, year1, year2, numyrs = yr, nsamps = 3, measure = 'Jbeta',
                               duration_group = paste0(yr, 'min3')), by = .(rarefyID)]
  trends = rbind(trends, temp[!is.na(disstrend), ])
  temp <- bt[!is.na(Horn), calctrendnsamps(Horn, year1, year2, numyrs = yr, nsamps = 3, measure = 'Horn',
                                           duration_group = paste0(yr, 'min3')), by = .(rarefyID)]
  trends = rbind(trends, temp[!is.na(disstrend), ])
}


# All years available using first year comparisons
temp <- bt[, calctrendy1(Jtu, year1, year2, measure = 'Jtu', duration_group = 'Ally1'), by = .(rarefyID)]
trends = rbind(trends, temp[!is.na(disstrend), ])
temp <- bt[, calctrendy1(Jbeta, year1, year2, measure = 'Jbeta', duration_group = 'Ally1'), by = .(rarefyID)]
trends = rbind(trends, temp[!is.na(disstrend), ])
temp <- bt[!is.na(Horn), calctrendy1(Horn, year1, year2, measure = 'Horn', duration_group = 'Ally1'), by = .(rarefyID)]
trends = rbind(trends, temp[!is.na(disstrend), ])

# Write out first year and last subset -----------------
saveRDS(trends, file = here('temp', 'trendstemp.rds'))



## Calcs  with all subsets ---------------
# Uses all sub
print('Calculating for temp/trendstempallmin3.rds')

yrslist <- 3:20
for(yr in yrslist){
  print(yr)
  temp <- bt[, calctrendnsampsall(Jtu, year1, year2, numyrs = yr, nsamps = 3, measure = 'Jtu', 
                                  duration_group = paste0(yr, 'min3')), by = .(rarefyID)]
  if(yr == min(yrslist)) trendsallmin3 = temp[!is.na(disstrend), ] # make a new dataset if first iteration through
  if(yr > min(yrslist)) trendsallmin3 = rbind(trendsallmin3,temp[!is.na(disstrend), ]) # otherwise append
  
  temp <- bt[, calctrendnsampsall(Jbeta, year1, year2, numyrs = yr, nsamps = 3, measure = 'Jbeta',
                                  duration_group = paste0(yr, 'min3')), by = .(rarefyID)]
  trendsallmin3 = rbind(trendsallmin3, temp[!is.na(disstrend), ])
  temp <- bt[!is.na(Horn), calctrendnsampsall(Horn, year1, year2, numyrs = yr, nsamps = 3, measure = 'Horn',
                                              duration_group = paste0(yr, 'min3')), by = .(rarefyID)]
  trendsallmin3 = rbind(trendsallmin3, temp[!is.na(disstrend), ])
}

# Write out with all subsets -----------------
saveRDS(trendsallmin3, file = here('temp', 'trendstempallmin3.rds'))




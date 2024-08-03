# utility functions

# require(Hmisc) # for wtd.var
require(scales) # for trans_new


# Modified from http://monkeysuncle.stanford.edu/?p=485
# Upper and lower are se (or sd)
error.bar <- function(x, y, upper, lower=upper, length=0.1, dir='y', ...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
    if(dir=='y') arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
    if(dir=='x') arrows(x-lower,y, x+upper, angle=90, code=3, length=length, ...)
}

# for vector y in which each entry is associated with start and end years (year1 and year2)
# pick the value of y that has the longest duration
picklongest <- function(year1, year2, y){
    dy <- year2 - year1
    keep <- which.max(dy)
    return(y[keep])
}

# for vector y in which each entry is associated with start and end years (year1 and year2)
# pick the value of y that has the earliest start and end
pickfirst <- function(year1, year2, y){
    keep1 <- which.min(year1)
    keep2 <- which.min(year2)
    keep <- intersect(keep1, keep2)
    if(length(keep)==1){
        return(y[keep])
    }
}

# sign of temperature change
signneg11 <- function(x){ # assign 0 a sign of 1 so that there are only 2 levels
    out <- sign(x)
    out[out == 0] <- 1
    return(out)
}


# transformation for 2 categories. Eq. 1 in Douma & Weedon 2019 MEE
transform01 <- function(x) (x * (length(x) - 1) + 0.5) / (length(x))

# function for returning the slope of a regression and catching errors related to NAs
lmNAcoef <- function(y, x){
    if(sum(!is.na(y)) > 1){ # if enough data points to calculate a slope
        b <- coef(lm(y ~ x))[2]
        return(b)
    } else {
        return(NA_real_)
    }
}

# plot pearson's correlation on the lower triangle of a pairs plot
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = 'pairwise.complete.obs')
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt) #, cex = cex.cor * r)
}

# function to scale a variable using the definitions in the scalingall object (assumes this is loaded)
scaleme <- function(x, nm){
    if(!exists('scalingall')) stop('Have to load scalingall')
    if(!(nm %in% scalingall[,var])) stop('nm not found in scalingall')
    if(scalingall[var==nm, log]){
        x.sc <- (log(x + scalingall[var==nm, plus]) - scalingall[var == nm, center]) / scalingall[var == nm, scale]  
    } else {
        x.sc <- (x  + scalingall[var==nm, plus] - scalingall[var == nm, center]) / scalingall[var == nm, scale]
    }
    return(x.sc)
}
unscaleme <- function(x.sc, nm){
    if(!exists('scalingall')) stop('Have to load scalingall')
    if(!(nm %in% scalingall[,var])) stop('nm not found in scalingall')
    if(scalingall[var==nm, log]){
        x <- exp(x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center]) - scalingall[var==nm, plus]
    } else {
        x <- x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center] - scalingall[var==nm, plus]
    }
    return(x)
}


# given a duration, make a Gaussian white noise timeseries and return the slope
calcslopeGauss <- function(dur){
    x <- 1:dur
    y <- rnorm(dur)
    return(coef(lm(y~x))[2])
}



# useful transformations for ggplot axes
signedsqrt = function(x) sign(x)*sqrt(abs(x))
signedsq = function(x) sign(x) * x^2
signedsqrttrans <- trans_new(name = 'signedsqrt', transform = signedsqrt, inverse = signedsq)


# from https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
    newplot <- myPlot +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
               color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = textSize), 
              legend.text  = element_text(size = textSize),
              legend.key.size = unit(spaceLegend, "lines"))
    return(newplot)
}

# binomial ci (95% by default)
binomci <- function(x, n){
    out <- as.numeric(binom.test(x, n)$conf.int)
    return(list(out[1], out[2]))
}

# downsample a data.table so that there are t rows for each rarefyID, where t is the number of unique year1 and year2 values
# (of the t*(t-1)/2 rows that currently exist for each rarefyID)
# inds is a vector of TRUE/FALSE that indicates which rows are being kept before calling this function
# seed is a random number seed, for reproducibility
# returns a vector of TRUE/FALSE to indicate which rows to keep after downsampling
downsampRarefyID <- function(dat, inds = rep(TRUE, nrow(dat)), seed=1){
    set.seed(seed)
    dat[, inds := inds]
    tslen <- dat[inds == TRUE, .(t = length(unique(c(year1, year2)))), by = rarefyID] # length of each rarefyID
    dat <- merge(dat, tslen, all.x = TRUE) # merge back. NAs for the rarefyIDs not included by inds
    dat[inds == TRUE , inds2 := sample(c(rep(TRUE, min(.N, unique(t))), rep(FALSE, max(0,.N-unique(t))))), by = rarefyID] # sample t rows (and make sure t <= existing number of rows)
    dat[is.na(inds2), inds2 := FALSE]
    return(dat$inds2)
}


# slopes and SEs from resampling
# n: number of resamples, other columns are the x, y, and se of y variables
# colnames refer to the output
slopesamp <- function(n, duration, Jtu, Jtu.se, colnames = c('slope', 'slope.se')){
    if(length(duration) != length(Jtu)) stop('duration and Jtu are not the same length')
    if(length(duration) != length(Jtu.se)) stop('duration and Jtu.se are not the same length')
    if(length(Jtu) != length(Jtu.se)) stop('Jtu and Jtu.se are not the same length')
    
    samp <- rep(NA, n) # will hold the slopes of the sampled data
    for(j in 1:n){
        y <- rnorm(length(Jtu), mean = Jtu, sd = Jtu.se) # one sample
        samp[j] <- coef(lm(y ~ duration))[2] # fit line, get slope
    }
    out <- c(coef(lm(Jtu ~ duration))[2], sd(samp))
    names(out) <- colnames
    return(as.list(out)) # coercing to list will allow the data.table aggregate used later to create 2 columns
}

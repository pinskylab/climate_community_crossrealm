# utility functions

require(Hmisc) # for wtd.var

# geometric mean. remove NAs and x<0
geomean = function(x){ 
    exp(sum(log(x[x > 0 & !is.na(x)])) / length(x[x > 0 & !is.na(x)]))
}

# geometric standard deviation
geosd <- function(x){ 
    exp(sd(log(x[x > 0 & !is.na(x)])))
}

# weighted mean, but returns regular mean if any w's are NA or if sum(w)==0
meanwt <- function(x, w){ 
    w <- w[!is.na(x)] # remove NAs
    x <- x[!is.na(x)]
    if(any(is.na(w))){
        return(mean(x))  
    } else {
        if(sum(w) == 0){
            return(mean(x))
        } else {
            return(weighted.mean(x, w))
        }
    }
}

# weighted sd
sdwt <- function(x, w){ 
    w <- w[!is.na(x)] # remove NAs
    x <- x[!is.na(x)]
    if(any(is.na(w))){
        return(sd(x))
    } else {
        if(sum(w) == 0){
            return(sd(x))
        } else {
            return(sqrt(wtd.var(x = x, weights = w, normwt = TRUE)))        
        }
    }
}

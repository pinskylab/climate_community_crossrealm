# utility functions

# require(Hmisc) # for wtd.var
require(scales) # for trans_new

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

# weighted sd (need Hmisc)
# sdwt <- function(x, w){ 
#     w <- w[!is.na(x)] # remove NAs
#     x <- x[!is.na(x)]
#     if(any(is.na(w))){
#         return(sd(x))
#     } else {
#         if(sum(w) == 0){
#             return(sd(x))
#         } else {
#             return(sqrt(wtd.var(x = x, weights = w, normwt = TRUE)))        
#         }
#     }
# }


# generate ellipse coefficients for the CI on the linear regression coefs, given a dataset with points on a line and standard error for y on the points
# assumes dat has x in col 1, y in col 2, and SE of y in col 3
# use the mean SE for now
# alpha is the alpha level for the CI
# adapted from https://stats.stackexchange.com/questions/70629/calculate-uncertainty-of-linear-regression-slope-based-on-data-uncertainty
# return coefs are for an allipse where x is the slope and y is the intercept
# developed 6 Jan 2022
ellipse.coefs <- function(dat, alpha = 0.05){
    a <- sum(dat[,1]^2)
    b <- nrow(dat)
    c <- 2*sum(dat[,1])
    d <- -2*sum(dat[,1]*dat[,2])
    e <- -2*sum(dat[,2])
    sigma <- mean(dat[,3])
    f <- sum(dat[,2]^2) - qchisq(alpha, nrow(dat), lower.tail = FALSE) * sigma^2
    out <- list(a=a,b=b,c=c,d=d,e=e,f=f)
}

# make ellipse points from coefficients
# https://stackoverflow.com/questions/41820683/how-to-plot-ellipse-given-a-general-equation-in-r
# a*x^2 + b*y^2 + c*x*y + d*x + e*y + f = 0
# in our application, x is the linear regression slope and y is the intercept
# developed 6 Jan 2022 to go with ellipse.coefs
get.ellipse <- function (a, b, c, d, e, f, n.points = 1000) {
    ## solve for centre
    A <- matrix(c(a, c / 2, c / 2, b), 2L)
    B <- c(-d / 2, -e / 2)
    mu <- solve(A, B)
    ## generate points on circle
    r <- sqrt(a * mu[1] ^ 2 + b * mu[2] ^ 2 + c * mu[1] * mu[2] - f)
    theta <- seq(0, 2 * pi, length = n.points)
    v <- rbind(r * cos(theta), r * sin(theta))
    ## transform for points on ellipse
    z <- backsolve(chol(A), v) + mu
    return(list(x = z[1,], y = z[2,]))
}


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

# produce ggplot-style colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
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

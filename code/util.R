# utility functions

# require(Hmisc) # for wtd.var

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

#Reproducible example of trying to calculate slope SE from a beta regression



########################
# Functions

library(glmmTMB)
library(ggplot2)
library(data.table)

# generate ellipse coefficients for the CI on the linear regression coefs, given a dataset
# assumes dat has x in col 1, y in col 2, and SE of y in col 3
# use the mean SE for now
# alpha is the alpha level for the CI
# adapted from https://stats.stackexchange.com/questions/70629/calculate-uncertainty-of-linear-regression-slope-based-on-data-uncertainty
# return coefs are for an allipse where x is the slope and y is the intercept
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

##################################################################################
## Example of trying to calculate a slope while considering SE in the response
# (not beta)
dat <- data.frame(x = 0:8, y = seq(0,16, length.out=9), y.se = 3)

mod <- lm(y ~ x, dat) # naive model, not considering error in y
summary(mod)
preds <- predict(mod, se.fit = TRUE)

plot(dat$x, dat$y, ylim=c(-7,24))
arrows(dat$x, dat$y-1.96*dat$y.se, dat$x, dat$y+1.96*dat$y.se, length=0)
polygon(c(dat$x, rev(dat$x)), c(preds$fit+preds$se.fit, rev(preds$fit-preds$se.fit)))

coefs <- ellipse.coefs(dat, alpha = 0.05)
xy <- get.ellipse(coefs$a, coefs$b, coefs$c, coefs$d, coefs$e, coefs$f)
plot(xy, xlab='slope', ylab='intercept') #
range(xy$x) # the 95% CI for the slope

####################################################
## An example with some error

# the data
set.seed(5)
dat <- data.frame(x = 0:8, y = seq(0,16, length.out=9)+rnorm(9, 0, 0.5), y.se = 3)

# fit a naive model, not considering error in y
mod <- lm(y ~ x, dat)
summary(mod)
preds <- predict(mod, se.fit = TRUE)

plot(dat$x, dat$y, ylim=c(-7,22))
arrows(dat$x, dat$y-1.96*dat$y.se, dat$x, dat$y+1.96*dat$y.se, length=0)

# plot the confidence interval on the linear regression
polygon(c(dat$x, rev(dat$x)), c(preds$fit+preds$se.fit, rev(preds$fit-preds$se.fit)), col = 'grey')

# get 95% CI for the slope parameter, assuming points fit a linear regression perfectly
coefs <- ellipse.coefs(dat, alpha = 0.05)
xy <- get.ellipse(coefs$a, coefs$b, coefs$c, coefs$d, coefs$e, coefs$f)
range(xy$x) # the 95% CI for the slope


######################################
## example of a beta regression

# set up the data. y is constrained 0-1. there is one covariate (z)
set.seed(5)
dat <- data.frame(x = rep(1:10, 2),
                  z = rep(c(0,2), c(10,10)))

dat$y_trans <- 0.1*dat$x + 0.1*dat$z*dat$x + rnorm(nrow(dat), 0, 0.3)
dat$y <- 1/(1+exp(-(dat$y_trans)))
plot(dat$x, dat$y, col =factor(dat$z))

# fit a beta regression where a affects the slope
mod <- glmmTMB(y ~ x + z:x, data = dat, family = beta_family(link='logit'))
summary(mod)

# predict from the model with SEs on the points
newdat <- data.table(expand.grid(z = seq(0,2, length.out=10), x = seq(1, 10, length.out=10)))
preds <- predict(mod, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$fit <- preds$fit
newdat$fit.se <- preds$se.fit

ggplot(newdat, aes(x, fit, color = z, group = z)) +
    geom_point() +
    geom_pointrange(aes(ymin=fit-fit.se, ymax=fit+fit.se))

# calculate slopes for each level of a
# including SE from the amount that points don't fit a straight line
slopes <- newdat[, .(slope = coef(lm(fit ~ x))[2], slope.se = sqrt(diag(vcov(lm(fit ~ x))))[2]), by = z]

ggplot(slopes, aes(z, slope)) +
    geom_line() +
    geom_pointrange(aes(ymin=slope-slope.se, ymax=slope+slope.se))

slopes$slope.se # SE is up to  0.0025


# 95% CI of the slope from SE on the individual points
coefs <- newdat[, ellipse.coefs(cbind(x, fit, fit.se), alpha = 0.05), by = z]
xy <- coefs[, get.ellipse(a, b, c, d, e, f), by = z]
xy[, .(ellipseCI = diff(range(x))/1.96/2), by = z] # the SE for the slope, converted from the 95% CI


# compare the ellipse CI and the lm CI
# ellipse CI is much higher when the points fit a line well, but about equivalent otherwise
merge(slopes[, .(z,slope.se)], xy[, .(ellipseCI = diff(range(x))/1.96), by = z])

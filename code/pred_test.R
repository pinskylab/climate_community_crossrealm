# Reproducible example of trying to calculate slope SE from a beta regression
library(glmmTMB)
library(ggplot2)
library(data.table)

# return x and y coords of an elipse given the conic parameters
# https://stackoverflow.com/questions/41820683/how-to-plot-ellipse-given-a-general-equation-in-r
# for Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
ellipse <- function(A, B, C, D, E, F){
    M0 <- matrix(c(F,D/2,E/2, D/2, A, B/2, E/2, B/2, C), nrow=3, byrow=TRUE)
    M <- matrix(c(A,B/2,B/2,C), nrow=2)
    lambda <- eigen(M)$values
    
    # assuming abs(lambda[1] - A) < abs(lambda[1] - C), if not, swap lambda[1] and lambda[2] in the following equations:
    if(abs(lambda[1] - A) > abs(lambda[2] - C) ){
        lambda <- lambda[2:1]
    }
    a <- sqrt(-det(M0)/(det(M)*lambda[1]))  
    b <- sqrt(-det(M0)/(det(M)*lambda[2]))
    xc <- (B*E-2*C*D)/(4*A*C-B^2)
    yc <- (B*D-2*A*E)/(4*A*C-B^2)
    phi <- pi/2 - atan((A-C)/B)*2
    
    t <- seq(0, 2*pi, 0.01) 
    x <- xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)
    y <- yc + a*cos(t)*cos(phi) + b*sin(t)*cos(phi)
    return(list(x=x,y=y))
}

# https://stackoverflow.com/questions/41820683/how-to-plot-ellipse-given-a-general-equation-in-r
# a*x^2 + b*y^2 + c*x*y + d*x + e*y + f = 0
plot.ellipse <- function (a, b, c, d, e, f, n.points = 1000) {
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
    ## plot points
    #return(z)
    plot(t(z), type = "l")
}

dat <- data.frame(x = rep(1:10, 2),
                  a = rep(c(0,2), c(10,10)))

dat$y_trans <- 0.1*dat$x + 0.1*dat$a*dat$x + rnorm(nrow(dat), 0, 0.1)
dat$y <- 1/(1+exp(-(dat$y_trans)))
plot(dat$x, dat$y, col =factor(dat$a))

mod <- glmmTMB(y ~ x + a:x, data = dat, family = beta_family(link='logit'))
summary(mod)


newdat <- data.table(expand.grid(a = seq(0,2, length.out=10), x = seq(1, 10, length.out=10)))
preds <- predict(mod, newdata = newdat, se.fit = TRUE, re.form = NA, allow.new.levels=TRUE, type = 'response')
newdat$fit <- preds$fit
newdat$fit.se <- preds$se.fit

slopes <- newdat[, .(slope = coef(lm(fit ~ x))[2]), by = a]

ggplot(newdat, aes(x, fit, color = a, group = a)) +
    geom_point() +
    geom_pointrange(aes(ymin=fit-fit.se, ymax=fit+fit.se))

ggplot(slopes, aes(a, slope)) +
    geom_line()



## Example of trying to calculate a slope while considering SE in the response
dat <- data.frame(x = 0:8, y = seq(0,16, length.out=9), y.se = 3)

mod <- lm(y ~ x, dat) # naive model, not considering error in y
summary(mod)
preds <- predict(mod, se.fit = TRUE)

plot(dat$x, dat$y, ylim=c(-10,30))
arrows(dat$x, dat$y-dat$y.se, dat$x, dat$y+dat$y.se, length=0)
polygon(c(dat$x, rev(dat$x)), c(preds$fit+preds$se.fit, rev(preds$fit-preds$se.fit)))

# https://stats.stackexchange.com/questions/70629/calculate-uncertainty-of-linear-regression-slope-based-on-data-uncertainty
chsq <- qchisq(0.05, nrow(dat), lower.tail = FALSE)
sigma <- mean(dat$y.se)


# plot equations from Ontiveros et al. 2021 Ecology
require(glmmTMB) # for fitting beta models

logit <- function(x) return(log(x/(1-x)))

# colonization and extinction rates
c = 2.04e-4; e = 4.1e-4; t = 1:5000 # island birds
c = 5.64e-3; e = 1.14e-2; t= 1:300 # arthropods
c = 7.47e-4; e = 2.84e-5; t= 1:3000 # forest birds

p = c/(e+c) # equilibrium proportion of species present
J = p*(c+e*exp(-(e+c)*t))/(p*(e+c*exp(-(e+c)*t)) + c*(1-exp(-(e+c)*t))) # Jaccard similarity
J2 = (c^2+e*c*exp(-(e+c)*t))/(2*e*c+c^2 - e*c*exp(-(e+c)*t)) # simplified equation for p at equilibrium
Jasym = p*c/(p*e+c) # asymptotic J
min(J)

plot(t,J); lines(t, J2, col='red')
plot(t,(1-J)/(1-Jasym)) # transform to 0-1 with Jasym
plot(t,(1-J)/(1-min(J))) # transform to 0-1 with observed min

plot(t, logit(1-J))
plot(t, logit((1-J)/(1-Jasym))) # logit of 0-1 transform with Jasym
plot(t, logit((1-J)/(1-min(J)))) # logit of 0-1 transform with min J
adj = 100; plot(t, logit((1-J+adj)/(1-Jasym+adj))) # logit of 0-1 Jasym transform with a baseline adjustment
adj = 100; plot(t, logit((1-J+adj)/(1-min(J)+adj))) # logit of 0-1 minJ transform with a baseline adjustment

adj = 100; plot(t, logit((1-J+adj)/(1+adj))) # logit of baseline adjustment (no 0-1 transform)

# other transformations
plot(t, exp(1-J))
plot(t, exp(J))
plot(t, (1-J)^10)
plot(t, -log(J))
plot(t, log(J))
plot(t, log(1/J))
plot(log(t),J) # log(t)
plot(log(t),1-J) # log(t)
plot(log(t),logit(1-J)) # log(t) and logit(J)

# fit a glm and related models
# fit to similarity
mod <-glmmTMB(J ~ t, family=binomial(link='log'), start = list(beta=c(-0.4, -1.6e-4)))
mod2 <-glmmTMB(J ~ t, family=ordbeta(link='logit'))
mod3 <-glmmTMB(J ~ t, family=ordbeta(link='probit'))
mod4 <-glmmTMB(J ~ t, family=ordbeta(link='cloglog'))
Jpred <- predict(mod, type = 'response')
Jpred2 <- predict(mod2, type = 'response')
Jpred3 <- predict(mod3, type = 'response')
Jpred4 <- predict(mod4, type = 'response')
plot(t, J, ylim = c(0,1)); lines(t,Jpred, col='blue'); lines(t, Jpred2, col='red'); lines(t, Jpred3, col='green'); lines(t, Jpred4, col='orange')
plot(t, log(J)); lines(t,log(Jpred), col='blue'); lines(t, log(Jpred2), col='red'); lines(t, log(Jpred3), col='green'); lines(t, log(Jpred4), col='orange') #log scale. linearizes the glm

# fit to dissimilarity
mod <-glmmTMB(I(1-J) ~ t, family=binomial(link='log'), start = list(beta=c(-0.4, -1.6e-4)))
mod2 <-glmmTMB(I(1-J) ~ t, family=ordbeta(link='logit'))
mod3 <-glmmTMB(I(1-J) ~ t, family=ordbeta(link='probit'))
mod4 <-glmmTMB(I(1-J) ~ t, family=ordbeta(link='cloglog'))
Jpred <- predict(mod, type = 'response')
Jpred2 <- predict(mod2, type = 'response')
Jpred3 <- predict(mod3, type = 'response')
Jpred4 <- predict(mod4, type = 'response')
plot(t,1-J, ylim = c(0,1)); lines(t,Jpred, col='blue'); lines(t, Jpred2, col='red'); lines(t, Jpred3, col='green'); lines(t, Jpred4, col='orange')

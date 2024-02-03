###############################
####   PERMUTATION TESTS   ####
####         ANOVA         ####
###############################

# One-way ANOVA
# (p=1, g=6)

head(chickwts)
attach(chickwts)
summary(chickwts)

g <- nlevels(feed)
n <- dim(chickwts)[1]


layout(cbind(1,2))
plot(feed, weight, xlab='treat',col=rainbow(g),main='Original Data')

# H0: tau1 = tau2 = tau3 = tau4 = tau5 = tau6 = 0
# the chickens belong to the same population

# H1: (H0)^c
# the chickens belong to several different population
# Parametric test:
fit <- aov(weight ~ feed)
summary(fit)

# Permutation test:
# Test statistic: F stat
T0 <- summary(fit)[[1]][1,4]
T0

# what happens if we permute the data?
permutazione <- sample(1:n)
weight_perm <- weight[permutazione]
fit_perm <- aov(weight_perm ~ feed)
summary(fit_perm)

plot(feed, weight_perm, xlab='treat',col=rainbow(g),main='Permuted Data')


# CMC to estimate the p-value
B <- 1000 # Number of permutations
T_stat <- numeric(B) 
n <- dim(chickwts)[1]

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  weight_perm <- weight[permutation]
  fit_perm <- aov(weight_perm ~ feed)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

layout(1)
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat))
abline(v=T0,col=3,lwd=4)

# p-value
p_val <- sum(T_stat>=T0)/B
p_val
# we reject the null hypothesis

detach(chickwts)
###################################
# MANOVA

data(iris)
attach(iris)
species.name <- factor(Species, labels=c('setosa','versicolor','virginica'))
iris4        <- iris[,1:4]
plot(iris4,col=species.name)

i1 <- which(species.name=='setosa')
i2 <- which(species.name=='versicolor')
i3 <- which(species.name=='virginica')
n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)
n  <- n1+n2+n3

g  <- length(levels(species.name))
p  <- 4
detach(iris)

# MANOVA
# parametric test:
fit <- manova(as.matrix(iris4) ~ species.name)
summary.manova(fit,test="Wilks")

# How to perform a permutation test in this case?
# Multivariate framework -> We permute the labels associated to each unit
permutation <- sample(1:n)
species.name.perm <- species.name[permutation]
plot(iris4,col=species.name.perm)

fit.perm <- manova(as.matrix(iris4) ~ species.name.perm)
summary.manova(fit.perm,test="Wilks")

# TEST
# Test statistics: Wilks Lambda
T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]
T0
# Note that wilk's lambda is significant for small values!
# It is sufficient to change its sign to use it in a permutation test

# Permutations
B <- 1000
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n)
  species.name.perm <- species.name[permutation]
  fit.perm <- manova(as.matrix(iris4) ~ species.name.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}

layout(1)
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat))
abline(v=T0,col=3,lwd=4)

# p-value
p_val <- sum(T_stat>=T0)/B
p_val


###################################################################################?
# Two-way ANOVA
# Variable:  Distance     [km/l]
# factor1:   Fuel station (0=Esso, 1=Shell)
# factor2:   Type of fuel (0=95,   1=98)

km          <- c(18.7, 16.8, 20.1, 22.4, 14.0, 15.2, 22.0, 23.3)
station     <- factor(c('Esso','Esso','Esso','Esso','Shell','Shell','Shell','Shell'))
fuel        <- factor(c('95','95','98','98','95','95','98','98'))
station_fuel<- factor(c('Esso95','Esso95','Esso98','Esso98','Shell95','Shell95','Shell98','Shell98'))

M             <- mean(km)
Mstation      <- tapply(km,      station, mean)
Mfuel         <- tapply(km,       fuel, mean)
Mstation_fuel <- tapply(km, station_fuel, mean)


plot(station_fuel, km, col=rainbow(5)[2:5], ylim=c(0,24))

# Parametric test:
summary.aov(aov(km ~ station + fuel + station:fuel))
# Without interaction
summary.aov(aov(km ~ station + fuel))
# Without station
summary.aov(aov(km ~ fuel))



# Permutation test
# The model we have to test is the following:
# km = mu + alpha*station + beta*fuel + gamma*station*fuel
# We have 3 different tests:
# 1. factor station   (H0: alpha=0)
# 2. factor fuel      (H0: beta=0)
# 3. interaction      (H0: gamma=0)

# We apply different permutations for developing the different tests!
# We start by testing the interaction:
# H0: gamma=0 against H_1: gamma!=0

# test statistic:
summary.aov(aov(km ~ station + fuel + station:fuel))
T0_station_fuel <- summary.aov(aov(km ~ station + fuel + station:fuel))[[1]][3,4]
T0_station_fuel

# permutation
# the idea is to permute the residuals under H0:
# km = mu + alpha*station + beta*fuel
# additive model
aov.H0station_fuel <- aov(km ~ station + fuel)
aov.H0station_fuel
residuals.H0station_fuel <- aov.H0station_fuel$residuals
n <- 8
permutation <- sample(1:n)
residuals.H0station_fuel <- residuals.H0station_fuel[permutation]
# permuted y values:
km.perm.H0station_fuel <- aov.H0station_fuel$fitted + residuals.H0station_fuel
summary.aov(aov(km.perm.H0station_fuel ~ station + fuel + station:fuel))
# How data has changed?
layout(rbind(1:2))
plot(station_fuel, km, col=rainbow(5)[2:5], ylim=c(0,24),main='Original data')
plot(station_fuel, km.perm.H0station_fuel, col=rainbow(5)[2:5], ylim=c(0,24),main='Permuted data')


# TEST of interaction
B <- 1000
T_station_fuel <- numeric(B)
for(perm in 1:B){
  permutation <- sample(n)
  residuals.H0station_fuel <- residuals.H0station_fuel[permutation]
  km.perm.H0station_fuel <- aov.H0station_fuel$fitted + residuals.H0station_fuel
  T_station_fuel[perm] <- summary.aov(aov(km.perm.H0station_fuel ~ station + fuel + station:fuel))[[1]][3,4]
}

# p-value
sum(T_station_fuel >= T0_station_fuel)/B
# The interaction is not significant. 
# We can remove it and perform a test for the two main effects

# TEST OF FACTOR STATION   (H0: alpha=0)
T0_station <- summary.aov(aov(km ~ station + fuel))[[1]][1,4]
# residuals under H0:
# km = mu + beta*fuel
aov.H0station <- aov(km ~ fuel)
residuals.H0station <- aov.H0station$residuals
# permuted y values:
km.perm.H0station <- aov.H0station$fitted + residuals.H0station[permutation]
summary.aov(aov(km.perm.H0station ~ station + fuel))
# How data has changed?
layout(rbind(1:2))
plot(station_fuel, km, col=rainbow(5)[2:5], ylim=c(0,24),main='Original data')
plot(station_fuel, km.perm.H0station, col=rainbow(5)[2:5], ylim=c(0,24),main='Permuted data')


# TEST OF FACTOR FUEL   (H0: beta=0)
T0_fuel <- summary.aov(aov(km ~ station + fuel))[[1]][2,4]
# residuals under H0:
# km = mu + alpha*station
aov.H0fuel <- aov(km ~ station)
residuals.H0fuel <- aov.H0fuel$residuals
# permuted y values:
km.perm.H0fuel <- aov.H0fuel$fitted + residuals.H0fuel[permutation]
summary.aov(aov(km.perm.H0fuel ~ station + fuel))
# How data has changed?
layout(rbind(1:2))
plot(station_fuel, km, col=rainbow(5)[2:5], ylim=c(0,24),main='Original data')
plot(station_fuel, km.perm.H0fuel, col=rainbow(5)[2:5], ylim=c(0,24),main='Permuted data')

# TEST OF FACTOR STATION AND TEST OF FACTOR FUEL
# p-values
B <- 1000
T_fuel <- T_station <- numeric(B)
for(perm in 1:B){
  permutation <- sample(n)
  
  km.perm.H0station <- aov.H0station$fitted + residuals.H0station[permutation]
  T_station[perm] <- summary.aov(aov(km.perm.H0station ~ station + fuel))[[1]][1,4]
  
  km.perm.H0fuel <- aov.H0fuel$fitted + residuals.H0fuel[permutation]
  T_fuel[perm] <- summary.aov(aov(km.perm.H0fuel ~ station + fuel))[[1]][2,4]
}

sum(T_station >= T0_station)/B
sum(T_fuel >= T0_fuel)/B

# Comparison with the parametric test
summary.aov(aov(km ~ station + fuel))

# We can remove also the factor station
# TEST ON THE FACTOR FUEL
T0_fuel <- summary.aov(aov(km ~  fuel))[[1]][1,4]
# residuals under H0
# km = mu
residuals.H0fuel <- km - M

# Note that in this case, permuting the residuals under H0 
# and permuting the data is exactly the same:
permutation <- sample(n)
km.perm.H0fuel <- M + residuals.H0fuel[permutation]
km.perm        <- km[permutation]

km.perm.H0fuel
km.perm

# CMC to estimate the p-value
B <- 1000
T_fuel <- numeric(B)
for(perm in 1:B){
  permutation <- sample(n)
  km.perm <- km[permutation]
  T_fuel[perm] <- summary.lm(aov(km.perm ~ fuel ))$f[1]
  
}
sum(T_fuel >= T0_fuel)/B
# We can reject H0 on the factor fuel

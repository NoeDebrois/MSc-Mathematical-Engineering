#################################
####    PERMUTATION TEST     ####
####     Linear models       ####
#################################

################################################################################
# Linear models ################################################################
# Example on simulated data ####################################################
################################################################################
set.seed(24021979)
n <- 50
# covariate values (we generate data here)
x1 <- runif(n,0,10)
x2 <- (1:n)/5
x3 <- rnorm(n,5,5)

# generating model (please notice the fact that b3 = 0 and we would like to test the nullity of the bi...)
b0 <- 2
b1 <- 4
b2 <- -2
b3 <- 0
Y <- b0 + b1*x1 + b2*x2 + b3*x3 + runif(n,-5,5)

plot(x1,Y,pch=16)
plot(x2,Y,pch=16)
plot(x3,Y,pch=16)

# parametric inference
result <- lm(Y ~ x1 + x2 + x3)
summary(result)
qqnorm(result$residuals)
shapiro.test(result$residuals)

################################################################################
# PERMUTATION INFERENCE : we want to perform different tests ###################
################################################################################

################################################################################
# Overall model (Global F-Test, see page 8/12 polycopié) #######################
# H0 : beta1 = beta2 = beta3 = 0 ###############################################
################################################################################
# test statistic (F Test see page 8/12 polycopié)
T0_glob <- summary(result)$f[1]
T0_glob

# permutations
permutazione <- sample(n)
Y.perm.glob <- Y[permutazione]

# in this case permuting the responses or the residuals is the same
res.H0glob <- Y - mean(Y) 
Y.perm.glob
mean(Y) + res.H0glob[permutazione]

layout(matrix(1:6,nrow=3,byrow=FALSE))
plot(x1,Y,main='Y vs x1 (original data)',pch=16)
plot(x2,Y,main='Y vs x2 (original data)',pch=16)
plot(x3,Y,main='Y vs x3 (original data)',pch=16)
plot(x1,Y.perm.glob,main='Y vs x1 (permuted data)',pch=16)
plot(x2,Y.perm.glob,main='Y vs x2 (permuted data)',pch=16)
plot(x3,Y.perm.glob,main='Y vs x3 (permuted data)',pch=16)

# OU BIEN ON PEUT FAIRE DES PARTIAL-TESTS POUR CHACUN DES bi :

################################################################################
# Test on variable x1 (partial T-Tests, see page 8/12 polycopié) ###############
# H0: beta1 = 0 ################################################################
################################################################################
# test statistic
summary(result)$coefficients
T0_x1 <- abs(summary(result)$coefficients[2,3])
T0_x1

# permutations
# residuals of the reduced model

# reduced model:
# Y = beta0 + beta2*x2 + beta3*x3
regr.H01 <- lm(Y ~ x2 + x3)
residui.H01 <- regr.H01$residuals
residui.H01.perm <- residui.H01[permutazione]

# permuted y:
Y.perm.H01 <- regr.H01$fitted + residui.H01.perm

layout(matrix(1:6,nrow=3,byrow=FALSE))
plot(x1,Y,main='Y vs x1 (original data)',pch=16)
plot(x2,Y,main='Y vs x2 (original data)',pch=16)
plot(x3,Y,main='Y vs x3 (original data)',pch=16)
plot(x1,Y.perm.H01,main='Y vs x1 (permuted data)',pch=16)
plot(x2,Y.perm.H01,main='Y vs x2 (permuted data)',pch=16)
plot(x3,Y.perm.H01,main='Y vs x3 (permuted data)',pch=16)

################################################################################
# Test on variable x2 (partial T-Tests, see page 8/12 polycopié) ###############
# H0: beta2 = 0 ################################################################
################################################################################
# test statistic
summary(result)$coefficients
T0_x2 <- abs(summary(result)$coefficients[3,3])
T0_x2

# permutations
# residuals of the reduced model

# reduced model:
# Y = beta0 + beta1*x1 + beta3*x3
regr.H02 <- lm(Y ~ x1 + x3)
residui.H02 <- regr.H02$residuals
residui.H02.perm <- residui.H02[permutazione]

# permuted y:
Y.perm.H02 <- regr.H02$fitted + residui.H02.perm

layout(matrix(1:6,nrow=3,byrow=FALSE))
plot(x1,Y,main='Y vs x1 (original data)',pch=16)
plot(x2,Y,main='Y vs x2 (original data)',pch=16)
plot(x3,Y,main='Y vs x3 (original data)',pch=16)
plot(x1,Y.perm.H02,main='Y vs x1 (permuted data)',pch=16)
plot(x2,Y.perm.H02,main='Y vs x2 (permuted data)',pch=16)
plot(x3,Y.perm.H02,main='Y vs x3 (permuted data)',pch=16)

################################################################################
# Test on variable x3 (partial T-Tests, see page 8/12 polycopié) ###############
# H0: beta3 = 0 ################################################################
################################################################################
# test statistic
summary(result)$coefficients
T0_x3 <- abs(summary(result)$coefficients[4,3])
T0_x3

# permutations
# residuals of the reduced model

# reduced model:
# Y = beta0 + beta1*x1 + beta2*x2
regr.H03 <- lm(Y ~ x1 + x2)
residui.H03 <- regr.H03$residuals
residui.H03.perm <- residui.H03[permutazione]

# permuted y:
Y.perm.H03 <- regr.H03$fitted + residui.H03.perm


layout(matrix(1:6,nrow=3,byrow=FALSE))
plot(x1,Y,main='Y vs x1 (original data)',pch=16)
plot(x2,Y,main='Y vs x2 (original data)',pch=16)
plot(x3,Y,main='Y vs x3 (original data)',pch=16)
plot(x1,Y.perm.H03,main='Y vs x1 (permuted data)',pch=16)
plot(x2,Y.perm.H03,main='Y vs x2 (permuted data)',pch=16)
plot(x3,Y.perm.H03,main='Y vs x3 (permuted data)',pch=16)

# MONTE CARLO :
# p-values of the tests
B <- 1000
T_H0glob <- T_H01 <- T_H02 <- T_H03 <- numeric(B)

for(perm in 1:B){
  permutazione <- sample(n)
  
  Y.perm.glob <- Y[permutazione]
  T_H0glob[perm] <- summary(lm(Y.perm.glob ~ x1 + x2 + x3))$f[1]
  
  residui.H01.perm <- residui.H01[permutazione]
  Y.perm.H01 <- regr.H01$fitted + residui.H01.perm
  T_H01[perm] <- abs(summary(lm(Y.perm.H01 ~ x1 + x2 + x3))$coefficients[2,3])
  
  residui.H02.perm <- residui.H02[permutazione]
  Y.perm.H02 <- regr.H02$fitted + residui.H02.perm
  T_H02[perm] <- abs(summary(lm(Y.perm.H02 ~ x1 + x2 + x3))$coefficients[3,3])
  
  residui.H03.perm <- residui.H03[permutazione]
  Y.perm.H03 <- regr.H03$fitted + residui.H03.perm
  T_H03[perm] <- abs(summary(lm(Y.perm.H03 ~ x1 + x2 + x3))$coefficients[4,3])
  
}

# We compute the p-values :
sum(T_H0glob>=T0_glob)/B
sum(T_H01>=T0_x1)/B
sum(T_H02>=T0_x2)/B
sum(T_H03>=T0_x3)/B
# Only H0 for b3 cannot be rejected.

# comparison with the parametric test
summary(result)

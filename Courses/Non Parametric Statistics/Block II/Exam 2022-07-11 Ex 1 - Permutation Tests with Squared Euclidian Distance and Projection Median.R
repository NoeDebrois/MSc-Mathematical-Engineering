################################################################################
################################ EXERCICE I ####################################
################################################################################
library(DepthProc)

percebes <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-07-11/percebes_1.rds")

# Dr. Andreas Qapos Ph.D. is a Galician mathematician/fisherman, who proves 
# theorems by day, and picks up goose barnacles, percebes in Spanish, an 
# extremely pricey seafood delicacy typical of Galicia, by night. Goose 
# barnacles live on rocks washed by the strong waves of the Atlantic Ocean.
# Dr. Qapos has identified two spots on the cliffs close to his hometown, 
# Pontevedra, identified (since they are very secret) as spot A and spot B.

# Dr. Qapos suspects that the goose barnacles picked out of spot B tend to be
# slightly shorter than the ones in spot A, but plumpier. For this reason he has
# collected some barnacles from spot A and some (less, access to the cliff is 
# terrible!) from spot B, and he measured their weight [g] and length [mm]. 
# You can find the data in percebes_1.rds.
# Now, assuming the tuple in each group (weight, length) being i.i.d:

# 1. To help Dr. Qapos in assessing his hypothesis, start by computing the 
# projection sample medians of the two groups.
p.median.a = depthMedian(percebes$spot.A,depth_params = list(method='Projection'))
p.median.b = depthMedian(percebes$spot.B,depth_params = list(method='Projection'))
p.median.a
p.median.b

# 2. To test the equality of the two theoretical projection medians, perform a 
# two-sample permutation test using as a test statistics the squared euclidean 
# distance between the two sample projection medians. Please describe briefly 
# the properties of permutation tests and present the empirical cumulative 
# distribution function of the permutational test statistic as well as p-value 
# for the test.

# Extract two subsets for testing
t1 <- percebes$spot.A 
t2 <- percebes$spot.B 

# Compute a proper test statistic (squared distance between two sample mean vectors)
n1 <- dim(t1)[1]
n2 <- dim(t2)[1]
n  <- n1 + n2
T20 <- as.numeric((p.median.a - p.median.b) %*% (p.median.a - p.median.b))

# Estimate the permutational distribution under H0
B <- 1000
T2 <- numeric(B)

for(perm in 1:B) {
  # Random permutation of indexes
  t_pooled <- rbind(t1, t2)
  permutation <- sample(n)
  t_perm <- t_pooled[permutation,]
  t1_perm <- t_perm[1:n1,]
  t2_perm <- t_perm[(n1 + 1):n,]
  
  # Evaluation of the test statistic on permuted data
  t1_perm.median.a = depthMedian(t1_perm,depth_params = list(method='Projection'))
  t2_perm.median.b = depthMedian(t2_perm,depth_params = list(method='Projection'))
  T2[perm] <- as.numeric((t1_perm.median.a - t2_perm.median.b) %*% (t1_perm.median.a - t2_perm.median.b))
}

# Plot the permutational distribution under H0
hist(T2, xlim=range(c(T2, T20)), breaks=1000)
abline(v=T20, col=3, lwd=4)

plot(ecdf(T2))
abline(v=T20, col=3, lwd=4)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
p_val
# We have to reject H0 so : the two theoretical projection medians are not equal

# 3. Is Dr. Qapos justified in risking his life to go for spot B, given that the
# barnacles are sold by their weight and not by their length? Justify your answer.
# Yes ! The projection medians are not equal and we see that the weights in spot B
# are bigger than the weights in spot A.




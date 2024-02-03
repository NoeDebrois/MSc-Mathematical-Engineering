#####################################
####      PERMUTATION TESTS      ####
##   TWO INDEPENDENT POPULATIONS   ##
#####################################
# We sample n1 data from a first population X1 
# and n2 data from a second population  X2
# Test:
# H0:  X1 =^d  X2
# H1:  X1 !=^d X2

# Clear workspace
rm(list=ls())

# Define parameters
n1 <- n2 <- 10  # Sample sizes for two populations
n <- n1 + n2    # Total sample size

# Set seed for reproducibility
set.seed(240279)

################################################################################
# Case 1: H0 FALSE - Populations with different means ##########################
################################################################################
# To understand the behaviour of the test, we sample from two populations 
# with different means
x1 <- runif(n1, 0, 4) # Generate data from first population
x2 <- runif(n2, 0, 4) + 3 # Generate data from second population with different means

x_pooled <- c(x1, x2) # Combine data from both populations

# Visualize original data
par(mfrow=c(1, 2))
boxplot(x1, x2, main='Original data')

# Apply one random permutation to see how data changes
permutation <- sample(1:n) 
x_perm <- x_pooled[permutation] 
x1_perm <- x_perm[1:n1]
x2_perm <- x_perm[(n1 + 1):n]

# Visualize permuted data
boxplot(x1_perm, x2_perm, main='Permuted data')

# Calculate absolute mean differences for original and permuted data
print("Absolute mean differences for original data :")
abs(mean(x1) - mean(x2))
print("Absolute mean differences for permuted data :")
abs(mean(x1_perm) - mean(x2_perm))
# The mean difference is lower.
# Was that a case ? We need to perform the test in order to be sure...

################################################################################
# Perform permutation test : ###################################################
################################################################################
T0 <- abs(mean(x1) - mean(x2))  # Test statistic for original data
print("Test statistic for original data :")
T0

# Cardinality of the permutational space
print("Cardinality of the permutational space :")
factorial(n)

# Number of distinct values of T* under H0
print("Number of distinct values of T* under H0 :")
factorial(n) / (2 * factorial(n1) * factorial(n2))

# Minimum achievable p-value
print("Minimum achievable p-value :")
1 / (factorial(n) / (2 * factorial(n1) * factorial(n2)))

# Monte Carlo simulation to estimate the p-value
B <- 100000            # Number of permutations
T_stat <- numeric(B)   # Vector to store test statistics

# Loop through permutations to calculate test statistics
for (perm in 1:B) {
  permutation <- sample(1:n)
  x_perm <- x_pooled[permutation]
  x1_perm <- x_perm[1:n1]
  x2_perm <- x_perm[(n1 + 1):n]
  T_stat[perm] <- abs(mean(x1_perm) - mean(x2_perm))
}

# Visualize permutational distribution of T
hist(T_stat, xlim=range(c(T_stat, T0)), breaks=30)
abline(v=T0, col=3, lwd=2)

# Plot empirical cumulative distribution function (ECDF)
plot(ecdf(T_stat))
abline(v=T0, col=3, lwd=2)

# Calculate p-value
p_val <- sum(T_stat >= T0) / B
print("The P-Value for this test is :")
p_val
# It's too small, we have to reject H0.


################################################################################
# Case 2: H0 TRUE - Populations with the same distribution #####################
################################################################################
x1 <- runif(n1, 0, 4) # Generate data from first population
x2 <- runif(n2, 0, 4) # Generate data from second population with the same distribution

x_pooled <- c(x1, x2) # Combine data from both populations

# Visualize original data
par(mfrow=c(1, 2))
boxplot(x1, x2, main='Original data')

# Apply one random permutation to see how data changes
permutation <- sample(1:n) 
x_perm <- x_pooled[permutation] 
x1_perm <- x_perm[1:n1]
x2_perm <- x_perm[(n1 + 1):n]

# Visualize permuted data
boxplot(x1_perm, x2_perm, main='Permuted data')

# Calculate absolute mean differences for original and permuted data
abs(mean(x1) - mean(x2))
abs(mean(x1_perm) - mean(x2_perm))
# The means are close.
# Is it significant ? We need to perform the test in order to be sure...

################################################################################
# Perform permutation test #####################################################
################################################################################
T0 <- abs(mean(x1) - mean(x2))  # Test statistic for original data

# Cardinality of the permutational space
factorial(n)

# Number of distinct values of T* under H0
factorial(n) / (2 * factorial(n1) * factorial(n2))

# Minimum achievable p-value
1 / (factorial(n) / (2 * factorial(n1) * factorial(n2)))

# Monte Carlo simulation to estimate the p-value
B <- 100000            # Number of permutations
T_stat <- numeric(B)   # Vector to store test statistics

# Loop through permutations to calculate test statistics
for (perm in 1:B) {
  permutation <- sample(1:n)
  x_perm <- x_pooled[permutation]
  x1_perm <- x_perm[1:n1]
  x2_perm <- x_perm[(n1 + 1):n]
  T_stat[perm] <- abs(mean(x1_perm) - mean(x2_perm))
}

# Visualize permutational distribution of T
hist(T_stat, xlim=range(c(T_stat, T0)), breaks=30)
abline(v=T0, col=3, lwd=2)

# Plot empirical cumulative distribution function (ECDF)
plot(ecdf(T_stat))
abline(v=T0, col=3, lwd=2)

# Calculate p-value
p_val <- sum(T_stat >= T0) / B
print("The P-Value for this test is :")
p_val
# It's quite big, we can't reject H0.

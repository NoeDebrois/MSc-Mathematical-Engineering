#####################################
####      PERMUTATION TESTS      ####
####       (MULTIVARIATE)        ####
#####################################

# In this file :
# 1 Two independent multivariate population test (page 1 cours Permutation tests)
# 2 One sample multivariate center of symmetry test (page 11 cours Permutation tests - 2ème partie)
# 3 Two paired multivariate population test (page 11 cours Permutation tests - 2ème partie)

################################################################################
# Example 1: Two (independent) multivariate population test ####################
# Hourly accesses to AreaC: working days vs week-end days ######################
################################################################################

# (paired dans le sens où l'on regarde chaque jour aux mêmes instants (toutes les 30'))

# Load data from CSV files
d1 <- read.csv('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/accessi-orari-areac-2016-09-12-00_00_00.csv', header=T)
d2 <- read.csv('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/accessi-orari-areac-2016-09-13-00_00_00.csv', header=T)
d3 <- read.csv('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/accessi-orari-areac-2016-09-14-00_00_00.csv', header=T)
d4 <- read.csv('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/accessi-orari-areac-2016-09-15-00_00_00.csv', header=T)
d5 <- read.csv('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/accessi-orari-areac-2016-09-16-00_00_00.csv', header=T)
d6 <- read.csv('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/accessi-orari-areac-2016-09-17-00_00_00.csv', header=T)
d7 <- read.csv('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/accessi-orari-areac-2016-09-18-00_00_00.csv', header=T)

# Combine data from all days into a single matrix
# C'est là qu'on crée le 'multivariate' (pour tout i, X_i \in R^7)
week <- rbind(d1[, 2], d2[, 2], d3[, 2], d4[, 2], d5[, 2], d6[, 2], d7[, 2])

# Plot the data for each day
matplot(seq(0, 47)/2, t(week), type='l', col=c(1, 1, 1, 1, 1, 2, 2), lty=1)

# Extract two subsets for testing
t1 <- week[1:5,] # lundi à vendredi
t2 <- week[6:7,] # samedi & dimanche

# Compute a proper test statistic (squared distance between two sample mean vectors)
t1.mean <- colMeans(t1)
t2.mean <- colMeans(t2)
n1 <- dim(t1)[1]
n2 <- dim(t2)[1]
n  <- n1 + n2
T20 <- as.numeric((t1.mean - t2.mean) %*% (t1.mean - t2.mean))

# Number of possible data point permutations and different values of the test statistic
factorial(7)
choose(7, 5)

# Estimate the permutational distribution under H0
B <- 100000
T2 <- numeric(B)

for(perm in 1:B) {
  # Random permutation of indexes
  t_pooled <- rbind(t1, t2)
  permutation <- sample(n)
  t_perm <- t_pooled[permutation,]
  t1_perm <- t_perm[1:n1,]
  t2_perm <- t_perm[(n1 + 1):n,]
  
  # Evaluation of the test statistic on permuted data
  t1.mean_perm <- colMeans(t1_perm)
  t2.mean_perm <- colMeans(t2_perm)
  T2[perm] <- as.numeric((t1.mean_perm - t2.mean_perm) %*% (t1.mean_perm - t2.mean_perm))
}

# Plot the permutational distribution under H0
hist(T2, xlim=range(c(T2, T20)), breaks=1000)
abline(v=T20, col=3, lwd=4)

plot(ecdf(T2))
abline(v=T20, col=3, lwd=4)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
p_val

################################################################################
# Example 2: Center of symmetry of one multivariate population #################
# Relative humidity in Milan during the summer months ##########################
################################################################################

# Read data
hum <- read.csv2('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/307_Umidita_relativa_2008_2014.csv', header=T)
hum <- hum[,3] # on garde que l'humidité moyenne
hum <- matrix(hum, ncol=12, byrow=T)[,6:9] # transforme le vecteur de longueur 84 en une matrice
# de 12 colonnes (chaque colonne est un mois), puis ne garde que les 4 colonnes 6 à 9.

# Boxplot and line plot of humidity data (POUR LES 4 MOIS QU'ON A GARDES)
boxplot(hum) 
matplot(t(hum), type='l', lty=1)

# Center of symmetry under H0
mu0 <- c(65, 65, 65, 65)

# Compute a proper test statistic (squared distance between sample mean vector and hypothesized center of symmetry)
x.mean <- colMeans(hum)
n <- dim(hum)[1]
p <- dim(hum)[2]
T20 <- as.numeric((x.mean - mu0) %*% (x.mean - mu0))

# Number of possible data point reflections and different values of the test statistic
2^7 # (cf le cours)
2^7/2

# Estimate the permutational distribution under H0
B <- 100000
T2 <- numeric(B) 
for(perm in 1:B) {
  # In this case, use changes of signs in place of permutations
  signs.perm <- rbinom(n, 1, 0.5)*2 - 1 # créer des ±1 aléatoirement
  hum_perm <- matrix(mu0, nrow=n, ncol=p, byrow=T) + (hum - matrix(mu0, nrow=n, ncol=p, byrow=T)) * matrix(signs.perm, nrow=n, ncol=p, byrow=FALSE)
  x.mean_perm <- colMeans(hum_perm)
  T2[perm] <- as.numeric((x.mean_perm - mu0) %*% (x.mean_perm - mu0))
}

# Plot the permutational distribution under H0
hist(T2, xlim=range(c(T2, T20)), breaks=100)
abline(v=T20, col=3, lwd=4)

plot(ecdf(T2))
abline(v=T20, col=3, lwd=4)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
p_val
# On peut rejeter "H0 : center of symmetry = C0".

################################################################################
# Example 3: Two paired multivariate population test ###########################
# (Same days for Milan and Barcelona) ##########################################
################################################################################
# The data set contains observations of temperature, humidity, and wind in
# Milan and Barcelona on 50 different days

# Read data
t1 <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/barcellona.txt', header=T)
t2 <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/02 Permutation Tests/milano.txt', header=T)

# 3D plot of differences between Milan and Barcelona
library(rgl)
open3d()
plot3d(t1 - t2, size=3, col='orange', aspect = F)
points3d(0, 0, 0, size=6)
rglwidget()

p  <- dim(t1)[2]
n1 <- dim(t1)[1]
n2 <- dim(t2)[1]
n <- n1 + n2

# Evaluate the test statistic
t1.mean <- colMeans(t1)
t2.mean <- colMeans(t2)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1) * t1.cov + (n2-1) * t2.cov) / (n1 + n2 - 2)
Spinv   <- solve(Sp)

delta.0 <- c(0, 0, 0)

diff <- t1 - t2
diff.mean <- colMeans(diff)
diff.cov <- cov(diff)
diff.invcov <- solve(diff.cov)

# Test statistic
T20 <- as.numeric(n1 * (diff.mean - delta.0) %*% diff.invcov %*% (diff.mean - delta.0))

# Number of possible data point reflections
2^50

# Estimate the permutational distribution under H0
B <- 10000
T2 <- numeric(B)

for(perm in 1:B) {
  # Random permutation
  signs.perm <- rbinom(n1, 1, 0.5)*2 - 1
  diff_perm <- diff * matrix(signs.perm, nrow=n1, ncol=p, byrow=FALSE)
  diff.mean_perm <- colMeans(diff_perm)
  diff.cov_perm <- cov(diff_perm)
  diff.invcov_perm <- solve(diff.cov_perm)
  
  T2[perm] <- as.numeric(n1 * (diff.mean_perm - delta.0) %*% diff.invcov_perm %*% (diff.mean_perm - delta.0))
}

# Plot the permutational distribution under H0
hist(T2, xlim=range(c(T2, T20)), breaks=100)
abline(v=T20, col=3, lwd=4)

plot(ecdf(T2))
abline(v=T20, col=3, lwd=4)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
p_val
# Again we can reject H0


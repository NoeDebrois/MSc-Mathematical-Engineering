################################################################################
################################ EXERCICE I ####################################
################################################################################
library(DepthProc)

milk_samples_1 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-01-21/milk_samples_1.Rds")

N = 382

# 1. Provide the Tukey median of the milk samples
tukey_depth = depth(u=milk_samples_1,method='Tukey')
# tukey_median <- milk_samples_1[which.max(tukey_depth),]
tukey_median <- depthMedian(milk_samples_1,depth_params = list(method='Tukey'))

# 2. Plot a bagplot matrix of the three collected variables, and to determine a 
# vector of row indexes identifying the milk samples that are outliers according
# to any panel of the bagplot matrix1. Report the plot of the bagplot matrix and
# briefly describe the process that led to the generation of the aforementioned 
# vector of row indexes.
library(aplpack)

bagplot_matrix <- aplpack::bagplot.pairs(milk_samples_1)

## Voir comment ils font en dim >2 les outlier detections

## Ma méthode, à la mano :

df_1 = milk_samples_1[1:2]
df_2 = milk_samples_1[2:3]
df_3 = milk_samples_1[c(1, 3)]

bagplot_1 = aplpack::bagplot(df_1,show.whiskers = F,main="Bagplot")
bagplot_2 = aplpack::bagplot(df_2,show.whiskers = F,main="Bagplot")
bagplot_3 = aplpack::bagplot(df_3,show.whiskers = F,main="Bagplot")

outlying_obs_1 <- bagplot_1$pxy.outlier
outlying_obs_2 <- bagplot_2$pxy.outlier
outlying_obs_3 <- bagplot_3$pxy.outlier

ind_outlying_obs_1 <- which(apply(df_1,1,function(x) all(x %in% outlying_obs_1))) 
ind_outlying_obs_1

ind_outlying_obs_2 <- which(apply(df_2,1,function(x) all(x %in% outlying_obs_2))) 
ind_outlying_obs_2

ind_outlying_obs_3 <- which(apply(df_3,1,function(x) all(x %in% outlying_obs_3))) 
ind_outlying_obs_3

outliers = unique(c(ind_outlying_obs_1, ind_outlying_obs_2, ind_outlying_obs_3))
outliers = sort(outliers)
outliers

# 3. Test whether the 341 milk samples, obtained by discarding the 41 units that
# were flagged as outliers by the bagplots, comply with the gold standard in terms
# of milk quality, for which κ-casein must be equal to 6 grams per liter, CMS to 174 nm
# and Milk pH to 7. 
# Perform a permutation test using as test statistic the squared
# euclidean distance between the sample mean and the gold standard. Specify assumptions,
# the null and the alternative hypthotesis you are testing, provide the histogram of the
# permuted distribution of the test statistic, its cumulative distribution function and
# the p-value, commenting the results.

# Multivariate permutation test for center of symmetry :
# H0 : C = gold_standard
# H1 : C != gold_standard 

df_without_outliers = milk_samples_1[-outliers, ]
gold_standard = c(6, 174, 7)

# Compute a proper test statistic (squared distance between sample mean vector and hypothesized center of symmetry)
x.mean <- colMeans(df_without_outliers)
n <- dim(df_without_outliers)[1]
p <- dim(df_without_outliers)[2]

#permutation test
# Computing a proper test statistic
# (i.e., L2 norm between the sample mean vector and the hypothesized center of symmetry)
T20 <- as.numeric((x.mean - gold_standard) %*% (x.mean - gold_standard))

# Estimate the permutational distribution under H0
B <- 1000
T2 <- numeric(B) 

set.seed(2022)
for(perm in 1:B) {
  # In this case, use changes of signs in place of permutations
  signs.perm <- rbinom(n, 1, 0.5)*2 - 1 # créer des ±1 aléatoirement
  df_without_outliers_perm = matrix(gold_standard, nrow=n, ncol=p, byrow=T) + (df_without_outliers - matrix(gold_standard, nrow=n, ncol=p, byrow=T)) * matrix(signs.perm, nrow=n, ncol=p, byrow=FALSE)
  x.mean_perm <- colMeans(df_without_outliers_perm)
  T2[perm] <- as.numeric((x.mean_perm - gold_standard) %*% (x.mean_perm - gold_standard))
}

# Plot the permutational distribution under H0
hist(T2, xlim=range(c(T2, T20)), breaks=100)
abline(v=T20, col=3, lwd=4)

plot(ecdf(T2))
abline(v=T20, col=3, lwd=4)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
p_val
# We cannot reject H0.

# 4. Perform the exact same test but this time consider the original N = 382 milk samples.
# Does your conclusion change? If so, what kind of estimator would you propose that does not
# imply the a-priori removal of a portion of the data?

# Compute a proper test statistic (squared distance between sample mean vector and hypothesized center of symmetry)
x.mean <- colMeans(milk_samples_1)
n <- dim(milk_samples_1)[1]
p <- dim(milk_samples_1)[2]
T20 <- as.numeric((x.mean - gold_standard) %*% (x.mean - gold_standard))

# Estimate the permutational distribution under H0
B <- 1000
T2 <- numeric(B) 
for(perm in 1:B) {
  # In this case, use changes of signs in place of permutations
  signs.perm <- rbinom(n, 1, 0.5)*2 - 1 # créer des ±1 aléatoirement
  milk_samples_1_perm = matrix(gold_standard, nrow=n, ncol=p, byrow=T) + (milk_samples_1 - matrix(gold_standard, nrow=n, ncol=p, byrow=T)) * matrix(signs.perm, nrow=n, ncol=p, byrow=FALSE)
  x.mean_perm <- colMeans(milk_samples_1_perm)
  T2[perm] <- as.numeric((x.mean_perm - gold_standard) %*% (x.mean_perm - gold_standard))
}

# Plot the permutational distribution under H0
hist(T2, xlim=range(c(T2, T20)), breaks=100)
abline(v=T20, col=3, lwd=4)

plot(ecdf(T2))
abline(v=T20, col=3, lwd=4)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
p_val
# We have to reject H0. The conclusion changes.
# # Possible solution:
# Use a different test statistic considering testing, for example,
# the equality of the Tukey median to the gold standard.

# 5. Provide a Full Conformal 1 − α = 90% prediction region for the κ-casein and
# Milk pH of a new milk sample, using the euclidean distance between the new data 
# point and the sample Tukey median of the augmented data set as non-conformity measure.
# After having discussed the theoretical properties of the prediction region, provide a plot of it.

# PROBLEME : A CHECKER !

data_predict = df_without_outliers[c(1, 3)]

plot(1, type="n", ylim=c(6, 7), xlab="Votre étiquette pour l'axe x", ylab="Votre étiquette pour l'axe y")
# Ajouter les points de votre jeu de données
plot(data_predict)

# Initialiser les paramètres et les grilles
alpha <- 0.1
x.obs <- data_predict
x1.new.grid <- seq(min(x.obs[,1]) - 0.25*diff(range(x.obs[,1])), max(x.obs[,1]) + 0.25*diff(range(x.obs[,1])), length = 20)
x2.new.grid <- seq(min(x.obs[,2]) - 0.25*diff(range(x.obs[,2])), max(x.obs[,2]) + 0.25*diff(range(x.obs[,2])), length = 20)
p.value <- matrix(nrow = length(x1.new.grid), ncol = length(x2.new.grid))

# Définir la fonction de prédiction non conforme (NC)
NC <- function(z.aug, i){
  sum((z.aug[i,] - colMeans(z.aug[-i,]))^2) # Euclidean distance
}

# Boucles pour calculer les p-values
for(k in 1:length(x1.new.grid)) {
  for(h in 1:length(x2.new.grid)) {
    # Ajouter une nouvelle observation aux données existantes
    x.obs.aug <- rbind(x.obs, c(x1.new.grid[k],x2.new.grid[h]))
    
    # Initialiser un vecteur pour stocker les scores
    scores <- numeric(dim(x.obs.aug)[1])
    
    # Calculer les scores non conformes pour chaque observation augmentée
    for (i in 1:dim(x.obs.aug)[1]) {
      scores[i] <- NC(x.obs.aug, i)
    }
    
    # Calculer la p-value pour la nouvelle observation
    p.value[k,h] <- sum(scores >= scores[dim(x.obs.aug)[1]])/(dim(x.obs.aug)[1])
    
    # Imprimer les indices de la grille pour suivre la progression
    print(c(k,h))
  }
}

# Tracer les p-values et la région de prédiction
image(x1.new.grid, x2.new.grid, p.value, zlim=c(0,1), asp=1)
points(data_predict, pch=16)
contour(x1.new.grid, x2.new.grid, p.value, levels = alpha, add=T)











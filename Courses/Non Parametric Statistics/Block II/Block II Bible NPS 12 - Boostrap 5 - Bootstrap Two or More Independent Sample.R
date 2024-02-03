################################################################################
### Bootstrap for Independent Samples ##########################################
################################################################################
# Ce code utilise la méthode bootstrap pour estimer la distribution de la ######
# différence entre les médianes de deux échantillons indépendants ##############
# (échantillons associés aux variables x1 et x2). Ensuite, il trace la #########
# fonction de répartition empirique de cette différence et calcule un ##########
# intervalle de confiance en utilisant la méthode des Reverse Percentile #######
# Intervals (RP intervals). L'intervalle de confiance est ensuite tracé ########
# sur la fonction de répartition empirique. ####################################
################################################################################
rm(list=ls())

# Import data
grades <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/03 Bootstrap/parziali.txt', header=T)
OE <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/03 Bootstrap/matricola.txt', header=T)

x1 <- grades[OE == 1,'PI']
x2 <- grades[OE == 0,'PI']

# Plot data
par(mfrow=c(1,2))
boxplot(x1, ylim=range(c(x1,x2)), main = 'Odd')
boxplot(x2, ylim=range(c(x1,x2)), main = 'Even')

################################################################################
### Computing the bootstrap distribution of the difference of the two sample ###
### medians ####################################################################
################################################################################
# Fixer la graine aléatoire pour la reproductibilité
set.seed(24021979)
# Nombre d'itérations Bootstrap
B <- 10000
# Initialiser des vecteurs pour stocker les estimations Bootstrap de la différence de médianes
x1.obs <- x1
x2.obs <- x2
diff.Q2.obs <- quantile(x1, 0.50) - quantile(x2, 0.50) # différence des médianes

T.boot.diff.Q2 <- numeric(B)

# Boucle sur les itérations Bootstrap
for(b in 1:B)
{
  # Ré-échantillonner les deux échantillons avec remplacement
  x1.b <- sample(x1.obs, replace = TRUE)
  x2.b <- sample(x2.obs, replace = TRUE)
  # Calculer la différence de médianes pour les données Bootstrap
  T.boot.diff.Q2[b] <- quantile(x1.b, 0.50) - quantile(x2.b, 0.50)
}

# Tracer la fonction de répartition empirique de la différence de médianes
plot(ecdf(T.boot.diff.Q2), main='Sample Median Odd - Sample Median Even')
abline(v = diff.Q2.obs, lty=2)

################################################################################
### Reverse Percentile Intervals (RP intervals) ################################
################################################################################
# Niveau de confiance
alpha <- 0.05

# Calculer les quantiles pour la différence de médianes
right.quantile <- quantile(T.boot.diff.Q2, 1 - alpha/2)
left.quantile  <- quantile(T.boot.diff.Q2, alpha/2)

# Valeur observée
diff.Q2.obs
# Calculer l'intervalle de confiance Reverse Percentile
CI.RP <- c(diff.Q2.obs - (right.quantile - diff.Q2.obs), diff.Q2.obs - (left.quantile - diff.Q2.obs))
CI.RP

# Tracer la fonction de répartition empirique avec l'intervalle de confiance
abline(v = CI.RP)


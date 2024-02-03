################################################################################
### Bootstrap Methods ##########################################################
### Page 5/6 du cours sur le Bootstrap #########################################
################################################################################
# Ce code génère des échantillons de distribution uniforme, puis applique ######
# les méthodes de bootstrap naive, smooth, et parametric pour estimer les ######
# fonctions de répartition empiriques des quartiles. Les résultats sont ########
# ensuite comparés aux quartiles de l'échantillon original. Les graphiques #####
# résultants permettent de visualiser la variabilité des estimations ###########
# bootstrap par rapport aux quartiles de l'échantillon original. ###############
################################################################################
# Effacer toutes les variables de l'environnement
rm(list=ls())

# Simulation: Unif(0,4)

################################################################################
### Naive Bootstrap ############################################################
################################################################################
# Fixer la graine aléatoire pour la reproductibilité
set.seed(24021979)
# Tailles d'échantillons
n.grid  <- 2^(1:8)

# Diviser la fenêtre graphique en 2 lignes et 4 colonnes
par(mfrow = c(2,4))
# Boucle sur différentes tailles d'échantillons
for(n in n.grid)
{
  # Générer un échantillon de distribution uniforme
  x <- runif(n, 0, 4)
  # Tracer la fonction de répartition empirique
  plot(ecdf(x), main=paste('n =', n), xlim = c(-1,5))
  # Tracer la vraie fonction de répartition (Unif(0,4))
  lines(seq(-1, 5, by=0.1), punif(seq(-1, 5, by=0.1),0,4), type='l', col='red', lty=2)
}

################################################################################
### Smooth Bootstrap ###########################################################
################################################################################
# Fixer une nouvelle graine aléatoire pour la reproductibilité
set.seed(24021979)
# Réinitialiser les tailles d'échantillons
n.grid  <- 2^(1:8)

# Diviser la fenêtre graphique en 2 lignes et 4 colonnes
par(mfrow = c(2,4))
# Boucle sur différentes tailles d'échantillons
for(n in n.grid)
{
  # Générer un échantillon de distribution uniforme
  x <- runif(n, 0, 4)
  # Fonction de lissage (kernel density estimation)
  smooth.ecdf <- function(x.0){sum(pnorm(x.0, x, 0.5))}
  # Tracer la fonction de répartition empirique lissée
  plot(seq(-1,5,by=0.1), as.numeric(apply(data.frame(matrix(seq(-1,5,by=0.1),1)), 2, smooth.ecdf))/n, main=paste('n =', n), xlim = c(-1,5), type='l', ylim=c(0,1))
  # Tracer la vraie fonction de répartition (Unif(0,4))
  lines(seq(-1, 5, by=0.1), punif(seq(-1, 5, by=0.1),0,4), type='l', col='red', lty=2)
}

################################################################################
### Parametric Bootstrap #######################################################
################################################################################
# Fixer une nouvelle graine aléatoire pour la reproductibilité
set.seed(24021979)
# Réinitialiser les tailles d'échantillons
n.grid  <- 2^(1:8)

# Diviser la fenêtre graphique en 2 lignes et 4 colonnes
par(mfrow = c(2,4))
# Boucle sur différentes tailles d'échantillons
for(n in n.grid)
{
  # Générer un échantillon de distribution uniforme
  x <- runif(n, 0, 4)
  # Fonction de lissage (kernel density estimation)
  smooth.ecdf <- function(x.0){sum(pnorm(x.0, x, 0.5))}
  # Déterminer les bornes de la distribution uniforme
  min.unif <- min(x)
  max.unif <- max(x)
  # Tracer la fonction de répartition empirique de la distribution uniforme
  plot(seq(-1, 5, by=0.1), punif(seq(-1, 5, by=0.1), min.unif, max.unif), main=paste('n =', n), xlim = c(-1,5), type='l', ylim=c(0,1))
  # Tracer la vraie fonction de répartition (Unif(0,4))
  lines(seq(-1, 5, by=0.1), punif(seq(-1, 5, by=0.1),0,4), type='l', col='red', lty=2)
}

################################################################################
################################################################################
### Naive, Smooth, and Parametric Bootstrap Estimation of quartiles ############
################################################################################
################################################################################
# Importer les données depuis un fichier
grades <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/03 Bootstrap/parziali.txt', header=T)
grades.PI <- grades[,'PI']
# Tracer une boîte à moustaches des données
boxplot(grades.PI)

# Calculer les quartiles de l'échantillon original
Q1 <- quantile(grades.PI, 0.25)
Q2 <- quantile(grades.PI, 0.50)
Q3 <- quantile(grades.PI, 0.75)

# Diviser la fenêtre graphique en 3 lignes et 1 colonne
par(mfrow=c(3,1))

################################################################################
### Naive Bootstrap ############################################################
################################################################################
# Fixer une nouvelle graine aléatoire pour la reproductibilité
set.seed(24021979)
# Nombre d'itérations Bootstrap
B <- 10000
# Échantillon observé
x.obs <- grades.PI
# Initialiser des vecteurs pour stocker les estimations Bootstrap des quartiles
T.boot.Q1 <- numeric(B)
T.boot.Q2 <- numeric(B)
T.boot.Q3 <- numeric(B)

# Boucle sur les itérations Bootstrap
for(b in 1:B)
{
  # Échantillon Bootstrap
  x.b <- sample(x.obs, replace = TRUE)
  # Estimation Bootstrap des quartiles
  T.boot.Q1[b] <- quantile(x.b, 0.25)
  T.boot.Q2[b] <- quantile(x.b, 0.50)
  T.boot.Q3[b] <- quantile(x.b, 0.75)
}

# Tracer les fonctions de répartition empiriques des quartiles
plot(ecdf(T.boot.Q1), main='Naive Bootstrap: Sample Q1, Q2, and Q3', col='green', xlim=c(15,26))
lines(ecdf(T.boot.Q2), col='black')
lines(ecdf(T.boot.Q3), col='red')

# Tracer les lignes verticales pour les quartiles de l'échantillon original
abline(v=c(Q1,Q2,Q3), col=c('green','black','red'), lty=2)

################################################################################
### Smooth Bootstrap ###########################################################
################################################################################
# Fixer une nouvelle graine aléatoire pour la reproductibilité
set.seed(24021979)
# Nombre d'itérations Bootstrap
B <- 10000
# Échantillon observé
x.obs <- grades.PI
# Initialiser des vecteurs pour stocker les estimations Bootstrap des quartiles
T.boot.Q1 <- numeric(B)
T.boot.Q2 <- numeric(B)
T.boot.Q3 <- numeric(B)

# Boucle sur les itérations Bootstrap
for(b in 1:B)
{
  # Échantillon Bootstrap avec ajout de bruit gaussien
  x.b <- sample(x.obs, replace = TRUE) + rnorm(length(x.obs), 0, 0.25)
  # Estimation Bootstrap des quartiles
  T.boot.Q1[b] <- quantile(x.b, 0.25)
  T.boot.Q2[b] <- quantile(x.b, 0.50)
  T.boot.Q3[b] <- quantile(x.b, 0.75)
}

# Tracer les fonctions de répartition empiriques des quartiles
plot(ecdf(T.boot.Q1), main='Smooth Bootstrap (Gaussian kernel): Sample Q1, Q2, and Q3', col='green', xlim=c(15,26))
lines(ecdf(T.boot.Q2), col='black')
lines(ecdf(T.boot.Q3), col='red')

# Tracer les lignes verticales pour les quartiles de l'échantillon original
abline(v=c(Q1,Q2,Q3), col=c('green','black','red'), lty=2)

################################################################################
### Parametric Bootstrap (assuming Gaussian data) ##############################
################################################################################
# Fixer une nouvelle graine aléatoire pour la reproductibilité
set.seed(24021979)
# Nombre d'itérations Bootstrap
B <- 10000
# Échantillon observé
x.obs <- grades.PI
# Initialiser des vecteurs pour stocker les estimations Bootstrap des quartiles
T.boot.Q1 <- numeric(B)
T.boot.Q2 <- numeric(B)
T.boot.Q3 <- numeric(B)

# Boucle sur les itérations Bootstrap
for(b in 1:B)
{
  # Échantillon Bootstrap en générant des valeurs aléatoires suivant une distribution normale
  x.b <- rnorm(length(x.obs), mean(x.obs), sd(x.obs))
  # Estimation Bootstrap des quartiles
  T.boot.Q1[b] <- quantile(x.b, 0.25)
  T.boot.Q2[b] <- quantile(x.b, 0.50)
  T.boot.Q3[b] <- quantile(x.b, 0.75)
}

# Tracer les fonctions de répartition empiriques des quartiles
plot(ecdf(T.boot.Q1), main='Parametric Bootstrap (Gaussian): Sample Q1, Q2, and Q3', col='green', xlim=c(15,26))
lines(ecdf(T.boot.Q2), col='black')
lines(ecdf(T.boot.Q3), col='red')

# Tracer les lignes verticales pour les quartiles de l'échantillon original
abline(v=c(Q1,Q2,Q3), col=c('green','black','red'), lty=2)


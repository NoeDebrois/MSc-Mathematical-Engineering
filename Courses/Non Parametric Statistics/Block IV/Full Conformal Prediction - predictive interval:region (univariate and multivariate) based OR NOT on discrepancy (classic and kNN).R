# In this file :
# - Univariate predictive interval
# - multivariate predictive region
# - univariate predictive interval NOT based on discrepancy (kNN)
# - univariate predictive region NOT based on discrepancy (kNN)

##################################################
### Example of Univariate Predictive Intervals ###
##################################################
# Charger la bibliothèque resample qui contient les fonctions nécessaires
library(resample)

# Lire les données depuis un fichier texte
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block IV/04 Conformal Prediction-2/parziali.txt')
P1 <- X$PI

# Analyser la distribution des données
hist(P1, breaks=15)  # Tracer l'histogramme des données
shapiro.test(P1)     # Effectuer le test de normalité de Shapiro sur les données :
# Hypothèse nulle (H0) : Les données suivent une distribution normale.
# Hypothèse alternative (H1) : Les données ne suivent pas une distribution normale.

# Paramètres pour la construction de l'intervalle de prédiction conformal
alpha <- 0.1  # Niveau de confiance
x.obs <- P1    # Données observées

# Définir une grille de valeurs pour les nouvelles observations
x.new.grid <- seq(min(x.obs) - 0.5*diff(range(P1)), max(x.obs) + 0.5*diff(range(P1)), length = 1000) 

# Initialiser un vecteur pour stocker les p-values
p.value <- numeric(length(x.new.grid))

# Fonction de prédiction non conforme (NC)
NC <- function(z.aug, i){
  abs(z.aug[i] - mean(z.aug[-i]))  # sample mean as a predictor (i.e., T intervals)
  # Autres exemples de prédicteurs peuvent être utilisés en commentant/décommentant les lignes suivantes
  #abs(z.aug[i] - median(z.aug[-i]))   # a robust predictor
  #abs(z.aug[i] - (mean(z.aug[-i])+10)) # a biased predictor
  #abs(z.aug[i] - 18)                  # a deterministic predictor
  #abs(z.aug[i] - random.number)       # a fully random predictor
}

# Boucle pour calculer les p-values pour chaque valeur de la grille
for(k in 1:length(x.new.grid)) {
  # Ajouter la valeur actuelle de la grille aux données observées augmentées
  x.obs.aug <- c(x.obs, x.new.grid[k])
  
  # Initialiser un vecteur pour stocker les scores
  scores <- numeric(length(x.obs.aug))
  
  # Calculer les scores pour chaque observation augmentée
  for (i in 1:length(x.obs.aug)) {
    scores[i] <- NC(x.obs.aug, i)
  }
  
  # Calculer la p-value pour la valeur actuelle de la grille
  p.value[k] <- sum(scores >= scores[length(x.obs.aug)]) / length(x.obs.aug)
}

# Tracer les p-values
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))  # Ligne horizontale à 0 et 1
abline(h=alpha, col='red', lty=2)  # Ligne horizontale à la valeur alpha en rouge (niveau de confiance)
points(x.obs, numeric(length(x.obs)), pch=3)  # Points pour les données observées

# Calculer l'intervalle de prédiction
PI.grid <- x.new.grid[which(p.value >= alpha)]
PI <- c(min(PI.grid), max(PI.grid))
PI
abline(v = PI, col='red')  # Ligne verticale pour l'intervalle de prédiction en rouge
points(PI.grid, numeric(length(PI.grid)), pch=16, col='red')  # Points pour les valeurs de la grille incluses dans l'intervalle


##################################################
### Example of Multivariate Predictive Regions ###
##################################################
# Lire les données à partir d'un fichier
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block IV/04 Conformal Prediction-2/parziali.txt')

# Tracer le nuage de points pour les variables (PI, PII)
plot(X, asp=1)

# Initialiser les paramètres et les grilles
alpha <- 0.1
x.obs <- X
x1.new.grid <- seq(min(x.obs[,1]) - 0.25*diff(range(x.obs[,1])), max(x.obs[,1]) + 0.25*diff(range(x.obs[,1])), length = 20)
x2.new.grid <- seq(min(x.obs[,2]) - 0.25*diff(range(x.obs[,2])), max(x.obs[,2]) + 0.25*diff(range(x.obs[,2])), length = 20)
p.value <- matrix(nrow = length(x1.new.grid), ncol = length(x2.new.grid))

# Définir la fonction de prédiction non conforme (NC)
NC <- function(z.aug, i){
  # Différentes options pour le calcul de la distance ou du score non conforme
  #sum((z.aug[i,] - colMeans(z.aug[-i,]))^2) # Euclidean distance
  #as.numeric( as.matrix(z.aug[i,] - colMeans(z.aug[-i,])) %*% 
  #              solve(cov(z.aug[-i,])) %*% 
  #              as.matrix(t(z.aug[i,] - colMeans(z.aug[-i,]))) ) # Mahalanobis Distance
  max((z.aug[i,] - colMeans(z.aug[-i,]))^2/colVars(z.aug[-i,])) # Standardized Maxi-distance
  #abs( z.aug[i,2] - sum(coefficients(lm(z.aug[-i,2]  ~ z.aug[-i,1]))*c(1, z.aug[i,1]))) # Regression
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
points(X, pch=16)
contour(x1.new.grid, x2.new.grid, p.value, levels = alpha, add=T)


##################################################
### Example of Univariate Predictive sets      ###
### not based on discrepancy (e.g., KNN)       ###
##################################################
# Lire les données à partir d'un fichier
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block IV/04 Conformal Prediction-2/parziali.txt')
P1 <- X$PI

# Tracer l'histogramme de la variable P1
hist(P1, breaks=15)

# Effectuer le test de normalité de Shapiro sur P1
shapiro.test(P1)

# Initialiser les paramètres pour la méthode conformal
alpha <- 0.1
x.obs <- P1
x.new.grid <- seq(min(x.obs) - 0.5*diff(range(P1)), max(x.obs) + 0.5*diff(range(P1)), length = 1000) 
p.value <- numeric(length(x.new.grid))
K <- 5  # Nombre de voisins pour le calcul de la distance
NC <- function(z.aug, i){
  distances2 <- ( as.matrix(dist(z.aug))[i,-i] )^2 # semble calculer les carrés des distances entre la i-ème 
  # observation augmentée (z.aug[i,]) et toutes les autres observations dans le jeu de données augmenté z.aug,
  # à l'exception de l'observation i-ème elle-même.
  
  mean(sort(distances2)[1:K]) # Moyenne des distances pour la liaison moyenne (average linkage)
  #min(sort(distances2)[1:K])  # Distance minimale pour la liaison simple (single linkage)
  #max(sort(distances2)[1:K])  # Distance maximale pour la liaison complète (complete linkage)
}

# Boucle pour calculer les p-values
for(k in 1:length(x.new.grid)) {
  x.obs.aug <- c(x.obs, x.new.grid[k])
  scores <- numeric(length(x.obs.aug))
  for (i in 1:length(x.obs.aug)) {
    scores[i] <- NC(x.obs.aug, i)
  }
  p.value[k] <- sum(scores >= scores[length(x.obs.aug)])/(length(x.obs.aug))
}

# Tracer les p-values
plot(x.new.grid, p.value, type='l', ylim=c(0,1))
abline(h=c(0,1))
abline(h=alpha, col='red', lty=2)
points(x.obs, numeric(length(x.obs)), pch=3)

# Calculer l'ensemble de prédiction
PI.grid <- x.new.grid[which(p.value >= alpha)]
points(PI.grid, numeric(length(PI.grid)), pch=16, col='red')


##################################################
### Example of Multivariate Predictive sets    ###
### not based on discrepancy (e.g., KNN)       ###
##################################################
# Lire les données à partir d'un fichier
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block IV/04 Conformal Prediction-2/parziali.txt')

# Tracer le nuage de points pour les variables (PI, PII)
plot(X, asp=1)

# Initialiser les paramètres pour la méthode conformale
alpha <- 0.5
x.obs <- X
x1.new.grid <- seq(min(x.obs[,1]) - 0.25*diff(range(x.obs[,1])), max(x.obs[,1]) + 0.25*diff(range(x.obs[,1])), length = 20)
x2.new.grid <- seq(min(x.obs[,2]) - 0.25*diff(range(x.obs[,2])), max(x.obs[,2]) + 0.25*diff(range(x.obs[,2])), length = 20)
p.value <- matrix(nrow = length(x1.new.grid), ncol = length(x2.new.grid))
K <- 5

# Définir la fonction de distance non conforme (NC)
NC <- function(z.aug, i){
  distances2 <- ( as.matrix(dist(z.aug))[i,-i] )^2
  mean(sort(distances2)[1:K]) # Liaison moyenne (average linkage)
  #min(sort(distances2)[1:K])  # Liaison simple (single linkage)
  #max(sort(distances2)[1:K])  # Liaison complète (complete linkage)
}

# Boucles pour calculer les p-values
for(k in 1:length(x1.new.grid)) {
  for(h in 1:length(x2.new.grid)) {
    # Ajouter une nouvelle observation aux données existantes
    x.obs.aug <- rbind(x.obs, c(x1.new.grid[k], x2.new.grid[h]))
    
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
points(X, pch=16)
contour(x1.new.grid, x2.new.grid, p.value, levels = alpha, add=T)


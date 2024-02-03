# Ce script R démontre l'utilisation de la méthodologie conformale pour construire des intervalles
# prédictifs univariés et des régions prédictives multivariées. Voici une explication du code :

## Univarié :
# - Charger la bibliothèque et les données :
#     La bibliothèque "resample" est chargée.
#     Les données sont lues à partir d'un fichier pour les analyses univariées.
# - Tracer l'histogramme et effectuer un test de normalité :
#     Un histogramme de la variable P1 est tracé.
#     Un test de normalité de Shapiro-Wilk est effectué sur P1.
# - Initialiser les paramètres pour la méthode univariée :
#     Le niveau de confiance alpha est défini.
#     Un ensemble d'entraînement est sélectionné à partir des données.
# - Définir la fonction de distance non conforme (NC) :
#     La fonction NC mesure la distance non conforme entre une nouvelle observation et
#     la moyenne de l'ensemble d'entraînement.
# - Boucle pour calculer les p-values et tracer l'intervalle de prédiction :
#     Une boucle itère sur une grille de nouvelles valeurs.
#     Les p-values sont calculées en utilisant la distance non conforme.
#     Les p-values sont tracées, et un intervalle de prédiction est calculé et également tracé.
    
## Multivarié :
# - Analyse multivariée :
#     Les mêmes étapes sont répétées pour une analyse multivariée avec les variables (PI, PII).
# - Définir la fonction de distance non conforme (NC) pour les variables multivariées :
#     La fonction NC est adaptée pour gérer les variables (PI, PII).
# - Boucles pour calculer les p-values et tracer la région de prédiction :
#     Des boucles itèrent sur les grilles de nouvelles valeurs pour les variables (PI, PII).
#     Les p-values sont calculées et utilisées pour tracer une région de prédiction multivariée.

##################################################
### Exemple d'Intervalles Prédictifs Univariés  ###
##################################################
library(resample)

# Charger la librairie 'resample' pour les fonctions conformes

# Lire les données à partir d'un fichier
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block IV/04 Conformal Prediction-2/parziali.txt')
P1 <- X$PI

# Tracer l'histogramme de la variable P1
hist(P1, breaks=15)

# Effectuer le test de normalité de Shapiro sur P1
shapiro.test(P1)

# Initialiser les paramètres pour la méthode conformale univariée
alpha <- 0.1
n <- length(P1)
training.prop <- 0.75
set.seed(240279)
#set.seed(240278)
training.id   <- sample(1:n, ceiling(n*training.prop), replace = F) 

training.set  <- P1[training.id]   # Ensemble d'entraînement approprié
x.obs         <- P1[-training.id]  # Ensemble de calibration

# Créer une grille de nouvelles valeurs pour P1
x.new.grid <- seq(min(x.obs) - 0.5*diff(range(P1)), max(x.obs) + 0.5*diff(range(P1)), length = 1000) 
p.value <- numeric(length(x.new.grid))

# Calculer la moyenne de l'ensemble d'entraînement
training.mean <- mean(training.set)

# Définir la fonction de distance non conforme (NC)
NC <- function(z.aug, i){
  abs(z.aug[i] - training.mean)  # Moyenne de l'ensemble d'entraînement comme prédicteur (intervalle T)
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
abline(v=mean(training.set), lwd=2)
points(x.obs, numeric(length(x.obs)), pch=3)

# Calculer et tracer l'intervalle de prédiction
PI.grid <- x.new.grid[which(p.value >= alpha)]
PI <- c(min(PI.grid), max(PI.grid))
PI
abline(v = PI, col='red')
points(PI.grid, numeric(length(PI.grid)), pch=16, col='red')


##################################################
### Exemple de Régions Prédictives Multivariées ###
##################################################

# Charger à nouveau les données
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block IV/04 Conformal Prediction-2/parziali.txt')

# Tracer le nuage de points pour les variables (PI, PII)
plot(X, asp=1)

# Initialiser les paramètres pour la méthode conformale multivariée
alpha <- 0.1
n <- dim(X)[1]
training.prop <- 0.75
set.seed(240279)
#set.seed(240278)
training.id   <- sample(1:n, ceiling(n*training.prop), replace = F) 

training.set  <- X[training.id,]    # Ensemble d'entraînement approprié
x.obs         <- X[-training.id,]  # Ensemble de calibration

# Créer des grilles de nouvelles valeurs pour les variables (PI, PII)
x1.new.grid <- seq(min(x.obs[,1]) - 0.25*diff(range(x.obs[,1])), max(x.obs[,1]) + 0.25*diff(range(x.obs[,1])), length = 30)
x2.new.grid <- seq(min(x.obs[,2]) - 0.25*diff(range(x.obs[,2])), max(x.obs[,2]) + 0.25*diff(range(x.obs[,2])), length = 30)
p.value <- matrix(nrow = length(x1.new.grid), ncol = length(x2.new.grid))

# Calculer la moyenne et la variance de l'ensemble d'entraînement
training.mean <- colMeans(training.set)
training.var <- colVars(training.set)

# Définir la fonction de distance non conforme (NC) pour les variables (PI, PII)
NC <- function(z.aug, i){
  max((z.aug[i,] - training.mean)^2/training.var)
}

# Boucles pour calculer les p-values
for(k in 1:length(x1.new.grid)) {
  for(h in 1:length(x2.new.grid)) {
    x.obs.aug <- rbind(x.obs, c(x1.new.grid[k],x2.new.grid[h]))
    scores <- numeric(dim(x.obs.aug)[1])
    for (i in 1:dim(x.obs.aug)[1]) {
      scores[i] <- NC(x.obs.aug, i)
    }
    p.value[k,h] <- sum(scores >= scores[dim(x.obs.aug)[1]])/(dim(x.obs.aug)[1])
    print(c(k,h))
  }
}

# Tracer les p-values et la région de prédiction
image(x1.new.grid, x2.new.grid, p.value, zlim=c(0,1), asp=1)
points(x.obs, pch=16)
contour(x1.new.grid, x2.new.grid, p.value, levels = alpha, add=T)
points(training.mean[1], training.mean[2], pch=3, cex=4)

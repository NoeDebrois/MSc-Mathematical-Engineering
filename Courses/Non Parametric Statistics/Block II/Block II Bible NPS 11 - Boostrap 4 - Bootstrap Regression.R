################################################################################
### Bootstrap Regression #######################################################
# Le code effectue une analyse de régression linéaire avec des méthodes ########
# de bootstrap pour estimer l'intervalle de confiance des paramètres de ########
# la régression, ainsi que l'intervalle de confiance de la moyenne #############
# conditionnelle à une valeur spécifique du prédicteur #########################
################################################################################
# Effacer toutes les variables de l'environnement
rm(list=ls())

# Importer les données depuis un fichier
grades <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/03 Bootstrap/parziali.txt', header=T)

# Variables de réponse et de prédicteur
response  <- grades[,'PII']
regressor <- grades[,'PI']

# Tracer les données
plot(regressor, response, asp=1)
grid()
abline(0,1, lty=3)

# Ajuster un modèle linéaire
fm <- lm(response ~ regressor)

# Tracer la droite de régression et les valeurs ajustées
abline(coefficients(fm), col='red')
points(regressor, fitted(fm), col='red', pch=16)

### Inférence paramétrique classique
summary(fm)

################################################################################
### Bootstrap Inference ########################################################
################################################################################
# Calculer les résidus et les valeurs ajustées
fitted.obs <- fitted(fm)
res.obs    <- residuals(fm)

# Coefficients estimés du modèle original
b0.obs <- coefficients(fm)[1] # ordonnée à l'origine de la droite de régression
b1.obs <- coefficients(fm)[2] # coefficient directeur de la droite de régression

# Fixer la graine aléatoire pour la reproductibilité
set.seed(24021979)
# Nombre d'itérations Bootstrap
B <- 10000
# Initialiser des vecteurs pour stocker les estimations Bootstrap des coefficients
T.boot.b0 <- numeric(B)
T.boot.b1 <- numeric(B)

# Boucle sur les itérations Bootstrap
for(b in 1:B)
{
  # Ré-échantillonner les résidus avec remplacement
  response.b <- fitted.obs + sample(res.obs, replace = TRUE)
  # Ajuster un modèle linéaire sur les données Bootstrap
  fm.b <- lm(response.b ~ regressor)
  # Stocker les coefficients estimés Bootstrap
  T.boot.b0[b] <- coefficients(fm.b)[1]
  T.boot.b1[b] <- coefficients(fm.b)[2]
}

# Diviser la fenêtre graphique en 2 lignes et 1 colonne
par(mfrow=c(2,1))
# Tracer les fonctions de répartition empiriques des coefficients Bootstrap
plot(ecdf(T.boot.b0), main='Intercept')
abline(v=b0.obs, lty=2)
plot(ecdf(T.boot.b1), main='Slope', col='red')
abline(v=b1.obs, lty=2, col='red')

# Estimer les écarts-types des coefficients Bootstrap (sans calcul analytique)
sd(T.boot.b0)
sd(T.boot.b1)
cov(T.boot.b0, T.boot.b1)

# Pour la construction des intervalles de confiance, voir le fichier 
# Bootstrap 2 sur les intervalles de confiance pour avoir plus de détails :

################################################################################
### Intervalle de confiance Reverse Percentile pour l'ordonnée à l'origine #####
################################################################################
# Diviser la fenêtre graphique en 3 lignes et 1 colonne
par(mfrow=c(3,1))

# Niveau de confiance
alpha <- 0.05

# Calculer les quantiles pour l'intercept
right.quantile.b0 <- quantile(T.boot.b0, 1 - alpha/2)
left.quantile.b0  <- quantile(T.boot.b0, alpha/2)

# Valeurs observées
b0.obs
right.quantile.b0 - b0.obs
left.quantile.b0  - b0.obs

# Calculer l'intervalle de confiance Reverse Percentile
CI.RP.b0 <- c(b0.obs - (right.quantile.b0 - b0.obs), b0.obs - (left.quantile.b0 - b0.obs))
CI.RP.b0

# Tracer la fonction de répartition empirique de l'intercept avec l'intervalle de confiance
plot(ecdf(T.boot.b0), main='Intercept')
abline(v = b0.obs, lty=2)
abline(v = CI.RP.b0)

################################################################################
### Intervalle de confiance Reverse Percentile pour le coefficient directeur ###
################################################################################
alpha <- 0.05

# Calculer les quantiles pour la pente
right.quantile.b1 <- quantile(T.boot.b1, 1 - alpha/2)
left.quantile.b1  <- quantile(T.boot.b1, alpha/2)

# Valeurs observées
b1.obs
right.quantile.b1 - b1.obs
left.quantile.b1  - b1.obs

# Calculer l'intervalle de confiance Reverse Percentile
CI.RP.b1 <- c(b1.obs - (right.quantile.b1 - b1.obs), b1.obs - (left.quantile.b1 - b1.obs))
CI.RP.b1

# Tracer la fonction de répartition empirique de la pente avec l'intervalle de confiance
plot(ecdf(T.boot.b1), main='Slope', col='red')
abline(v = b1.obs, lty=2, col='red')
abline(v = CI.RP.b1, col='red')

################################################################################
### Intervalles de confiance Reverse Percentile ################################
### pour la moyenne conditionnelle à x0 = 24 ###################################
################################################################################
alpha <- 0.05

# Valeur spécifique du prédicteur
x0 <- 24
# Moyenne conditionnelle observée à x0 = 24
mean.x0.obs <- b0.obs + b1.obs * x0 # on remplace x par x0=24 et en moyenne epsilon = 0 donc c'est facile...
# Estimations Bootstrap de la moyenne conditionnelle à x0 = 24
T.boot.mean.x0 <- T.boot.b0 + T.boot.b1 * x0 # idem pour le bootstrap

# Calculer les quantiles pour la moyenne conditionnelle
right.quantile.mean.x0 <- quantile(T.boot.mean.x0, 1 - alpha/2)
left.quantile.mean.x0  <- quantile(T.boot.mean.x0, alpha/2)

# Valeurs observées
mean.x0.obs
right.quantile.mean.x0 - mean.x0.obs
left.quantile.mean.x0  - mean.x0.obs

# Calculer l'intervalle de confiance Reverse Percentile pour la moyenne conditionnelle
CI.RP.mean.x0 <- c(mean.x0.obs - (right.quantile.mean.x0 - mean.x0.obs), mean.x0.obs - (left.quantile.mean.x0 - mean.x0.obs))
CI.RP.mean.x0

# Tracer la fonction de répartition empirique de la moyenne conditionnelle avec l'intervalle de confiance
plot(ecdf(T.boot.mean.x0), main='Conditional mean at x0 = 24', col='green')
abline(v = mean.x0.obs, lty=2, col='green')
abline(v = CI.RP.mean.x0, col='green')


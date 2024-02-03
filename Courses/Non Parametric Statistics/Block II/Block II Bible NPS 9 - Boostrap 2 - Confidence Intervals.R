################################################################################
### Bootstrap Confidence Intervals #############################################
################################################################################
# Ce code en R concerne la construction d'intervalles de confiance à l'aide ####
# de la méthode de bootstrap pour deux estimations différentes (95e percentile #
# et probabilité d'obtenir une note supérieure ou égale à 18). #################
# En résumé, ce code importe des données, effectue des calculs de percentile ###
# et de probabilité, puis utilise la méthode de bootstrap pour estimer les #####
# distributions bootstrap de deux estimateurs. Ensuite, il construit des #######
# intervalles de confiance en utilisant différentes méthodes (REVERSE ##########
# PERCENTILES, PERCENTILES, et T-BOOTSTRAP *) et les trace sur des graphiques ##
# avec les fonctions de répartition empiriques des distributions bootstrap. ####
# L'objectif est de visualiser les intervalles de confiance pour évaluer la ####
# variabilité et la précision des estimations. #################################
# * : Voir pages 3 et 4 du cours sur le Boostrap pour les 3 types de CI. #######
################################################################################
# Supprimer toutes les variables de l'environnement
rm(list=ls())

# Importer les données depuis un fichier
grades <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/03 Bootstrap/parziali.txt', header=T)
grades.PI <- grades[,'PI']
boxplot(grades.PI)

# Calculer le 95e percentile et la probabilité d'obtenir une note supérieure ou égale à 18
Q95 <- quantile(grades.PI, 0.95)
P18 <- sum(grades.PI >= 18)/length(grades.PI)

# Calculer la distribution bootstrap des estimateurs
set.seed(24021979)
B <- 10000
x.obs <- grades.PI
T.boot.Q95 <- numeric(B)
T.boot.P18 <- numeric(B)

for(b in 1:B)
{
  x.b <- sample(x.obs, replace = TRUE)
  T.boot.Q95[b] <- quantile(x.b, 0.95)
  T.boot.P18[b] <- sum(x.b >= 18)/length(x.b)
}

# Tracer les fonctions de répartition empiriques des distributions bootstrap
par(mfrow=c(1,2))
plot(ecdf(T.boot.Q95), main='95th Percentile')
abline(v=Q95, lty=2)
plot(ecdf(T.boot.P18), main='P(grade >= 18)')
abline(v=P18, lty=2)

################################################################################
### Intervalle de confiance ####################################################
################################################################################
# Niveau de confiance
alpha <- 0.05

# Choix de l'estimateur et de la distribution bootstrap correspondante
T.boot <- T.boot.Q95
T.obs  <- Q95

# Tracer la fonction de répartition empirique de la distribution bootstrap
plot(ecdf(T.boot))
abline(v=T.obs, lty=2)

################################################################################
### Intervalle de confiance par la méthode des inverse percentiles #############
### (Reverse Percentile Intervals) #############################################
################################################################################
right.quantile <- quantile(T.boot, 1 - alpha/2)
left.quantile  <- quantile(T.boot, alpha/2)

# Valeurs
T.obs
right.quantile - T.obs
left.quantile  - T.obs

# Calculer l'intervalle
CI.RP <- c(T.obs - (right.quantile - T.obs), T.obs - (left.quantile - T.obs))
CI.RP

# Tracer l'intervalle sur le graphique
abline(v = CI.RP, col='red')

################################################################################
### Intervalle de confiance par la méthode des percentiles #####################
### (Percentile Intervals) #####################################################
################################################################################
right.quantile <- quantile(T.boot, 1 - alpha/2)
left.quantile  <- quantile(T.boot, alpha/2)

# Valeurs
CI.P <- c(left.quantile, right.quantile)
CI.P

# Tracer l'intervalle sur le graphique
abline(v = CI.P)

################################################################################
### Intervalle de confiance t-Bootstrap ########################################
### (Bootstrap T-Intervals) ####################################################
################################################################################
set.seed(24021979)
B <- 10000
x.obs <- grades.PI
T.boot.stud.P18 <- numeric(B)

# Calculer la distribution t-bootstrap de l'estimateur
for(b in 1:B)
{
  x.b <- sample(x.obs, replace = TRUE)
  P18.boot <- sum(x.b >= 18)/length(x.b)
  T.boot.stud.P18[b] <- ( P18.boot - P18 ) / sqrt(P18.boot * (1 - P18.boot) / length(x.obs))
}

# Valeurs
right.t.quantile <- quantile(T.boot.stud.P18, 1 - alpha/2)
left.t.quantile  <- quantile(T.boot.stud.P18, alpha/2)
right.t.quantile
left.t.quantile

# Calculer l'intervalle
CI.t <- c(P18 - right.t.quantile * sqrt(P18 * (1 - P18) / length(x.obs)), P18 - left.t.quantile  * sqrt(P18 * (1 - P18) / length(x.obs)))
CI.t

# Tracer les fonctions de répartition empiriques des distributions bootstrap
abline(v=T.obs, lty=2)
abline(v = CI.t, col='green')


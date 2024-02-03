################################################################################
### Real world vs Boostrap World ###############################################
################################################################################
################################################################################
### Real world vs Bootstrap World ###############################################
################################################################################
# Supprimer toutes les variables de l'environnement
rm(list=ls()) # clear all

################################################################################
### Simulation : Loi des grands nombres ########################################
### Cette section utilise la loi des grands nombres pour illustrer comment #####
### la distribution empirique converge vers la distribution théorique ##########
### lorsque la taille de l'échantillon augmente ################################
################################################################################
# Fixer la graine pour la reproductibilité
set.seed(24021979)

# Tailles d'échantillons à considérer
n.grid  <- 2^(1:8) # de 2 à 2^8

# Diviser la fenêtre graphique en 2 lignes et 4 colonnes
par(mfrow = c(2,4))

# Boucle sur différentes tailles d'échantillons
for(n in n.grid)
{
  # Générer un échantillon aléatoire
  x <- runif(n, 0, 4)
  
  # Tracer la fonction de répartition empirique
  plot(ecdf(x), main=paste('n =', n), xlim = c(0,4))
  
  # Tracer la fonction de répartition de la distribution uniforme
  lines(seq(-1, 5, by=0.1), punif(seq(-1, 5, by=0.1),0,4), type='l', col='red', lty=2)
}

################################################################################
### Distribution réelle et distribution Bootstrap de la moyenne ################
### Ces sections génèrent des échantillons aléatoires, calculent la ############
### statistique d'intérêt, la moyenne, pour la distribution réelle et ##########
### utilisent la méthode de bootstrap pour estimer la distribution bootstrap ###
### de ces statistiques ########################################################
################################################################################
set.seed(24021979)
n.grid  <- 2^(1:8) # Crée une séquence de tailles d'échantillons en utilisant des puissances de 2, allant de 2 à 2^8 (2, 4, 8, ..., 256)
M <- 10000 # Définit le nombre de répétitions pour la simulation de la distribution REELLE. Chaque répétition générera un échantillon de taille n et calculera la moyenne
B <- 1000 # Définit le nombre de répétitions pour la méthode de BOOTSTRAP. Chaque répétition effectuera un échantillonnage avec remplacement à partir de l'échantillon observé pour estimer la distribution bootstrap de la moyenne
x.obs.long <- runif(max(n.grid), 0, 4) # Génère un échantillon long (taille maximale parmi les tailles d'échantillons définies) à partir d'une distribution uniforme entre 0 et 4
# Cet échantillon "x.obs.long" sera utilisé comme population à partir de laquelle les échantillons seront tirés pour la méthode bootstrap

par(mfrow = c(2,4))

# Boucle sur différentes tailles d'échantillons
for(n in n.grid) # itère sur différentes tailles d'échantillons
{
  ### Distribution réelle de la moyenne ########################################
  T.real <- numeric(M) # Initialise un vecteur vide pour stocker les moyennes des échantillons de la distribution réelle
  for(m in 1:M) # génère M échantillons de taille n à partir de la distribution uniforme entre 0 et 4 et calcule la moyenne de chaque échantillon, stockant ces moyennes dans le vecteur T.real
  {
    x <- runif(n, 0, 4)
    T.real[m] <- mean(x)
  }
  
  # Échantillon observé pour la méthode Bootstrap
  x.obs <- x.obs.long[1:n] # Sélectionne les n premiers éléments de l'échantillon long pour être utilisés comme échantillon observé dans la méthode bootstrap
  
  ### Distribution Bootstrap de la moyenne #####################################
  T.boot <- numeric(B) # Initialise un vecteur vide pour stocker les moyennes bootstrap
  for(b in 1:B) # Effectue la méthode de bootstrap : elle effectue un échantillonnage avec remplacement à partir de l'échantillon observé (x.obs), calcule la moyenne de chaque échantillon bootstrap et stocke ces moyennes dans le vecteur T.boot
  {
    x.b <- sample(x.obs, replace = TRUE)
    T.boot[b] <- mean(x.b)
  }
  
  # Tracer la fonction de répartition empirique avec la commande "ecdf(...)" :
  plot(ecdf(T.real), main=paste('Sample mean: n =', n), col='red') # Tracer la fonction de répartition empirique des moyennes réelles avec un titre indiquant la taille de l'échantillon
  lines(ecdf(T.boot), col='black') # Ajoute la fonction de répartition empirique des moyennes bootstrap au graphique
}

################################################################################
### Les sections suivantes suivent une structure similaire pour l'écart-type ###
### et la médiane ##############################################################
################################################################################
# Real distribution and Bootstrap distribution of the Sample Standard Deviation
set.seed(24021979)
n.grid  <- 2^(1:8)
M <- 10000
B <- 1000
x.obs.long <- runif(max(n.grid),0,4)

par(mfrow = c(2,4))
for(n in n.grid)
{
  T.real <- numeric(M)
  for(m in 1:M)
  {
    x <- runif(n,0,4)
    T.real[m] <- sd(x)
  }
  
  x.obs <- x.obs.long[1:n]
  T.boot <- numeric(B)
  for(b in 1:B)
  {
    x.b <- sample(x.obs, replace = T)
    T.boot[b] <- sd(x.b)
  }
  
  plot(ecdf(T.real), main=paste('Sample SD: n =', n), xlim = c(0,3), col='red')
  lines(ecdf(T.boot), col='black')
}

# Real distribution and Bootstrap distribution of the Sample Median
set.seed(24021979)
n.grid  <- 2^(1:8)
M <- 10000
B <- 1000
x.obs.long <- runif(max(n.grid),0,4)

par(mfrow = c(2,4))
for(n in n.grid)
{
  T.real <- numeric(M)
  for(m in 1:M)
  {
    x <- runif(n,0,4)
    T.real[m] <- median(x)
  }
  
  x.obs <- x.obs.long[1:n]
  T.boot <- numeric(B)
  for(b in 1:B)
  {
    x.b <- sample(x.obs, replace = T)
    T.boot[b] <- median(x.b)
  }
  
  plot(ecdf(T.real), main=paste('Sample Median: n =', n), col='red')
  lines(ecdf(T.boot), col='black')
}

################################################################################
### Bootstrap Estimation du biais, de la variance et de la MSE #################
################################################################################
# Importer les données depuis un fichier
grades <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/03 Bootstrap/parziali.txt', header=TRUE)
grades.PI <- grades[,'PI']

# Tracer une boîte à moustaches pour les données
boxplot(grades.PI)

# Évaluer la "qualité statistique" de trois quartiles échantillonnés
Q1 <- quantile(grades.PI, 0.25)
Q2 <- quantile(grades.PI, 0.50)
Q3 <- quantile(grades.PI, 0.75)

# Calculer la distribution Bootstrap de trois quartiles échantillonnés
set.seed(24021979)
B <- 10000
x.obs <- grades.PI
T.boot.Q1 <- numeric(B)
T.boot.Q2 <- numeric(B)
T.boot.Q3 <- numeric(B)

# Boucle pour la méthode Bootstrap
for(b in 1:B)
{
  x.b <- sample(x.obs, replace = TRUE)
  T.boot.Q1[b] <- quantile(x.b, 0.25)
  T.boot.Q2[b] <- quantile(x.b, 0.50)
  T.boot.Q3[b] <- quantile(x.b, 0.75)
}

# Tracer la fonction de répartition empirique des quartiles Bootstrap
plot(ecdf(T.boot.Q1), main='Sample Q1, Q2, and Q3', col='green', xlim=range(c(T.boot.Q1, T.boot.Q2, T.boot.Q3)))
lines(ecdf(T.boot.Q2), col='black')
lines(ecdf(T.boot.Q3), col='red')

# Tracer une ligne verticale pour les quartiles réels
abline(v=c(Q1, Q2, Q3), col=c('green','black','red'), lty=2)

# Estimation Bootstrap de l'écart-type des trois quartiles échantillonnés
sd(T.boot.Q1)
sd(T.boot.Q2)
sd(T.boot.Q3)

# Estimation Bootstrap du biais des trois quartiles échantillonnés
mean(T.boot.Q1) - Q1
mean(T.boot.Q2) - Q2
mean(T.boot.Q3) - Q3

# Estimation Bootstrap de la racine carrée de la MSE des trois quartiles échantillonnés
sqrt(sd(T.boot.Q1)^2 + (mean(T.boot.Q1) - Q1)^2)
sqrt(sd(T.boot.Q2)^2 + (mean(T.boot.Q2) - Q2)^2)
sqrt(sd(T.boot.Q3)^2 + (mean(T.boot.Q3) - Q3)^2)


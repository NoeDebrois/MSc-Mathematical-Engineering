#########################################################
## NON-PARAMETRIC-STATISTICS - BLOCK II.1 - RANK TESTS ##
## Based on the course of Pr Simone Vantini #############
## Modified by me, MSc Student Noé Debrois ##############
## December 2023 ######################################## 
#########################################################
## UTILISATION : ########################################
## Run selected lines : shift+enter #####################
## Clear all plots : shift+option+back ##################
## Clear workbench : shift+back #########################
#########################################################
## ATTENTION : dans ce doc, #############################
## X correspond aux notes des élèves ; ##################
## G correspond à leur ID (odd or even : 1 or 0). #######
#########################################################
## REGLE POUR REJETER OU NON H0 : VOIR A LA FIN DU DOC ##
#########################################################

# DATA
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/01 Sign and Rank Tests/parziali.txt')
G <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/01 Sign and Rank Tests/matricola.txt')

#########################################################
### TWO SAMPLE INDEPENDENT RANK SUM TEST (MW U-TEST) ####
### TWO-SIDED MANN-WHITNEY U-TEST #######################
#########################################################
# H0: P(PI.odd > PI.even)  = 0.5 ########################
# H1: P(PI.odd > PI.even) != 0.5 ########################
#########################################################
PI1 <- X$PI[G$OE == 1]
PI2 <- X$PI[G$OE == 0]
PI   <- X$PI
n1 <- length(PI1)
n2 <- length(PI2)
n  <- length(PI)

# Visualisation de la structure de dépendance
diff <- expand.grid(PI1, PI2)[,2] - expand.grid(PI1, PI2)[,1]
image(1:n1, 1:n2, matrix(sign(diff), n1, n2))

# Calcul des rangs pour l'ensemble des PI
ranks.PI <- rank(PI)

# Calcul des rangs et des statistiques U pour chaque échantillon
R1 <- sum(ranks.PI[G$OE == 1])
U1 <- R1 - n1*(n1+1)/2  # Nombre de victoires du 1er échantillon

R2 <- sum(ranks.PI[G$OE == 0])
U2 <- R2 - n2*(n2+1)/2  # Nombre de victoires du 2e échantillon

n1*n2 # Nombre total de combats
U1 - n1*n2/2  # Déséquilibre par rapport à la moyenne sous H0
U2 - n1*n2/2  # Déséquilibre par rapport à la moyenne sous H0

# Calcul de la p-value par simulation Monte Carlo
set.seed(24021979)
B <- 100000
U1.sim <- numeric(B)
U2.sim <- numeric(B)
for (k in 1:B)
{
  ranks.temp <- sample(1:n)
  R1.temp <- sum(ranks.temp[1:n1])
  R2.temp <- sum(ranks.temp[(n1+1):(n1+n2)])
  U1.temp <- R1.temp - n1*(n1+1)/2
  U2.temp <- R2.temp - n2*(n2+1)/2
  U1.sim[k] <- U1.temp
  U2.sim[k] <- U2.temp
}

# Tracé des histogrammes des simulations
hist(U1.sim, breaks = 50)
abline(v = c(U1, U2), col='red')
abline(v = n1*n2/2, lwd=3)

hist(U2.sim, breaks = 50)
abline(v = c(U1, U2), col='red')
abline(v = n1*n2/2, lwd=3)

# Sélection du maximum entre U1 et U2
U.star <- max(U1, U2)

# Calcul de la p-value
p.value <- 2 * sum(U1.sim >= U.star)/B # 2* car two sided
p.value

# Transformation non linéaire monotone
PI.exp  <- -exp(-PI)
PI1.exp <- -exp(-PI1)
PI2.exp <- -exp(-PI2)

# Tracé des graphiques avant et après transformation
plot(PI, PI.exp, col=G$OE+1)
boxplot(PI)
boxplot(PI.exp)

# Calcul des rangs et des statistiques U après transformation
ranks.PI.exp <- rank(PI.exp)
R1.exp <- sum(ranks.PI.exp[G$OE == 1])
U1.exp <- R1.exp - n1*(n1+1)/2  

# Calcul de la p-value après transformation
U.star.exp <- n1*n2/2 + abs(U1.exp - n1*n2/2)
p.value.exp <- 2 * sum(U1.sim >= U.star.exp)/B

# Affichage des p-values
phrase <- "The p-value of this TWO SAMPLE TWO-SIDED INDEPENDENT RANK SUM MANN-WHITNEY U-TEST is : "
resultat <- paste0(phrase, p.value) # concaténation
cat(resultat) # affichage
phrase <- "The p-value of this TWO SAMPLE TWO-SIDED INDEPENDENT RANK SUM MANN-WHITNEY U-TEST, after monotonic transformation is : "
resultat <- paste0(phrase, p.value.exp) # concaténation
cat(resultat) # affichagep.value.exp
print("The MANN-WHITNEY U TEST is invariant WRT monotonic transformation of data, since the two p-values are the same.")

# Test de t pour comparaison
t.test(PI1,     PI2    )$p.value
t.test(PI1.exp, PI2.exp)$p.value
print("The CLASSICAL T-TEST is NOT invariant WRT monotonic transformation of data, since the two p-values are different.")


#########################################################################################
## REGLE POUR REJETER OU NON H0 #########################################################
#########################################################################################
# La règle générale pour accepter ou rejeter une hypothèse testée statistiquement repose sur la 
# comparaison de la valeur p (p-value) avec un niveau de signification prédéfini, généralement noté
# alpha (α). Voici la règle générale :
# 1. **Si la valeur p est inférieure à alpha (p < α) :** Rejeter l'hypothèse nulle.
#   - Cela suggère que les résultats de l'expérience sont statistiquement significatifs, et on a des
#     preuves pour rejeter l'hypothèse nulle en faveur de l'hypothèse alternative.
# 2. **Si la valeur p est supérieure ou égale à alpha (p ≥ α) :** Ne pas rejeter l'hypothèse nulle.
#   - Cela suggère que les résultats de l'expérience ne fournissent pas suffisamment de preuves pour
#     rejeter l'hypothèse nulle.
#
# En d'autres termes, le choix de la valeur alpha est crucial. 
# Il représente le seuil au-delà duquel vous considérez les résultats comme suffisamment extrêmes
# pour rejeter l'hypothèse nulle. Une valeur couramment utilisée pour alpha est 0.05, ce qui signifie
# que vous êtes prêt à rejeter l'hypothèse nulle si la probabilité d'observer les résultats obtenus 
# (ou quelque chose de plus extrême) sous l'hypothèse nulle est inférieure à 5%.
#
# Il est important de noter que la non-rejection de l'hypothèse nulle ne prouve pas que l'hypothèse 
# nulle est vraie; cela signifie simplement qu'il n'y a pas suffisamment de preuves pour la rejeter
# avec le niveau de signification choisi.
#
# En résumé, la décision de rejeter ou ne pas rejeter l'hypothèse nulle dépend de la comparaison de 
# la valeur p avec le niveau de signification alpha.
#########################################################################################
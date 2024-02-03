#########################################################
## NON-PARAMETRIC-STATISTICS - BLOCK II.1 - SIGN TESTS ##
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
## PI correspond à Y dans le cours ; ####################
## PII correspond à X dans le cours. ####################
#########################################################
## REGLE POUR REJETER OU NON H0 : VOIR A LA FIN DU DOC ##
#########################################################

# DATA :
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/01 Sign and Rank Tests/parziali.txt')

################################################
### TWO SAMPLE PAIRED RIGHT-SIDED SIGN TEST ####
################################################
# H0: P(PII > PI) = 0.5 (i.e.MED(PII-PI)=0)) ###
# H1: P(PII > PI) > 0.5 (i.e.MED(PII-PI)>0)) ###
################################################
differences <- X$PII - X$PI # Calcul des différences entre deux colonnes
boxplot(differences) # Création d'un boxplot pour visualiser les différences

n <- length(differences) # Nombre d'observations
signs <- sign(differences) # Calcul des signes des différences

W <- sum(signs==1) # Calcul de la statistique de test W (nombre de signes positifs)
#W <- (n+sum(signs))/2 # Deuxième manière de calculer W (cf cours)

plot(0:n, dbinom(0:n, n, 0.5)) # Tracé de la distribution binomiale
abline(v = W, col='red') # Ajout d'une ligne rouge à la position W
points(0:n, dbinom(0:n, n, 0.5), col= (0:n >= W) + 1, pch=16) # Points pour indiquer les probabilités

p.value <- 1 - pbinom(W-1, n, 0.5) # Calcul de la p-value (probabilité d'observer W signes positifs ou plus).
# Calcul de la pvalue en utilisant la fonction de distribution cumulative binomiale (cdf). 
# Cette p-value représente la probabilité d'observer W signes positifs ou plus dans
# une distribution binomiale avec une probabilité de succès de 0.5. En d'autres termes :
# pval = P(having >W HEADS in a coin-toss-game) = P(X_binom > W) = 1 - P(X_binom <= W) = 1 - cdf(X_binom)
# with X_binom following a binomial of parameter (n, p=0.5).

# Print of the result
phrase <- "The p-value of this TWO SAMPLE PAIRED RIGHT-SIDED SIGN TEST is : "
resultat <- paste0(phrase, p.value) # concaténation
cat(resultat) # affichage
################################################
################################################

################################################
### TWO SAMPLE PAIRED TWO-SIDED SIGN TEST ######
################################################
# H0: P(PII > PI)  = 0.5 (i.e.MED(PII-PI)=0)) ##
# H1: P(PII > PI) != 0.5 (i.e.MED(PII-PI)!=0)) #
################################################
differences <- X$PII - X$PI # Calcul des différences entre deux colonnes
boxplot(differences) # Création d'un boxplot pour visualiser les différences

n <- length(differences) # Nombre d'observations
signs <- sign(differences) # Calcul des signes des différences

# Calcul du nombre de signes positifs (W) et du maximum entre W et n-W (W.max) :
W <- sum(signs ==  1) # W : nombre de positifs
W.max <- max(W, n-W)  # n-W : nombre de négatifs

plot(0:n, dbinom(0:n, n, 0.5)) # Tracé de la distribution binomiale avec une probabilité de succès de 0.5
abline(v = c(W, n-W), col='red') # Ajout de lignes rouges aux positions W et n-W
points(0:n, dbinom(0:n, n, 0.5), col= (0:n >= max(W,n-W) | 0:n <= min(W,n-W)) + 1, pch=16) # Points pour indiquer les probabilités

# Calcul de la valeur p en utilisant une distribution binomiale à deux côtés
p.value <- 2*(1 - pbinom(W.max-1, n, 0.5) ) # *2 car on prend les deux côtés (TWO-SIDED) !
# Le facteur 2 est utilisé car il s'agit d'un test bilatéral (à deux côtés)
# P( B(n,0.5) >= 51 OR B(n,0.5) <= 8) = 2* P( B(n,0.5) >= 51 ) = 2*P( B(n,0.5) <= 8)

# Print of the result
phrase <- "The p-value of this TWO SAMPLE PAIRED TWO-SIDED SIGN TEST is : "
resultat <- paste0(phrase, p.value) # concaténation
cat(resultat) # affichage
################################################
################################################

################################################
### TWO SAMPLE PAIRED TWO-SIDED SIGN TEST ######
### WITH 1 BIG OUTLIER #########################
################################################
# H0: P(PII > PI)  = 0.5 (i.e.MED(PII-PI)=0)) ##
# H1: P(PII > PI) != 0.5 (i.e.MED(PII-PI)!=0)) #
################################################
differences <- X$PII - X$PI
differences.out <- differences
#differences.out[1] <- differences.out[1] - 1000 # sign consistent outlier
differences.out[1] <- differences.out[1] + 1000 # sign non-consistent outlier
boxplot(differences)
boxplot(differences.out)

n <- length(differences)
signs <- sign(differences)
signs.out <- sign(differences.out)

W <- sum(signs ==  1)
W.max <- max(W, n-W) 

W.out <- sum(signs.out ==  1)
W.max.out <- max(W.out, n-W.out) 

p.value <- 2*(1 - pbinom(W.max-1, n, 0.5) )
p.value.out <- 2*(1 - pbinom(W.max.out-1, n, 0.5) )

# Comparison with the parametric t-test
p.value
p.value.out

t.test(differences)$p.value
t.test(differences.out)$p.value

W.max
W.max.out

t.test(differences)$statistic
t.test(differences.out)$statistic

################################################
### ONE SAMPLE PAIRED TW0-SIDED SIGN TEST ######
### FOR THE MEDIAN #############################
################################################
# H0: P(PII > C0) = 0.5 (i.e. MED(PII)=C0)) ####
# H1: P(PII > C0) != 0.5 (i.e. MED(PII)!=C0)) ##
################################################
C_0 = 25

PII <- X$PII
differences <- PII - C_0
boxplot(differences)

n <- length(differences) # Nombre d'observations
signs <- sign(differences)

# Calcul du nombre de signes positifs (W) et du maximum entre W et n-W (W.max) :
W <- sum(signs ==  1) # W : nombre de positifs
W.max <- max(W, n-W)  # n-W : nombre de négatifs

plot(0:n, dbinom(0:n, n, 0.5)) # Tracé de la distribution binomiale avec une probabilité de succès de 0.5
abline(v = c(W, n-W), col='red') # Ajout de lignes rouges aux positions W et n-W
points(0:n, dbinom(0:n, n, 0.5), col= (0:n >= max(W,n-W) | 0:n <= min(W,n-W)) + 1, pch=16) # Points pour indiquer les probabilités

# Calcul de la valeur p en utilisant une distribution binomiale à deux côtés
p.value <- 2*(1 - pbinom(W.max-1, n, 0.5) ) # *2 car on prend les deux côtés (TWO-SIDED) !
# Le facteur 2 est utilisé car il s'agit d'un test bilatéral (à deux côtés)
# P( B(n,0.5) >= 51 OR B(n,0.5) <= 8) = 2* P( B(n,0.5) >= 51 ) = 2*P( B(n,0.5) <= 8)

# Print of the result
phrase <- "The p-value of this ONE SAMPLE PAIRED TWO-SIDED SIGN TEST FOR THE MEDIAN C_0 is : "
resultat <- paste0(phrase, p.value) # concaténation
cat(resultat) # affichage

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
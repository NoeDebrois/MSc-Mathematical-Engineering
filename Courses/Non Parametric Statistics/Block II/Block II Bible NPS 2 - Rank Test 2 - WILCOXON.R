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
## PII : correspond à X dans le cours ; #################
## PI : correspond à Y dans le cours. ###################
#########################################################
## REGLE POUR REJETER OU NON H0 : VOIR A LA FIN DU DOC ##
#########################################################

# DATA
X <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/01 Sign and Rank Tests/parziali.txt')
G <- read.table('/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/Block II/01 Sign and Rank Tests/matricola.txt')

#########################################################
## TWO SAMPLE (PAIRED) SIGNED-RANK TEST #################
## WILCOXON TEST (for odd students) #####################
## LEFT, RIGHT AND TWO-SIDED ############################
#########################################################
# H0: P(PII > PI | odd)  = 0.5 ##########################
# H1: P(PII > PI | odd) != 0.5 ##########################
#########################################################
# Extraction des différences entre les échantillons appariés
differences <- X$PII[G$OE == 1] - X$PI[G$OE == 1]

# Création d'un boxplot pour visualiser les différences
boxplot(differences)

# Nombre d'observations
n <- length(differences)

# Calcul des rangs pour les différences absolues
ranks <- rank(abs(differences))

# Calcul des statistiques de test W+
W.plus <- sum(ranks[differences > 0])
W.minus <- sum(ranks[differences < 0])

# Calcul de la statistique de test W
W <- W.plus - W.minus

# Simulation Monte Carlo de la p-value
set.seed(24021979)
B <- 1000000
W.sim <- numeric(B)
for (k in 1:B)
{
  ranks.temp <- sample(1:n)
  signs.temp <- 2*rbinom(n, 1, 0.5) - 1 # génère des -1 et des +1 (les signes) sous H0 (loi binomiale (n,0.5)) (car 2*B(0.5)-1 \in {-1, 1})
  W.temp <- sum(signs.temp * ranks.temp) # on peut calculer directement W (W+-W-) puisqu'on a déjà les +1 et -1 des signes... (ce qui revient au même que de faire les [positifs] - [les négatifs] avec les indicatrices...)
  W.sim[k] <- W.temp
}

# Tracé de l'histogramme des simulations avec une ligne rouge pour la statistique de test observée
hist(W.sim, xlim=c(-n*(n+1)/2, n*(n+1)/2), breaks = 50)
abline(v = W, col='red')
abline(v = 0, lwd=3)

# Test à deux côtés
p.value_two_sided <- 2 * sum(W.sim >= abs(W))/B

# Test unilatéral droit
p.value_right_sided <- sum(W.sim >= W)/B

# Test unilatéral gauche
p.value_left_sided <- sum(W.sim <= W)/B

# Affichage des p-values
p.value_two_sided
p.value_right_sided
p.value_left_sided

# Affichage des p-values
phrase <- "The p-value of this TWO SAMPLE (PAIRED) SIGNED-RANK WILCOXON TEST in TWO-SIDED MODE is : "
resultat <- paste0(phrase, p.value_two_sided) # concaténation
cat(resultat) # affichage
phrase <- "The p-value of this TWO SAMPLE (PAIRED) SIGNED-RANK WILCOXON TEST in RIGHT-SIDED MODE is : "
resultat <- paste0(phrase, p.value_right_sided) # concaténation
cat(resultat) # affichage
phrase <- "The p-value of this TWO SAMPLE (PAIRED) SIGNED-RANK WILCOXON TEST in LEFT-SIDED MODE is : "
resultat <- paste0(phrase, p.value_left_sided) # concaténation
cat(resultat) # affichage

# La différence entre un test statistique à deux côtés (bilatéral) et un test unilatéral droit ou gauche réside
# dans la manière dont l'hypothèse alternative est formulée et, par conséquent, dans la façon dont la région de
# rejet est déterminée.
# - Test statistique à deux côtés (bilateral) :
# Hypothèse nulle (H0) : Aucun effet ou différence.
# Hypothèse alternative (H1) : Il existe un effet ou une différence, mais la direction n'est pas spécifiée.
# Région de rejet : Les extrêmes des deux côtés de la distribution.
# Exemple : Une différence significative dans les deux directions, inférieure ou supérieure à zéro.
# - Test unilatéral droit (one-tailed, right) :
# Hypothèse nulle (H0) : Aucun effet ou différence.
# Hypothèse alternative (H1) : Il existe un effet ou une différence, et la direction est spécifiée comme étant vers la droite (positif).
# Région de rejet : Seulement du côté droit de la distribution.
# Exemple : On teste si quelque chose est significativement plus grand que zéro.
# - Test unilatéral gauche (one-tailed, left) :
# Hypothèse nulle (H0) : Aucun effet ou différence.
# Hypothèse alternative (H1) : Il existe un effet ou une différence, et la direction est spécifiée comme étant vers la gauche (négatif).
# Région de rejet : Seulement du côté gauche de la distribution.
# Exemple : On teste si quelque chose est significativement plus petit que zéro.

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

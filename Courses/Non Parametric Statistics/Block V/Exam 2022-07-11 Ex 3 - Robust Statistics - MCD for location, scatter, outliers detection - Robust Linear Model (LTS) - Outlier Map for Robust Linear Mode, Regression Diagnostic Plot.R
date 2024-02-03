################################################################################
################################## EXERCISE III ################################
################################################################################
# Dr. Andrea Qapos, Ph.D would like to be able to assess the length of his 
# barnacles from the cozyness of his seaside mansion. To do so he requires very 
# good quality glass. Luckily, a former collegue of his is Dr. Mateo De la 
# Fuente, Ph.D, a Mexican former academic statistician, who now turned to the 
# glass-making industry. To this extent, he collects measurements of the 
# presence of chemical constituents in 76 pieces of glass lenses he has produced
# : his aim is to evaluate how these compounds affect the refractive index (RI) 
# of the glass pieces. In details, the considered chemical constituents (unit 
# measurement weight percent in corresponding oxide) are : sodium oxide (Na20), 
# magnesium oxide (MgO), aluminum oxide (Al2O3), silcon oxide (SiO2) potassium 
# oxide (K2O) and calcium oxide (CaO). The resulting samples are contained in 
# the glass_3.rds file. Mateo De la Fuente wants to reduce his scraps to zero,
# he thus prefers to employ robust methods to analyse his data as some defective
# parts may have been produced. He therefore asks you to:

glass_3 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-07-11/glass_3.rds")

# 1. Compute the Minimum Covariance Determinant estimator for the chemical 
# constituents (i.e., all the variables in the dataset but the refractive index 
# RI). Consider 1000 subsets for initializing the algorithm and set alpha equal 
# to 0.5. Report the reweighed MCD estimates of location and scatter. Define a 
# vector ind_out_MCD_rw of row indexes identifying the glass samples that are 
# outliers according to the final MCD estimates and report it.

chem = glass_3[,-1]
# Get the number of rows in the dataset
N <- nrow(chem)

# Set a seed for reproducibility
set.seed(2022)

# Fit the Minimum Covariance Determinant (MCD) estimator
fit_MCD <- covMcd(x = chem, alpha = 0.5, nsamp = 1000)

# Report the reweighed MCD estimate of location
fit_MCD$center # reweighed i.e NOT the raw one...

# Report the reweighed MCD estimate of scatter (covariance matrix)
fit_MCD$cov # reweighed i.e NOT the raw one...

# Identify the glass samples that are outliers according to the final MCD estimates :
# Attention il existe plusieurs critères pour identifier les outliers ! (cf Ex 3 11/02/2022)
# Ici : certains ne sont pas sélectionnés par les 'best' mais collent quand même
# à la meilleure approximation : on ne veut donc pas les décrire comme outliers !
# ICI : "vrais outliers". A FAIRE SI PAS CLAIR !
ind_out_MCD_rw <- (1:N)[fit_MCD$mcd.wt==0]
ind_out_MCD_rw

# 2. Build a robust linear model to regress the refractive index (RI) on the 
# sodium oxide (Na20), calcium oxide (CaO) and magnesium oxide (MgO) using a 
# Least Trimmed Squares (LTS) approach, setting the hyperparameter α = 0.75.
# Report the table of robustly estimated coefficients. Plot the resulting 
# outlier map and report the row indexes of bad leverage points and vertical 
# outliers present in the dataset according to the diagnostic plot.

# Fit a robust linear model using Least Trimmed Squares (LTS)
fit_lts <- ltsReg(RI ~ Na2O + CaO + MgO, alpha = 0.75, data = glass_3)
summary(fit_lts)

fit_lts$coefficients # table of robustly estimated coefficients

plot(fit_lts, which="rdiag") # resulting outlier map
# 37 and 59 BAD LEVERAGE POINTS
# 55 VERTICAL OUTLIER

# 3. Simplify the robust linear model of the previous exercise regressing the 
# refractive index (RI) on the sodium oxide (Na20) only. Using a bootstrap 
# approach, provide the reverse percentile confidence intervals for the 
# corresponding mean value of the RI. Plot the result.

# Ajustement du modèle LTS pour la régression robuste de RI sur Na2O
fit_lts <- ltsReg(RI ~ Na2O, alpha = 0.75, data = glass_3)
fitted.obs <- fit_lts$fitted.values  # Valeurs ajustées
res.obs <- fit_lts$residuals  # Résidus du modèle

# Création d'une grille de valeurs de Na2O pour la prédiction
Na2O_grid = seq(range(glass_3$Na2O)[1],
                range(glass_3$Na2O)[2],
                length.out = 100)

# Prédiction des valeurs à partir du modèle LTS
preds = c(cbind(1, Na2O_grid) %*% fit_lts$coefficients)

# Tracé des points de données et de la ligne de régression LTS
with(glass_3,
     plot(Na2O, RI))
lines(Na2O_grid, preds, col="blue")

# Initialisation de la boucle bootstrap
set.seed(2022)
B <- 1000
preds_boot_container <- matrix(ncol = length(Na2O_grid), nrow = B)
pb = progress::progress_bar$new(total = B,
                                format = " Processing [:bar] :percent eta: :eta")
# Boucle bootstrap
for(b in 1:B) {
  # Nouvelles valeurs de RI générées par échantillonnage des résidus
  RI_boot <- fitted.obs + sample(res.obs, replace = T)
  
  # Réajustement du modèle LTS sur les nouvelles données bootstrap
  fit_lts_boot <- ltsReg(RI_boot ~ glass_3$Na2O, alpha = 0.75)
  
  # Prédiction des valeurs à partir du modèle bootstrap
  preds_boot = c(cbind(1, Na2O_grid) %*% fit_lts_boot$coefficients)
  preds_boot_container[b,] = preds_boot
  pb$tick()
}

# Calcul des intervalles de confiance inversés (reverse percentile)
alpha <- 0.05
right.quantile.preds <- apply(preds_boot_container, 2, quantile, probs = 1 - alpha/2)
left.quantile.preds  <- apply(preds_boot_container, 2, quantile, probs = alpha/2)

# Calcul des intervalles de confiance inversés pour les prédictions
CI.preds <- list(low = preds - (right.quantile.preds - preds),
                 up = preds - (left.quantile.preds - preds))

# Tracé des intervalles de confiance inversés sur le graphique
matlines(
  Na2O_grid ,
  cbind(CI.preds$up, CI.preds$low) ,
  lwd = 1,
  col = " blue",
  lty = 3
)












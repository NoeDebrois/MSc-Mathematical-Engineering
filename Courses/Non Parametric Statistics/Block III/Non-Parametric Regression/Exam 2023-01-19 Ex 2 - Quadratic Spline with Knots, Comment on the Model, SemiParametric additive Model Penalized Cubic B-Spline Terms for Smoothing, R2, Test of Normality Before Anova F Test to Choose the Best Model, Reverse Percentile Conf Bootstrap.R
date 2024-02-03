################################################################################
################################ EXERCICE II ###################################
################################################################################
library(splines)
library(mgcv)
df_2 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2023-01-19/df_2.rds")

# Dr. Matteus Fontansen is now interested in building some regression models to
# predict the penguins body mass (measured in grams) as a function of species, 
# flipper and bill lengths: the related data are contained in the df_2.rds file.
# He therefore asks you to:

# 1. Build a quadratic spline model to regress body mass as a function of the
# flipper length only. Set two knots at 190 and 210 mm for the explanatory 
# variable. Provide a plot of the regression line with standard errors for the 
# prediction, a table summarizing the coefficients and comment the results.
knots <- c(190,210)
model_quad_splines <- lm(body_mass_g ~ bs(flipper_length_mm, knots = knots, degree = 2), data = df_2)

new_data <- with(df_2, data.frame(flipper_length_mm = seq(range(flipper_length_mm)[1], range(flipper_length_mm)[2], by = 1)))
preds = predict(model_quad_splines, new_data,se=T)

se.bands=cbind(preds$fit + 2 * preds$se.fit, preds$fit - 2 * preds$se.fit)

with(df_2, plot(flipper_length_mm, body_mass_g, xlim=range(new_data$flipper_length_mm) ,cex =.5, col = "darkgrey"))
lines(new_data$flipper_length_mm, preds$fit, lwd =2, col = "blue")
matlines(new_data$flipper_length_mm, se.bands ,lwd =1, col =" blue",lty =3)

knots_pred=predict(model_quad_splines,list(flipper_length_mm=knots))
points(knots,knots_pred, col='blue',pch=19)
abline(v = knots, lty=3, col="blue")

summary(model_quad_splines) # Commentaires : 
# Ce tableau est une sortie de régression linéaire multiple qui utilise des 
# splines cubiques pour modéliser la relation entre la variable de réponse 
# (body_mass_g) et la variable explicative (flipper_length_mm). 
# Voici comment interpréter les différentes parties du tableau :
  
#  Call: Indique la formule du modèle. Dans ce cas, le modèle de régression 
# linéaire multiple est spécifié avec la fonction de base spline cubique (bs)
# pour la variable flipper_length_mm, avec un degré 2.

#  Residuals:
# -Min: La valeur minimale des résidus (différence entre les valeurs observées 
# et les valeurs prédites par le modèle).
# -1Q: Le premier quartile des résidus.
# -Median: La médiane des résidus.
# -3Q: Le troisième quartile des résidus.
# -Max: La valeur maximale des résidus.

#  Coefficients:
# -Intercept: Représente l'ordonnée à l'origine du modèle, c'est-à-dire la 
# valeur attendue de la variable de réponse lorsque toutes les variables 
# explicatives sont égales à zéro.
# -bs(flipper_length_mm, knots = knots, degree = 2)1 à 4: Ce sont les 
# coefficients associés aux termes de base spline cubique pour flipper_length_mm.
# Chaque terme représente une partie différente de la fonction spline.

#  Estimate: Les valeurs estimées des coefficients.
#  Std. Error: Les erreurs standard associées à chaque coefficient.
#  t value: La statistique de test t pour tester l'hypothèse nulle que le 
# coefficient est égal à zéro.
#  Pr(>|t|): La valeur p associée à la statistique de test t. Elle indique la 
# probabilité d'observer une statistique de test aussi extrême que celle observée,
# sous l'hypothèse nulle.
#  Signif. codes: Indique le niveau de signification des coefficients. Les 
# étoiles (*) sont souvent utilisées pour indiquer les niveaux de significativité,
# où plus d'étoiles indiquent une plus grande significativité.
#  Residual standard error: C'est l'estimation de l'écart-type des résidus. Il 
# mesure la dispersion des résidus autour de la ligne de régression.
#  Multiple R-squared / Adjusted R-squared: Ces mesures indiquent la proportion
# de la variance de la variable de réponse expliquée par le modèle. Adjusted 
# R-squared prend en compte le nombre de variables dans le modèle.
#  F-statistic / p-value: Le F-statistic teste l'hypothèse nulle selon laquelle 
# tous les coefficients de régression sont égaux à zéro (aucun effet global). 
# La p-value associée indique la probabilité d'observer un tel F-statistic sous
# l'hypothèse nulle.

# Dans l'ensemble, ce modèle semble expliquer une proportion significative de la
# variance de la variable de réponse, et certains des termes spline cubique pour
# flipper_length_mm sont statistiquement significatifs.

# 2. Build a semiparametric additive model for regressing body mass on the 
# flipper length, bill length and species, using penalized cubic b-spline terms
# for smoothing the continuous predictors. After having written in proper 
# mathematical terms the additive model you have estimated, report the adjusted
# R2 and the p-values of the tests.

# Our model :
# body massi = β0+f(flipper lengthi)+f(bill lengthi)+β1speciesi_Chinstrap+β2speciesi_Gentoo+εi, 
# i = 1, . . . , N

# Fit a GAM model
fit_gam <- gam(body_mass_g ~ s(flipper_length_mm, bs = "cr") + s(bill_length_mm, bs = "cr") + species, 
               data = df_2)

# Summarize the results of the GAM fit
table_fit_gam <- summary(fit_gam)

# Extract and print the R-squared value
(r_2_squared <- table_fit_gam$r.sq)

table_fit_gam$p.pv # Extracts the p-values associated with the estimated smooth terms in the GAM.
# A low p-value (typically below a chosen significance level, e.g., 0.05)
# suggests that the corresponding smooth term is statistically significant, and you may reject
# the null hypothesis that the term has no effect.

table_fit_gam$s.table # You can inspect s_table to obtain detailed information about each smooth term in the GAM.

# 3. Build a reduced version of the semiparametric additive model of the previous
# exercise by letting the covariate bill_length_mm enter linearly in the model specification.
# Use an appropriate Anova F test statistic to assess whether a smooth function 
# is needed for bill_length_mm, specifying the null and the alternative hypothesis 
# you are testing and report the resulting p-value. Comment on the results.
fit_gam_reduced <- gam(body_mass_g ~ s(flipper_length_mm, bs = "cr") + bill_length_mm + species, 
               data = df_2)

# Test of normality of the data :
shapiro.test(fit_gam$residuals)
shapiro.test(fit_gam_reduced$residuals)
# pval > 0.05 for both.
# Residuals for both models can be assumed to be normally distributed at significance level alpha=0.05,
# so I DO NOT NEED to use a permutational approach
# (I could have done it of course, but it is a lot of work...).

# So we wimply do an ANOVA F Test :
anova(fit_gam_reduced, fit_gam, test = "F")
# At significance level alpha=0.05, we cannot reject H0 ("there is no significant
# improvement between the two models"). So I can be satisfied with the reduced model.

# 4. By using a bootstrap approach, provide a reverse percentile confidence 
# interval for the bill_length_mm coefficient of the reduced semiparametric 
# additive model employed in the previous exercise.

# Obtain observed fitted values and residuals
fitted.obs <- fit_gam_reduced$fitted.values
res.obs <- fit_gam_reduced$residuals

summary(fit_gam_reduced)$p.table
bill_length_mm_obs <- summary(fit_gam_reduced)$p.table[2, 1]
# So, bill_length_mm_obs contains the observed value of the estimate for the 
# linear term of bill_length_mm_obs in the semi-parametric model.

# Bootstrap settings
B <- 1000  # Number of bootstrap samples
T.boot.bill_length_mm <- numeric(B)

# Bootstrap procedure
set.seed(2022)
for (b in 1:B) {
  # Generate a bootstrap sample by resampling residuals
  body_mass_boot <- fitted.obs + sample(res.obs, replace = TRUE)
  
  # Fit a semi-parametric model to the bootstrap sample
  fit_semi_boot <- gam(body_mass_boot ~ s(flipper_length_mm, bs = "cr") + bill_length_mm + species, 
                       data = df_2)
  fit_semi_boot_table <- summary(fit_semi_boot)
  
  # Record the parametric effect for wavenumber 300cm−1
  T.boot.bill_length_mm[b] = fit_semi_boot_table$p.table[2, 1]
}

# Calculate reverse percentile confidence interval
alpha <- 0.05
right.quantile <- quantile(T.boot.bill_length_mm, 1 - alpha/2)
left.quantile <- quantile(T.boot.bill_length_mm, alpha/2)
CI.RP.bill_length <- c(
  bill_length_mm_obs - (right.quantile - bill_length_mm_obs),
  bill_length_mm_obs,
  bill_length_mm_obs - (left.quantile - bill_length_mm_obs)
)
names(CI.RP.bill_length) <- c('lwr', 'pointwise', 'upr')

CI.RP.bill_length




























# Load the splines package
library(splines) # Q1
library(mgcv) # Q2
#############
## EX II : ##
#############

# SEE : Lab 12 - Generalized Additive Models (GAM), for more !!

# Read the RDS file
milk_samples_2 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/SAVE COURS/PREVIOUS EXAMS-20231223/Data/2022-01-21/milk_samples_2.Rds")

##

# 1) Build a natural cubic spline model to regress Milk pH on the milk absorbance at
# wavenumber 70 cmˆ{-1}. Set knots at the three quartiles of the explanatory variable
# and boundary knots at the 5th and 95th percentiles. Provide a plot of the regression
# line with standard errors for the prediction, a table summarizing the coefficients and comment the results.

# Calculate knots and boundary knots
knots <- quantile(milk_samples_2$wave_70, probs = c(0.25, 0.5, 0.75))
boundary_knots <- quantile(milk_samples_2$wave_70, probs = c(0.05, 0.95))

# Fit a natural cubic spline model
fit_spline <- lm(Native_pH ~ ns(wave_70, knots = knots, Boundary.knots = boundary_knots),
                 data = milk_samples_2)

# Display coefficients summary in a table
knitr::kable(broom::tidy(summary(fit_spline)))

# Generate sequence of values for predictions
new_data_seq <- seq(min(milk_samples_2$wave_70), max(milk_samples_2$wave_70), length.out = 100)

# Predict values and standard errors
# So, preds$fit will contain the predicted values of the response variable (Native_pH) 
# for each value in new_data_seq, and preds$se.fit will contain the standard errors of these predictions.
preds <- predict(fit_spline, newdata = list(wave_70 = new_data_seq), se = TRUE) # se = TRUE : standard error computed
se.bands <- cbind(preds$fit + 2 * preds$se.fit, preds$fit - 2 * preds$se.fit)
# So, se.bands is a matrix where the first column represents the upper bounds of 
# the standard error bands, and the second column represents the lower bounds. 
# This matrix is used to plot the standard error bands around the predicted regression line in the subsequent plot.

# Plot the data, regression line, and standard error bands
plot(y = milk_samples_2$Native_pH, x = milk_samples_2$wave_70, cex = 0.5, col = "darkgrey ")
lines(new_data_seq, preds$fit, lwd = 2, col = "blue")
matlines(new_data_seq, se.bands, lwd = 1, col = "blue", lty = 3)

# Plot knots and boundary knots
knots_pred <- predict(fit_spline, list(wave_70 = knots))
points(knots, knots_pred, col = 'blue', pch = 19)
boundary_pred <- predict(fit_spline, list(wave_70 = boundary_knots))
points(boundary_knots, boundary_pred, col = 'red', pch = 19)

# Add vertical lines for knots and boundary knots
abline(v = knots, lty = 3, col = "blue")
abline(v = boundary_knots, lty = 3, col = "red")

## 

# 2) Build an additive model for regressing Milk pH on the milk absorbance at wavenumbers 70cm−1 and 
# 300cm−1 and their interaction, using penalized cubic b-spline terms for smoothing. After having written 
# in proper mathematical terms the additive model you have estimated, report the adjusted R2 and the p-values
# of the tests.

# Our additive model : pH_i = β_0 + f(wave300_i) + f(wave70_i) + f(wave70_i ∗ wave300_i) + ε_i i=1,...,N.
# Fit a GAM model
fit_gam <- gam(Native_pH ~ s(wave_70, bs = "cr") + s(wave_300, bs = "cr") + s(I(wave_300 * wave_70), bs = "cr"), 
               data = milk_samples_2)
# Native_pH ~ s(wave_70, bs = "cr") + s(wave_300, bs = "cr") + s(I(wave_300 * wave_70), bs = "cr")
# : Specifies the model formula. It models the response variable Native_pH as a smooth function of 
# wave_70, wave_300, and the interaction term wave_300 * wave_70, USING REGRESSION SPLINES BASIS (bs = "cr").

# Summarize the results of the GAM fit
table_fit_gam <- summary(fit_gam)

# Extract and print the R-squared value
(r_2_squared <- table_fit_gam$r.sq)

table_fit_gam$p.pv # Extracts the p-values associated with the estimated smooth terms in the GAM.
# A low p-value (typically below a chosen significance level, e.g., 0.05)
#suggests that the corresponding smooth term is statistically significant, and you may reject
# the null hypothesis that the term has no effect.

table_fit_gam$s.table # You can inspect s_table to obtain detailed information about each smooth term in the GAM.

##

# 3) Build an additive model for regressing Milk pH on the milk absorbance at wavenumbers 70cm−1
# and 300cm−1 with no interaction, using again cubic b-spline terms. Assuming the residuals to be normal,
# use the Anova F test statistic to assess the significance of the interaction term, specifying the null
# and the alternative hypothesis you are testing and report the resulting p-value. Comment on the results.

# Our new additive model : pH_i = β_0 + f(wave300_i) + f(wave70_i) + ε_i i=1,...,N.

## Fit the GAM model
fit_gam_reduced <- gam(Native_pH ~ s(wave_70, bs = "cr") + s(wave_300, bs = "cr"), data = milk_samples_2)

## Performing Anova F Test:
# You use the anova function to compare the full GAM (fit_gam) with the reduced GAM (fit_gam_reduced).
# The test is specified as an F test using test = "F".
# The null hypothesis H0 is that the interaction term (s(I(wave_300 * wave_70), bs = "cr")) is not significant, 
# meaning it does not contribute significantly to explaining the variance in Native_pH.
# The alternative hypothesis H1 is that at least one of the terms in the interaction is significant.
anova_test <- anova(fit_gam_reduced, fit_gam, test = "F")
# 1) Residual Degrees of Freedom (Resid. Df):
#    Model 1 has 371.08 residual degrees of freedom.
#    Model 2 has 368.37 residual degrees of freedom.
# 2) Residual Deviance (Resid. Dev):
#    Model 1 has a residual deviance of 2.9385.
#    Model 2 has a slightly lower residual deviance of 2.8638.
# 3) Degrees of Freedom Change (Df):
#    The change in degrees of freedom is calculated as Df(Model 2) - Df(Model 1).
#    In this case, it's 2.7137.
# 4) Deviance Change (Deviance):
#    The change in deviance is calculated as Deviance(Model 1) - Deviance(Model 2).
#    In this case, it's 0.074648.
# 5) F-statistic (F):
#    The F-statistic is calculated as Deviance Change / Df.
#    Here, it's 0.074648 / 2.7137 = 3.5616.
# 6) p-value (Pr(>F)):
#    The p-value associated with the F-statistic tests the null hypothesis that the additional terms (interaction term) do not significantly improve the model fit.
#    In this case, the p-value is 0.01767.
#
# Interpretation:
# The p-value (Pr(>F)) is less than 0.05, indicating that the interaction term significantly 
# improves the model fit.
# The asterisk (*) next to the p-value suggests that the improvement is statistically significant
# at the 0.05 significance level.
# Full model is thus preferred

##

# 3) Build a semi-parametric model for regressing Milk pH on the milk absorbance at wavenumbers 70cm−1
# and 300cm−1, assuming the effect of wavenumber 300cm−1 to be linear on the response and no interaction.
# Using a bootstrap approach, provide a reverse percentile confidence interval for the parametric effect.

# Fit a semi-parametric model
fit_semi <- gam(Native_pH ~ s(wave_70, bs = "cr") + wave_300, data = milk_samples_2)

# Obtain observed fitted values and residuals
fitted.obs <- fit_semi$fitted.values
res.obs <- fit_semi$residuals

summary(fit_semi)$p.table
wave300_obs <- summary(fit_semi)$p.table[2, 1]

# Bootstrap settings
B <- 1000  # Number of bootstrap samples
T.boot.wave300 <- numeric(B)

# Bootstrap procedure
set.seed(2022)
for (b in 1:B) {
  # Generate a bootstrap sample by resampling residuals
  Native_pH_boot <- fitted.obs + sample(res.obs, replace = TRUE)
  
  # Fit a semi-parametric model to the bootstrap sample
  fit_semi_boot <- gam(Native_pH_boot ~ s(milk_samples_2$wave_70, bs = "cr") + milk_samples_2$wave_300)
  fit_semi_boot_table <- summary(fit_semi_boot)
  
  # Record the parametric effect for wavenumber 300cm−1
  T.boot.wave300[b] = fit_semi_boot_table$p.table[2, 1]
}

# Calculate reverse percentile confidence interval
alpha <- 0.05
right.quantile.wave300 <- quantile(T.boot.wave300, 1 - alpha/2)
left.quantile.wave300 <- quantile(T.boot.wave300, alpha/2)
CI.RP.wave300 <- c(
  wave300_obs - (right.quantile.wave300 - wave300_obs),
  wave300_obs - (left.quantile.wave300 - wave300_obs)
)
names(CI.RP.wave300) <- c('lwr', 'upr')

CI.RP.wave300






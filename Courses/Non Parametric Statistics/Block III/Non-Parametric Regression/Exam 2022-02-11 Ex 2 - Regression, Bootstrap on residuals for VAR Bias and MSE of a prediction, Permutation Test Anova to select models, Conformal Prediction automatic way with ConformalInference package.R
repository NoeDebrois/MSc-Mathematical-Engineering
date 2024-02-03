################################################################################
################################ EXERCISE II ###################################
################################################################################
# Matthew O’Fountain knows that spectroscopy is the state-of-the-art technology to
# employ when it comes to evaluate milk quality. Motivated by this he is interested
# in building a nonparametric model to predict κ-casein by means of the absorbance
# values at wavenumbers 280 and 700 cm−1, contained in the milk_samples_2.Rds file.

##

# 1. Build a smoothing spline model to regress κ-casein on the milk absorbance at 
# wavenumber 280 cm−1, selecting λ by means of Generalized Cross Validation. Provide
# a plot of the regression line and a point-wise estimate for κ-casein when the milk
# absorbance at wavenumber 280 cm−1 is equal to 1.05. By using a bootstrap approach
# on the residuals, calculate the bias variance and MSE of such prediction
# (fix the λ value to the one obtained via Generalized Cross Validation).

milk_samples_2 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-02-11/milk_samples_2.Rds")

## SMOOTHING SPLINE MODEL WITH GCV (FOR LAMBDA) AND POINTWISE ESTIMATE (PREDICTION) :
### My proposition :
# data_1 = milk_samples_2[-3]
# fit_smooth_spline_GCV <- with(data_1, smooth.spline(wave_280,kappa_casein,cv = FALSE))
# with(data_1, plot(wave_280, kappa_casein, cex =.5, col =" darkgrey "))
# lines(fit_smooth_spline_GCV,col="blue",lwd=2, lty=2)
# legend(20000, 30, legend=c("GCV"),
#        col=c("blue"), lty=1:2, cex=0.8)
# lambda_GCV = fit_smooth_spline_GCV$lambda
# prediction_at_1_05 = predict(fit_smooth_spline_GCV,1.05) # predict for a certain value of absorbance

### Correction : (more compact)
plot(milk_samples_2$wave_280, milk_samples_2$kappa_casein)
fit_smooth <- smooth.spline(x = milk_samples_2$wave_280, 
                            y = milk_samples_2$kappa_casein,
                            cv = FALSE) # FALSE pour utiliser GCV pour lambda
plot(milk_samples_2$wave_280, milk_samples_2$kappa_casein)
lines(fit_smooth, col="red")

pred <- predict(fit_smooth,x=1.05)
y_hat <- pred$y
y_hat

lambda_GCV = fit_smooth$lambda
lambda_GCV

## BOOTSTRAP APPROACH ON RESIDUALS TO CALCULATE BIAS, VAR, AND MSE OF SUCH PREDICTION (1.05's prediction) :
### Correction :
fitted=predict(fit_smooth,milk_samples_2$wave_280)$y # get all the predictions
residuals=milk_samples_2$kappa_casein-fitted # residuals = vraies_valeurs - predictions

N <- nrow(milk_samples_2)
B <- 1000
boot_d=numeric(B)
set.seed(2022)
for(sample in 1:B){
  kappa_casein_boot=fitted+sample(residuals,N,replace=T) # tirage bootstrap avec les residuals
  new_model = smooth.spline(x = milk_samples_2$wave_280, # fit un nouveau modèle avec en x le wave_280
                            y = kappa_casein_boot, # fit un nouveau modèle avec en y le kappa_casein_boot "en cours"
                            lambda = lambda_GCV) # fit un nouveau modèle avec lambda_GCV
  boot_d[sample]=predict(new_model, 1.05)$y # nouvelle prédiction à la valeur x = 1.05
}

(variance_105 = var(boot_d))
(bias_105  = mean(boot_d)-y_hat)
(MSE_105 = variance_105 + bias_105*bias_105)

##

# 2. Build an additive model for regressing κ-casein on the milk absorbance at 
# wavenumbers 280 cm−1 and 700 cm−1, using degree 2 b-spline bases with just one 
# knot at the median as univariate smoother for the two predictors. Report the 
# summary table including a comment on statistical significance. Provide an 
# histogram of the residuals of the model.

model_quad_splines <- lm(kappa_casein ~ 
                           bs(wave_280,degree = 2, knots = median(milk_samples_2$wave_280)) # préciser le degré, et le noeud pour chaque covariate 
                           + bs(wave_700,degree = 2, knots = median(milk_samples_2$wave_700)), # préciser le degré, et le noeud pour chaque covariate 
                           data = milk_samples_2) # data
summary(model_quad_splines) 
### Comment on statistical significance :
# The 280's coefficients represent the estimated effects of the basis splines for the 
# variable wave_280. The asterisks in the p-values indicate statistical significance.
# The second and third basis functions are highly significant (p-values < 0.05), 
# suggesting that they contribute significantly to the model.
# The 700's coefficients represent the estimated effects of the basis splines for the 
# variable wave_700. None of the coefficients for wave_700 are statistically significant
# (all p-values > 0.05), suggesting that these terms may not be contributing significantly
# to the model.

hist(model_quad_splines$residuals) # histogramme des residuals du modèle

##

# 3. Build a reduced version of the previous model considering only the contribution
# of the milk absorbance at wavenumber 280 cm−1 for explaining κ-casein. Employ a 
# permutational Anova (using the F value as test statistic) to validate which model
# should be preferred, specifying the null and the alternative hypothesis you are 
# testing and report the resulting p-value. Comment on the results.

model_quad_splines_reduced <- lm(kappa_casein ~
                                   bs(wave_280, degree = 2, df = 3),
                                   data = milk_samples_2)
fitted.obs <- model_quad_splines_reduced$fitted.values
res.obs <- model_quad_splines_reduced$residuals

T_0 <- anova(model_quad_splines_reduced,model_quad_splines)[2,5]
T_0

# Estimating the permutational distribution under H0 : 
B <- 1000
T2 <- numeric(B)
set.seed(2022)

for (perm in 1:B) {
  res_reduced_perm <- res.obs[sample(1:N)]
  y_perm <- fitted.obs + res_reduced_perm
  # Création d'un nouveau dataset avec les réponses permutées :
  milk_samples_2_perm <- milk_samples_2
  milk_samples_2_perm$kappa_casein <- y_perm
  # Fit du modèle classique et du modèle réduit :
  model_quad_splines_perm <- lm(kappa_casein ~
                                  bs(wave_280, degree = 2) + bs(wave_700, degree = 2),
                                  data = milk_samples_2_perm)
  model_quad_splines_reduced_perm <- lm(kappa_casein ~ bs(wave_280, degree = 2),
                                        data = milk_samples_2_perm)
  
  T2[perm] <- anova(model_quad_splines_reduced_perm,model_quad_splines_perm)[2,5]
}

hist(T2, xlim = range(c(T2, T_0)))
abline(v = T_0, col = 3, lwd = 4)

plot(ecdf(T2))
abline(v = T_0, col = 3, lwd = 4)

p_val <- sum(T2 >= T_0) / B
p_val
# H0 : the 700 cm^-1 is not significant ;
# H1 : the 700 cm^-1 is significant.
# p_val > 0.05 so :
# Cannot reject H0 : so the reduced model is better
# CF : ~/Desktop/Politecnico/Non Parametric Statistics/BIBLE/Block II/Block II Bible NPS 6 - Permutation Tests 3 - (M)ANOVA.R
# (CHEMICAL CONFIRMATION: wave 700 is associated with the presence of water in the sample,
# and it is known to be not informative of milk quality)

# 4. Compute the prediction bands for the regression model selected according to
# the test performed in the previous exercise, using a full conformal approach and
# setting α = 0.05 as the miscoverage level

# cf : "Tutorial on ConformalInference package" from /Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/BIBLE/Block IV/Lab 14 - Conformal Prediction.pdf

library(conformalInference)

# According to what we did just before, we keep the reduced model (only 280cm^-1) :

# Generating a Grid of 100 `wave_280` Values (between min and max) :
wave_280_grid = seq(range(milk_samples_2$wave_280)[1],
                    range(milk_samples_2$wave_280)[2],
                    length.out = 100)

# Predictions with the Reduced Model :
preds = predict(model_quad_splines_reduced, # modèle réduit (choisit précédemment)
                list(wave_280 = wave_280_grid), # prédictions sur la grille wave_280_grid
                se = T) # standard errpr

# Scatter Plot of Data:
with(
  milk_samples_2,
  plot(
    wave_280,
    kappa_casein,
    xlim = range(wave_280_grid),
    cex = 0.5,
    col = "darkgrey")
)

# Conformal Prediction Interval Calculation:
# Step 1: Training Functions for Linear Model (Including Intercept):
lm_train = lm.funs(intercept = T)$train.fun
lm_predict = lm.funs(intercept = T)$predict.fun

# Step 2: Design Matrix for the Original Data (`wave_280`):
# - Creates a matrix of basis splines for the `wave_280` variable with a quadratic spline (degree = 2) (and 3 degrees of freedom).
design_matrix = bs(milk_samples_2$wave_280, degree = 2) # df = 3

# Step 3: Design Matrix for the Prediction Grid (`wave_280_grid`):
# - Creates a matrix of basis splines for the prediction grid with the same specifications as the original data.
pred_grid = matrix(bs(wave_280_grid, degree = 2), nrow = length(wave_280_grid)) # df = 3 

# Step 4: Conformal Prediction Interval Computation:
c_preds = conformal.pred(
  x = design_matrix,                  # Design matrix for the original data
  y = milk_samples_2$kappa_casein,    # Response variable
  pred_grid,                          # Design matrix for the prediction grid
  alpha = 0.05,                       # Significance level for the prediction interval
  verbose = F,                        # Disable verbosity (optional)
  train.fun = lm_train,               # Training function for the linear model
  predict.fun = lm_predict,           # Prediction function for the linear model
  num.grid.pts = 200                  # Number of grid points for visualization
)

# Plotting Conformal Prediction Interval:
lines(
  wave_280_grid,
  c_preds$pred,
  lwd = 2,
  col = "red",
  lty = 3
)

matlines(
  wave_280_grid,
  cbind(c_preds$up, c_preds$lo),
  lwd = 1,
  col = "blue",
  lty = 3
)















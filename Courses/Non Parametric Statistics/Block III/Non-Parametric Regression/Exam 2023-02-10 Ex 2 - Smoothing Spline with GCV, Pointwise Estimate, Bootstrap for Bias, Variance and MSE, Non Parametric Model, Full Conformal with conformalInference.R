################################################################################
################################# EXERCICE II ##################################
################################################################################
df_2 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2023-02-10/df_2.rds")
# Professor Franziska Iünz is now interested in knowing more about the 
# relationship between the amount of money (in k$) that shall be requested as a 
# function of the duration of funding (measured in month). To do so, she has at
# her disposal a dataset on 67 proposals, previously funded, contained in the 
# df_2.Rds file. Help Professor Iünz by solving the following tasks.
library(ISLR2)
library(car)
# 1. Build a smoothing spline model to regress the money response on the 
# funding.month variable, selecting λ by means of Generalized Cross Validation
# (GCV). Provide a plot of the regression line and report the best λ value
# identified via GCV (round it at the third decimal digit).
N <- nrow(df_2)
fit_smooth <- smooth.spline(x = df_2$funding.month,
                            y = df_2$money,
                            cv = FALSE) # To have GCV !
plot(df_2$funding.month, df_2$money)
lines(fit_smooth, col="red")

lambda_GCV = fit_smooth$lambda
round(lambda_GCV, 3)

# 2. Ideally, Professor Iünz would like to have a research proposal lasting for 
# five years. Using the model estimated in the previous exercise, provide a 
# point-wise estimate for the amount of money (in k$) she shall request. 
# In addition, by using a bootstrap approach on the residuals, estimate the bias,
# variance and MSE of such prediction (fix the λ value to the one obtained via
# Generalized Cross Validation).
y_hat = predict(fit_smooth, 60)$y # Point-wise estimate for the amount of money (in k$) she shall request.
y_hat

fitted=predict(fit_smooth,df_2$funding.month)$y # get all the predictions
residuals=df_2$money-fitted # residuals = vraies_valeurs - predictions

B <- 1000
boot_d=numeric(B)
set.seed(2022)
for(sample in 1:B){
  money_boot=fitted+sample(residuals,N,replace=T) # tirage bootstrap avec les residuals
  new_model = smooth.spline(x = df_2$funding.month, # fit un nouveau modèle avec en x le funding.month
                            y = money_boot, # fit un nouveau modèle avec en y money
                            lambda = lambda_GCV) # fit un nouveau modèle avec lambda_GCV
  boot_d[sample]=predict(new_model, 60)$y # nouvelle prédiction à la valeur x = 60 (mois)
}

(variance_105 = var(boot_d))
(bias_105  = mean(boot_d)-y_hat)
(MSE_105 = variance_105 + bias_105*bias_105)

# 3. Build a nonparametric model to regress the money response on the 
# funding.month variable via a step regression procedure, breaking the range of 
# funding.month into 4 bins, using the quartiles of the covariate as breaking 
# points. Provide a plot of the regression line and report the summary table 
# including a comment on statistical significance.
funding.month_quart <- quantile(df_2$funding.month) # to get the quartiles
m_cut=lm(money ~ cut(funding.month, breaks=c(-Inf,funding.month_quart[c(-1,-5)],Inf)), data = df_2)

month_grid = with(df_2, seq(range(funding.month)[1], range(funding.month)[2],by=1)) 
preds=predict(m_cut, list(funding.month=month_grid), se=T)

# plot with the standard error bands :
se.bands = cbind(preds$fit + 2 * preds$se.fit, preds$fit - 2 * preds$se.fit)
with(df_2, plot(funding.month, money, xlim=range(month_grid), cex =.5, col = "darkgrey")) 
lines(month_grid, preds$fit, lwd = 2, col = "blue", type="s")
matlines(month_grid, se.bands, lwd =1, col ="blue", lty = 3, type="s")

# Comment :
coef(m_cut)
summary(m_cut)
# All coefs are highly significant we look at the pvalues

# 4. Compute the prediction bands for the regression model of the previous 
# exercise, using a full conformal approach and setting α = 0.1 as the 
# miscoverage level.

# Load the 'conformalInference' package for conformal prediction
library(conformalInference)

# Generate a grid of values for 'funding.month'
funding.month_grid = seq(range(df_2$funding.month)[1],
                                range(df_2$funding.month)[2],
                                length.out = 100)

# Plot the original data points with reduced transparency
with(df_2,
     plot(
       funding.month,
       money,
       xlim = range(funding.month),
       cex = .5,
       col = "darkgrey",
       ylim=c(200,2500)
     )
)

# Extract training and prediction functions from linear regression model
lm_train = lm.funs(intercept = T)$train.fun
lm_predict = lm.funs(intercept = T)$predict.fun

# Create a design matrix for the model
design_matrix = model.matrix(m_cut)[,-1]

pred_cut = cut(month_grid,breaks=c(-Inf,funding.month_quart[c(-1,-5)],Inf))
pred_grid = model.matrix(lm(rep(1,length(pred_cut))~pred_cut))[,-1]

# Perform conformal prediction on the design matrix
c_preds = conformal.pred(
  x = design_matrix,
  y = df_2$money,
  pred_grid,
  alpha = 0.1,
  verbose = F,
  train.fun = lm_train,
  predict.fun = lm_predict,
  num.grid.pts = 200
)

# Plot the central predicted values in red with dashed line
lines(
  x = month_grid,
  y = c_preds$pred,
  lwd = 2,
  col = "red",
  lty = 3,
  type = "s"
)

# Plot the upper and lower prediction bands in blue with dashed line
matlines(
  month_grid ,
  cbind(c_preds$up, c_preds$lo) ,
  lwd = 1,
  col = " blue",
  lty = 3,type="s"
)








################################################################################
################################## EXERCISE II #################################
################################################################################
library(conformalInference)
percebes_2 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-07-11/percebes_2.Rds")

# Dr. Andreas Qapos, Ph.D, after having deducted from what you have told him 
# after Exercise 1 that the relationship between length and weight of the 
# barnacles is not a linear one, wants to understand what is the best time to
# harvest his precious goose barnacles. To do so, he has collected 999 barnacles
# from his best spot, and wants to build a prediction model for the weight of 
# the barnacle (which can be assessed only after picking it up. . . ) as a
# function of the length of the barnacle (which can be assessed using a special 
# barnacle caliber). You can find the data in percebes_2.rds, Now :

# 1. Help Dr. Qapos in building his model: specifically, he would like to use a
# smoothing spline model of order 4, with lambda selected via cross-validation. 
# Please report the lambda value with 3 decimal digits, as well as a plot of the
# regression line, alongside the pointwise prediction of the weight of a barnacle
# whose length is 33 mm, which from a visual inspection of the data seems the 
# optimal value for harvesting.
plot(percebes_2$length, percebes_2$weight)
fit_smooth <- smooth.spline(x = percebes_2$length, 
                            y = percebes_2$weight,
                            cv = TRUE) # TRUE pour utiliser CV pour lambda

pred <- predict(fit_smooth,x=33)
y_hat <- pred$y
y_hat

plot(percebes_2$length, percebes_2$weight)
lines(fit_smooth, col="red")
points(33, y_hat, col='blue') # quite an optimal value for harvesting !

lambda_GCV = fit_smooth$lambda
lambda_GCV

# 2. Assess the uncertainty of this prediction by calculating, via residual 
# bootstrapping, its bias, its variance and its MSE.
fitted=predict(fit_smooth,percebes_2$length)$y # get all the predictions
residuals=percebes_2$weight-fitted # residuals = vraies_valeurs - predictions

N <- nrow(percebes_2)
B <- 1000
boot_d=numeric(B)
set.seed(2022)
for(sample in 1:B){
  weight_boot=fitted+sample(residuals,N,replace=T) # tirage bootstrap avec les residuals
  new_model = smooth.spline(x = percebes_2$length, # fit un nouveau modèle avec en x le wave_280
                            y = weight_boot, # fit un nouveau modèle avec en y le kappa_casein_boot "en cours"
                            lambda = lambda_GCV) # fit un nouveau modèle avec lambda_GCV
  boot_d[sample]=predict(new_model, 33)$y # nouvelle prédiction à la valeur x = 1.05
}

(variance_105 = var(boot_d))
(bias_105  = mean(boot_d)-y_hat)
(MSE_105 = variance_105 + bias_105*bias_105)

# 3. Like every good statistician, Dr. Qapos is not happy at all with a pointwise
# prediction. After having defined the distributional hypotheses behind this approach,
# and stated its coverage properties, build a prediction interval for the weight
# of a barnacle 33mm long using a full conformal prediction approach, using as a
# non-conformity score the absolute value of the regression residuals.

# We use smoothing spline, so it is non-linear. In a “nonlm” case, custom 
# functions for the point predictor need to be created, In particular, we need 
# to explicitly define the prediction model objects, which will be used to 
# output predictions.

# Now define the function for training the point predictor and the function 
# actually used for prediction.
train_ss <- function(x, y, out=NULL){
  smooth.spline(x, y, lambda=lambda_GCV)
}
predict_ss <- function(obj, new_x){
  predict(obj, new_x)$y
}

# These ad-hoc functions can be given as input to the conformal.pred function:
set.seed(2022)
pred = conformal.pred(percebes_2$length,
                      percebes_2$weight,
                      33,
                      train.fun=train_ss,
                      predict.fun=predict_ss,
                      alpha=0.05)
data.frame(lwr=pred$lo,
           pred=pred$pred,
           upr=pred$up)




















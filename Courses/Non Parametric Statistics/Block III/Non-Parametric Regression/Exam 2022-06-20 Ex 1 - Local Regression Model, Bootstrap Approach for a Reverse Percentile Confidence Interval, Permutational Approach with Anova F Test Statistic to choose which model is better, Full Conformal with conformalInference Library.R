################################################################################
################################## EXERCISE I ##################################
################################################################################
# Mr András Kaponzyi, an Hungarian data scientist, is about to become a father of
# his first daughter. During the gynecological examination of the second quarter,
# the fetus weighted 353 grams. Mr András Kaponzyi is worried that the newborn
# will be very chubby, he is therefore interested in estimating the expected
# birth weight of his daughter. To do so, he has collected 240 measurements of
# fetus’weights (in grams) at 20 weeks pregnancy (baby_weight_20_weeks), and the
# respective weight in grams at birth (baby_weight_birth). The resulting samples
# are contained in the df_1.Rds file. Assume that:
# baby_weight_birth = f (baby_weight_20_weeks) + ε

library(np)

df_1 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-06-20/df_1.Rds")
 
# 1. Fit a local regression model with a Gaussian kernel and bandwidth equal to
# 10 grams to predict the birth weight as a function of the weight at 20 weeks
# pregnancy. Provide a plot of the regression curve and compute the point-wise
# prediction (round it at the second decimal digit) for the birth weight of
# Kaponzyi’s daughter. Using a bootstrap approach, provide a reverse percentile
# confidence interval (α = 0.05) for the expected weight of Kaponzyi’s daughter.

m_loc = npreg(baby_weight_birth ~ baby_weight_20_weeks,
              ckertype = 'gaussian', # kernel type
              bws = 10, # bandwidth 
              residuals=TRUE, # TO GET THE RESIDUALS (later useful for the bootstrap !)
              data = df_1)

measures=data.frame(baby_weight_20_weeks = with(df_1, seq(
  range(baby_weight_20_weeks)[1],
  range(baby_weight_20_weeks)[2],
  by=1)))

preds=predict(m_loc, newdata=measures, se=T)
se.bands=cbind(preds$fit + 2 * preds$se.fit, preds$fit - 2 * preds$se.fit)

# Plot :
with(df_1, 
     plot(baby_weight_20_weeks,
          baby_weight_birth,
          xlim = range(df_1$baby_weight_20_weeks),
          cex = .5,
          col = " darkgrey ",
          main = 'Local Averaging - bws10 - Gaussian kernel'
          )
     )
lines(measures$baby_weight_20_weeks,preds$fit ,lwd =2, col =" blue")
matlines(measures$baby_weight_20_weeks,se.bands ,lwd =1, col =" blue",lty =3)

# Point-wise prediction :
pred_kaponzyi <- predict(m_loc,newdata=list(baby_weight_20_weeks=353))
round(pred_kaponzyi,2)
points(353, pred_kaponzyi, col = "red", pch = 16) # add it to the plot, in red

# Reverse percentile confidence interval, by bootstrap :
B = 1000
res.obs <- m_loc$resid # because at the line 28 we mentioned : "residuals=TRUE"
fitted.obs <- m_loc$mean

boot_pred_kaponzyi = numeric(B)
set.seed(2022)
for (b in 1:B) {
  baby_weight_birth_boot <- fitted.obs + sample(res.obs, replace = T)
  df_boot <- data.frame(baby_weight_birth_boot,
                        baby_weight_20_weeks = df_1$baby_weight_20_weeks)
  fit_kernel_boot <- npreg(baby_weight_birth_boot ~ baby_weight_20_weeks, 
                           ckertype = 'gaussian', bws = 10, data = df_boot)
  boot_pred_kaponzyi[b] = predict(fit_kernel_boot, newdata = list(baby_weight_20_weeks = 353))
}

alpha <- 0.05

right.quantile.pred_kaponzyi <- quantile(boot_pred_kaponzyi, 1 - alpha / 2)
leftquantile.pred_kaponzyi <- quantile(boot_pred_kaponzyi, alpha / 2)

CI.RP.pred_kaponzyi <- c(pred_kaponzyi - (right.quantile.pred_kaponzyi - pred_kaponzyi),
                         pred_kaponzyi - (leftquantile.pred_kaponzyi - pred_kaponzyi))
names(CI.RP.pred_kaponzyi) = c('lwr', 'upr')
CI.RP.pred_kaponzyi

# 2. Fit a polynomial regression model of degree 5 to predict the birth weight as
# a function of the weight at 20 weeks pregnancy. Provide a plot of the 
# regression line and a table with the estimated model coefficients. By using a
# permutational approach, use the Anova F test statistic to check whether a
# polynomial regression model of degree 6 should be preferred. Specify the null
# and the alternative hypothesis you are testing and report the resulting
# p-value. Comment on the result.
fit_poly5 <- lm(baby_weight_birth ~ poly(baby_weight_20_weeks,
                                         degree = 5,
                                         raw = F),
                data = df_1)
preds_poly5 = predict(fit_poly5,
                      list(baby_weight_20_weeks = measures$baby_weight_20_weeks),
                      se = T)

se.bands_poly5 = cbind(preds_poly5$fit + 2 * preds_poly5$se.fit ,
                       preds_poly5$fit - 2 * preds_poly5$se.fit)
with( df_1,
      plot(
        baby_weight_20_weeks ,
        baby_weight_birth ,
        xlim = range(measures$baby_weight_20_weeks) ,
        cex = .5,
        col = " darkgrey ",
        main = 'Degree 5 Poly fit'
      ) )

lines(measures$baby_weight_20_weeks,
      preds_poly5$fit ,
      lwd = 2,
      col = " blue")

matlines(
  measures$baby_weight_20_weeks ,
  se.bands_poly5 ,
  lwd = 1,
  col = " blue",
  lty = 3)

knitr::kable(broom::tidy(summary(fit_poly5)))

# Permutational anova :
# The anova F-test is used to assess whether the more complex model provides a
# significantly better fit than the simpler model. If the F-value is
# significantly larger than expected by chance, it suggests that the more
# complex model is justified.

fit_poly5_reduced <- fit_poly5
fit_poly6_full <- lm(baby_weight_birth ~ poly(baby_weight_20_weeks,
                                              degree = 6,
                                              raw = F),
                     data = df_1)

fitted.obs <- fit_poly5_reduced$fitted.values
res.obs <- fit_poly5_reduced$residuals

T_0 <- anova(fit_poly5_reduced,fit_poly6_full)[2,5]

# Estimating the permutational distribution under H0 : 
B <- 1000
T2 <- numeric(B)
set.seed(2022)
N = nrow(df_1)

for (perm in 1:B) {
  res_reduced_perm <- res.obs[sample(1:N)]
  y_perm <- fitted.obs + res_reduced_perm # Under H0 : no problem
  # Création d'un nouveau dataset avec les réponses permutées :
  df_1_perm <- df_1
  df_1_perm$baby_weight_birth <- y_perm
  # Fit du modèle classique et du modèle réduit :
  fit_poly6_full_perm <- lm(baby_weight_birth ~ poly(baby_weight_20_weeks,
                                                     degree = 6,
                                                     raw = F),
                            data = df_1_perm)
  fit_poly5_reduced_perm <- lm(baby_weight_birth ~ poly(baby_weight_20_weeks,
                                                        degree = 5,
                                                        raw = F),
                               data = df_1_perm)
  
  T2[perm] <- anova(fit_poly5_reduced_perm, fit_poly6_full_perm)[2, 5]
}

hist(T2, xlim = range(c(T2, T_0)))
abline(v = T_0, col = 3, lwd = 4)

plot(ecdf(T2))
abline(v = T_0, col = 3, lwd = 4)

p_val <- sum(T2 >= T_0) / B
p_val
# H0 : The null hypothesis assumes that there is no significant
# difference in the model fit between the polynomial regression model of degree
# 6 and the alternative model of degree 5. 
# H1 : The alternative hypothesis suggests that there is a significant 
# improvement in the model fit when using the polynomial regression model of
# degree 6 compared to the alternative model.
# p_val > 0.05 so : cannot reject H0 so the reduced model is better (well, no
# significant improvement in the full model so "autant choisir le modèle réduit").

# 3. Compute the prediction bands for the regression model selected according to 
# the test performed in the previous exercise, using a full conformal approach
# and setting α = 0.1 as the miscoverage level.

# cf : "Tutorial on ConformalInference package" from /Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/BIBLE/Block IV/Lab 14 - Conformal Prediction.pdf

# Load the 'conformalInference' package for conformal prediction
library(conformalInference)

# Generate a grid of values for 'baby_weight_20_weeks'
baby_weight_20_weeks_grid = seq(range(df_1$baby_weight_20_weeks)[1],
                                range(df_1$baby_weight_20_weeks)[2],
                                length.out = 100)

# Plot the original data points with reduced transparency
with(df_1,
     plot(
       baby_weight_20_weeks,
       baby_weight_birth,
       xlim = range(baby_weight_20_weeks_grid),
       cex = .5,
       col = "darkgrey"
     )
)

# Extract training and prediction functions from linear regression model
lm_train = lm.funs(intercept = T)$train.fun
lm_predict = lm.funs(intercept = T)$predict.fun

# Create a design matrix for polynomial regression of degree 5
design_matrix = matrix(poly(df_1$baby_weight_20_weeks, degree = 5), ncol = 5)

# Generate a matrix for polynomial regression prediction on the grid
pred_grid = matrix(poly(
  baby_weight_20_weeks_grid,
  degree = 5,
  coefs = attr(poly(df_1$baby_weight_20_weeks, degree = 5), "coefs")
), ncol = 5)

# Perform conformal prediction on the design matrix
c_preds = conformal.pred(
  design_matrix,
  df_1$baby_weight_birth,
  pred_grid,
  alpha = 0.1,
  verbose = F,
  train.fun = lm_train,
  predict.fun = lm_predict,
  num.grid.pts = 200
)

# Plot the central predicted values in red with dashed line
lines(
  x = baby_weight_20_weeks_grid,
  y = c_preds$pred,
  lwd = 2,
  col = "red",
  lty = 3
)

# Plot the upper and lower prediction bands in blue with dashed line
matlines(
  baby_weight_20_weeks_grid,
  cbind(c_preds$up, c_preds$lo),
  lwd = 1,
  col = "blue",
  lty = 3
)


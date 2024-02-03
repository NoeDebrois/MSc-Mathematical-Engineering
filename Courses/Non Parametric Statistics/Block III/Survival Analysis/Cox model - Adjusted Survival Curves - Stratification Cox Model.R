library(survival)
library(tidyverse)
library(survminer)

# Cox model example :
# Cox model is used to assess how the covariates affect the hazard function.

# Ovarian Cancer Survival Data
tibble::glimpse(ovarian)

# fustat: censoring status (0 = censored, 1 = observed)
# rx: treatment group (1 = cyclophosphamide alone,
# 2 = cyclophosphamide plus adriamycin)

ovarian$rx <- factor(ovarian$rx)
ovarian$resid.ds <- factor(ovarian$resid.ds)
ovarian$ecog.ps <- factor(ovarian$ecog.ps)

fit_ovarian <- coxph(Surv(futime, fustat) ~ age + rx,
                 data=ovarian)

summary(fit_ovarian)

################################################################################
################################################################################
# Adjusted Survival Curves #####################################################
################################################################################
################################################################################
# (!!) cf SLIDE 11/24 du cours (!!)
# Essentially, this part of the code generates adjusted survival curves for 
# hypothetical individuals with specific characteristics, as predicted by the Cox
# model fitted earlier. The survival curves show how the survival probability 
# changes over time for individuals with the specified characteristics.

# This data frame is passed to survfit() via the newdata argument:
new_ovarian_age <- data.frame(age = c(50,60,70),
                          resid.ds = factor(rep(2,3)),
                          rx = factor(c(1,1,1)),
                          ecog.ps = factor(rep(1,3)))

fit_rx_age <- survfit(fit_ovarian, newdata = new_ovarian_age)
# Use survfit to calculate adjusted survival curves:
# The survfit function takes the fitted Cox model (fit_ovarian) and uses it to 
# calculate survival probabilities over time for the new individuals specified in
# the new_ovarian_age data frame. The result (fit_rx_age) is an object that 
# contains information about the adjusted survival curves for the specified 
# characteristics (age, treatment group, etc.).


# Estimated adjusted survival curves:
plot(fit_rx_age, col=c("dodgerblue2","navy","darkmagenta"), lwd=2, lty=1,
     xlab='Time [days]', ylab='Survival Probability',
     main='Adjusted Survival Probability Plot')
legend('topright', c("Age = 50", "Age = 60", "Age = 70"),
       lty=c(1,1,1), lwd=c(2,2,2), col=c("dodgerblue2","navy","darkmagenta"))

################################################################################
################################################################################
# Goodness of Fit ##############################################################
################################################################################
################################################################################
# (!!) cf SLIDE 13/24 du cours (!!)

# Martingale Residuals (!!) Cours 14/24 (!!)
plot(
  predict(fit_ovarian),
  residuals(fit_ovarian, type = 'martingale'),
  xlab = 'Fitted values',
  ylab = 'Martingale residuals', 
  main = 'Residual Plot',
  las = 1
)
# Add a line for residual=0
abline(h=0, col='red')
# Fit a smoother for the points
lines(smooth.spline(predict(fit_ovarian), residuals(fit_ovarian, type='martingale')), col='blue')

# Alternatively Martingale Residuals (!!) Cours 14/24 (!!)
ggcoxdiagnostics(fit_ovarian, type = "martingale")

# Deviance residuals (!!) Cours 15/24 (!!)
ggcoxdiagnostics(fit_ovarian, type = "deviance")

# Schoenfeld residuals (!!) Cours 16/24 (!!)
cox_test_ovarian <- cox.zph(fit_ovarian)
plot(cox_test_ovarian)

# Log-negative-log-KM (!!) Cours 17/24 (!!)
rx_km <- survfit(Surv(futime, fustat) ~ rx, data = ovarian) # Only for categorical covariates !
ggsurvplot(rx_km,fun = "cloglog") # Are the curves parallel ? if YES : PH assumption is satisfied
cox_test_ovarian

################################################################################
################################################################################
# Stratified Cox PH Model ######################################################
################################################################################
################################################################################
fit_ovarian_full <- coxph(Surv(futime, fustat) ~ ., data=ovarian)

# Hint: use the ggforest function to visualize Hr and its CIs
ggforest(fit_ovarian_full, data=ovarian)

coxtest_HR_full <- cox.zph(fit_ovarian_full)

ecog.ps_km <- survfit(Surv(futime, fustat) ~ ecog.ps, data = ovarian)

ggsurvplot(ecog.ps_km,fun = "cloglog")

# Stratification not really necessary, but let us try to stratify by ecog.ps anyway

fit_ovarian_strata <-
  coxph(Surv(futime, fustat) ~ age + rx + resid.ds + strata(ecog.ps),
        data = ovarian)

summary(fit_ovarian_strata)


cox.zph(fit_ovarian_strata)

################################################################################
################################################################################
# Optional : manual calculation std schoenfeld residuals #######################
################################################################################
################################################################################
sch1 <- residuals(fit_ovarian,"scaledsch")
sresid <- residuals(fit_ovarian, "schoenfeld")

varnames <- names(fit_ovarian$coefficients)
nvar <- length(varnames)
ndead <- length(sresid)/nvar
sch2 <- sresid %*% fit_ovarian$var*ndead + rep(fit_ovarian$coefficients, each = nrow(sresid))

sch1
sch2

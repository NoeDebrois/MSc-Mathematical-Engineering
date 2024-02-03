################################################################################
################################ EXERCICE I ####################################
################################################################################

                      #################################
###################### !!!! CORRECTION AT THE END !!!! #########################
                      #################################

library(DepthProc)
library(progress)
library(survival)
milk_samples_1 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-02-11/milk_samples_1.Rds")
scaled = scale(milk_samples_1[,1:3])
scaled_milk_samples_1 = milk_samples_1
scaled_milk_samples_1[,1:3] = scaled
# An Irish farmer, Matthew O’Fountain, owns N = 382 cows and he is the best 
# milk-maker in Ireland. Nevertheless, Matthew O’Fountain is still not satisfied
# with this result, and he aims at becoming the best milk-maker in the whole world.
# In order to do so he must expand his market presence in other countries, and he 
# would like to know how to preserve milk quality once it is shipped from Ireland.
# Specifically, he is testing which type of pasteurization (Pasteurized or 
# Ultra-Pasteurized) ensures longer shelf life. He thus conducts the following 
# experiment: he collects N = 382 milk samples, one for each cow, and he monitors
# for T = 100 days whether the milk gets spoiled or not. The variable time 
# indicates at which day the sample quality deteriorated (spoiled=2) or if it was
# still fresh after 100 days (spoiled=1). Together with the pasteurization type, 
# he monitors three milk quality traits, namely Milk pH, Casein Micelle Size 
# (CMS), expressed in nm, and κ-casein (grams per liter). The resulting samples 
# are contained in the milk_samples_1.Rds file.

# 1. First off, Matthew O’Fountain is interested in knowing whether the type of 
# pasteurization alters the milk quality traits (Milk pH, Casein Micelle Size and
# κ-casein). By employing a permutation test on the standardized data1, and using
# as test statistic the maximum absolute difference between the sample 
# multivariate Tukey medians of the pasteurization types, check whether the milk
# samples differ in median in the two groups. Plot the permutational cumulative
# distribution function of the test statistic, report the p-value of the test and
# comment it.

# Create d1 with pasteurization_type = 1
t1 <- subset(scaled_milk_samples_1, pasteurization_type == "Pasteurized")
t1 <- t1[,1:3]

# Create d2 with pasteurization_type = 2
t2 <- subset(scaled_milk_samples_1, pasteurization_type == "Ultra-Pasteurized")
t2 <- t2[,1:3]

# Compute a proper test statistic (squared distance between two sample mean vectors)
t1.tukey.median <- depthMedian(t1,depth_params = list(method='Tukey'))
t2.tukey.median <- depthMedian(t2,depth_params = list(method='Tukey'))
n1 <- dim(t1)[1]
n2 <- dim(t2)[1]
n  <- n1 + n2
T20 <- max(abs(as.numeric((t1.tukey.median - t2.tukey.median))))

# Estimate the permutational distribution under H0
B <- 1000

# Create a progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = B
)

T2 <- numeric(B)

for(perm in 1:B) {
  # Random permutation of indexes
  t_pooled <- rbind(t1, t2)
  permutation <- sample(n)
  t_perm <- t_pooled[permutation,]
  t1_perm <- t_perm[1:n1,]
  t2_perm <- t_perm[(n1 + 1):n,]
  
  # Evaluation of the test statistic on permuted data
  t1.tukey.median_perm <- depthMedian(t1_perm,depth_params = list(method='Tukey'))
  t2.tukey.median_perm <- depthMedian(t2_perm,depth_params = list(method='Tukey'))
  T2[perm] <- max(abs(as.numeric((t1.tukey.median_perm - t2.tukey.median_perm))))
  pb$tick()
}

# Plot the permutational distribution under H0
hist(T2, xlim=range(c(T2, T20)), breaks=1000)
abline(v=T20, col=3, lwd=4)

plot(ecdf(T2))
abline(v=T20, col=3, lwd=4)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
p_val # We cannot reject H0 : so the milk samples do not differ in median in the two groups.

# 2. Compute the Kaplan-Meier estimation of the survival curves for the two 
# pasteurization types and plot it. Report median survival times and test if the
# time-to-event distributions of the two behavioral groups are equal via a 
# Log-rank test. Report the p-value and comment the result.

## Group 'Pasteurized' :
# Fitting the Kaplan-Meier Survival Curve
fit_KM_pasteurized <- survfit(Surv(time, spoiled) ~ 1, data = subset(milk_samples_1, pasteurization_type == "Pasteurized"))
# Plotting the Kaplan-Meier Survival Curve
plot(fit_KM_pasteurized, ylab = "Survival", xlab = "Time (in days)", mark.time = TRUE)
# Customizing the Kaplan-Meier Survival Curve Plot
ggsurvplot(fit_KM_pasteurized)
# Customizing the Kaplan-Meier Survival Curve Plot with additional features
ggsurvplot(fit_KM_pasteurized, surv.median.line = "hv", risk.table = TRUE, ggtheme = theme_bw())
# Median Survival Time Calculation
survminer::surv_median(fit_KM_pasteurized)
median_St_pasteurized <- fit_KM_pasteurized$time[fit_KM_pasteurized$surv <= 0.5][1]
median_St_pasteurized

## Group 'Ultra-Pasteurized' :
# Fitting the Kaplan-Meier Survival Curve
fit_KM_ultra_pasteurized <- survfit(Surv(time, spoiled) ~ 1, data = subset(milk_samples_1, pasteurization_type == "Ultra-Pasteurized"))
# Plotting the Kaplan-Meier Survival Curve
plot(fit_KM_ultra_pasteurized, ylab = "Survival", xlab = "Time (in days)", mark.time = TRUE)
# Customizing the Kaplan-Meier Survival Curve Plot
ggsurvplot(fit_KM_ultra_pasteurized)
# Customizing the Kaplan-Meier Survival Curve Plot with additional features
ggsurvplot(fit_KM_ultra_pasteurized, surv.median.line = "hv", risk.table = TRUE, ggtheme = theme_bw())
# Median Survival Time Calculation
survminer::surv_median(fit_KM_ultra_pasteurized)
median_St_ultra_pasteurized <- fit_KM_ultra_pasteurized$time[fit_KM_ultra_pasteurized$surv <= 0.5][1]
median_St_ultra_pasteurized

## Log-Rank Test :
fit_KM_by_pasteurization_type <- survfit(Surv(time, spoiled) ~ pasteurization_type, data = milk_samples_1)
fit_KM_by_pasteurization_type
plot(fit_KM_by_pasteurization_type, ylab = "Survival", xlab = "Time (in days)", col = 1:2, mark.time = TRUE)
ggsurvplot(fit_KM_by_pasteurization_type,risk.table.col = "strata",risk.table = T,surv.median.line = "hv", conf.int = TRUE) 

fit_log_rank <- survdiff(Surv(time, spoiled) ~ pasteurization_type, data = milk_samples_1)
fit_log_rank 
fit_log_rank$pvalue # We reject H0 ("S1(.)=S2(.)"). So S1(.)!=S2(.). Kaplan-Meier curves are statistically different.

# 3. Fit a suitable Cox model for long-term survival as a function of all the
# available covariates. Interpret the estimated coefficients for covariates Milk
# pH and pasteurization type, including a comment on statistical significance

fit_cox <- coxph(Surv(time, spoiled) ~ ., data = milk_samples_1)
summary(fit_cox)
# - Using p-values : p_val_pH is 0.968 so that we conclude that Milk PH is not statistically significant in predicting the survival time.
# p_val_Pasteurization_type is less than 2e-16 so that the pasteurization type is highly significant in predicting the survival time.
# - Using exp(coef) : exp(pH) = 1.01 so native PH has quite no effect in the hazard.
# exp(Pasteurization) = 8.68 so holding other variables constant, the hazard of experiencing the event is 8.685 times higher
# for "Pasteurized" milk samples compared to "Ultra-Pasteurized" ones.

# 4. Using the previously estimated Cox model, provide an estimate of the median
# survival time, under both pasteurization types, for the gold standard sample in
# terms of milk quality, for which κ-casein must be equal to 6 grams per liter, 
# CMS to 174 nm and Milk pH to 7.

# This data frame is passed to survfit() via the newdata argument:
new_data <- data.frame(kappa_casein = c(6,6),
                       Casein_micelle_size = c(174, 174),
                       Native_pH = c(7,7),
                       pasteurization_type = c("Pasteurized", "Ultra-Pasteurized"))

fit_new_data <- survfit(fit_cox, newdata = new_data)
fit_new_data

################################################################################
# CORRECTION ###################################################################
################################################################################

# 1
milk_samples_1 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-02-11/milk_samples_1.Rds")
milk_traits = scale(milk_samples_1[, 1:3])
n <- nrow(milk_traits)
table(milk_samples_1$pasteurization_type)

n1 = table(milk_samples_1$pasteurization_type)[1]
n2 = table(milk_samples_1$pasteurization_type)[2]

groups <- list()
groups$Pasteurized <- milk_traits[milk_samples_1$pasteurization_type=="Pasteurized",]
groups$`Ultra-Pasteurized` <- milk_traits[milk_samples_1$pasteurization_type=="Ultra-Pasteurized",]

median_pasturized = depthMedian(groups$Pasteurized, depth_params = list(method ='Tukey'))
median_ultra_pasturized = depthMedian(groups$`Ultra-Pasteurized`, depth_params = list(method = 'Tukey'))

t_stat = max(abs(median_ultra_pasturized - median_pasturized))

B = 1000
T_dist = numeric(B)
set.seed(2022)
pb = progress::progress_bar$new(total = B,format = " Processing [:bar] :percent eta: :eta")

for (index in 1:B) {
  perm = sample(1:n)
  milk_traits.p = milk_traits[perm, ]
  
  median_pasturized.p = depthMedian(milk_traits.p[1:n1, ], depth_params = list(method ='Tukey'))
  median_ultra_pasturized.p = depthMedian(milk_traits.p[(n1 + 1):n, ], depth_params = list(method ='Tukey'))
  
  T_dist[index] = max(abs(median_ultra_pasturized.p - median_pasturized.p))
  pb$tick()}

hist(T_dist, xlim = range(c(T_dist, t_stat)))
abline(v = t_stat, col = 3, lwd = 4)

plot(ecdf(T_dist))
abline(v = t_stat, col = 3, lwd = 4)

p_val <- sum(T_dist >= t_stat) / B
p_val
# p_val > 0.05 so : type of pasteurization does not alter the quality of the milk. 

# 2
fit <- survfit(Surv(time, spoiled == 2) ~ pasteurization_type, data = milk_samples_1)
ggsurvplot(
  fit,
  risk.table = TRUE,
  # Add risk table
  risk.table.col = "strata",
  # Change risk table color by groups
  surv.median.line = "hv",
  # Specify median survival
  ggtheme = theme_bw(),
  # Change ggplot2 theme
)

surv_median(fit)

log_rank_test <- survdiff(Surv(time, spoiled == 2) ~ pasteurization_type, data = milk_samples_1)
log_rank_test
# Kaplan-Meier curves are statistically different.

# 3
fit_cox <- coxph(Surv(time, spoiled) ~ ., data = milk_samples_1)
summary(fit_cox)
# Milk pH is not significant (see p_val).
# The HR for pasteurization_typePasteurized is exp(coef) = 8.685.
# Holding the other covariates constant, using a Pasteurized pasteurization_type
# increases the hazard by a factor of 8.67.

# 4 
standard_cow <- c(6,174,7)
new_df <- data.frame("kappa_casein" = standard_cow[1],
                     "Casein_micelle_size" = standard_cow[2],
                     "Native_pH" = standard_cow[3],
                     pasteurization_type=c("Pasteurized","Ultra-Pasteurized"))
fit_new <- survfit(fit_cox, newdata = new_df)
fit_new



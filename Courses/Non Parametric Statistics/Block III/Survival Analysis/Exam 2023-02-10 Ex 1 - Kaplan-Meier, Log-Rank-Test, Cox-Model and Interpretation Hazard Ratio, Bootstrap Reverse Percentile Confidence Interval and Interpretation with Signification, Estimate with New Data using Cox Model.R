################################################################################
################################# EXERCICE I ###################################
################################################################################
# Professor Franziska Iünz, a German biostatistician, is about to apply for a 
# grant. Despite her colleagues being quite positive about her success, she is 
# keen in knowing more about the grants submission process and the subsequent 
# publication of disruptive research. To this extent, she has collected data 
# about 77 previously funded proposals, contained in the df_1.Rds file. 
# In details, the following information is recorded:
# • mainpaper.event: main paper published (1=yes, 0=censored)
# • time.from.funding: time in years from funding until main paper was published (or censored)
# • members: number of investigators
# • funding.years: duration of funding in years
# • requested_money: a factor indicating whether SMALL FUNDING or LARGE FUNDING was requested
# Help Professor Iünz by solving the following tasks.
df_1 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2023-02-10/df_1.rds")
library(survival)
library(survminer)
# 1. Compute the Kaplan-Meier curves of the main paper publication event for the
# two funding types (according to the requested_money variable) and plot them. 
# Report the median times for publishing the main paper and test if the time-to-
# event distributions of the two groups are equal via a Log-rank test. Report 
# the p-value and comment the result.

# dans Surv : time, variable_censor :
fit <- survfit(Surv(time.from.funding, mainpaper.event) ~ requested_money, data = df_1)
ggsurvplot(
  fit,
  risk.table = FALSE,
  risk.table.col = "strata",
  surv.median.line = "hv",
  ggtheme = theme_bw(),
)

fit_info = survminer::surv_median(fit)
fit_info
fit_info[1,2] # median for LARGE FUNDING
fit_info[2,2] # median for SMALL FUNDING

## Log-Rank Test :
fit_log_rank <- survdiff(Surv(time.from.funding, mainpaper.event) ~ requested_money, data = df_1)
fit_log_rank 
fit_log_rank$pvalue # p_val > 0.05. We cannot reject H0 ("S1(.)=S2(.)"). So S1(.)=S2(.).
# Time-to-event distributions of the two groups have no significant difference.

# 2. Fit a suitable Cox model for time from funding until publication as a 
# function of all the available covariates. Interpret the estimated coefficient 
# for the funding.years covariate, including a comment on statistical significance.
fit_cox <- coxph(Surv(time.from.funding, mainpaper.event) ~ ., data = df_1) # . pour 'all the available covariates'
summary(fit_cox)
# - Using p-values : p_val_funding.years = 0.00767 < 0.05 so that we conclude 
# that funding.years is statistically significant in predicting the time (before publication) from funding.
# - Using exp(coef) : exp(funding.years) = 0.6383 so funding.years has an effect in the hazard.
# Longer study length increases the time to the main paper (the hazard of experiencing the event 
# of publication fo the paper is smaller (by a factor 0.6383)).

# 3. Employing a naive bootstrap approach, provide a reverse percentile 
# confidence interval for the funding.years coefficient of the Cox model 
# constructed in the previous exercise. Briefly describe the procedure and 
# comment on the result.
N <- nrow(df_1)

funding_coef <- fit_cox$coefficients[2]

B <- 1000
T.boot_funding = numeric(B)
set.seed(2022)
for(b in 1:B) {
  boot_id <- sample(x = 1:N,size = N,replace = TRUE)
  df_boot <- df_1[boot_id,]
  fit_cox_boot <- coxph(Surv(time.from.funding, mainpaper.event) ~ ., data = df_boot)
  T.boot_funding[b] = fit_cox_boot$coefficients[2]
}

alpha <- 0.05

right.quantile.funding <- quantile(T.boot_funding, 1 - alpha/2)
left.quantile.funding  <- quantile(T.boot_funding, alpha/2)

CI.RP.funding <-
  c(
    funding_coef - (right.quantile.funding - funding_coef),
    funding_coef,
    funding_coef - (left.quantile.funding- funding_coef))

names(CI.RP.funding)=c('lwr','pointwise','upr')
CI.RP.funding
# According to naive bootstrap, the coef is not significant at 5%
# Explanation :
# The confidence interval for the coefficient includes zero. When the confidence 
# interval contains zero, it generally indicates that the coefficient is not 
# statistically significant at the chosen significance level (e.g., 5%).
# In other words, if the coefficient were truly zero (no effect), you would 
# expect to obtain similar values through resampling, and the resulting 
# confidence interval would likely span zero. Therefore, the correction is 
# suggesting that, based on the naive bootstrap results, there is not enough
# evidence to reject the null hypothesis that the coefficient is equal to zero 
# at the 5% significance level.

# 4. Using the previously estimated Cox model, provide an estimate of the median
# time to publish the main paper for a 2-year long grant with 5 members involved
# and with LARGE FUNDING requested.

# This data frame is passed to survfit() via the newdata argument:
new_data <- data.frame(mainpaper.event = 1,
                       members = 5,
                       funding.years = 2,
                       requested_money = "LARGE FUNDING")

fit_new_data <- survfit(fit_cox, newdata = new_data)
fit_new_data
# 5.48





















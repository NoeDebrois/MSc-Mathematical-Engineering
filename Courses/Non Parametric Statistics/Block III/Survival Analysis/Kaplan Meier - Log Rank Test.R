# Loading Required Libraries
library(survival)
library(tidyverse)
library(survminer)

# The dataset being referred to as aml is part of the survival package in R. 
# This package is commonly used for survival analysis, which involves the 
# analysis of time-to-event data. The aml dataset specifically contains 
# information related to acute myelogenous leukemia (AML), a type of cancer that
# affects the blood and bone marrow.

# In the context of survival analysis, the dataset includes variables such as:
# - time : This variable represents the time until death occurs ;
# - status : This variable indicates whether the death has occurred (e.g., 1 for
# death, 0 for censored data) ;
# - x : "Maintained" or "Notmaintained" the chemiotherapy. 

# NB : interpretation of the Hazard Ratio in this case at the end. 

################################################################################
################################################################################
# KAPLAN-MEIER #################################################################
################################################################################
################################################################################

# Creating a Tibble from the 'aml' dataset in the 'survival' package
df <- as_tibble(survival::aml)
# Displaying the first 3 rows of the data
n_patients <- nrow(df)
head(df, 3)

# Defining the Survival Object
# Surv(df$time, df$status)

# Fitting the Kaplan-Meier Survival Curve ######################################
fit_KM <- survfit(Surv(time, status) ~ 1, data = df)

# Plotting the Kaplan-Meier Survival Curve #####################################
plot(fit_KM, ylab = "Survival", xlab = "Time (in days)", mark.time = TRUE)

# Customizing the Kaplan-Meier Survival Curve Plot
ggsurvplot(fit_KM)

# Customizing the Kaplan-Meier Survival Curve Plot with additional features ####
ggsurvplot(fit_KM, surv.median.line = "hv", risk.table = TRUE, ggtheme = theme_bw())

# Median Survival Time Calculation #############################################
survminer::surv_median(fit_KM)
median_St <- fit_KM$time[fit_KM$surv <= 0.5][1]
median_St

################################################################################
# LOG-CONFIDENCE-INTERVAL AT 95% (W/ TSIATIS FORMULA) & STD ERROR : ############
# (it's the default return when conf.type is not specified in "survfit") #######
################################################################################
# Standard Error with Tsiatis Formula
# USING VARIANCE OF ln(S_hat(t)) (!! PAGE 21/36 COURS !!) :
fit_KM$std.err # qui donne bien la même chose que :
(se_log_S <- sqrt(cumsum(fit_KM$n.event / (fit_KM$n.risk * (fit_KM$n.risk - fit_KM$n.event)))))

# Confidence Interval upper limit with Tsiatis Formula (!! PAGE 23/36 COURS !!)
pmin(1, fit_KM$surv * exp(qnorm(0.975) * se_log_S)) # qui donne la même chose que :
fit_KM$upper
# Checking if the calculated upper confidence limit matches survminer's result
all(dplyr::near(pmin(1, fit_KM$surv * exp(qnorm(0.975) * se_log_S)), fit_KM$upper))

# Confidence Interval Lower Limit with Tsiatis Formula (!! PAGE 23/36 COURS !!)
pmax(0, fit_KM$surv * exp(-qnorm(0.975) * se_log_S)) # qui donne la même chose que :
fit_KM$lower
# Checking if the calculated lower confidence limit matches survminer's result
all(dplyr::near(pmax(0, fit_KM$surv * exp(-qnorm(0.975) * se_log_S)), fit_KM$lower))

################################################################################
# GREENWOOD FORMULA : PLAIN-CONFIDENCE-INTERVAL AND STANDARD ERROR #############
# This time we specify 'conf.type = "plain"' ###################################
################################################################################
# Standard Error with Greenwood Formula
# USING THE STANDARD ERROR FORMULA (!! PAGE 22/36 COURS !!) :
model_gr <- survfit(Surv(time, status) ~ 1, data = df, conf.type = "plain")
se_S <- sqrt(model_gr$surv ^ 2 * cumsum(fit_KM$n.event / (fit_KM$n.risk * (fit_KM$n.risk - fit_KM$n.event))))

# Confidence Interval Lower Limit with Greenwood Formula (!! PAGE 23/36 COURS !!)
model_gr$lower
pmax(0, fit_KM$surv - qnorm(0.975) * se_S)

# Confidence Interval Upper Limit with Greenwood Formula (!! PAGE 23/36 COURS !!)
model_gr$upper
pmin(1, fit_KM$surv + qnorm(0.975) * se_S)

################################################################################
################################################################################
# LOG RANK TEST ################################################################
################################################################################
################################################################################

fit_KM_by_x <- survfit(Surv(time, status) ~ x, data = df)

plot(fit_KM_by_x, ylab = "Survival", xlab = "Time (in days)", col = 1:2, mark.time = TRUE)
ggsurvplot(fit_KM_by_x,risk.table.col = "strata",risk.table = T,surv.median.line = "hv", conf.int = TRUE) 

fit_log_rank <- survdiff(Surv(time, status) ~ x, data = df)
fit_log_rank

# x = maintained
n_j <- 23
d_j <- 2
n_kj <- 11

(E_X <- d_j*n_kj/n_j) # number of expected events in group k at tj* (!! PAGE 29/36 COURS !!)

# x = not maintained
n_j <- 23
d_j <- 2
n_kj <- 12

(E_X <- d_j*n_kj/n_j) # number of expected events in group k at tj* (!! PAGE 29/36 COURS !!)

################################################################################
################################################################################
# HAZARD RATIO #################################################################
################################################################################
################################################################################
# (!! PAGE 34/36 COURS !!)
hazard_ratio <-
  (fit_log_rank$obs[1] / fit_log_rank$exp[1]) /
  (fit_log_rank$obs[2] / fit_log_rank$exp[2])
hazard_ratio

# HR = 0.435 < 1 indicating that the risk of deaths for those patients
# who have maintained the chemotherapy is 0.435 times the risk of those who have
# not maintained it : maintaining the chemotherapy is a protective factor.


################################################################################
################################################################################
# OPTIONAL: MANUAL DERIVATIONS #################################################
################################################################################
################################################################################

# Manual calculation survival

df_KM <- df %>% 
  mutate(status=if_else(status==1,"death", "censored")) %>% 
  arrange(time)

df_KM_enriched <- df_KM %>% 
  count(time,status,name = "n.event")

n_risk_manual <- numeric(n_distinct(df_KM_enriched$time))
n_censored_manual <- df_KM %>% 
  group_by(time) %>% 
  summarise(n_censored_manual=sum(status=="censored")) %>% 
  pull(n_censored_manual)

n_event_manual <- df_KM %>% 
  group_by(time) %>% 
  summarise(n_event_manual=sum(status=="death")) %>% 
  pull(n_event_manual)

n_risk_manual[1] <- nrow(df_KM)
time_manual <- unique(df_KM_enriched$time)

for(i in 2:length(n_risk_manual)){
  previous_t <- df_KM_enriched %>% 
    filter(time==time_manual[i-1])
  n_death_previous_t <- pull(.data = filter(previous_t,status=="death"), n.event)
  n_death_previous_t <- ifelse(length(n_death_previous_t)==0,0,n_death_previous_t)
  n_censored_previous_t <- pull(.data = filter(previous_t,status=="censored"), n.event)
  n_censored_previous_t <- ifelse(length(n_censored_previous_t)==0,0,n_censored_previous_t)
  n_risk_manual[i] <- n_risk_manual[i-1]-n_death_previous_t-n_censored_previous_t
}

df_KM_fit <- distinct(df_KM, time) %>%
  mutate(
    n_risk = n_risk_manual,
    n_censored = n_censored_manual,
    n_event = n_event_manual,
    h_est = (n_event ) / n_risk,
    S_t_KM = cumprod(1 - h_est), # Kaplan-Meier estimator
    S_t_NA = exp(-cumsum(h_est)), # Nelson-Aalen estimator
    cumhaz=cumsum(h_est) # or -log(S_t_KM)
  )


dplyr::near(df_KM_fit$S_t_KM,fit_KM$surv) # survfit returns KM estimator 


plot(df_KM_fit$time,df_KM_fit$S_t_KM,type="s") 
lines(df_KM_fit$time,df_KM_fit$S_t_NA,type="s",col="red") # very similar estimates

# Log rank test manual

tidy_km <- broom::tidy(fit_KM_by_x)

log_rank_calculation_long <- tidy_km %>% 
  select(time,strata,n.risk,n.event,n.censor) %>% 
  complete(time, strata,fill = list(n.event=0,n.censor=0)) %>% 
  group_by(strata) %>% 
  fill(n.risk,.direction = "updown") %>% 
  group_by(time, strata) %>% 
  summarise(d_kj=sum(n.event), n_kj=sum(n.risk)) %>% 
  ungroup() %>% 
  add_count(time, wt = d_kj,name = "d_j") %>% 
  add_count(time, wt = n_kj,name = "n_j") %>% 
  group_by(time,strata)


log_rank_approx_manual <- log_rank_calculation_long %>% 
  summarise(e_kj=d_j*n_kj/n_j, d_kj=d_kj) |> 
  group_by(strata) %>% 
  summarise(Observed=sum(d_kj), Expected=sum(e_kj)) %>% 
  mutate(`(O-E)^2/E`=(Observed-Expected)^2/Expected)

log_rank_manual <- log_rank_calculation_long %>% 
  summarise(E_dkj=d_j*n_kj/n_j,V_dkj=(n_kj*(n_j-n_kj)*d_j*(n_j-d_j))/(n_j^2*(n_j-1))) |> 
  group_by(strata) |> 
  summarise(E_k=sum(E_dkj),V_k=sum(V_dkj))

chisq_manual <- (fit_log_rank$obs-fit_log_rank$exp)^2/fit_log_rank$var[1,1]
pchisq(chisq_manual,df = 1,lower.tail = FALSE)

log_rank_manual |> 
  left_join(log_rank_approx_manual) |> 
  mutate(`(O-E)^2/V`=(Observed-Expected)^2/V_k)


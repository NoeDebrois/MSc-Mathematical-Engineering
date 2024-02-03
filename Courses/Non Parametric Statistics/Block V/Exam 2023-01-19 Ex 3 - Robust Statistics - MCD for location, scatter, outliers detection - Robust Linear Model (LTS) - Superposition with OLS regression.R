################################################################################
################################ EXERCICE III ##################################
################################################################################
# Dr. Matteus Fontansen is afraid some of the measurements he has collected may 
# have been wrongly recorded. To this extent, he is interested in performing 
# some statistical analyses using robust methods.

df_3 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2023-01-19/df_3.rds")

# 1. Focusing on the Gentoo species, compute the Minimum Covariance Determinant 
# estimator for the flipper_length_mm and bill_length_mm variables contained in 
# the df_3.rds dataset. Consider 1000 subsets for initializing the algorithm and
# set the sample size of H, the subset over which the determinant is minimized, 
# equal to 100. Report the (reweighted) MCD estimates of location and scatter. 
# Define a vector ind_out_MCD of row indexes identifying the samples (within the
# Gentoo subpopulation) that are outliers according to the MCD call and report it.

# Extract the relevant columns (body_mass_g, flipper_length_mm, bill_length_mm) 
# from the dataset
X_mcd <- df_3[df_3$species == "Gentoo", 2:3]

# Get the number of rows in the dataset
N <- nrow(X_mcd)

# Set a seed for reproducibility
set.seed(2022)

# Fit the Minimum Covariance Determinant (MCD) estimator
fit_MCD <- covMcd(x = X_mcd,
                  alpha = 100 / N, # or alpha = (N - 19) / N since N = 119
                  nsamp = 1000)

# Report the raw MCD estimate of location
fit_MCD$center

# Report the raw MCD estimate of scatter (covariance matrix)
fit_MCD$cov

# Identify outliers using the MCD estimator and create a vector of row indexes
# Attention il existe plusieurs critères pour identifier les outliers ! (cf Ex 3 11/07/2022)
# Ici : on veut décrire comme outliers ceux qui ne sont pas sélectionnés par les "best".
# ICI : "ceux qui ne sont pas les best" : pas les "vrais" outliers.
# FAIRE L'AUTRE SI PAS CLAIR !
ind_out_MCD <- setdiff(1:N, fit_MCD$best)

# Report the vector of row indexes identifying the outliers
ind_out_MCD

# 2. Dr. Matteus Fontansen believes some Gentoo penguins have been wrongly 
# labeled as to belong to the Adelie and Chinstrap species. Employing the 
# (reweighted) MCD estimates obtained in the previous exercise, check if some
# Adelie and Chinstrap penguins have been wrongly labeled by computing robust 
# squared Mahalanobis distances using χ2_{2,0.975} as cut-off value, with χ2_{p,α} 
# denoting the α-quantile of a χ2 distribution with p degrees of freedom.

X_no_gentoo <- df_3[df_3$species != "Gentoo", 2:3]
ind_wrongly_labeled_obs <-
  which(
    mahalanobis(
      x = X_no_gentoo,
      center = fit_MCD$center,
      cov = fit_MCD$cov
    ) <= qchisq(p = .975, df = 2)
  )
length(ind_wrongly_labeled_obs)
# 11 obs may have been wrongly labeled!

plot(df_3$flipper_length_mm, df_3$bill_length_mm, type = "n")
points(X_no_gentoo$flipper_length_mm,
       X_no_gentoo$bill_length_mm,
       cex = .2) # plot the no-gentoo df
points(X_no_gentoo$flipper_length_mm[ind_wrongly_labeled_obs],
       X_no_gentoo$bill_length_mm[ind_wrongly_labeled_obs]) # plot the wrongly-labeled no-gentoo df
X_gentoo = X_mcd
points(
  X_gentoo$flipper_length_mm,
  X_gentoo$bill_length_mm,
  cex = .2,
  col = "blue"
) # plot the gentoo df

# 3. Since the species variable could be unreliable, Dr. Matteus Fontansen asks 
# you to build a robust linear model for the entire dataset, to regress body mass
# on the bill length using a Least Trimmed Squares (LTS) approach, setting the 
# hyperparameter α = 0.75. Provide a plot of the regression line, flagging the 
# units (i.e., color them in red in the scatterplot) whose squared residuals were
# not minimized in the LTS call. Superimpose the fit we would obtain if we were
# to use OLS and comment accordingly.

# Fit a robust linear model using Least Trimmed Squares (LTS)
fit_lts <- ltsReg(body_mass_g ~ bill_length_mm, alpha = 0.75, data = df_3)

# Scatterplot with color-coded units based on whether their squared residuals were minimized
with(df_3, plot(bill_length_mm, body_mass_g, col = ifelse(1:nrow(df_3) %in% fit_lts$best, "black", "red")))

# Add the LTS regression line to the scatterplot
abline(fit_lts)

# Superposition with the fit we would obtain if we were to use OLS
abline(lm(body_mass_g ~ bill_length_mm, data = df_3), col="blue")
# Penguins with long bills but low body mass can bias the estimates 
# (we know they actually are Chinstrap)









################################################################################
################################ EXERCISE III ##################################
################################################################################
library(robustbase)
library(MASS)

milk_samples_3 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-02-11/milk_samples_3.Rds")

# Matthew O’Fountain has recently become passionate about robust statistics: he
# therefore would like to exploit these modern statistical methods to further 
# analyze his milk samples.

# 1. Compute the Minimum Covariance Determinant estimator for the Milk pH, Casein
# Micelle Size (CMS), expressed in nm, and κ-casein (grams per liter) variables
# contained in the milk_samples_3.Rds dataset. Consider 1000 subsets for
# initializing the algorithm and set the sample size of H, the subset over which
# the determinant is minimized, equal to 341. Report the raw MCD estimates of
# location and scatter. Define a vector ind_out_MCD of row indexes identifying
# the milk samples that are outliers according to the MCD call and report it.

# Extract the relevant columns (Milk pH, Casein Micelle Size, κ-casein) from the dataset
X_mcd <- milk_samples_3[, 1:3]

# Get the number of rows in the dataset
N <- nrow(milk_samples_3)

# Set a seed for reproducibility
set.seed(2022)

# Fit the Minimum Covariance Determinant (MCD) estimator
fit_MCD <- covMcd(x = X_mcd, alpha = 341 / N, nsamp = 1000) # N - 41 = 341... donc alpha = 341 / N = (N-41) / N

# Report the raw MCD estimate of location
fit_MCD$raw.center

# Report the raw MCD estimate of scatter (covariance matrix)
fit_MCD$raw.cov

# Identify outliers using the MCD estimator and create a vector of row indexes
# Attention il existe plusieurs critères pour identifier les outliers ! (cf Ex 3 11/07/2022)
# Ici : on veut décrire comme outliers ceux qui ne sont pas sélectionnés par les "best".
# ICI : "ceux qui ne sont pas les best" : pas les "vrais" outliers.
# FAIRE L'AUTRE SI PAS CLAIR !
ind_out_MCD <- setdiff(1:N, fit_MCD$best)

# Report the vector of row indexes identifying the outliers
ind_out_MCD

# 2. Build a robust linear model to regress κ-casein on the milk absorbance at
# wavenumber 280 cm−1 using a Least Trimmed Squares (LTS) approach, setting the
# hyperparameter α = 0.75. Provide a plot of the regression line, flagging the 
# units (i.e., color them in red in the scatterplot) whose squared residuals were
# not minimized in the LTS call.

# Fit a robust linear model using Least Trimmed Squares (LTS)
fit_lts <- ltsReg(kappa_casein ~ wave_280, alpha = 0.75, data = milk_samples_3)

# Scatterplot with color-coded units based on whether their squared residuals were minimized
with(milk_samples_3, plot(wave_280, kappa_casein, col = ifelse(1:nrow(milk_samples_3) %in% fit_lts$best, "black", "red")))

# Add the regression line to the scatterplot
abline(fit_lts)


# 3. Provide the outlier map for the robust linear model estimated in the
# previous exercise. Are bad leverage points present in the dataset according
# to the diagnostic plot?

# Plot diagnostics for the LTS model (provided by robustbase)
# plot(fit_lts) # to have all plots !
plot(fit_lts, which="rdiag")
# It seems that no bad leverage points are present, only vertical outliers and good leverage points.












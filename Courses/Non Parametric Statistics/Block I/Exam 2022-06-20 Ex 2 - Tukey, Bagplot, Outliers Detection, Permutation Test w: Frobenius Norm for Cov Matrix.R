################################################################################
################################ EXERCISE II ###################################
################################################################################
# Since becoming a parent is not stressful enough, Mr András Kaponzyi is also
# considering buying a house in Budapest. In order to do so, he has collected
# information on 505 properties on sale in the capital, recording the price 
# (in mln Hungarian forint), the Dimension in m^2 and the area in which the 
# property is located (either Buda or Pest). The resulting samples are contained
# in the df_2.Rds file. He now requires you to:
library(DepthProc)
library(aplpack)

df_2 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-06-20/df_2.rds")

N = nrow(df_2)

# 1. Provide the Tukey median of the properties located in Pest.
tukey_median <- depthMedian(df_2[,1:2], depth_params = list(method='Tukey'))
tukey_median

# 2. The real estate agent helping Mr András Kaponzyi with the house hunting
# told him that the Pest area is outliers-free, while some peculiar properties
# may be found in the Buda area. Therefore, for the properties in the Buda area,
# plot a bagplot of the two collected variables and determine a vector of row 
# indexes identifying the samples that are outliers according to such procedure.
# What characteristics do the outlying properties showcase?

buda <- df_2[df_2$Area=="Buda",1:2]
bagplot_buda <- bagplot(buda)

# Extract outliers :
outlying_obs_buda <- bagplot_buda$pxy.outlier
outlying_obs_buda # values. Now we want a vector of row indexes...
ind_outlying_obs_buda <- which(apply(buda,1,function(x) all(x %in% outlying_obs_buda))) 
ind_outlying_obs_buda
unname(ind_outlying_obs_buda) # Juste pour l'affichage des bons index !
# They are very expensive and very big !

# 3. After having eliminated the properties identified as outliers in the Buda
# area, use a permutation test to check whether the variance covariance matrix
# of Price and Dimension differ in the two areas (Buda and Pest). To do so, 
# consider the following procedure: firstly center the Buda (no outliers) and 
# Pest observations by subtracting their respective sample means. Then, to 
# perform the test, use as test statistic the squared Frobenius distance
# between the sample covariance matrices of the centered data. Specify
# assumptions, the null and the alternative hypothesis you are testing, provide
# the histogram of the permuted distribution of the test statistic, its 
# cumulative distribution function and the p-value, commenting the results.

buda_no_out = buda[-ind_outlying_obs_buda,]
buda_no_out_centered = scale(buda_no_out, center = TRUE, scale = FALSE)
pest_no_out = df_2[df_2$Area=="Pest",1:2]
pest_no_out_centered = scale(pest_no_out, center = TRUE, scale = FALSE)

# Permutation test
S_buda <- cov(buda_no_out_centered)
S_pest <- cov(pest_no_out_centered)
T20 <- norm(S_buda-S_pest,"F")^2

x_pooled = rbind(buda_no_out_centered, pest_no_out_centered) 

# Estimating the permutational distribution under H0
B <- 1000
T2 <- numeric(B)
n = nrow(x_pooled)
n1 = nrow(buda_no_out_centered)
set.seed(2022)
for (perm in 1:B) {
  permutation <- sample(1:n)
  x_perm <- x_pooled[permutation,]
  x1_perm <- x_perm[1:n1,]
  x2_perm <- x_perm[(n1 + 1):n,]
  S_1 <- cov(x1_perm)
  S_2 <- cov(x2_perm)
  T2[perm] <- norm(S_1-S_2,"F")^2
}

# Visualize permutational distribution of T
hist(T2, xlim=range(c(T2, T20)), breaks=30)
abline(v=T20, col=3, lwd=2)

# Plot empirical cumulative distribution function (ECDF)
plot(ecdf(T2))
abline(v=T20, col=3, lwd=2)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
print("The P-Value for this test is :")
p_val
# H0 : the variance covariance matrix of Price and Dimension don't differ in the
# two areas (S_buda = S_pest).
# H1 : the variance covariance matrix of Price and Dimension differ in the two 
# areas (S_buda != S_pest).
# It's quite big, we can't reject H0. So The variance covariance matrix of Price
# and Dimension don't differ in the two areas.











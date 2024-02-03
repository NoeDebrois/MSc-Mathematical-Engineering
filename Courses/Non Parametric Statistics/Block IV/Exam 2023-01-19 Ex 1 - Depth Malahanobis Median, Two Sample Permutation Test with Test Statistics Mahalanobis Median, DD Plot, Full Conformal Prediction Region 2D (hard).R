################################################################################
################################# EXERCISE I ###################################
################################################################################
df_1 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2023-01-19/df_1.Rds")

# Dr. Matteus Fontansen is a Norwegian ecologist leading an integrative study of
# the population structure of Pygoscelis penguins along the western Antarctic 
# Peninsula. To this aim, he has collected information about 333 penguins from 
# three different species, namely Adélie, Chinstrap and Gentoo. 
# Particularly, he is interested in knowing how the flipper length 
# (flipper_length_mm) and the bill length (bill_length_mm), both measured in 
# millimeters, vary across species. The resulting samples are contained in the 
# df_1.Rds file. Help Dr. Matteus Fontansen by solving the following tasks:

# 1. Provide the Mahalanobis medians for the three penguin species and 
# superimpose them to the scatterplot of flipper length vs bill length.
df_split <- split(x = df_1[, 1:2], f = df_1$species)
(maha_medians <- lapply(df_split, function(x)
  depthMedian(x, depth_params = list(method = 'Mahalanobis'))))

plot(df_1[, 1:2], col = df_1$species, cex = .5) 
for (i in 1:length(maha_medians)) {
  points(
    x = maha_medians[[i]][1],
    y = maha_medians[[i]][2],
    col = i,
    cex = 2,
    pch = "x"
  ) }

# 2. Test the equality of the theoretical Mahalanobis medians for the Adelie and
# Chinstrap species by performing a two-sample permutation test using as a test
# statistics the squared euclidean distance between the two sample Mahalanobis 
# medians. Please describe briefly the properties of permutation tests and 
# present the empirical cumulative distribution function of the permutational 
# test statistic as well as the p-value for the test.
# Extract two subsets for testing
t1 <- df_1[df_1$species == "Adelie",]
t2 <- df_1[df_1$species == "Chinstrap",] 

# Compute a proper test statistic (squared distance between two sample mean vectors)
n1 <- dim(t1)[1]
n2 <- dim(t2)[1]
n  <- n1 + n2

# Test statistics
t1_maha_median <- depthMedian(t1[,1:2], depth_params = list(method='Mahalanobis'))
t2_maha_median <- depthMedian(t2[,1:2], depth_params = list(method='Mahalanobis'))
T20 <- as.numeric((t1_maha_median - t2_maha_median) %*% (t1_maha_median - t2_maha_median))

# Estimate the permutational distribution under H0
B <- 1000
T2 <- numeric(B)

for(perm in 1:B) {
  # Random permutation of indexes
  t_pooled <- rbind(t1, t2)
  permutation <- sample(n)
  t_perm <- t_pooled[permutation,]
  t1_perm <- t_perm[1:n1,]
  t2_perm <- t_perm[(n1 + 1):n,]
  
  # Evaluation of the test statistic on permuted data
  t1_perm.maha.median = depthMedian(t1_perm[,1:2],depth_params = list(method='Mahalanobis'))
  t2_perm.maha.median = depthMedian(t2_perm[,1:2],depth_params = list(method='Mahalanobis'))
  T2[perm] <- as.numeric((t1_perm.maha.median - t2_perm.maha.median) %*% (t1_perm.maha.median - t2_perm.maha.median))
}

# Plot the permutational distribution under H0
hist(T2, xlim=range(c(T2, T20)), breaks=1000)
abline(v=T20, col=3, lwd=4)

plot(ecdf(T2))
abline(v=T20, col=3, lwd=4)

# Calculate p-value
p_val <- sum(T2 >= T20) / B
p_val 
# p_val < 0.05.
# We have to reject H0 so : the theoretical Mahalanobis medians for the Adelie and
# Chinstrap species are not equal.

# 3. Employing the Mahalanobis depth, build a DD-plot for the empirical 
# distributions of the Adelie vs Chinstrap species. Can you conclude that the 
# two empirical distributions are identical?

Adelie <- df_1[df_1$species == "Adelie",]
Chinstrap <- df_1[df_1$species == "Chinstrap",]
ddPlot(x = Adelie[,1:2],y = Chinstrap[,1:2],depth_params = list(method='Mahalanobis'))
# The two empirical distributions are clearly not identical

# 4. Provide a Full Conformal 1 − α = 90% prediction region of flipper and bill
# lengths for a new Gentoo penguin, using the squared euclidean distance 
# between the new data point and the sample Mahalanobis median of the augmented 
# data set as non-conformity measure. After having discussed the theoretical
# properties of the prediction region, provide a plot of it.

# Extract Gentoo penguin data from the data frame
data_predict = df_split$Gentoo

# Set the number of grid points for the prediction region
n_grid = 20
# Set the grid factor for extending the range of the grid
grid_factor = 0.25
# Set the significance level (1 - alpha) for the prediction region
alpha = .1
# Get the number of observations in the data
n = nrow(data_predict)

# Calculate the range of the x and y coordinates in the data
range_x = range(data_predict[, 1])[2] - range(data_predict[, 1])[1]
range_y = range(data_predict[, 2])[2] - range(data_predict[, 2])[1]

# Create grid points for testing the prediction region
test_grid_x = seq(
  min(data_predict[, 1]) - grid_factor * range_x,
  max(data_predict[, 1]) + grid_factor * range_x,
  length.out = n_grid
)

test_grid_y = seq(
  min(data_predict[, 2]) - grid_factor * range_y,
  max(data_predict[, 2]) + grid_factor * range_y,
  length.out = n_grid
)

# Create a grid of test points for the prediction region
xy_surface = expand.grid(test_grid_x, test_grid_y)
colnames(xy_surface) = colnames(data_predict)

# Define a function to compute the non-conformity measure for each test point
wrapper_multi_conf = function(test_point) {
  newdata = rbind(test_point, data_predict)
  newmedian = depthMedian(newdata, depth_params = list(method = 'Mahalanobis'))
  depth_surface_vec = rowSums(t(t(newdata) - newmedian) ^ 2)
  sum(depth_surface_vec[-1] >= depth_surface_vec[1]) / (n + 1)
}

# Apply the non-conformity measure function to each test point in the grid
pval_surf = pbapply::pbapply(xy_surface, 1, wrapper_multi_conf)
data_plot = cbind(pval_surf, xy_surface)

# Identify the points outside the prediction region based on the significance level
p_set = xy_surface[pval_surf > alpha, ]
poly_points = p_set[chull(p_set), ]

# Plot the prediction region
ggplot() +
  geom_raster(data = data_plot, aes(flipper_length_mm, bill_length_mm, fill = pval_surf)) +
  geom_point(data = data.frame(data_predict), aes(flipper_length_mm, bill_length_mm)) +
  geom_polygon(
    data = poly_points,
    aes(flipper_length_mm, bill_length_mm),
    color = 'red',
    size = 1,
    alpha = 0.01
  )



























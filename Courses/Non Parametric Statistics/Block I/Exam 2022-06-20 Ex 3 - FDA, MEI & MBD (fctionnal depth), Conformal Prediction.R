################################################################################
################################ EXERCISE III ##################################
################################################################################

# SAME EXERCISE AS 21/01/2022 EX3 (but different data) !!

# Mr András Kaponzyi is now interested in deeper understanding how the height of
# his soon-to-be-born daughter might look like in the years to come. To do so, 
# he asks you to analyze (a simplified version of) the famous Berkeley Growth 
# Study dataset. The df_3.Rds file contains a list of two components.
df_3 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-06-20/df_3.rds")
library(roahd)

# • height_values : a 54 by 18 numeric matrix giving the heights in centimeters 
# of 54 girls.
# • age_range : a numeric vector of length 18 giving the ages at which the
# heights were measured.

# 1. By suitably defining a functional data object, plot the resulting height 
# for the N = 54 girls contained in the sample. Then, compute the median age 
# using the modified band depth and superimpose it to the previous plot.

## PLOT THE RESULTING SPECTRA :
# The way roahd package represents functional objects is by providing an evenly
# spaced grid I = [t0, t1, ... , tP−1] (tj − tj−1 = h > 0, ∀j = 1, ... , P − 1).
grid <- df_3$age_range # !! we need only the values of the grid ! directly in age_range
f_data <- fData(grid, df_3$height_values) # functional data (using roahd package)
plot(f_data) # it works !
class(f_data) # print the class : it's functional !

## COMPUTE THE MEDIAN SPECTRUM USING MBD :
median_curve <- median_fData(fData = f_data, type = "MBD") # still an fData object

## SUPERPOSITION :
lines(grid,median_curve$values)

# 2. Produce a functional boxplot and use it to identify the outlying curves. 
# Are there amplitude outliers? 
# Without resorting to the outliergram function, manually build a scatterplot of 
# Modified Epigraph Index (MEI) vs Modified Band Depth (MBD), and flag the curves
# (i.e., color them in red in the scatterplot) whose di ≥ Q_3(d) + 2*IQR(d), with
# and :
# di = a_0 + a_1*MEI(f) + a_2*N^2*MEI2(f) − MBD(f), i=1,...,N 
# a_0 = a_2 = −2/(N(N−1)), 
# a_1 = 2*(N+1)/(N−1).
N <- f_data$N

## FUNCTIONAL BOXPLOT AND IDENTIFICATION OF THE OUTLYING CURVES :
fb_plot <- fbplot(f_data, main="Magnitude outliers")
fb_plot # prints all the info
outliers = c(fb_plot$ID_outliers)
outliers # yes there is one amplitude outlier !

## SCATTERPLOT OF MEI VS MBD AND FLAG OF THE CURVES BASED ON THE CONDITION :
# outgr_plot <- outliergram(f_data,Fvalue = 2) # we don't have the right to use 
# it here !
MEI_f <- MEI(f_data)
MBD_f <- MBD(f_data)

a_0 <- a_2 <- -2/(N*(N-1))
a_1 <- 2*(N+1)/(N-1)

d_manual <- a_0+a_1*MEI_f+a_2*f_data$N^2*MEI_f^2-MBD_f

ID_outliers_manual <- which(d_manual > quantile(d_manual,probs = .75) + 2 * IQR(x = d_manual))
ID_outliers_manual

plot(MEI_f, MBD_f, col=ifelse(1:N %in% ID_outliers_manual,"red","black"), pch=19, main="Outliers manual")
# The automatic way (but forbidden here ha ha) :
outgr_plot <- outliergram(f_data,Fvalue = 2) # check that we get the same outliers ! (yes)

# 3. Build a split conformal prediction band for a height curve using the 
# previously computed functional median as point-wise predictor, the sup norm of
# absolute value of the functional residuals from the functional median as a NCM
# and the identity function as the modulation function. Consider a 70-30 percent
# split for calibration and validation sets, respectively.

# Définir le niveau de confiance alpha
alpha = 0.1

# Nombre total d'échantillons
n = N

# Fixer la graine aléatoire pour la reproductibilité
set.seed(2022)

# Sélection aléatoire de 2/3 des échantillons pour le training
i1 = sample(1:n, 2*n/3)

# Définir l'ensemble d'apprentissage (training set) et l'ensemble de test (complémentaire)
t_set = df_3$height_values[i1,]
c_set = df_3$height_values[-i1,]

# Calculer la médiane fonctionnelle sur l'ensemble d'apprentissage
# (attention : la faire sur l'ensemble d'apprentissage ou bien sur tout df_3$$height_values ?)
mu <- median_fData(fData = fData(grid, t_set), type = "MBD")
mu = mu$values

# Calculer les résidus fonctionnels entre l'ensemble de test et la médiane fonctionnelle
res=abs(sweep(x = c_set,MARGIN = 2,STATS = t(mu),"-"))
# La fonction sweep est utilisée pour effectuer des opérations du type X_{ij}−STATS_j
# où : X est la matrice, STATS est un vecteur statistique et j est une "dimension"
# spécifiée par l'argument MARGIN (colonne, ou ligne...).

# Calculer la Non Conformity Measure (NCM) pour chaque variable
ncm = apply(res, 2, max) # appliquer le max sur chaque colonne (1 : ligne, 2 : colonne) 

# Trier les NCM
ncm_sort = c(sort(ncm), Inf)

# Calculer la distance pour obtenir la largeur de la bande de prédiction conforme
d = ncm_sort[ceiling((n/3 + 1) * (1 - alpha))]

# Tracer la bande de prédiction en utilisant matplot
matplot(cbind(t(mu), t(mu) + d, t(mu) - d), type = 'l')












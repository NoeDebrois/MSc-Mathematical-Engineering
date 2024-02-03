################################################################################
############################### EXERCICE III ###################################
################################################################################
library(roahd)

milk_samples_3 <- readRDS("/Users/noedebrois/Desktop/Politecnico/Non Parametric Statistics/LECTURES/PREVIOUS EXAMS-20231223/Data/2022-01-21/milk_samples_3.Rds")

# Andrew O’Cappor has recently become passionate about functional data analysis (FDA):
# Dr. Alexander House is then ready to unleash the power of the MilkoScan FT6000 milk
# analyzer, producing MIRS transmittance spectra in the entire mid-infrared light region
# (wavenumbers 70 and 300 cm−1 of the previous exercise were just two points of this region).
# The MIRS for the N = 382 milk samples are contained in the milk_samples_3.Rds file.
# Andrew O’Cappor would like to exploit these modern statistical methods to further analyze
# his milk samples. To this extent, he requires you to:

# 1. By suitably defining a functional data object, plot the resulting spectra for the N = 382
# milk samples. Then, compute the median spectrum using the modified band depth and superimpose
# it to the previous plot.

## PLOT THE RESULTING SPECTRA :
# The way roahd package represents functional objects is by providing an evenly
# spaced grid I = [t0, t1, ... , tP−1] (tj − tj−1 = h > 0, ∀j = 1, ... , P − 1).
grid <- 1:ncol(milk_samples_3)
f_data <- fData(grid, milk_samples_3) # functional data (using roahd package)
plot(f_data) # it works !
class(f_data) # print the class : it's functional !

## COMPUTE THE MEDIAN SPECTRUM USING MBD :
median_curve <- median_fData(fData = f_data, type = "MBD") # still an fData object

## SUPERPOSITION :
lines(grid,median_curve$values)

# 2. Produce a functional boxplot and use it to identify the outlying curves. 
# Are there amplitude outliers? 
# Without resorting to the outliergram function, manually build a scatterplot of 
# Modified Epigraph Index (MEI) vs Modified Band Depth (MBD), and flag the spectra
# (i.e., color them in red in the scatterplot) whose di ≥ Q_3(d) + 2*IQR(d), with
# and :
# di = a_0 + a_1*MEI(f) + a_2*N^2*MEI2(f) − MBD(f), i=1,...,N 
# a_0 = a_2 = −2/(N(N−1)), 
# a_1 = 2*(N+1)/(N−1).

## FUNCTIONAL BOXPLOT AND IDENTIFICATION OF THE OUTLYING CURVES :
fbplot <- fbplot(f_data, main="Magnitude outliers") # Prints all the info...
outliers = c(fbplot$ID_outliers) 
outliers # yes there are !

## SCATTERPLOT OF MEI VS MBD AND FLAG SOME SPECTRA BASED ON THE CONDITION :
modified_band_depth <- MBD(Data = f_data) # MBD
modified_epigraph_index <- MEI(Data = f_data) # MEI 
plot(modified_epigraph_index, modified_band_depth) # MBD vs MEI

N = length(milk_samples_3[,1])
a02 = -2 / (N*(N-1))
a1 = 2*(N+1)/(N-1)
di = a02 + a1*modified_epigraph_index + a02*N^2*modified_epigraph_index^2 - modified_band_depth
# flag = di(which(di >= )) # Je ne sais pas ce que sont Q3 et IQR... pas de correction...

# 3. Build a split conformal prediction band for a new milk spectrum using the 
# previously computed functional median as pointwise predictor, the sup norm of 
# absolute value of the functional residuals from the functional median as a NCM
# and the identity function as the modulation function

# Définir le niveau de confiance alpha
alpha = 0.1

# Nombre total d'échantillons
n = nrow(milk_samples_3)

# Fixer la graine aléatoire pour la reproductibilité
set.seed(2022)

# Sélection aléatoire de la moitié des échantillons
i1 = sample(1:n, n/2)

# Définir l'ensemble d'apprentissage (training set) et l'ensemble de test (complémentaire)
t_set = milk_samples_3[i1,]
c_set = milk_samples_3[-i1,]

# Calculer la médiane fonctionnelle sur l'ensemble d'apprentissage
mu <- median_fData(fData = fData(grid, t_set), type = "MBD")
mu = mu$values

# Calculer les résidus fonctionnels entre l'ensemble de test et la médiane fonctionnelle
res = c_set - t(mu) # t(mu) = transposée de mu

# Calculer la Non Conformity Measure (NCM) pour chaque variable
ncm = apply(res, 2, max)

# Trier les NCM
ncm_sort = c(sort(ncm), Inf)

# Calculer la distance pour obtenir la largeur de la bande de prédiction conforme
d = ncm_sort[ceiling((n/2 + 1) * (1 - alpha))]

# Tracer la bande de prédiction en utilisant matplot
matplot(cbind(t(mu), t(mu) + d, t(mu) - d), type = 'l')

















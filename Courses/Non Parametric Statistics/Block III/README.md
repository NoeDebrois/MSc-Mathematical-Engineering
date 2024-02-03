# Block III : Regression :

The method is always to : build a model on the dataset, then build new data, then predict on these new data, then compute the se bands, then plot : the points, the regression line and the se bands. If you want you also have to precise some knots.

## LAB 11 : polynomial and step regression (2 first ways to move beyond linearity)
- Polynomial Regression : the simplest way to extend linear regression to a nonlinear one it to employ polynomial functions. 
    - In R : “m_quad=lm(wage ~ age + I(age^2))”
- Step Regression : 
    - Using polynomial functions of the features as predictors in a linear model imposes a global structure on the non-linear function of X. 
    - We can instead use step regression in order to avoid imposing such a global structure. The idea is to break the range of X into intervals, and treat the fact that x is in that interval as a factor. 
    - In R (with the “cut”) : 
        - If you just want to cut in  (e.g) 4 equal intervals : “m_cut=lm(prestige ~ cut(income,breaks=4),data = Prestige)”
        - If you want to precise your interval cuts :“m_cut=lm(prestige ~ cut(income,breaks = c(-Inf,10000,max(income))), data=Prestige)”
        - See p59/171 of “Global_Labs.pdf”

## LAB 11 bis : local regression and splines (2 other methods that provide non-linear fitting)
- Local regression : it involves computing the fit at a target point x0 using only the nearby training observations. Local regression is sometimes referred to as a memory-based procedure because, like nearest-neighbors, we need all the training data each time we wish to compute a prediction. 
    - Try different kernels : uniform, gaussian…
    - Change the bandwidth…
- Splines : cf “BSpline, Bootstrap Var Bias MSE, Smoothing Spline, CV.pdf”
    - Piecewise Polynomial Functions : Splines are piecewise polynomial functions which are constrained to join smoothly at knots.
    - “splines” package 
    - Why natural splines ? To avoid weird behavior at the boundaries of the domain. Then smoothing splines which are just Natural Splines with a curvature penalty in the objective function. 
    - Smoothing splines : With smoothing splines we look for a function f( ⋅ ) that makes RSS = ∑n (y − f(x ))^2 small, but that is also smooth but penalyzing it with the integral of its second derivatives squared.
        - This is operationally done by performing a penalized regression over the natural spline basis, placing knots at all the inputs. This means that, remarkably, the problem defined on an infinite-dimensional function space has a finite-dimensional, unique minimizer which is a natural cubic spline! 
        - This is automatically performed in R through the smooth.spline function. 
        - Generally, one wants to optimize the LOOCV (Leave One Out Cross Validation) or the GCV (Generalized Cross Validation) error

## LAB 12 : generalised additive models (GAM)
- Multivariate linear model :
    - Without interaction :
    - With interaction :
- GAM : cf “GAM, Reduced Model, Linear Model, AVG.pdf”
    - GAM with cubic smoothing spline :
    - GAM with natural spline smoothing :
- Semiparametric model :
- F-Test for comparing two models :
- GAM : “GAM, Reduced Model, Linear Model, AVG.pdf”
    - With interaction : 
    - With thin plate splines : 
- Nonparametric inference for NP regression :

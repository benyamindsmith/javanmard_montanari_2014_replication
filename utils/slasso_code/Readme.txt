This folder contains the source code of the method proposed in " Confidence Intervals and Hypothesis Testing for High-Dimensional Regression", by Adel Javanmard and Andrea Montanari.

*Description:
Script 'lasso_inference'' contains the main function "SSLasso" for computing confidence intervals and p-values.  Script 'lasso_test_synth.R' contains a synthetic example.

*Usage:
SSLasso(X, y, alpha = 0.05, lambda = NULL, mu = NULL, intercept = TRUE,
        resol = 1.3, maxiter = 50, threshold = 1e-2, verbose = TRUE)


* Arguments:
 X     :  design matrix
 y     :  response
 alpha :  significance level
 lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
 mu    :  Infinity-norm constraint on M (if null, searches)
 resol :  step parameter for the function that computes M
 maxiter: iteration parameter for computing M
 threshold : tolerance criterion for computing M
 verbose : verbose?


* Returns:
 noise.sd : Estimate of the noise standard deviation
 norm0    : Estimate of the number of significant coefficients
 coef     : Lasso estimated coefficients
 unb.coef : Unbiased coefficient estimates
 low.lim  : Lower limits of confidence intervals
 up.lim   : upper limit of confidence intervals
 pvals    : p-values for the coefficients	


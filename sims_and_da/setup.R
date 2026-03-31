circulant_symmetric_mat <- function(p) {
  # 1. Initialize a vector of zeros of length p
  v <- rep(0, p)
  
  # 2. Set the diagonal entry (j = k)
  v[1] <- 1
  
  # 3. Set the first set of off-diagonals: k in {j+1, ..., j+5}
  # For the first row (j=1), this is indices 2, 3, 4, 5, 6
  v[2:6] <- 0.1
  
  # 4. Set the wrap-around off-diagonals: k in {j+p-5, ..., j+p-1}
  # For the first row (j=1), this is indices p-4, p-3, p-2, p-1, p
  v[(p-4):p] <- 0.1
  
  # 5. Use the toeplitz function to expand this row into a symmetric matrix
  Sigma <- toeplitz(v)
  
  return(Sigma)
}

compute_metrics <- function(thetahat, ci_lower, ci_upper, pvals, theta0, S, alpha) {
  p <- length(theta0)
  Sc <- setdiff(1:p, S) 
  s0 <- length(S)
  
  # Interval Lengths (Eq 57 & 58)
  lengths <- ci_upper - ci_lower
  avg_len <- mean(lengths)
  avg_len_S <- mean(lengths[S])
  avg_len_Sc <- mean(lengths[Sc])
  
  # Coverage Flags
  covered <- (theta0 >= ci_lower) & (theta0 <= ci_upper)
  
  # Average Coverage (Eq 59, 60, 61)
  dCov <- mean(covered)
  dCov_S <- mean(covered[S])
  dCov_Sc <- mean(covered[Sc])
  
  # Error Rates
  rejected <- (pvals < alpha) 
  TPR <- sum(rejected[S]) / s0              
  FPR <- sum(rejected[Sc]) / (p - s0)       
  
  # Estimation Accuracy (Single numbers for easy averaging)
  if (!is.null(thetahat)) {
    # We store the MEAN bias and MEAN squared error of the vector
    vec_bias <- mean(thetahat - theta0)
    vec_mse  <- mean((thetahat - theta0)^2)
  } else {
    vec_bias <- NA
    vec_mse  <- NA
  }
  
  # Return ONLY single numeric values
  return(c(
    AvgLength = avg_len,
    AvgLength_S = avg_len_S,
    AvgLength_Sc = avg_len_Sc,
    dCov = dCov,
    dCov_S = dCov_S,
    dCov_Sc = dCov_Sc,
    TPR = TPR,
    FPR = FPR,
    MeanBias = vec_bias,
    MeanMSE = vec_mse
  ))
}


# Helper to build CI dataframe
get_ci_data <- function(method_name, indices, bhat, low, up) {
  data.frame(
    Gene = gene_ids[indices], 
    Index = factor(indices, levels = top_genes),
    Estimate = if(is.null(bhat)) NA else bhat[indices],
    Low = low[indices],
    Up = up[indices],
    Method = method_name
  )
}


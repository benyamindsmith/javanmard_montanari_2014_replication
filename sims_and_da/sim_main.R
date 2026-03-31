library(mvnfast)
library(hdi)
library(parallel)
library(doParallel)
library(foreach)

source("lasso_inference.R")
source("setup.R")

# 1. SETUP CLUSTER
cores_to_use <- 10
cl <- makeCluster(cores_to_use, outfile = "") 
registerDoParallel(cl)

# Due to resource constrain, the
# first 6 of the 12 configurations are considered
configurations <- list(
  c(n=1000, p=600, s0=10, b=0.5),
  c(n=1000, p=600, s0=10, b=0.25),
  c(n=1000, p=600, s0=10, b=0.1),
  c(n=1000, p=600, s0=30, b=0.5),
  c(n=1000, p=600, s0=30, b=0.25),
  c(n=1000, p=600, s0=30, b=0.1)
)

config_no <- length(configurations)
res <- vector("list", config_no) 
n_sims <- 20
alpha <- 0.05
# Critical value for 95% CI construction for the Projection methods
z_crit <- qnorm(1 - alpha/2)

for(i in 7:config_no){
  dat <- configurations[[i]] 
  n <- dat["n"]; p <- dat["p"]; s0 <- dat["s0"]; b <- dat["b"]
  
  Sigma <- circulant_symmetric_mat(p)
  S <- sample(1:p, size = s0, replace = FALSE)
  theta0 <- rep(0, p) 
  theta0[S] <- b 
  
  cat(sprintf("\n--- Starting Config %d (n=%d, p=%d) ---\n", i, n, p))
  
  sim_results <- foreach(sim = 1:n_sims, 
                         .packages = c("mvnfast", "hdi"),
                         .export = c("SSLasso", "compute_metrics")) %dopar% {
                           
                           source("lasso_inference.R", local = TRUE)
                           source("setup.R", local = TRUE)
                           set.seed(i * 1000 + sim)
                           
                           X <- rmvn(n, mu=rep(0, p), sigma=Sigma)
                           W <- rnorm(n) 
                           y <- as.numeric(X %*% theta0 + W) 
                           
                           # --- FIT MODELS ---
                           fit_ss <- SSLasso(X, y, alpha=alpha)
                           fit_ms <- multi.split(X, y)
                           fit_rp <- ridge.proj(X, y)
                           fit_lp <- lasso.proj(X, y)
                           
                           # --- MAP OUTPUTS TO COMPUTE_METRICS ---
                           
                           # 1. SS Lasso (Uses unb.coef, low.lim, up.lim, pvals)
                           m_ss <- compute_metrics(fit_ss$unb.coef, fit_ss$low.lim, fit_ss$up.lim, fit_ss$pvals, theta0, S, alpha)
                           
                           # 2. Multi-Split (Uses estimate [if exists], lci, uci, pval)
                           # Note: multi.split often doesn't return a point estimate. We pass NULL if missing.
                           m_ms <- compute_metrics(NULL, fit_ms$lci, fit_ms$uci, fit_ms$pval.corr, theta0, S, alpha)
                           
                           # 3. Ridge Proj (Uses bhat, pval, and SE to build intervals)
                           ci_low_rp <- fit_rp$bhat - z_crit * fit_rp$se
                           ci_up_rp  <- fit_rp$bhat + z_crit * fit_rp$se
                           m_rp <- compute_metrics(fit_rp$bhat, ci_low_rp, ci_up_rp, fit_rp$pval, theta0, S, alpha)
                           
                           # 4. Lasso Proj (Uses bhat, pval, and SE to build intervals)
                           ci_low_lp <- fit_lp$bhat - z_crit * fit_lp$se
                           ci_up_lp  <- fit_lp$bhat + z_crit * fit_lp$se
                           m_lp <- compute_metrics(fit_lp$bhat, ci_low_lp, ci_up_lp, fit_lp$pval, theta0, S, alpha)
                           
                           list(sslasso = m_ss, multi = m_ms, ridge = m_rp, lasso_proj = m_lp)
                         }
  
  # --- AGGREGATE ---
  # These are now all named numeric vectors, so colMeans works perfectly.
  res[[i]] <- list(
    params = list(n = n, p = p, s0 = s0, b = b),
    metrics = list(
      sslasso    = colMeans(do.call(rbind, lapply(sim_results, `[[`, "sslasso")), na.rm=TRUE),
      multi_split= colMeans(do.call(rbind, lapply(sim_results, `[[`, "multi")), na.rm=TRUE),
      ridge_proj = colMeans(do.call(rbind, lapply(sim_results, `[[`, "ridge")), na.rm=TRUE),
      lasso_proj = colMeans(do.call(rbind, lapply(sim_results, `[[`, "lasso_proj")), na.rm=TRUE)
    )
  )
  # If you want to have a checkpoint uncomment the line below.
  # saveRDS(res, "res_backup.rds")
}

stopCluster(cl)
saveRDS(res, "res.rds")
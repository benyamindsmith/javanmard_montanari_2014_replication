# Create Latex Tables for the Report
library(xtable)
res <- readRDS("./res.rds")

# Construct appropriate Data Frames
# Table 1 Dataframe
tab1 <- do.call(rbind, lapply(res, function(x) {
  # Extract params (handling them as vectors)
  p_str <- paste0("(", x$params[["n"]], ",", x$params[["p"]], ",", 
                  x$params[["s0"]], ",", x$params[["b"]], ")")
  
  # Extract metrics using [ ] since sslasso is a named vector
  res <- x$metrics$sslasso
  
  data.frame(
    configuration = p_str,
    ell     = res["AvgLength"],
    ell_S   = res["AvgLength_S"],
    ell_Sc  = res["AvgLength_Sc"],
    dCov    = res["dCov"],
    dCov_S  = res["dCov_S"],
    dCov_Sc = res["dCov_Sc"],
    row.names = NULL # Prevents messy row names in the final table
  )
}))
# Table 2 Dataframe
tab2 <- do.call(rbind, lapply(res, function(x) {
  # Extract params (handling them as vectors)
  p_str <- paste0("(", x$params[["n"]], ",", x$params[["p"]], ",", 
                  x$params[["s0"]], ",", x$params[["b"]], ")")
  
  # Extract metrics using [ ] since sslasso is a named vector
  sslasso <- x$metrics$sslasso
  multi_split <- x$metrics$multi_split
  ridge_proj <- x$metrics$ridge_proj
  lasso_proj <- x$metrics$lasso_proj
  data.frame(
    configuration = p_str,
    sslasso_fp = sslasso["FPR"],
    sslasso_tp = sslasso["TPR"],
    multi_split_fp= multi_split["FPR"],
    multi_split_tp= multi_split["TPR"],
    ridge_proj_fp= ridge_proj["FPR"],
    ridge_proj_tp= ridge_proj["TPR"],
    lasso_proj_fp= lasso_proj["FPR"],
    lasso_proj_tp= lasso_proj["TPR"])
}))



# Latex Output- Table 1

colnames(tab1) <- c("\\diagbox{Configuration}{Measure}", 
                    "$\\ell$", "$\\ell_S$", "$\\ell_{S^c}$", 
                    "$\\widehat{Cov}$", "$\\widehat{Cov}_S$", "$\\widehat{Cov}_{S^c}$")

xtab1 <- xtable(tab1, align = "ll|c|c|c|c|c|c|", 
                caption = "Simulation results for synthetic data. Results correspond to 95\\% confidence intervals.",
                digits =4)

print(xtab1, 
      include.rownames = FALSE, 
      sanitize.text.function = function(x) x,
      hline.after = c(-1, 0, nrow(tab1)),
      type = "latex")

colnames(tab2) <- c("Configuration", "FP", "TP", "FP", "TP", "FP", "TP", "FP", "TP")


# Latex Output- Table 2
header_list <- list()
header_list$pos <- list(-1)
header_list$command <- c(
  paste0(
    "\\hline & \\multicolumn{2}{c|}{", " \\makecell{Javanmard \\\\ \\& Montanari}} ",
    "& \\multicolumn{2}{c|}{", " \\makecell{Multisample- \\\\ splitting}} ",
    "& \\multicolumn{2}{c|}{", " \\makecell{Ridge projection  \\\\ estimator}} ",
    "& \\multicolumn{2}{c|}{", " \\makecell{LASSO projection \\\\  estimator}} \\\\ \\cline{2-9} "
  )
)
xtab2 <- xtable(tab2, align = "ll|cc|cc|cc|cc|",
                caption = "Simulation results for FP and TP rates at significance level $\\alpha = 0.05$.",
                digits = 4)

print(xtab2, 
      include.rownames = FALSE,
      add.to.row = header_list,
      hline.after = c(0, nrow(tab2)),
      sanitize.text.function = function(x) x,
      type = "latex")

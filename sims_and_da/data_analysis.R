library(mvnfast)
library(hdi)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 
source("lasso_inference.R")
source("setup.R")

data("riboflavin")
alpha <- 0.05
z_crit <- qnorm(1 - alpha/2)

X<- riboflavin$x
y<- riboflavin$y
fit.multi <- multi.split(X,y)
fit.ridge.proj <- ridge.proj(X, y)
fit.sslasso <- SSLasso(X,y)
fit.lasso.proj <- lasso.proj(X,y)

# Save Model Objects
# data_analyis_mods <- list(
#   fit.multi,
#   fit.ridge.proj,
#   fit.sslasso,
#   fit.lasso.proj
# )
# 
# saveRDS(data_analyis_mods, "data_analysis_mods.rds")

# Setup IDs and Thresholds
gene_ids <- colnames(X)
p <- ncol(X)
alpha <- 0.05
z_crit <- qnorm(1 - alpha/2)
bonferroni_thresh <- -log10(alpha / p)


# 1. MANHATTAN PLOT DATA
pvals_df <- data.frame(
  Gene_Index = 1:p,
  Gene_ID = gene_ids,
  `Multi-Split` = fit.multi$pval.corr,
  `Ridge-Proj`  = fit.ridge.proj$pval,
  `SSLasso`     = fit.sslasso$pvals,
  `Lasso-Proj`  = fit.lasso.proj$pval,
  check.names = FALSE
)

pvals_long <- pvals_df %>%
  pivot_longer(cols = -c(Gene_Index, Gene_ID), names_to = "Method", values_to = "P_Value") %>%
  mutate(LogP = -log10(P_Value)) 

# 2. FOREST PLOT DATA
# Identify top 10 genes by significance across all methods
top_genes <- pvals_df %>%
  rowwise() %>%
  mutate(Min_P = min(c_across(c(`Multi-Split`, `Ridge-Proj`, `SSLasso`, `Lasso-Proj`)), na.rm = TRUE)) %>%
  ungroup() %>% 
  arrange(Min_P) %>% 
  slice_head(n = 10) %>% 
  pull(Gene_Index)

# Calculate CIs for Projection methods
ci_low_rp <- fit.ridge.proj$bhat - z_crit * fit.ridge.proj$se
ci_up_rp  <- fit.ridge.proj$bhat + z_crit * fit.ridge.proj$se
ci_low_lp <- fit.lasso.proj$bhat - z_crit * fit.lasso.proj$se
ci_up_lp  <- fit.lasso.proj$bhat + z_crit * fit.lasso.proj$se

ci_df <- bind_rows(
  get_ci_data("SSLasso", top_genes, fit.sslasso$unb.coef, fit.sslasso$low.lim, fit.sslasso$up.lim),
  get_ci_data("Multi-Split", top_genes, NULL, fit.multi$lci, fit.multi$uci), 
  get_ci_data("Ridge-Proj", top_genes, fit.ridge.proj$bhat, ci_low_rp, ci_up_rp),
  get_ci_data("Lasso-Proj", top_genes, fit.lasso.proj$bhat, ci_low_lp, ci_up_lp)
)

# 3. APPLY ACADEMIC LABELS
method_levels <- c("Multi-Split", "Ridge-Proj", "Lasso-Proj", "SSLasso")
method_labels <- c(
  "Multisample-splitting", 
  "Ridge projection estimator", 
  "Lasso projection estimator", 
  "Javanmard & Montanari (2014)"
)

pvals_long$Method <- factor(pvals_long$Method, levels = method_levels, labels = method_labels)
ci_df$Method <- factor(ci_df$Method, levels = method_levels, labels = method_labels)
ci_df$Gene <- factor(ci_df$Gene, levels = gene_ids[top_genes])

# 4. PLOTTING
p_manhattan <- ggplot(pvals_long, aes(x = Gene_Index, y = LogP, color = Method)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = bonferroni_thresh, linetype = "dashed", color = "red") +
  guides(color = "none") + 
  facet_wrap(~ Method, ncol = 1) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "gray95")) +
  labs(title = "Global P-Value Comparison", x = "Gene Index", y = expression(-log[10](p-value)))

p_forest <- ggplot(ci_df, aes(y = Method, color = Method)) +
  geom_errorbar(aes(xmin = Low, xmax = Up), width = 0.3, linewidth = 0.8) +
  geom_point(aes(x = Estimate), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ Gene, ncol = 2, scales = "free_x") + 
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 11) + 
  labs(title = "Model Comparison per Gene", subtitle = "Top 10 discoveries", x = "Coefficient Estimate", y = NULL) +
  theme(
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold", margin = margin(t = 2, b = 2)),
    panel.spacing = unit(0.4, "lines"),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )

# Final Assembly
final_plot <- (p_manhattan | p_forest) + 
  plot_layout(widths = c(1, 1.2), guides = "collect") & 
  theme(
    legend.position = "bottom", 
    legend.box.margin = margin(t = 10),
    legend.text = element_text(size = 9)
  )

print(final_plot)
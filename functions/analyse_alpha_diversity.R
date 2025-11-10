analyse_alpha_diversity <- function(data, response_var, pairwise = FALSE, file_prefix = "alpha") {
  
  mycols <- RColorBrewer::brewer.pal(n = 6, name = "Set2")
  
  # 1. Plot boxplot
  p <- ggplot(data, aes_string(x = "Tissue_health", y = response_var)) +
    geom_boxplot(aes(fill = Site)) +
    theme_minimal() +
    geom_jitter(width = 0.05, alpha = 0.6) +
    scale_fill_manual(values = mycols) +
    scale_color_manual(values = mycols) +
    theme(
      axis.text = element_text(size = 10, colour = "black"),
      axis.line = element_line(color = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  print(p)
  
  # 2. Fit linear mixed model
  formula <- as.formula(paste0(response_var, " ~ Tissue_health + Site + ReadDepth_z + (1 | Site_Tree)"))
  model <- lmer(formula, data = data)
  print(summary(model))
  
  # 3. Type II ANOVA
  anova_res <- Anova(model, type = "II")
  print(anova_res)
  
  # 4. Model diagnostics
  check_model(model)
  sim_res <- simulateResiduals(fittedModel = model)
  plot(sim_res)
  
  # 5. Marginal and conditional R²
  r2_res <- r2(model)
  print(r2_res)
  
  # 6. Pairwise comparisons (optional)
  pairwise_res <- NULL
  if (pairwise) {
    pairwise_res <- emmeans(model, pairwise ~ Tissue_health, adjust = "tukey")
    print(pairwise_res)
  }
  
  # Return results as a list
  return(list(
    plot = p,
    model = model,
    anova = anova_res,
    r2 = r2_res,
    pairwise = pairwise_res
  ))
}

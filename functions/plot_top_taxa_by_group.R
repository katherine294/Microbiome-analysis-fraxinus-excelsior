plot_top_taxa_by_group <- function(ps_obj, 
                                   tax_level = "Family", 
                                   group = "Tissue_health", 
                                   top_n = 20,
                                   mycols = c("brown3", "steelblue", "grey50"),
                                   out_dir = "/rds/homes/k/kgh742/psf_wgs_project/02.MicrobiomeAnalysis/ITS_Figures/composition_barplots/test_scripts/") {
  
  # Ensure output directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message(paste0("Processing taxonomic level: ", tax_level))
  
  # Clean taxonomy
  ps_obj <- ps_obj %>%
    tax_fix(min_length = 5,
            unknowns = c("Unknown", "NA", "Incertae Sedis", "Unclassified"),
            sep = " ",
            anon_unique = TRUE,
            suffix_rank = "classified")
  
  # Agglomerate to taxonomic level
  ps_obj <- tax_glom(ps_obj, taxrank = tax_level, NArm = FALSE)
  
  # Transform to compositional data
  ps_obj <- microbiome::transform(ps_obj, "compositional")
  
  # Get group abundances
  grp_abund <- get_group_abundances(ps_obj, level = tax_level, group = group)
  
  # Clean OTU/taxon names
  grp_abund$OTUID <- gsub("(k__|p__|c__|o__|f__|g__|s__)", "", grp_abund$OTUID)
  grp_abund$OTUID <- ifelse(grp_abund$OTUID == "" | grp_abund$OTUID == " ", "Unclassified", grp_abund$OTUID)
  
  # Identify top taxa across all groups
  top_taxa <- grp_abund %>%
    group_by(OTUID) %>%
    summarise(overall_mean = mean(mean_abundance)) %>%
    arrange(desc(overall_mean)) %>%
    slice_head(n = top_n)
  
  grp_abund_top <- grp_abund %>%
    filter(OTUID %in% top_taxa$OTUID)
  
  # --- Plot ---
  p <- grp_abund_top %>%
    ggplot(aes(x = reorder(OTUID, mean_abundance),
               y = mean_abundance,
               fill = .data[[group]])) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(group, values = mycols) +
    coord_flip() +
    labs(
      x = tax_level,
      y = "Mean Relative Abundance",
      title = paste("Top", top_n, tax_level, "by", group)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "top"
    )
  
  # --- File names ---
  plot_file <- file.path(out_dir, paste0("barplot_test_", tolower(tax_level), "_top", top_n, ".pdf"))
  csv_file  <- file.path(out_dir, paste0(tolower(tax_level), "_taxa_summary_top", top_n, ".csv"))
  
  # --- Save outputs ---
  ggsave(plot_file, plot = p, width = 7, height = 5, device = "pdf")
  message(paste("Saved PDF:", plot_file))
  
  tx.sum <- taxa_summary(ps_obj, tax_level)
  write_csv(tx.sum, csv_file)
  message(paste("Saved CSV:", csv_file))
  
  return(p)
}

plot_taxonomic_barplot <- function(ps_obj, taxrank, n_taxa = 20, file_suffix = NULL, 
                                   width = 10, height = 8, unknowns = c("Incertae Sedis", "NA", "Unknown family"),
                                   facet_var = "Tissue_health", sample_order = ordered_ids) {
  
  # Clean taxonomy labels
  ps_obj <- tax_fix(ps_obj,
                    min_length = 5,
                    unknowns = unknowns,
                    sep = " ",
                    anon_unique = TRUE,
                    suffix_rank = "classified")
  
  # Agglomerate to desired taxonomic level if not Phylum
  if (taxrank != "Phylum") {
    ps_obj <- tax_glom(ps_obj, taxrank = taxrank, NArm = FALSE)
  }
  
  # Transform to compositional abundance
  ps_obj <- microbiome::transform(ps_obj, "compositional")
  
  # Clean taxonomy prefixes and unclassifieds
  tax_tab <- as.data.frame(tax_table(ps_obj))
  tax_tab[] <- lapply(tax_tab, function(x) gsub("^(k__|p__|c__|o__|f__|g__|s__)", "", x))
  tax_tab[] <- lapply(tax_tab, function(x) ifelse(x == "" | is.na(x), "Unclassified", x))
  tax_table(ps_obj) <- as.matrix(tax_tab)
  
  # Create file name suffix if none provided
  if (is.null(file_suffix)) file_suffix <- tolower(taxrank)
  
  # Create output paths
  pdf_path <- sprintf("/rds/homes/k/kgh742/psf_wgs_project/02.MicrobiomeAnalysis/ITS_Figures/composition_barplots/test_scripts/barplot_test_%s.pdf", file_suffix)
  csv_path <- sprintf("/rds/homes/k/kgh742/psf_wgs_project/02.MicrobiomeAnalysis/ITS_Figures/composition_barplots/test_scripts/%s_taxa_summary_test.csv", file_suffix)
  
  # Open PDF
  pdf(pdf_path, width = width, height = height)
  
  # Generate plot
  p <- ps_obj %>%
    comp_barplot(
      tax_level = taxrank,
      n_taxa = n_taxa,
      other_name = "Other",
      order_with_all_taxa = TRUE,
      merge_other = TRUE,
      bar_width = 0.7,
      bar_outline_colour = "grey5",
      facet_by = facet_var,
      sample_order = sample_order
    ) +
    coord_flip() +
    guides(fill = guide_legend(ncol = 1)) +
    theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white", colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
    )
  
  print(p)
  dev.off()
  
  # Save summary CSV
  tx_sum <- taxa_summary(ps_obj, taxrank)
  write_csv(tx_sum, csv_path)
  
  message(sprintf("Saved %s-level plot and data to:\n%s\n%s", taxrank, pdf_path, csv_path))
}

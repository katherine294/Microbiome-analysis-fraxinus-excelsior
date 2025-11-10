# Function to convert ANOVA or emmeans contrasts to dataframe
results_to_df <- function(res_list, index_name, emmeans = FALSE) {
  if (emmeans) {
    if (is.null(res_list$pairwise)) return(NULL)  # Skip if no pairwise results
    df <- as.data.frame(res_list$pairwise$contrasts) %>%
      mutate(Index = index_name)
  } else {
    df <- as.data.frame(res_list$anova) %>%
      mutate(Index = index_name)
  }
  df$Effect <- rownames(df)
  return(df)
}

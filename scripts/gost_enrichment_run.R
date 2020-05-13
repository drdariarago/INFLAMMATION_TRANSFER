library(tidyverse)
library(magrittr)
library(gprofiler2)

# Import list of ranked gene names
ranked_gene_list <- 
  read_rds(path = snakemake@input[["ranked_genes"]])

# Run GO enrichment
gost_result_list <-
  map(
    .x = ranked_gene_list,
    .f = ~ gprofiler2::gost(
      query = .x,
      ordered_query = TRUE, 
      multi_query = FALSE,
      organism = "mmusculus",
      sources = c("GO"), 
      domain_scope = "custom", 
      custom_bg = .x,
      correction_method = "gSCS",
      user_threshold = snakemake@params[["max_fdr"]],
      significant = FALSE,
      measure_underrepresentation = FALSE
    )
  )

write_rds(x = gost_result_list, path = snakemake@output[["raw_results"]])
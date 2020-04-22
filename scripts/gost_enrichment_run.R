library(tidyverse)
library(magrittr)
library(gprofiler2)

# Import background set of genes expressed in each tissue
tissue_bg <-
  read_rds(path = snakemake@input[["expressed_genes"]])

# Import list of ranked gene names
ranked_gene_list <- 
  read_rds(path = snakemake@input[["ranked_genes"]])

# Run GO enrichment
gost_result_list <-
  map2(
    .x = ranked_gene_list,
    .y = tissue_bg,
    .f = ~ gprofiler2::gost(
      query = .x,
      ordered_query = TRUE, 
      multi_query = TRUE,
      organism = "mmusculus",
      sources = c("GO"), 
      domain_scope = "custom", 
      custom_bg = .y,
      correction_method = "gSCS",
      user_threshold = snakemake@params[["max_fdr"]],
      significant = TRUE,
      measure_underrepresentation = FALSE
    )
  )

write_rds(x = gost_result_list, path = snakemake@output[["raw_results"]])
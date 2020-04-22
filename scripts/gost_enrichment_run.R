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
      user_threshold = snakemake@params[["min_fdr"]],
      significant = TRUE,
      measure_underrepresentation = FALSE
    )
  )

write_rds(x = gost_result_list, path = snakemake@output[["raw_results"]])

## Convert to summary result table

gost_results_table <- 
  gost_result_list %>% 
  map_df(.f = ~ .x$result, .id = "contrast") %>% 
  dplyr::select(contrast, term_id, term_name, significant, source, term_size, query_sizes, intersection_sizes) %>% 
  mutate_if(
    .predicate = is.list, 
    .funs = ~ map_chr(.x = ., .f = paste0, collapse = ":" )
  ) %>% 
  separate(
    col = significant,
    sep = ":",
    convert = TRUE,
    into = paste0("signif_", c(2,5,12,24), "_hrs")
  ) %>%  separate(
    col = query_sizes,
    sep = ":",
    convert = TRUE,
    into = paste0("query_size_", c(2,5,12,24), "_hrs")
  ) %>%  separate(
    col = intersection_sizes,
    sep = ":",
    convert = TRUE,
    into = paste0("intersection_size_", c(2,5,12,24), "_hrs")
  )

write_csv(x = gost_results_table, path = snakemake@output[["results_table"]])

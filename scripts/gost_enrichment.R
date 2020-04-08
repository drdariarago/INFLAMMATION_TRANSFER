library(tidyverse)
library(magrittr)
library(gprofiler2)

# Import fold change values
gene_list <-
  snakemake@input[[1]] %>% 
  readRDS() %>% 
  map(.x = ., .f = ~ filter(.x, treatment == "exposureLPS")) %>%
  map(.x = ., .f = ~ dplyr::select(.x, -treatment))  %>% 
  map(.x = ., .f = ~ mutate(.data = .x, gene_id = gsub('\\..*', '', gene_id))) %>% 
  map(.x = ., .f = ~ filter(.data = .x, AveExpr > snakemake@params[["expression_threshold"]]) )

# Create list of genes that are expressed in each tissue
tissue_bg <-
  map(
    .x = gene_list,
    .f = 
      ~ .x %>% 
      pull(gene_id) %>% 
      unique
  )

# Create list of queries sorted by LogFC for each timepoint
ranked_gene_list <-
  gene_list %>% 
  map(.f = ~ dplyr::arrange(.x, desc(LogFC))) %>% 
  map(.f = ~ 
        mutate(.data = .x,
               timepoint = fct_relevel(
                 timepoint,
                 c("timepoint2", "timepoint5", "timepoint12", "timepoint24") 
               )
        )
  ) %>% 
  map(.f = ~ split(x = .x, f = .x$timepoint)) %>% 
  modify_depth(.depth = 2, .f = ~ pull(.data = .x, gene_id))

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
      sources = c("GO", "TF"), 
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
  go_results %>% 
  dplyr::select(contrast, term_id, significant, source, term_name, query_sizes, intersection_sizes) %>% 
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
    into = paste0("total_annotated_", c(2,5,12,24), "_hrs")
  ) %>%  separate(
    col = intersection_sizes,
    sep = ":",
    convert = TRUE,
    into = paste0("significant_annotated_", c(2,5,12,24), "_hrs")
  )

write_csv(x = gost_results_table, path = snakemake@output[["results_table"]])
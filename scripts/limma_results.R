## Limma results

library(tidyverse)
library(stageR)
library(limma)
library(magrittr)

limma_results <-
  snakemake@input[["results"]] %>% 
  read_rds()

LOG_THRESHOLD = snakemake@params[['fold_change_threshold']] %>% log2
ALPHA = snakemake@params[['alpha']]
COEFS = 5:8

# Summarize results and filter by log fold change

limma_results_summary <-
  limma_results %>% 
  map(
    ~ topTable(
      .x, 
      sort.by = 'none', 
      number = Inf, 
      adjust.method = 'none',
      lfc = LOG_THRESHOLD, 
      coef = COEFS
    )
  )

# Extract global p values for initial screening of genes with at least 1 significant contrast

limma_global_pvals <-
  limma_results_summary %>%
  map(.f = ~ setNames(.x$P.Value, nm = .x$gene_id))

# Subset only gene results that pass initial filtering criteria
limma_filtered_results <-
  limma_results %>% 
  map2(
    .x = .,
    .y = limma_global_pvals,
    ~ .x$p.value %>% 
      magrittr::extract(, COEFS) %>% 
      magrittr::extract(row.names(.) %in% names(.y),)
  )

# Create stager object 
stager_input <-
  map2(
    .x = limma_global_pvals,
    .y = limma_filtered_results, 
    .f = ~ stageR(
      pScreen = .x,
      pConfirmation = .y,
      pScreenAdjusted = FALSE
    )
  )

stager_result <- 
  stager_input %>% 
  map(
    stageWiseAdjustment,
    method = 'holm',
    alpha = ALPHA
  )

stager_binary_result <- 
  stager_result %>% 
  map(
    getResults
  ) %>% 
  map(
    as.tibble,
    rownames = "gene_id"
  ) %>% 
  map(
    filter,
    padjScreen == 1
  ) %>% 
  map(
    dplyr::select, 
    -starts_with("padj")
  ) %>% 
  map(
    ~ set_colnames(
      x = .x,
      value = names(.x) %>% str_extract(pattern = "(^timepoint[0-9]{1,2}|gene_id)")
    )
  )

# Convert ENSEMBL IDs to gene symbols
gene_annotation <- read_rds(snakemake@input[["gene_data"]])

annotated_results <-
  stager_binary_result %>% 
  map(
    mutate,
    gene_id = str_extract(string = gene_id, pattern = "^[[:alnum:]]*")
  ) %>% 
  map(
    .f = ~ right_join(x = gene_annotation, y = .x, by = c("ensembl_gene_id" = "gene_id"))
  ) %>% 
  map(
    .f = ~ cbind(.x[,-3], .x[,3])
  )

write_rds(x = annotated_results, path = snakemake@output[["r_summaries"]])

map2(
  .x = annotated_results, 
  .y = snakemake@output[["csv_summaries"]],
  .f = ~ write_csv(
    x = .x,
    path = .y
  )
)

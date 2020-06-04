## Format fold change tables to generate global report

library(tidyverse)
library(magrittr)

# Load logFC tibbles for all tissues and subset only to responses
fold_change_list <-
  snakemake@input[["fold_change_summaries"]] %>% 
  setNames(
    object = ., 
    nm = 
      str_extract(string = ., pattern = "limma_[a-z,_]*") %>% 
      str_remove(string = ., pattern = "limma_")
  ) %>% 
  map(.f = read_rds) %>% 
  map(
    .f = ~ filter(.x, exposure == "response")
  )

# Reshape all non-placenta tissues to gene by time responses
response_matrix_list <-
  fold_change_list %>% 
  map(
    .f = ~ select(.x, ensembl_gene_id, timepoint, logFC)
  ) %>% 
  magrittr::extract(-3) %>% 
  map(
    .f = ~ pivot_wider(data = .x, names_from = timepoint, values_from = logFC)
  ) %>% 
  map(
    .f = ~ column_to_rownames(.x, "ensembl_gene_id")
  ) %>% 
  map(
    .f = as.matrix
  )

# Extract shared and differences in placental responses
shared_response_placenta <-
  fold_change_list[["placentas"]] %>% 
  filter(maternal == "shared") %>% 
  select(ensembl_gene_id, timepoint, logFC) %>% 
  pivot_wider(names_from = timepoint, values_from = logFC) %>% 
  column_to_rownames("ensembl_gene_id") %>% 
  as.matrix()

differences_response_placenta <-
  fold_change_list[["placentas"]] %>% 
  filter(maternal == "difference") %>% 
  select(ensembl_gene_id, timepoint, logFC) %>% 
  pivot_wider(names_from = timepoint, values_from = logFC) %>% 
  column_to_rownames("ensembl_gene_id") %>% 
  as.matrix()

# Add net maternal and net fetal response to the list of matrices
response_matrix_list[["maternal_placenta"]] <-
  shared_response_placenta + differences_response_placenta

response_matrix_list[["fetal_placenta"]] <-
  shared_response_placenta - differences_response_placenta

# Save as rds
write_rds(response_matrix_list, snakemake@output[[1]])

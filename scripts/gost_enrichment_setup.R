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

write_rds(x = tissue_bg, path = snakemake@output[["expressed_genes"]])

# Create list of queries sorted by LogFC for each timepoint
ranked_gene_list_descending <-
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

write_rds(x = ranked_gene_list_descending, path = snakemake@output[["genes_descending_lfc"]])

ranked_gene_list_ascending <-
  gene_list %>% 
  map(.f = ~ dplyr::arrange(.x, LogFC)) %>% 
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

write_rds(x = ranked_gene_list_ascending, path = snakemake@output[["genes_ascending_lfc"]])
## Limma results

library(tidyverse)
library(stageR)
library(limma)

limma_results <-
  here::here("results/limma/fitted_models.Rdata") %>% 
  read_rds()


## Filter out genes with less than 0.5 change across all contrasts
## Pick p-values

limma_global_pvals <-
  limma_results %>% 
  map(
    ~ topTable(.x, sort.by = 'none', number = Inf) %$%
      P.Value %>% 
      set_names(x = ., nm = row.names(.x$p.value))
  )

stager_input <-
  map2(
    .x = limma_global_pvals,
    .y = limma_results, 
    .f = ~ stageR(
      pScreen = .x,
      pConfirmation = .y$p.value[,5:8],
      pScreenAdjusted = FALSE
    )
  )

stager_result <- 
  stager_input %>% 
  map(
    stageWiseAdjustment,
    method = 'none',
    alpha = 0.05
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
    ~ set_names(
      x = .x,
      nm = names(.x) %>% str_extract(pattern = "(^timepoint[0-9]{1,2}|gene_id)")
    )
  )

# Convert ENSEMBL IDs to gene symbols

library(biomaRt)

gene_annotation <-
  getBM(
    mart = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = 'useast'),
    filters = "ensembl_gene_id", 
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values = map(
      .x = stager_binary_result, 
      .f = ~ .x$gene_id %>%
        str_extract(string = ., pattern = "^[[:alnum:]]*")
      ) %>% unlist %>% unique 
  )

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

write_rds(x = annotated_results, path = here::here("results/limma_results/annotated_results.Rdata"))

map2(
  .x = annotated_results, 
  .y = names(annotated_results),
  .f = ~ write_csv(
    x = .x,
    path = paste0("results/limma_results/significant_contrasts_", .y, ".csv") 
  )
)

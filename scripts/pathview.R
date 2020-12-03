### Visualize pathways with pathview

library(pathview)
library(tidyverse)
library(magrittr)

gene_data <-
  here::here("results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>% 
  read_csv() %>% 
  mutate(
    timepoint = factor(x = timepoint, levels = paste0("timepoint", c(2,5,12,24)))
  )

phospho_data <-
  here::here("results/phospho_import/fold_changes.csv") %>% 
  read_csv() %>% 
  mutate(
    timepoint = factor(x = timepoint, levels = paste0("timepoint", c(2,5,12,24)))
  )

pathways <- "mmu04010"
  # c("mmu04064", "mmu04657", "mmu04062", "mmu04514",  "mmu04668", "mmu03320", 
    # "mmu04020", "mmu04926", "mmu04920", "mmu04650", "mmu04620", "mmu04668", 
    # "mmu04330", "mmu04010")

tissues <- 
  c("fetal_liver", "placentas")

# Convert to gene by timepoint matrix for multi-state plotting

gene_sets <-
  map(
    .x = tissues, 
    .f = ~ 
      filter(
        .data = gene_data, 
        tissue == .x,
        exposure == "response",
        # q_value < 0.05
      )
  ) %>%
  setNames(
    object = .,
    nm = tissues
  ) %>% 
  modify(
    .x = .,
    .f = ~ pivot_wider(
      data = .x,
      id_cols = ensembl_gene_id, names_from = timepoint, values_from = logFC
    ) 
  ) %>% 
  modify(
    .x = .,
    .f = ~ .x %>% column_to_rownames("ensembl_gene_id") %>% as.matrix
  )

## Plot RNA data

imap(
  .x = gene_sets,
  .f = ~ pathview(
    gene.data = .x,
    out.suffix = .y,
    pathway.id =  pathways,
    species = 'mmu',
    gene.idtype = 'ENSEMBL',
    kegg.native = TRUE,
    multi.state = TRUE,
    limit = list(
      gene = c(-1.5,1.5),
      cpd = 1
    )
  )
)

# convert phospho data to protein to treatment matrix
phospho_sets <-
  tissues %>% 
  purrr::set_names(.) %>% 
  map(
    .f = ~ filter(.data = phospho_data, tissue == .x)
  ) %>% 
  map(
    .f = ~ pivot_wider(
      data = .x,
      id_cols = ensembl_gene_id, 
      names_from = timepoint, 
      values_from = logFC,
      values_fn = function(x){max(abs(x)) * sign(x[which.max(abs(x))]) }
    )
  ) %>% 
  modify(
    .f = ~ filter( .x, is.na(ensembl_gene_id) == FALSE )
  ) %>% 
  modify(
    .f = ~ column_to_rownames( .x, "ensembl_gene_id" ) %>% as.matrix
  )
  


imap(
  .x = phospho_sets,
  .f = ~ pathview(
    gene.data = .x,
    out.suffix = .y,
    pathway.id =  pathways,
    species = 'mmu',
    gene.idtype = 'ENSEMBL',
    kegg.native = TRUE,
    multi.state = TRUE,
    limit = list(
      gene = c(-1.5,1.5),
      cpd = 1
    )
  )
)

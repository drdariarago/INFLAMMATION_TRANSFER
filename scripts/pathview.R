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

pathways <-
  c("mmu04064", "mmu04657")
# "mmu04062", "mmu04514",  "mmu04668", "mmu03320", "mmu04020", "mmu04926", "mmu04920",

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
        q_value < 0.05,
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

imap(
  .x = gene_sets,
  .f = ~ pathview(
    gene.data = .x,
    out.suffix = .y,
    pathway.id =  pathways,
    species = 'mmu',
    gene.idtype = 'ENSEMBL',
    kegg.native = TRUE,
    multi.state = TRUE
  )
)

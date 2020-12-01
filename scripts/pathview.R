### Visualize pathways with pathview

library(pathview)
library(tidyverse)
library(magrittr)

gene_data <-
  here::here("results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>% 
  read_csv()

pathways <-
  c("mmu04064", "mmu04657")
# "mmu04062", "mmu04514",  "mmu04668", "mmu03320", "mmu04020", "mmu04926", "mmu04920",

contrasts <-
  expand_grid(
    tissue = c("fetal_liver", "placentas"),
    timepoint = paste0("timepoint", c(2,5,12,24))
  )

tissues <-
  contrasts %>% 
  pull(tissue)

timepoints <-
  contrasts %>% 
  pull(timepoint)

gene_sets <-
  map2(
    .x = tissues, 
    .y = timepoints, 
    .f = ~ 
      filter(
        .data = gene_data, 
        tissue == .x,
        timepoint == .y,
        exposure == "response",
        q_value < 0.05,
      )
  ) %>%
  setNames(
    object = .,
    nm = paste(tissues, timepoints, sep = "_")
  ) %>% 
  modify(
    .x = .,
    .f = ~ .x %$% 
      setNames(
        object = logFC,
        nm = ensembl_gene_id
      )
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
  )
)

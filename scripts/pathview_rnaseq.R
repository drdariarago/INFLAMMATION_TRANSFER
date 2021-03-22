### Visualize pathways with pathview for RNAseq data

library(pathview)
library(tidyverse)
library(magrittr)

gene_data <-
  snakemake@input[['rna']] %>% 
  read_csv() %>% 
  mutate(
    timepoint = factor(x = timepoint, levels = paste0("timepoint", c(2,5,12,24)))
  )

pathways <- snakemake@params[['pathways']]
tissues <- snakemake@params[['tissues']]

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

# Move RNA results to appropriate directory and clean up leftover files
dir.create( path = here::here(snakemake@output[['rna']]))

list.files(path = here::here(""), pattern = "mmu.*multi.png") %>% 
  file.rename(
    from = .,
    to = here::here( snakemake@output[['rna']], .)
  )

list.files(path = here::here(""), pattern = "mmu.*") %>% 
  file.remove(.)

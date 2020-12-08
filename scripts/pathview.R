### Visualize pathways with pathview

library(pathview)
library(tidyverse)
library(magrittr)

gene_data <-
  snakemake@input[['rna']] %>% 
  read_csv() %>% 
  mutate(
    timepoint = factor(x = timepoint, levels = paste0("timepoint", c(2,5,12,24)))
  )

phospho_data <-
  snakemake@input[['phospho']] %>% 
  read_csv() %>% 
  mutate(
    timepoint = factor(x = timepoint, levels = paste0("timepoint", c(2,5,12,24)))
  )

pathways <- snakemake@params[['pathways']]
tissues <- snakemake@params[['tissues']]

snakemake@output[['phospho']] %>% 
  dirname(.) %>% 
  here::here() %>% 
  dir.create( path = .)

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

# Move phospho results to appropriate directory and clean up leftover files

dir.create( path = here::here(snakemake@output[['phospho']]))

list.files(path = here::here(""), pattern = "mmu.*multi.png") %>% 
  file.rename(
    from = .,
    to = here::here( snakemake@output[['phospho']], .)
  )

list.files(path = here::here(""), pattern = "mmu.*") %>% 
  file.remove(.)

# Merge phospho and gene data 
# Current approach: select the strongest if two genes have both data

merged_data <-
  bind_rows(
    gene_data %>% filter(exposure == "response") %>% select(tissue, timepoint, ensembl_gene_id, logFC), 
    phospho_data %>% select(tissue, timepoint, ensembl_gene_id, logFC)
  ) %>% 
    group_by(tissue, timepoint, ensembl_gene_id) %>% 
    summarise(
      logFC = max(abs(logFC)) * sign(logFC[which.max(abs(logFC))])
    )

merged_sets <-
  tissues %>% 
  purrr::set_names(.) %>% 
  map(
    .f = ~ filter(.data = merged_data, tissue == .x)
  ) %>% 
  map(
    .f = ~ pivot_wider(
      data = .x,
      id_cols = ensembl_gene_id, 
      names_from = timepoint, 
      values_from = logFC
    )
  ) %>% 
  modify(
    .f = ~ filter( .x, is.na(ensembl_gene_id) == FALSE )
  ) %>% 
  modify(
    .f = ~ column_to_rownames( .x, "ensembl_gene_id" ) %>% as.matrix
  )

imap(
  .x = merged_sets,
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

# Move merged results to appropriate directory and clean up leftover files
dir.create( path = here::here(snakemake@output[['merged']]))

list.files(path = here::here(""), pattern = "mmu.*multi.png") %>% 
  file.rename(
    from = .,
    to = here::here( snakemake@output[['merged']], .)
  )

list.files(path = here::here(""), pattern = "mmu.*") %>% 
  file.remove(.)

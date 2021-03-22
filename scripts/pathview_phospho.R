### Visualize pathways with pathview

library(pathview)
library(tidyverse)
library(magrittr)

pathways <- snakemake@params[['pathways']]
tissues <- snakemake@params[['tissues']]

##### Plot gene-wise summaries ####

phospho_genes <-
  snakemake@input[['phospho_genes']] %>%
  read_rds() 


# convert phospho data to protein to treatment matrix
phospho_gene_sets <-
  tissues %>% 
  purrr::set_names(.) %>% 
  map(
    .f = ~ filter(.data = phospho_genes, tissue == .x)
  ) %>% 
  map(
    .f = ~ pivot_wider(
      data = .x,
      id_cols = ensembl_gene_id, 
      names_from = timepoint, 
      values_from = m_log_fc,
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
  .x = phospho_gene_sets,
  .f = ~ pathview(
    gene.data = .x,
    out.suffix = .y,
    pathway.id =  pathways,
    species = 'mmu',
    gene.idtype = 'ENSEMBL',
    kegg.native = TRUE,
    multi.state = TRUE,
    limit = list(
      cpd = 1
    )
  )
)

# Move phospho results to appropriate directory and clean up leftover files

dir.create( path = here::here(snakemake@output[['phospho_genes']]))

list.files(path = here::here(""), pattern = "mmu.*multi.png") %>% 
  file.rename(
    from = .,
    to = here::here( snakemake@output[['phospho_genes']], .)
  )

list.files(path = here::here(""), pattern = "mmu.*") %>%
  file.remove(.)

##### Now for site-specific data  ####

phospho_sites <-
  snakemake@input[['phospho_sites']] %>%
  read_rds() 

# convert phospho data to protein to treatment matrix
phospho_site_sets <-
  tissues %>% 
  purrr::set_names(.) %>% 
  map(
    .f = ~ filter(.data = phospho_sites, tissue == .x)
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
  .x = phospho_site_sets,
  .f = ~ pathview(
    gene.data = .x,
    out.suffix = .y,
    pathway.id =  pathways,
    species = 'mmu',
    gene.idtype = 'ENSEMBL',
    kegg.native = TRUE,
    multi.state = TRUE,
    limit = list(
      cpd = 1
    )
  )
)

# Move phospho results to appropriate directory and clean up leftover files

dir.create( path = here::here(snakemake@output[['phospho_sites']]))

list.files(path = here::here(""), pattern = "mmu.*multi.png") %>% 
  file.rename(
    from = .,
    to = here::here( snakemake@output[['phospho_sites']], .)
  )

list.files(path = here::here(""), pattern = "mmu.*") %>%
  file.remove(.)

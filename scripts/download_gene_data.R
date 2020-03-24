## Download gene names and annotation

library(biomaRt)
library(magrittr)
library(stringr)

gene_names <- 
  snakemake@input[[1]] %>% 
  read.csv() %$% 
  X

gene_annotation <-
  getBM(
    mart = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = 'useast'),
    filters = "ensembl_gene_id", 
    attributes = c("ensembl_gene_id", "mgi_symbol", "description"),
    values =  str_extract(string = gene_names, pattern = "^[[:alnum:]]*")
  )

readr::write_rds(x = gene_annotation, path = snakemake@output[[1]])

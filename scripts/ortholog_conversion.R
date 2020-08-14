#### Load libraries and human/mouse marts ####

library(biomaRt)
library(tidyverse)
library(magrittr)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# add mouse ensembl IDs, and, for safety, gene symbols to the human ones from the untangel perl script
human.start <- 
  snakemake@input[[1]] %>%
  read_tsv()

# this is a bit off: you someteims get multiple answers for te same gene - we can only really deal with 1...
#solutoin may be to just skip those mwith multiples beacsue we dont know anythnig. 

# current strategy: use perl to sort this out - bascially remove rows where one shows up all teh time
# so we use R only to get the unique genes, then perl again. A bit dirty but...

unique_human_ensembl <- 
  human.start %$%
  c(secreted.ensembl.gene.human, interacting.receptor.ensembl.gene.human) %>% 
  unique()

all_genes_of_interest <- 
  getLDS(
    filters = "ensembl_gene_id", 
    values = unique_human_ensembl, 
    mart = human, 
    attributes = c("hgnc_symbol", "ensembl_gene_id"), 
    martL = mouse, 
    attributesL = c( "mgi_symbol", "ensembl_gene_id")
  ) %>% 
  set_names(
    c("human.symbol", "human_ensembl", "mouse_symbol", "mouse_ensembl")
  )

# and then just print that
all_genes_of_interest

write_tsv(all_genes_of_interest, snakemake@output[[1]])
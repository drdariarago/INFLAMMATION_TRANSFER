## Import string data in easy to parse format
library(tidyverse)
library(magrittr)
library(biomaRt)

# Import STRING data

string_data <-
  list.files(
    path = here::here("data/string_plots"), 
    pattern = ".*string_interactions.tsv", 
    full.names = T
  ) %>% 
  set_names( str_extract(., "[:alnum:]*_?[2,5]?_string") ) %>% 
  map(read_tsv) %>% 
  map(janitor::clean_names) %>%
  map(
    ~dplyr::select(.x, 
                   mgi_prot_node1 = number_node1, mgi_prot_node2 = node2, 
                   ensembl_prot_node1 = node1_string_id, ensembl_prot_node2 = node2_string_id,
                   experimental = experimentally_determined_interaction, 
                   text = automated_textmining, score = combined_score
    )
  ) %>% 
  map(
    ~ mutate( .x,
      across( contains("ensembl"), ~str_extract(.x, "ENSMUSP[:digit:]*") )
    )
  )

protein_ids <-
  string_data %>% 
  map( ~ c( pull(.x, ensembl_prot_node1), pull(.x, ensembl_prot_node2) ) ) %>% 
  reduce( .f = union)

## Get ensembl gene IDs corresponding to protein IDs

gene_annotation <-
  getBM(
    mart = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl"),
    filters = "ensembl_peptide_id",
    attributes = c("ensembl_gene_id", "ensembl_peptide_id", "mgi_symbol", "description"),
    values =  protein_ids
  )

# Merge with protein IDs
protein_to_gene_list <-
  string_data %>% 
  map(
    ~ left_join(
      x = .x,
      y = gene_annotation[,-4],
      by = c( "ensembl_prot_node1" = "ensembl_peptide_id")
    )
  ) %>% 
  map(
    ~ rename(.x,  "ensembl_gene_node1" = "ensembl_gene_id", "mgi_gene_node1" =  "mgi_symbol") 
  ) %>% 
  map(
    ~ left_join(
      x = .x,
      y = gene_annotation[,-4],
      by = c( "ensembl_prot_node2" = "ensembl_peptide_id")
    ) 
  ) %>% 
  map( 
    ~ rename(.x, "ensembl_gene_node2" = "ensembl_gene_id", "mgi_gene_node2" =  "mgi_symbol") 
  ) %>% 
  map( 
    ~ dplyr::select( .data = .x, 
      contains("mgi_gene"), contains("ensembl"), experimental, text, score
    )
  )

write_rds(protein_to_gene_list, "results/string_import/protein_to_gene_list.rds")

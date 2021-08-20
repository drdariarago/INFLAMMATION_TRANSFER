## Save table of raw gene TPM and mean values per sample type for exploration

library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)

experiment_data <-
  snakemake@input[['expression_results']] %>%
  readRDS() 

## Format as sample by gene tibble
tpm_data_frame <-
  experiment_data %>% 
  assay("abundance") %>% 
  t() %>% 
  as.data.frame %>% 
  rownames_to_column("sample_id")

## Join with sample metadata
sample_metadata <-
  colData(experiment_data) %>% 
  as.data.frame() %>% 
  mutate(
    tissue2 = paste(maternal_fetal, tissue, sep = "_"),
    tissue = str_replace(tissue2, "(maternal_placenta|fetal_placenta)", "placentas")
  ) %>% 
  select(
    sample_id = names, tissue, tissue2, exposure, timepoint
  )
  
full_data <-
  right_join(
    x = sample_metadata, y = tpm_data_frame,
    by = "sample_id"
  )

## Filter out TiO2 samples
filtered_data <-
  full_data %>% 
  filter(exposure != "TiO2")

## Replace ENSEMBL IDs with MGI symbols

# Load MGI symbols 
mgi_data <-
  snakemake@input[['mgi_symbols']] %>% 
  read_rds()

# Make longer tibble of gene counts
long_data <-
  filtered_data %>% 
  as_tibble() %>% 
  pivot_longer(
    cols = starts_with("ENSMUSG"), 
    names_to = "ensembl_gene_id", 
    values_to = "TPM"
  )

# Join with MGI symbols
data_with_symbols <-
  long_data %>% 
  mutate(
    ensembl_gene_id = gsub(
      pattern = "\\..*", 
      x = ensembl_gene_id, 
      replacement = ""
    )
  ) %>% 
  left_join(
    x = ., y = mgi_data[,c("ensembl_gene_id", "mgi_symbol")],
    by = "ensembl_gene_id"
  )

write_rds(x = data_with_symbols, file = snakemake@output[['tpm_results']])

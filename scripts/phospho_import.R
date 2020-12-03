### Import phosphoproteomics data

library(tidyverse)
library(magrittr)
library(readxl)
library(biomaRt)

## Read files in

results_path <- 
  here::here('data/proteomics/Quant_TMM_EdgeR_4_experiments.xlsx')

phospho_results <-
  results_path %>% 
  excel_sheets() %>% 
  setNames(object = ., nm = .) %>% 
  map_dfr(
    .id = "contrast_id",
    .x = .,
    .f = ~ read_xlsx(
      results_path,
      sheet = .x
      ) %>% 
      rename_with(
        .cols = matches(".*[2,5,12,24]H.*"), 
        .fn = ~ gsub(.x, pattern = "[2,5,12,24]H_", replacement = "")
      )
  )

# Import conversion table from UNIPROT IDs to gene symbols and ENSEMBL IDs

protein_id <-
  getBM(
    mart = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl"),
    filters = "uniprot_gn_id", 
    attributes = c("ensembl_gene_id", "mgi_symbol", "uniprot_gn_id"),
    values =  phospho_results$Uniprot %>% unique
  )

# Export fold change by tissue:timepoint for graphviz plots

phospho_results %>% 
  dplyr::select(Uniprot, logFC, FDR, contrast_id) %>% 
  right_join(protein_id, ., by = c("uniprot_gn_id" = "Uniprot")) %>% 
  mutate(
    tissue = str_extract(string = contrast_id, pattern = "^[:alpha:]*") ,
    tissue = factor(
      x = tissue,
      levels = c("Liver", "Placenta"), 
      labels = c("fetal_liver", "placentas")
    ) ,
    timepoint = str_extract(string = contrast_id, pattern = "[0-9]") ,
    timepoint = paste0("timepoint", timepoint) %>% as.factor ,
    contrast_id = NULL,
    uniprot_id = uniprot_gn_id
  ) %>% 
  write_csv(x = ., file = here::here('results/phospho_import/fold_changes.csv'))

## Plot phospho raw data
phospho_results %>% 
  filter(
    grepl(pattern = 'Liver', contrast_id)
  ) %>% 
  magrittr::extract(,2:9) %>%
  as.matrix() %>% 
  apply(., c(1,2), log10) %>% 
  pheatmap::pheatmap(
    mat = .,
    scale = 'column'
  )

## Import raw data

data_path <-
  here::here('data/proteomics/All_tables_4_experiments.xlsx')

phospho_data <-
  data_path %>% 
  excel_sheets() %>% 
  setNames(object = ., nm = .) %>% 
  map_dfr(
    .id = "sample_id",
    .x = .,
    .f = ~ read_xlsx(
      data_path,
      sheet = .x) %>% 
      rename_with(
        .cols = matches(".*[2,5,12,24]H.*"), 
        .fn = ~ gsub(.x, pattern = "[2,5,12,24]H_", replacement = "")
      )
  )

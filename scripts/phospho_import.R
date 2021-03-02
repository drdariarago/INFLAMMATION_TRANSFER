### Import phosphoproteomics data

library(tidyverse)
library(magrittr)
library(readxl)
library(biomaRt)

## Read files in

results_path <- 
  snakemake@input[['results']]

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
      dplyr::select(
        !starts_with("RI")
      )
  ) %>% 
  mutate(
    tissue = str_extract(string = contrast_id, pattern = "^[[:alpha:]]*") %>% 
      factor(
        levels = c("Liver", "Placenta"), 
        labels = c("fetal_liver", "placentas")
      ),
    timepoint = str_extract(string = contrast_id, pattern = "[2,5]") %>% as.factor(),
    uniprot = str_extract(string = Acc, pattern = "[[:upper:],[:digit:]]{6}")
  ) 

# Import conversion table from UNIPROT IDs to gene symbols and ENSEMBL IDs

protein_id <-
  getBM(
    mart = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl"),
    filters = "uniprot_gn_id", 
    attributes = c("ensembl_gene_id", "mgi_symbol", "uniprot_gn_id"),
    values =  phospho_results$uniprot %>% unique
  )

# Export fold change by tissue:timepoint for graphviz plots

annotated_phospho_results <-
  phospho_results %>% 
  dplyr::select(Acc, uniprot, tissue, timepoint, logFC, p_value = PValue, FDR, m_ctrl = ave_Control, m_LPS = ave_LPS) %>% 
  left_join(protein_id, ., by = c("uniprot_gn_id" = "uniprot"))

annotated_phospho_results %>% 
  dplyr::select(-p_value) %>% 
  write_rds(x = ., file = snakemake@output[['sitewise_results']])

# Export values binned by gene, averaging p-values by Fisher's method and then re-running FDR correction

pdf(file = snakemake@output[["genewise_fdr_plot"]])

genewise_annotated_phospho_results <-
  annotated_phospho_results %>% 
  dplyr::select(-Acc) %>% 
  group_by(tissue, timepoint, ensembl_gene_id, mgi_symbol) %>% 
  summarise(
    m_ctrl = mean(m_ctrl),
    m_lps = mean(m_LPS),
    m_log_fc = log2(m_lps/m_ctrl),
    n_sites = n(),
    n_sign_sites = sum(p_value < 0.05),
    p_value = ifelse(
      n() > 1,
      metap::sumlog(p_value) %$% p,
      p_value
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    q_value = fdrtool::fdrtool(x = p_value, statistic = "pvalue") %$% qval
  )

dev.off()

genewise_annotated_phospho_results %>% 

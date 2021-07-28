## Select genes of interest from phospho data and compare with RNAseq expression over time
library(tidyverse)
library(magrittr)

## Import phospho and rnaseq data

phospho_data <-
  here::here("results/phospho_import/genewise_results.rds") %>% 
  read_rds() %>% 
  select(
    tissue, timepoint, ensembl_gene_id, mgi_symbol, phospho_log_fc = min_log_fc, phospho_q_value = q_value_min
  )

rnaseq_data <-
  here::here("results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>% 
  read_csv() %>% 
  select(
    tissue, timepoint, ensembl_gene_id, mgi_symbol, rna_log_fc = logFC, rna_q_value = q_value, exposure
  ) %>% 
  pivot_wider(names_from = exposure, values_from = c(rna_log_fc, rna_q_value) ) %>% 
  select( -rna_q_value_baseline, 
          rna_baseline  = rna_log_fc_baseline, 
          rna_log_fc = rna_log_fc_response, 
          rna_q_value = rna_q_value_response 
  ) %>% 
  mutate(
    timepoint = factor(timepoint, levels = paste0( "timepoint", c(2,5,12,24) ) )
  )

# Merge and compare to each other

matched_data <-
  rnaseq_data %>% 
  filter(
    tissue %in% c("fetal_liver", "placentas")
  ) %>% 
  mutate(
    timepoint = str_extract(string = timepoint, pattern = "[0-9]{1,2}$")
  ) %>% 
  right_join(
    x = .,
    y = phospho_data,
    by = c("ensembl_gene_id" = "ensembl_gene_id", "tissue" = "tissue", "mgi_symbol" = "mgi_symbol", "timepoint" = "timepoint")
  )

matched_data %>% 
  ggplot(
    aes( x = rna_log_fc, y = phospho_log_fc, label = mgi_symbol)
  ) +
  geom_hex() +
  geom_point(
    data = filter(
      matched_data, 
      phospho_q_value < 0.05, 
      rna_q_value < 0.05, 
      abs(phospho_log_fc) > 0.25,
      abs(rna_log_fc) > 0.25
    ),
    col = 'hotpink', alpha = 0.5 
  ) +
  ggrepel::geom_text_repel(
    data = filter(
      matched_data, 
      phospho_q_value < 0.05, 
      rna_q_value < 0.05, 
      abs(phospho_log_fc) > 0.25,
      abs(rna_log_fc) > 0.25
    ),
    col = 'darkorange', alpha = 1,
    max.overlaps = 20
    ) +
  facet_grid(tissue ~ timepoint) +
  scale_fill_viridis_b(trans = "log") +
  coord_cartesian(xlim = c(-1, 0.8))

# Load interesting liver genes and check RNA expression over time

liver_genes <-
  c("DDX21", "EP300", "HIF3A", "KI67", "SLBP", "CCNE2", "RBM8A", "SRRM1", "WBP4")
  
rnaseq_data %>%
  filter(
    tolower(mgi_symbol) %in% tolower(liver_genes)
  ) %>% 
  ggplot(
    aes( x = timepoint, y = rna_log_fc, group = ensembl_gene_id, col = rna_q_value < 0.05, label = mgi_symbol )
  ) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(
    data = filter(rnaseq_data, 
                  tolower(mgi_symbol) %in% tolower(liver_genes), 
                  rna_q_value < 0.05, 
                  timepoint == "timepoint5" 
    )
  ) +
  facet_wrap(~ tissue) 

# Import interaction pairs from string plots

string_data <-
  here::here("results/string_import/protein_to_gene_list.rds") %>% 
  read_rds()

# Check expression over time of liver genes
matched_data %>% 
  filter(ensembl_gene_id %in% liver_string$ensembl_node1)


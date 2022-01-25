### Import phosphoproteomics data

library(tidyverse)
library(magrittr)
library(readxl)
library(biomaRt)

## Read files in

results_path <- 
  snakemake@input[['results']]
  # here::here("data/proteomics/20210202_second_run/Analysed_data_TMM_EdgeR_4_experiments__Jan_2021.xlsx")

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
    mart = useEnsembl(
      biomart = "ensembl", 
      dataset = "mmusculus_gene_ensembl",
      mirror = "www"
      ),
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
  write_rds(x = ., file = snakemake@output[['sitewise_results']])

# Export values binned by gene, averaging p-values by Fisher's method and then re-running FDR correction

pdf(file = snakemake@output[["genewise_fdr_plot"]])

genewise_annotated_phospho_results <-
  annotated_phospho_results %>% 
  dplyr::select(-Acc) %>% 
  group_by(tissue, timepoint, ensembl_gene_id, mgi_symbol) %>% 
  arrange(.by_group = TRUE, p_value) %>% 
  summarise(
    site_mad = ifelse(
      n() > 1,
      mad(log10(m_LPS/m_ctrl)),
      NA
    ),
    min_log_fc = log2(first(m_LPS)/first(m_ctrl)),
    min_average = (first(m_LPS) + first(m_ctrl))/2,
    p_value_min = first(p_value),
    n_sites = n(),
    n_sign_sites = sum(p_value < 0.05),
    p_value_avg = ifelse(
      n() > 1,
      metap::sumlog(p_value) %$% p,
      p_value
    )
  ) %>% 
  ungroup() %>% 
  mutate(
    q_value_avg = fdrtool::fdrtool(x = p_value_avg, statistic = "pvalue") %$% qval,
    q_value_min = fdrtool::fdrtool(x = p_value_min, statistic = "pvalue") %$% qval
  )

dev.off()

# Save lists of signed averaged and min q-values for enrichment testing



# Check average distance from the mean fold-change in phosphorilation within each gene

genewise_annotated_phospho_results %>% 
  ggplot( aes (x = n_sites, y = site_mad)) +
  geom_hex() + 
  scale_x_log10() + 
  geom_smooth() + 
  viridis::scale_fill_viridis(trans = 'log') +
  facet_wrap( ~ n_sign_sites > 0)

genewise_annotated_phospho_results %>% 
  ggplot( aes (x = n_sign_sites > 0, y = site_mad, group = n_sign_sites > 0)) +
  geom_boxplot(varwidth = T, notch = T) + 
  geom_smooth() + 
  viridis::scale_fill_viridis(trans = 'log') + 
  facet_wrap(~ cut_interval(log10(n_sites), length = .5)) +
  ggtitle("Mean Absolute Deviation in logFC of sites within genes, split by significance \n 
          facets are by log10 number of sites within a gene")

genewise_annotated_phospho_results %>% 
  write_rds(x = ., file = snakemake@output[["genewise_results"]])

### Explore limma results

library(tidyverse)

results <-
  here::here("results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>% 
  read_csv()

wide_results <-
  results %>% 
  select(-contrast) %>% 
  pivot_wider(
    id_cols = c(tissue, ensembl_gene_id, mgi_symbol, timepoint), 
    names_from = exposure, values_from = c(q_value, logFC)
  ) %>% 
  select(-q_value_baseline) %>% 
  mutate(
    timepoint = factor(x = timepoint, levels = paste0("timepoint", c(2,5,12,24), sep = ""))
  )

## Explore fetal liver

wide_results %>% 
  filter(
    tissue == "fetal_liver",
    q_value_response < 0.05,
    logFC_response > 1,
    logFC_baseline > 1
  ) 


# Filter top 10 overexpressed at any stage and plot over time
top_10_table <-
  wide_results %>%
  filter(
    tissue %in% c("fetal_liver", "placentas")
  ) %>% 
  group_by(tissue, timepoint) %>% 
  filter(
    q_value_response < 0.05,
    logFC_baseline > 0.5
  ) %>% 
  arrange(-logFC_response) %>% 
  slice_head(
    n = 10
  ) 

top_10_genes <-
  top_10_table %>% 
  select(
    sign_tissue = tissue, sign_timepoint = timepoint, ensembl_gene_id
  )

top_10_tibble <-
  inner_join(
    x = top_10_genes,
    y = wide_results
  ) %>% 
  filter(
    tissue == sign_tissue
  )


top_10_tibble %>% 
  filter (sign_tissue == "placentas") %>% 
  ggplot(
    aes(x = timepoint, y = logFC_response, group = mgi_symbol, col = mgi_symbol)
  ) + 
  geom_point(aes(shape = q_value_response < 0.05)) +
  facet_wrap(~ sign_timepoint) +
  geom_line()

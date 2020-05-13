## Reformat GO enrichment results to human-readable and plot friendly format
library(tidyverse)
library(magrittr)

gost_results_table <-
  snakemake@input[[1]] %>%
  read_rds() %>%
  map_df(.f = ~ .x$result, .id = "contrast"
  ) %>% 
  dplyr::select(contrast, term_id, term_name, significant, source, term_size, intersection_size
  ) %>% 
  group_by(term_id) %>% 
  filter(any(significant == TRUE)) %>% 
  mutate_at(
    .vars = vars(starts_with("intersection")),
    .funs = ~ round(.x * 100 / term_size, digits = 1)
  ) %>% 
  dplyr::select(
    contrast, source, term_id, term_name, term_size, significant, intersection_percent = intersection_size
  ) 

write_csv(x = gost_results_table, path = snakemake@output[['go_long_format']])

gost_results_table_wider <-
  gost_results_table %>% 
  pivot_wider(
    names_from = contrast, values_from = c(significant, intersection_percent)
  ) 

write_csv(x = gost_results_table_wider, path = snakemake@output[['go_wide_format']])  
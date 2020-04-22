## Reformat GO enrichment results to human-readable and plot friendly format
library(tidyverse)
library(magrittr)

gost_results_table <- 
  snakemake@input[[1]] %>% 
  read_rds() %>% 
  map_df(.f = ~ .x$result, .id = "contrast") %>% 
  dplyr::select(contrast, term_id, term_name, significant, source, term_size, query_sizes, intersection_sizes) %>% 
  mutate_if(
    .predicate = is.list, 
    .funs = ~ map_chr(.x = ., .f = paste0, collapse = ":" )
  ) %>% 
  separate(
    col = significant,
    sep = ":",
    convert = TRUE,
    into = paste0("signif_", c(2,5,12,24), "_hrs")
  ) %>%  separate(
    col = query_sizes,
    sep = ":",
    convert = TRUE,
    into = paste0("query_size_", c(2,5,12,24), "_hrs")
  ) %>%  separate(
    col = intersection_sizes,
    sep = ":",
    convert = TRUE,
    into = paste0("intersection_size_", c(2,5,12,24), "_hrs")
  ) %>% mutate_at(
    .vars = vars(starts_with("intersection")),
    .funs = ~ round(.x * 100 / term_size, digits = 1)
  ) %>% dplyr::select(
    contrast, source, term_id, term_name, term_size, starts_with("signif"), starts_with("intersect")
  ) %>% pivot_longer(
    cols = starts_with("signif"), names_to = "timepoint", 
    names_prefix = "signif_", values_to = "significant"
  ) %>% pivot_longer(
    cols = starts_with("intersect"), names_to = "timepoint_2",
    names_prefix = "intersection_size_", values_to = "percent_enrichment"
  ) %>% dplyr::filter(
    timepoint == timepoint_2
  ) %>% dplyr::select(
    -timepoint_2
  ) %>% mutate(
    timepoint = factor(x = timepoint, levels = paste(c(2,5,12,24), "hrs", sep = "_"))
  )

write_csv(x = gost_results_table, path = snakemake@output[[1]])

library(here)
library(tidyverse)
library(magrittr)
library(gprofiler2)

# Import fold change values
go_results <-
  readRDS(here::here("results/gost_enrichment/enrichment_results.Rdata")) %>% 
  map_df(.f = ~ .x$result, .id = "contrast")

go_results_table <- 
  go_results %>% 
  dplyr::select(contrast, term_id, significant, source, term_name, query_sizes, intersection_sizes) %>% 
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
    into = paste0("total_annotated_", c(2,5,12,24), "_hrs")
  ) %>%  separate(
    col = intersection_sizes,
    sep = ":",
    convert = TRUE,
    into = paste0("significant_annotated_", c(2,5,12,24), "_hrs")
  )

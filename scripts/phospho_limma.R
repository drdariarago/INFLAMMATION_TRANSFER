### Apply standard Limma/Voom analysis to phospho data

library(readxl)
library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)

## Import results as single data.frame

results_path <-
  here::here('data/proteomics/20210202_second_run/MaxQuant_table_all4exp_Jan_2021.xlsx')

results_list <- 
  results_path %>% 
  readxl::excel_sheets() %>% 
  purrr::set_names() %>%
  map(.f = readxl::read_xlsx, path = results_path)

all_sites <-
  results_list %>% 
  map( .f = ~ .x %$% LP_position ) %>% 
  reduce( .f = union) %>% 
  data.frame( LP_position = .)

results_matrix <-
  results_list %>% 
  map(
    ~ select(.x, LP_position, starts_with(match = "Reporter") ) %>% 
      set_colnames( gsub( x = colnames(.), pattern = "Reporter Intensity_", replacement = "")) %>% 
      column_to_rownames( var = "LP_position")
  ) %>% 
  imap(
    ~ str_extract(string = .y, pattern = "(P[L,l]acenta|Liver)") %>% 
      tolower %>% 
      { set_colnames(.x, paste( ., colnames(.x), sep = "_" )) }
  ) %>% 
  map_dfc(
    ~ rownames_to_column(.x, var = "LP_position") %>% 
      left_join( all_sites, .) %>% 
      column_to_rownames( var = "LP_position")
  ) %>% 
  mutate( across( .cols = everything(), .fns = ~ replace( .x, is.na(.x), 0 ) ) ) %>% 
  as.matrix()

sample_info <-
  results_matrix %>%
  colnames %>%
  str_extract_all( pattern = "[[:alnum:]]{1,}", simplify = T ) %>%
  set_colnames( c( "tissue", "treatment", "time", "replicate" ) ) %>%
  set_rownames( colnames( results_matrix ) )

site_info <-
  results_matrix %>% 
  rownames %>% 
  str_extract_all( pattern = "[[:alnum:]]{1,}", simplify = T ) %>% 
  magrittr::extract(, 2:5) %>% 
  set_colnames( c( "uniprot_id", "mgi_id", "species", "position")) %>% 
  set_rownames( results_matrix %>% rownames )

## Convert to limma object

experiment_data <-
  DGEList(
    counts = results_matrix, 
    samples = sample_info, 
    genes = site_info, 
    group = 
      sample_info %>% as_tibble %>% 
      transmute( paste( tissue, treatment, time, sep = "_" ) ) %>%
      pull
  ) %>% 
  calcNormFactors(method = "TMM")

design <-
  # Set independent tissue baselines, 
  # then check for shared response to treatment
  # then check if tissue baseline changes over time
  # then check if responses change over time
  model.matrix(
    object = formula( ~ 0 + tissue / ( treatment * time ) ), 
    data = sample_info %>% as.data.frame
  )

## Filter low counts

results_plot_data <-
  results_matrix %>% 
  as_tibble( rownames = "site_id" ) %>%
  pivot_longer( 
    cols = contains("H"), 
    names_to = c( "tissue", "treatment", "time", "replicate"), 
    names_sep = "_", 
    values_to = "counts"
  ) 

results_plot_data %>% 
  ggplot(
    aes( x = counts + 100 , col = tissue, lty = treatment)
  ) +
  geom_density() + 
  scale_x_log10() +
  geom_vline( xintercept = 8E4, col = "hotpink", lty = 3 ) +
  geom_vline( xintercept = 2E3, col = "blue", lty = 3 )

results_plot_data %>% 
  group_by( site_id, tissue ) %>% 
  summarise( MAD = mad(counts) )

## Voom

experiment_data %>% 
  filterByExpr(
    y = ., 
    design = design, 
    min.count = 64
  ) %>% 
  experiment_data[., , keep.lib.sizes = T] %>% 
  voom(
    counts = .,
    design = design,
    lib.size = .$lib.size,
    plot = TRUE
  )

## Limma 

## fdr
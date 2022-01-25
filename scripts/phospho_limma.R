### Apply standard Limma/Voom analysis to phospho data

library(readxl)
library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)
library(patchwork)

## Import results as single data.frame

input_path <-
  snakemake@input[['data']]

results_path <-
  snakemake@params[['output_path']]

results_list <- 
  input_path %>% 
  readxl::excel_sheets() %>% 
  str_subset( pattern = 
                ifelse(
                  grepl(pattern = "liver", x = results_path,),
                  "Liver", "P[L,l]acenta")
  ) %>% 
  purrr::set_names() %>%
  map(.f = readxl::read_xlsx, path = input_path)

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
  set_colnames( c( "treatment", "time", "replicate" ) ) %>%
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
      transmute( paste( treatment, time, sep = "_" ) ) %>%
      pull
  ) %>% 
  calcNormFactors(method = "TMM")

design <-
  # check if tissue baseline changes over time
  # then check if responses are significant at each timepoint
  model.matrix(
    object = formula( ~ time / treatment ), 
    data = sample_info %>% as.data.frame
  )

## Filter low counts

results_plot_data <-
  results_matrix %>% 
  as_tibble( rownames = "site_id" ) %>%
  pivot_longer( 
    cols = contains("H"), 
    names_to = c( "treatment", "time", "replicate"), 
    names_sep = "_", 
    values_to = "counts"
  ) 

results_plot_data %>% 
  ggplot(
    aes( x = counts + 100 , lty = treatment )
  ) +
  geom_density() + 
  scale_x_log10() +
  geom_vline( xintercept = 8E4, col = "hotpink", lty = 3 ) +
  geom_vline( xintercept = 1E3, col = "blue", lty = 3 )

results_summary <-
  results_plot_data %>% 
  group_by(site_id) %>% 
  summarise( 
    MAD = mad(counts),
    average = mean(counts) 
  )

results_summary %>% 
  ggplot(
    aes( x = average, y = MAD)
  ) +
  geom_hex() + 
  geom_smooth(method = 'lm') +
  scale_x_log10() +
  scale_y_log10() +
  viridis::scale_fill_viridis( trans = 'log', breaks = c( 2, 10, 100, 500) )

## Voom

pdf( file = here::here(results_path, "voom_plot.pdf") )

filtered_experiment_data <-
  experiment_data %>% 
  filterByExpr(
    y = ., 
    design = 
      model.matrix(
        object = formula( ~ 1 ),
        data = sample_info %>% as.data.frame
      ),
    min.count = 5E3,
    min.prop = .5
  ) %>% 
  experiment_data[., , keep.lib.sizes = T] %>% 
  voom(
    counts = .,
    design =  design, 
    lib.size = .$lib.size,
    plot = TRUE
  )

dev.off()

## Limma 

linear_model <-
  lmFit(
    object = filtered_experiment_data, 
    design = filtered_experiment_data$design
  ) %>% 
  eBayes()

write_rds(x = linear_model, file = here::here( results_path, "linear_model.rds" ) )

## Volcano and MA plots

volcano_2H <-
  linear_model %>% 
  topTable(coef = "time2H:treatmentLPS", number = Inf, adjust.method = 'none', sort.by = 'P') %>% 
  as_tibble() %>% 
  ggplot(
    aes( x = logFC, y = P.Value)
  ) + 
  geom_hex() + 
  geom_hline( yintercept = 0.05) +
  scale_y_log10() +
  viridis::scale_fill_viridis( trans = 'log', breaks = 2^c( 0:5 * 2) )

volcano_5H <-
  linear_model %>% 
  topTable(coef = "time5H:treatmentLPS", number = Inf, adjust.method = 'none', sort.by = 'P') %>% 
  as_tibble() %>% 
  ggplot(
    aes( x = logFC, y = P.Value)
  ) + 
  geom_hex() + 
  geom_hline( yintercept = 0.05) +
  scale_y_log10() +
  viridis::scale_fill_viridis( trans = 'log', breaks = 2^c( 0:5 * 2) )

MA_2H <- 
  linear_model %>% 
  topTable(coef = "time2H:treatmentLPS", number = Inf, adjust.method = 'none', sort.by = 'P') %>% 
  as_tibble() %>% 
  ggplot(
    aes( x = AveExpr, y = logFC)
  ) + 
  geom_hex() +
  viridis::scale_fill_viridis( trans = 'log', breaks = 2^c( 0:5 * 2) )

MA_5H <- 
  linear_model %>% 
  topTable(coef = "time5H:treatmentLPS", number = Inf, adjust.method = 'none', sort.by = 'P') %>% 
  as_tibble() %>% 
  ggplot(
    aes( x = AveExpr, y = logFC)
  ) + 
  geom_hex() +
  viridis::scale_fill_viridis( trans = 'log', breaks = 2^c( 0:5 * 2) )

(volcano_2H + MA_2H) / (volcano_5H + MA_5H) +
  plot_annotation( title = glue::glue("Volcano and MA plots for {snakemake@params[['tissue']]}") )

ggsave( filename = here::here(results_path, "volcano_MA_plots.pdf"), units = 'mm', width = 297, height = 210 )

## Calculate fdr for treatment::time

q_values_meta <-
  linear_model %$%
  p.value %>% 
  as_tibble(x = ., rownames = 'site_id') %>% 
  select( site_id, contains('treatment') ) %>% 
  pivot_longer(data = ., cols = -contains('site'), names_to = 'contrast', values_to = 'p_value') %>% 
  mutate( gene_id = str_extract( string = site_id, pattern = "[:alnum:]*(?=_MOUSE)" ) ) %>% 
  group_by( gene_id, contrast ) %>% 
  summarise( meta_p = ifelse( n() > 1, metap::sumlog( p_value )$p, first( p_value) ) ) %>% 
  ungroup %>% 
  mutate( q_value = fdrtool::fdrtool(x = meta_p, statistic = "pvalue") %$% lfdr )

q_values_base <-
  linear_model %$%
  p.value %>% 
  as_tibble(x = ., rownames = 'site_id') %>% 
  select( site_id, contains('treatment') ) %>% 
  pivot_longer(data = ., cols = -contains('site'), names_to = 'contrast', values_to = 'p_value') %>% 
  mutate( gene_id = str_extract( string = site_id, pattern = "[:alnum:]*(?=_MOUSE)" ) ) %>% 
  mutate( q_value = fdrtool::fdrtool(x = p_value, statistic = "pvalue") %$% lfdr )

q_values_meta %>% 
  group_by( contrast ) %>% 
  summarise(
    sign_hits = sum(q_value < 0.05)
  )

# Merge with logFC values

filtered_q_values<-
  q_values_meta %>% 
  select( -meta_p ) %>% 
  pivot_wider(
    names_from = contrast, values_from = q_value
  ) %>% 
  select(
    gene_id = gene_id,
    response_2H = `time2H:treatmentLPS`,
    response_5H = `time5H:treatmentLPS`
  ) %>% 
  filter( response_2H < 0.05 | response_5H < 0.05 )

filtered_results_table <-
  topTable(linear_model, sort.by = "none", number = Inf) %>% 
  as_tibble() %>% 
  select( uniprot_id, mgi_id, position, AveExpr, contains("treatment")) %>% 
  inner_join(., filtered_q_values, by = c("mgi_id" = "gene_id")) %>% 
  arrange( mgi_id, position) %>% 
  select(
    uniprot_id, mgi_id, 
    p_pos = position, 
    average_expression = AveExpr,
    response_2h_q_value = response_2H,
    response_2h_log_fc = time2H.treatmentLPS,
    response_5h_q_value = response_5H,
    response_5h_log_fc = time5H.treatmentLPS,
  ) %>% 
  arrange( -min(response_2h_q_value, response_5h_q_value) )

write_csv(
  x = filtered_results_table, 
  file = here::here(results_path, "significant_results.csv")
  )


full_join(
  filtered_results_table %>% 
    select( mgi_id, response_2h_q_value ) %>% 
    filter( response_2h_q_value < 0.05 ) %>% 
    arrange( response_2h_q_value ) %>% 
    distinct(), 
  filtered_results_table %>% 
    select( mgi_id, response_5h_q_value ) %>% 
    filter( response_5h_q_value < 0.05 ) %>% 
    arrange( response_5h_q_value ) %>% 
    distinct() 
) %>% 
  write_csv( x = ., file = snakemake@output[['significant_genes']] )
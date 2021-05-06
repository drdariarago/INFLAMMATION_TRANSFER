## Compare placental response to LPS-induced inflammation

library(tidyverse)
library(magrittr)

q_value_threshold = 0.05
log_fc_threshold = 0.25
a4_short_long = c(210, 297)

## Import DE genes from inflammation study and from our study
de_inflammation_results <-
  "results/process_lien_results/fold_change_summary.rds" %>% 
  here::here() %>% 
  read_rds() 

de_transfer_results <-
  "results/limma_compile_results/limma_results_no_maternal_contrasts.csv" %>% 
  here::here() %>% 
  read_csv() %>%
  mutate(
    timepoint = factor(x = timepoint, levels = paste0("timepoint", c(2,5,12,24)))
  ) %>% 
  {list(
    placenta = filter(
      .data = .,
      tissue == "placentas"
    ),
    lungs = filter(
      .data = ., 
      tissue == "maternal_lung"
    )
  )}
  

# Perform rank-correlation between the two studies
# Filter only genes in both studies with FDR >= 0.05

shared_genes <-
  map(
    .x = de_transfer_results, 
    .f = ~ intersect(
      .x$mgi_symbol, 
      de_inflammation_results$mgi_symbol
    ) 
  )

filtered_de_transfer <-
  de_transfer_results %>% 
  map2(.x = ., .y = shared_genes, .f = ~ filter(.data = .x, mgi_symbol %in% .y)) %>% 
  map(.x = ., .f = select, ensembl_gene_id, mgi_symbol, q_value, log_fc = logFC, exposure, timepoint) %>% 
  map(.x = ., .f = pivot_wider, names_from = exposure, values_from = c(q_value, log_fc)) %>%
  map(.x = ., .f = select, 
      ensembl_gene_id, mgi_symbol, timepoint,
      ave_expr = log_fc_baseline, q_value = q_value_response, log_fc = log_fc_response)

filtered_de_inflammation <-
  de_inflammation_results %>% 
  select(
    mgi_symbol, q_value, ave_expr, log_fc
  ) 

shared_de <-
  map(
    .x = filtered_de_transfer,
    .f = left_join,
    y = filtered_de_inflammation,
    by = "mgi_symbol",
    suffix = c("_transfer", "_inflamed")
  )

# Filter only genes with significant logFC in at least one condition
filtered_shared_de <-
  map(
    .x = shared_de,
    .f = filter,
    (q_value_transfer < q_value_threshold) | (q_value_inflamed < q_value_threshold),
    (abs(log_fc_transfer) > log_fc_threshold) | (abs(log_fc_inflamed) > log_fc_threshold)
  ) %>% 
  map(
    .x = .,
    .f = mutate,
    de_class = case_when(
      (q_value_transfer < q_value_threshold) & (q_value_inflamed < q_value_threshold) ~ "both",
      (q_value_transfer < q_value_threshold) ~ "transfer_only",
      (q_value_inflamed < q_value_threshold) ~ "inflamed_only",
      TRUE ~ 'not_significant'),
    alpha = (log_fc_inflamed / 5) ^2 +  log_fc_transfer^2  %>% sqrt
  )
  
write_rds(
  x = filtered_shared_de, 
  file = "results/inflammation_comparison/filtered_shared_genes_list.rds" %>% here::here()
  )

# Calculate and compare baseline expression in inflamed and indirect placenta
filtered_shared_de[[1]] %>% 
  ggplot(
    aes(
      y = ave_expr_inflamed, x = ave_expr_transfer, 
      label = mgi_symbol
    )
  ) +
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, col = 'hotpink') +
  geom_smooth(method = 'lm', col = 'magenta')  +
  facet_grid( de_class ~ timepoint) + 
  viridis::scale_fill_viridis(trans = 'log', breaks = c(2,4,8,16,32,64,128), name = "number\nof genes") +
  theme(legend.position = 'right') +
  ggtitle(label = "Baseline expression of indirect and direct placenta inflammation at each timepoint")

ggsave(
  filename = "results/inflammation_comparison/placenta_baseline_plot.pdf" %>% here::here(), 
  width = a4_short_long[2], height = a4_short_long[1],
  units = 'mm'
)

## Same but with fold change expression values

filtered_shared_de[[1]] %>% 
  ggplot(
    aes(
      y = log_fc_inflamed, x = log_fc_transfer, 
      label = mgi_symbol
    )
  ) +
  geom_hex() +
  viridis::scale_fill_viridis(trans = 'log', breaks = c( 2^(2*0:5) )) +
  geom_abline(slope = 1, intercept = 0, col = 'hotpink') +
  geom_smooth(method = 'lm', col = 'magenta')  +
  facet_grid( de_class ~ timepoint) + 
  scale_y_continuous(limits = c(-3, 7.5)) +
  theme(legend.position = 'right') +
  ggrepel::geom_text_repel(
    data = filtered_shared_de[[1]] %>% 
      filter( 
        ( abs(log_fc_transfer) > 0.3 & abs(log_fc_inflamed) > 1.5 & de_class == 'both' ) |
          ( abs(log_fc_transfer) > 0.7 & de_class == 'transfer_only' ) |
          ( abs(log_fc_inflamed) > 4 & de_class == 'inflamed_only' )
      ), 
    col = 'darkorange'
  ) +
  ggtitle(label = "logFoldChange of indirect and direct placenta inflammation at each timepoint")

ggsave(
  filename = "results/inflammation_comparison/placenta_response_plot.pdf" %>% here::here(),
  width = a4_short_long[2], height = a4_short_long[1],
  units = 'mm'
)

## Compare baseline expression of inflamed placenta and inflamed lung

# Calculate and compare baseline expression in inflamed and indirect placenta
filtered_shared_de[[2]] %>% 
  ggplot(
    aes(
      y = ave_expr_inflamed, x = ave_expr_transfer, 
      label = mgi_symbol
    )
  ) +
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, col = 'hotpink') +
  geom_smooth(method = 'lm', col = 'magenta')  +
  facet_grid( de_class ~ timepoint) + 
  viridis::scale_fill_viridis(trans = 'log', breaks = c(2^c(1:7)), name = "number\nof genes") +
  theme(legend.position = 'right') +
  ggtitle(label = "Baseline expression of inflamed placenta and lung at each timepoint")

ggsave(
  filename = "results/inflammation_comparison/lung_baseline_plot.pdf" %>% here::here(),
  width = a4_short_long[2], height = a4_short_long[1],
  units = 'mm'
)

## Compare expression of inflamed placenta and inflamed lung

filtered_shared_de[[2]] %>% 
  ggplot(
    aes(
      y = log_fc_inflamed, x = log_fc_transfer, 
      label = mgi_symbol
    )
  ) +
  geom_hex()+
  geom_abline(slope = 1, intercept = 0, col = 'hotpink') +
  geom_smooth(method = 'lm', col = 'magenta')  +
  facet_grid( de_class ~ timepoint) + 
  theme(legend.position = 'right') +
  viridis::scale_fill_viridis(trans = 'log', breaks = c(3^c(1:8)), name = "number\nof genes") +
  ggrepel::geom_text_repel(
    data = filtered_shared_de[[2]] %>%
      filter(
        ( abs(log_fc_transfer) > 2.5 & abs(log_fc_inflamed) > 2.5 & de_class == 'both' ) |
          ( abs(log_fc_transfer) > 2 & de_class == 'transfer_only' ) |
          ( abs(log_fc_inflamed) > 2 & de_class == 'inflamed_only' )
      ),
    col = 'darkorange'
  ) +
  ggtitle(label = "Log Fold Change of indirect placenta and lung inflammation at 2 and 5 hpe")

ggsave(
  filename = "results/inflammation_comparison/lung_response_plot.pdf", 
  width = a4_short_long[2], height = a4_short_long[1], 
  units = 'mm'
)

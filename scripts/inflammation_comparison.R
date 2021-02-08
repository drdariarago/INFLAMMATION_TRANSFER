## Compare placental response to LPS-induced inflammation

library(tidyverse)
library(magrittr)

q_value_threshold = 0.05
log_fc_threshold = 0

divergence = function(x,y){
  difference = (x-y)^2
  mean = (x^2+y^2)
  sqrt(difference / mean)
}

## Import DE genes from inflammation study and from our study
de_inflammation_results <-
  here::here('results/process_lien_results/fold_change_summary.rds') %>% 
  read_rds()

de_transfer_results <-
  read_csv(file = "results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>% 
  filter(
    tissue == "placentas",
    timepoint == "timepoint5"
  )

# Perform rank-correlation between the two studies
# Filter only genes in both studies with FDR >= 0.05

shared_genes <-
  intersect(
  de_inflammation_results$mgi_symbol, 
  de_transfer_results$mgi_symbol
)

filtered_de_transfer <-
  de_transfer_results %>% 
  filter(mgi_symbol %in% shared_genes) %>% 
  select(
    ensembl_gene_id, mgi_symbol, q_value, log_fc = logFC, exposure
  ) %>% 
  pivot_wider(names_from = exposure, values_from = c(q_value, log_fc)) %>% 
  select(
    ensembl_gene_id, mgi_symbol, ave_expr = log_fc_baseline, q_value = q_value_response, log_fc = log_fc_response
  )

filtered_de_inflammation <-
  de_inflammation_results %>% 
  filter(mgi_symbol %in% shared_genes) %>% 
  select(
    mgi_symbol, q_value, ave_expr, log_fc
  ) 

shared_de <-
  inner_join(
    x = filtered_de_transfer,
    y = filtered_de_inflammation,
    by = "mgi_symbol",
    suffix = c("_transfer", "_inflamed")
  ) %>% 
  mutate(
    divergence = divergence(log_fc_inflamed, log_fc_transfer)
      # abs(log_fc_inflamed - log_fc_transfer) / (abs(log_fc_inflamed) + abs(log_fc_transfer))
  )

# Filter only genes with significant logFC in at least one condition
filtered_shared_de <-
  shared_de %>% 
  filter(
    (q_value_transfer < 0.05) | (q_value_inflamed < 0.05),
    (abs(log_fc_transfer) > 0.4 | abs(log_fc_inflamed) > 0.4)
  ) %>% 
  mutate(
    de_class = case_when(
      (q_value_transfer < 0.05) & (q_value_inflamed < 0.05) ~ "both",
      (q_value_transfer < 0.05) ~ "transfer_only",
      (q_value_inflamed < 0.05) ~ "inflamed_only",
      TRUE ~ 'not_significant'
    )
  )

# Calculate correlation
filtered_shared_de %$%
  cor.test(log_fc_inflamed, log_fc_transfer)

filtered_shared_de %$%
  cor.test(log_fc_inflamed, log_fc_transfer, method = 'spearman')

# Calculate and compare fold change and rank values
filtered_shared_de %>% 
  ggplot(
    aes(
      y = log_fc_transfer, x = log_fc_inflamed, 
      col = de_class
    )
  ) +
  geom_point(
    aes(alpha = 
          sqrt( (log_fc_inflamed / 5) ^2 + log_fc_transfer^2 ) 
    )) +
  geom_abline(slope = 1, intercept = 0, col = 'hotpink') +
  geom_smooth(method = 'lm')  +
  ggrepel::geom_text_repel(
    data = filtered_shared_de %>%
      filter(
      abs(log_fc_transfer) > 0.45,
      abs(log_fc_inflamed) > 0.45,
      divergence < 0.4
      ),
    aes(label = mgi_symbol)
  ) +
  ggrepel::geom_text_repel(
    data = filtered_shared_de %>%
      filter(
        abs(log_fc_inflamed) > 2.5
        ),
    aes(label = mgi_symbol)
  ) +
  ggrepel::geom_text_repel(
    data = filtered_shared_de %>%
      filter(
        abs(log_fc_transfer) > 0.7,
        abs(log_fc_inflamed) < 1,
        divergence > 0.4
      ),
    aes(label = mgi_symbol)
  ) +
  ggtitle(label = "Log Fold Change of indirect and direct placenta inflammation at 5 hpe")

filtered_shared_de %>% 
  # mutate(
  #   log_fc_inflamed = log_fc_inflamed %>% scale(center = F, scale = T) / 2,
  #   log_fc_transfer = log_fc_transfer %>% scale(center = F, scale = T)
  # ) %>% 
  ggplot(
    aes(
      x = log_fc_inflamed,
      y = log_fc_transfer,
      label = mgi_symbol
    )
  ) +
  geom_hex() +
  viridis::scale_fill_viridis(trans = 'log', breaks = c(1,10,100,500)) +
  facet_wrap( ~ de_class) +
  ggrepel::geom_text_repel(
    data = . %>% filter( 
      (log_fc_inflamed > 2.5) & de_class != 'transfer_only' | 
        (abs(log_fc_transfer) > 0.7) & de_class != 'inflamed_only' |
        ( (log_fc_transfer > 0.5) & (de_class == "both") )
      ),
    col = 'hotpink'
  ) +
  geom_abline(slope = 1, intercept = 0, color = 'blue')

table(filtered_shared_de$de_class)

filtered_shared_de %>% 
  filter(
    q_value_transfer < 0.05,
    q_value_inflamed < 0.05
  ) %>% 
  ggplot(
    aes(y = log_fc_transfer, x = log_fc_inflamed, col = de_class)
  ) +
  geom_point(alpha = .8) +
  geom_abline(slope = 1, intercept = 0, col = 'hotpink') +
  # ggrepel::geom_text_repel(
  #   data = filtered_shared_de %>%
  #     filter(
  #       de_class != "not_significant"
  #     ),
  #   aes(label = mgi_symbol, col = de_class) 
  # ) +
  ggtitle(label = "Log Fold Change of indirect and direct placenta inflammation at 5 hpe")

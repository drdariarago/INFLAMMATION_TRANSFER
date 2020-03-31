## Plot contrasts between maternal and fetal placentas

library(tidyverse)
library(ggrepel)

placentas_results <- 
  read_rds(path = here::here("results/limma_placentas/limma_placenta_results.Rdata"))

gene_names <- 
  read_rds(here::here("results/download_gene_data/gene_names.Rdata")) %>% 
  dplyr::select(.data = ., -3)

placentas_results_formatted <-
  placentas_results %>% 
  pivot_wider(names_from = maternal, values_from = c(expression, pval, qval)
  ) %>% 
  mutate(
    joint_qval = paste0(as.numeric(qval_maternal < 0.05), as.numeric(qval_fetal < 0.05)) %>% 
      factor(x = ., 
             levels = c("00", "01", "10", "11"), 
             labels = c("NS", "fetal_only", "convergent", "divergent")
             ),
    timepoint = factor(
      x = timepoint, 
      levels = c(2,5,12,24), 
      labels = paste( c(2,5,12,24), "hours post exposure")),
    gene_id = gsub('\\..*', '', gene_id)
  ) %>% 
  right_join(
    x = gene_names, 
    y = ., 
    by = c("ensembl_gene_id" = "gene_id")
  )

placenta_plot <- ggplot(
  data = placentas_results_formatted, 
  mapping = aes(
    x = expression_maternal,
    y = expression_fetal + expression_maternal,
    label = mgi_symbol)
) + 
  geom_bin2d(data = subset(placentas_results_formatted, joint_qval == "NS"),
             bins = 50) +
  scale_fill_gradient(name = "count", trans = "log", low = 'white', high  = 'pink', 
                      breaks = c(0,1,10,50,100,500,1000,3000)) +
  geom_vline(xintercept = 0, col = 'yellow') + 
  geom_abline(intercept = 0, slope = 1, col = 'yellow') +
  geom_point(data = subset(placentas_results_formatted, joint_qval != "NS"), 
             aes(col = joint_qval), alpha = 0.8) +
  scale_color_brewer(type = 'qual', palette = 2) +
  ylim(-4,4) + xlim(-2.5,3.7) +
  facet_wrap(~timepoint) +
  xlab(label = "Fold change between control and LPS in decidua") +
  ylab(label = "Fold change between control and LPS in placenta") +
  ggtitle("Expression responses in maternal placenta and decidua")

placenta_plot

# Save plot with labels for the 2 classes of discordant q-values
ggsave(plot = placenta_plot, 
       filename = here::here("results/MA_placentas/unlabelled_fold_changes.pdf"), 
       width = 13, height = 7.5)

placenta_plot + geom_label_repel(data = subset(placentas_results_formatted, joint_qval == "divergent")) 
ggsave(filename = here::here("results/MA_placentas/divergent_fold_changes.pdf"), 
       width = 13, height = 7.5)

placenta_plot + geom_label_repel(data = subset(placentas_results_formatted, joint_qval == "fetal_only")) 
ggsave(filename = here::here("results/MA_placentas/fetal_only_fold_changes.pdf"), 
       width = 13, height = 7.5)

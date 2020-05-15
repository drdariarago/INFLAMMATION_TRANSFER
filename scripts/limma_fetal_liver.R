## Compare expression between exposed and control fetal liver

library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)

experiment_data <-
  snakemake@input[['expression_results']] %>%
  readRDS() %>%
  {DGEList(
    counts  = assays(.) %>% extract2("counts"),
    samples = colData(.),
    genes   = rowData(.),
    group   =
      colData(.) %>%
      as_tibble %>%
      transmute(group = paste(tissue, timepoint, maternal_fetal, exposure, sep = ":")) %>%
      pull(group)
  )} %>% 
  .[,.$samples$exposure != "TiO2"] %>% 
  calcNormFactors(method = "TMM")

# Subset only the fetal liver
liver_data <- 
  experiment_data[, experiment_data$samples$tissue == "liver" & experiment_data$samples$maternal_fetal == "fetal"]

liver_data$samples <-
  liver_data$samples %>% 
  mutate(
    exposure = exposure %>% droplevels() %>% factor(x = ., levels = c('LPS', 'ctr')),
    timepoint = factor(x = liver_data$samples$timepoint, levels = c(2,5,12,24))
  )

# Design matrix
# Set intercept on controls, calculate coefs for each treatment separately at each stage
# Calculate response for maternal and compare with foetal

design <- 
  model.matrix(
    object = formula( ~ 0 + timepoint / exposure),
    contrasts.arg = list( 
      timepoint = "contr.sum",
      exposure = "contr.sum"
    ),
    data = liver_data$samples
  )

design_check <-
  design %>% 
  set_rownames(liver_data$samples[ ,1])

print(x = "Performing the following contrasts:")
colnames(design_check)

write_csv(x = as.data.frame(design_check), path = snakemake@output[['factor_design_matrix']])

# Filter low expression genes and apply voom transformation

filtered_data <-
  filterByExpr(
    y = liver_data, 
    design = design, 
    min.count = snakemake@params[['min_counts']],
    min.total.count = 1
  ) %>% 
  liver_data[., , keep.lib.sizes = T] %>% 
  voom(
    counts = .,
    design = design,
    lib.size = .$lib.size,
    plot = TRUE
  )

linear_model <-
  lmFit(
    object = filtered_data, 
    design = filtered_data$design
  ) %>% 
  eBayes()

write_rds(linear_model, snakemake@output[['linear_models']])

# Check distribution of q-values across different types of factors

pdf(file = snakemake@output[['fdr_plot']], width = 8.3, height = 11.7)

q_values <- 
  linear_model %$% 
  p.value %>% 
  as_tibble(x = ., rownames = 'gene_id') %>% 
  pivot_longer(data = ., cols = -contains('gene'), names_to = 'contrast', values_to = 'p_value') %>% 
  mutate(
    q_value = fdrtool::fdrtool(x = p_value, statistic = 'pvalue') %$% qval,
    contrast = factor(
      x = contrast, 
      levels = c(
        paste("timepoint", c(2,5,12,24), sep = ""),
        paste("timepoint", c(2,5,12,24), ":exposure1", sep = "")
      )),
    exposure = factor(x = grepl(pattern = "exposure", x = contrast), labels = c('other', 'exposure'))
  )

dev.off()

# Plot q-value distribution

q_distribution <-
  q_values %>% 
  ggplot(data = ., mapping = aes(x = q_value, fill = contrast)) +
  geom_histogram(position = 'dodge', binwidth = 0.08) +
  facet_wrap( ~ exposure) +
  scale_fill_manual(values = viridis::viridis(16)) +
  ggtitle(label = 'q-value distribution of different fetal liver contrasts')

q_distribution

ggsave(filename = snakemake@output[['q_values_plot']], 
       width = 297, height = 210, units = 'mm')

# Zoom to focus on interaction terms

q_distribution +
  coord_cartesian(ylim = c(0,6000)) 

ggsave(filename = snakemake@output[['q_values_plot_zoomed']], 
       width = 297, height = 210, units = 'mm')

# Merge with fold change values and save as table
fold_change <- 
  linear_model %$%
  coefficients %>% 
  as_tibble(x = ., rownames = 'gene_id') %>% 
  pivot_longer(data = ., cols = -contains('gene'), names_to = 'contrast', values_to = 'logFC') %>% 
  mutate(
    contrast = factor(
      x = contrast, 
      levels = c(
        paste("timepoint", c(2,5,12,24), sep = ""),
        paste("timepoint", c(2,5,12,24), ":exposure1", sep = "")
      )
    ),
    timepoint = factor(x = str_extract(string = contrast, pattern = 'timepoint[0-9]{1,2}'), levels = paste('timepoint', c(2,5,12,24), sep = "")),
    exposure = factor(x = grepl(pattern = 'exposure', x = contrast), levels = c(F,T), labels = c('baseline', 'response'))
  )

# Collate summary table with q-values, fold change and gene symbols

gene_annotation <- 
  read_rds(snakemake@input[['gene_data']])

result_summary_table <-
  full_join(
    x = q_values %>% dplyr::select(c('gene_id', 'contrast', 'q_value')),
    y = fold_change,
    by = c('gene_id', 'contrast')
  ) %>%
  mutate(gene_id = gsub(pattern = "\\..*", x = gene_id, replacement = "")) %>% 
  right_join(
    x = gene_annotation %>% dplyr::select(c('ensembl_gene_id', 'mgi_symbol')),
    y = .,
    by = c("ensembl_gene_id" = "gene_id")
  )


## Volcano plot for each contrast type
ggplot(result_summary_table, 
       aes(x = logFC, y = q_value, col = abs(logFC) > snakemake@params[['fold_change_threshold']])
) +
  facet_grid(timepoint ~ exposure) +
  geom_point(alpha = 0.3) + 
  scale_y_continuous(trans = 'log10', oob = scales::squish) +
  scale_x_continuous(breaks = seq(-5,5,2), minor_breaks = seq(-4,4,2)) +
  coord_cartesian(xlim = c(-5,5), ylim = c(1E-13,1)) +
  scale_color_brewer(type = 'qual', direction = -1) +
  geom_hline(aes(yintercept = snakemake@params[['alpha']])) +
  labs(col = 'Absolute LogFC over 0.5') + 
  ggtitle(label = 'q-value vs fold change across different contrasts')

ggsave(filename = snakemake@output[['volcano_plots']], width = 11.7, height = 8.3, units = "in", device = 'png', dpi = 'retina')


# Save gene-wise summaries 

write_csv(x = result_summary_table, path = snakemake@output[['summary_csv']])
write_rds(x = result_summary_table, path = snakemake@output[['summary_rds']])

# Create ranked gene lists for GO enrichment analysis

fold_change_grouped_tibble <-
  fold_change %>% 
  filter(grepl(pattern = 'exposure', x = contrast)) %>% 
  mutate(gene_id = gsub(pattern = "\\..*", x = gene_id, replacement = "")) %>% 
  group_by(contrast)

# Sort genes from highest to lowest FC to detect upregulation in pathways
fold_change_upregulated_genes <-
  fold_change_grouped_tibble %>% 
  arrange(-logFC, .by_group = TRUE) %>% 
  group_map(
    .f = ~ pull(.x, gene_id)
  ) %>% 
  set_names(
    group_keys(fold_change_grouped_tibble) %>% pull(contrast)
  )

write_rds(x = fold_change_upregulated_genes, path = snakemake@output[['ranked_genes_upregulated']])

# Sort genes from lowest to highest FC to detect downregulation in pathways
fold_change_downregulated_genes <-
  fold_change_grouped_tibble %>% 
  arrange(logFC, .by_group = TRUE) %>% 
  group_map(
    .f = ~ pull(.x, gene_id)
  ) %>% 
  set_names(
    group_keys(fold_change_grouped_tibble) %>% pull(contrast)
  )

write_rds(x = fold_change_downregulated_genes, path = snakemake@output[['ranked_genes_downregulated']])

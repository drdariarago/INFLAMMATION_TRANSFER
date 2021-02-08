### Process direct inflammation data, from
## Lien YC, Zhang Z, Barila G, Green-Brown A et al. 
## Intrauterine Inflammation Alters the Transcriptome and Metabolome in Placenta. 
## Front Physiol 2020;11:592689. PMID: 33250783

library(tidyverse)
library(magrittr)
library(limma)
library(edgeR)

count_data_names <-
  snakemake@input[[1]] %>%
  gzfile() %>%
  read.table(nrows = 1, stringsAsFactors = F) %>% 
  as.character() %>% 
  c("X1", .)

count_data <-
  snakemake@input[[1]] %>%
  gzfile() %>%
  read_tsv(col_names = count_data_names, skip = 1)

metadata <-
  count_data_names %>% 
  str_subset(string = ., pattern = "X[0-9]{4}[A-Z]+\\.[L,R][1-3][M,F]") %>% 
  as_tibble() %>% 
  set_colnames("sample_id") %>% 
  mutate(
    family = str_extract(string = sample_id, pattern = "X[0-9]{4}"),
    treatment = str_extract(string = sample_id, pattern = "[0-9]{4}[A-Z]{1,3}") %>% 
      str_extract(., pattern = "[A-Z]{1,3}"),
    sex = str_extract(string = sample_id, pattern =  "[M,F]$"),
    side = str_extract(string = sample_id, pattern = "[L,R][0-9]") %>% 
      str_extract(., pattern = "[L,R]"),
    replicate = str_extract(string = sample_id, pattern = "[R,L][1-3][M,F]") %>% 
      str_extract(., pattern = "[1-3]")
  )

## Convert to DGE list for edgeR/Limma processing

experiment_data <- 
  DGEList(
    counts = select(count_data, matches("X[0-9]{4}[A-Z]+\\.[L,R][1-3][M,F]")),
    samples = metadata,
    genes = select(count_data, c("Symbol", "Type_of_gene", "Description"))
  )

design <- 
  model.matrix(
    object = formula( ~ treatment + sex),
    data = experiment_data$samples
  )

## Filter low expression gene and convert to pseudocounts

filtered_data <-
  filterByExpr(
    y = experiment_data, 
    design = design, 
    min.count = snakemake@params[['min_counts']]
  ) %>% 
  experiment_data[., , keep.lib.sizes = T] %>% 
  voom(
    counts = .,
    design = design,
    lib.size = .$lib.size,
    plot = TRUE
  )

## Find outlier sample using PCA
pca_results <-
  filtered_data$E %>% 
  t %>% 
  prcomp()

pca_results$x %>% 
  as_tibble(rownames = 'sample_id') %>% 
  ggplot(data = ., mapping = aes(x = PC1, y = PC2, label = sample_id)) +
  geom_point() +
  ggrepel::geom_text_repel()

pdf(file = snakemake@output[["hclust"]])
filtered_data$E %>% 
  t() %>% 
  dist() %>% 
  hclust(method = 'average') %>% 
  plot()
dev.off()

## Removing sample X2864C.L2M

experiment_data_2 <-
  row.names(experiment_data$samples) %>% 
  setdiff(., "X2864C.L2M") %>% 
  experiment_data[,.] 

design_2 <-
  model.matrix(
    object = formula( ~ treatment + sex),
    data = experiment_data_2$samples
  )

filtered_experiment_data_2 <-
  experiment_data_2 %>%  
  filterByExpr(
    y = ., 
    design = design_2, 
    min.count = snakemake@params[['min_counts']]
  ) %>% 
  experiment_data_2[., , keep.lib.sizes = T] 

# Normalize
pdf(snakemake@output[["voom"]])
vst_data <-
  voom(
    counts = filtered_experiment_data_2,
    design = design_2,
    lib.size = filtered_experiment_data_2$lib.size,
    plot = TRUE
  )
dev.off()

linear_model <-
  vst_data %>% 
  lmFit(
    object = ., 
    design = .$design
  ) %>% 
  eBayes()

linear_model_results <-
  topTable(linear_model, number = Inf, coef = "treatmentLPS") %>% 
  as_tibble() %>% 
  mutate(
    q_value = fdrtool::fdrtool(x = P.Value, statistic = 'pvalue') %$% qval
  ) %>% 
  select(
    mgi_symbol = Symbol, 
    gene_type = Type_of_gene, 
    q_value, 
    ave_expr = AveExpr,
    log_fc = logFC,
    description = Description
  )

write_rds(x = linear_model_results, file = snakemake@output[['rds_results']])

linear_model_results %>% 
  ggplot(
    aes(
      x = q_value, 
      y = log_fc
    )
  ) +
  geom_hex() + 
  scale_x_log10() + 
  viridis::scale_fill_viridis(trans = 'log', breaks = c(1,10,100,1000)) + 
  ggtitle(label = "Volcano plot of Placenta LPS responses")

ggsave(filename = snakemake@output[['volcano_plot']])

linear_model_results %>% 
  ggplot(
    aes(
      x = ave_expr, 
      y = log_fc
    )
  ) +
  geom_hex() + 
  viridis::scale_fill_viridis(trans = 'log', breaks = c(5,50,500)) + 
  ggtitle(label = "MA plot of Placenta LPS responses")

ggsave(filename = snakemake@output[['MA_plot']])

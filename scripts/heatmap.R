library(here)
library(magrittr)
library(SummarizedExperiment)
library(DESeq2)

# Loading vsd object and integrating summarized experiment, to obtain all variables
variance_stabilized_counts <-
  readRDS(here('results/tximeta/vsd.Rdata')) %>%
  DESeqDataSet(se = ., design = ~ tissue + maternal_fetal + exposure + timepoint) %>%
  vst(object = ., blind = T)

# Subset matrix with only
lung_annotated_counts <- 
  variance_stabilized_counts[, variance_stabilized_counts$tissue == "lung"] 

lung_count_matrix <- 
  assay(lung_annotated_counts)

lung_annotation <- 
  colData(lung_annotated_counts) %>% 
  magrittr::extract(, c(5,7)) %>% 
  as.data.frame()

# Create matrix with top 1000 variance genes
top_genes <- 
  lung_count_matrix %>%
  rowVars() %>% 
  order(decreasing = T) %>% 
  magrittr::extract(seq_len(5000))

# Making heatmap of top 1000 variance genes
library(RColorBrewer)
colors <- 
  colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)


pheatmap::pheatmap(
  lung_count_matrix[top_genes,] %>% t(),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  scale = "row",
  cutree_rows = 10,
  cluster_rows = T,
  cluster_cols = T,
  display_numbers = F,
  cex= 0.9,
  show_colnames = F,
  annotation_row = lung_annotation,
  col = colors
)

library(here)
library(magrittr)
library(SummarizedExperiment)
library(DESeq2)

# Loading vsd object and integrating summarized experiment, to obtain all variables
variance_stabilized_counts <-
  readRDS(here('results/tximeta/vsd.Rdata')) %>%
  DESeqDataSet(se = ., design = ~ tissue + maternal_fetal + exposure + timepoint) %>%
  vst(object = ., blind = T)

# Remove TiO2 samples from analysis
annotated_counts = 
 variance_stabilized_counts[,!grepl("TiO2",variance_stabilized_counts$exposure)]

#annotated_counts <-
 # variance_stabilized_counts$exposure %in% c("LPS", "ctr")
#annotated_counts <-
# variance_stabilized_counts$exposure != "TiO2"


## MATERNAL LUNGS
  # Subset matrix with only maternal lung samples
lung_annotated_counts <- 
  annotated_counts[, annotated_counts$tissue == "lung"] 

# Extract only the genes by samples count matrix 
lung_count_matrix <- 
  assay(lung_annotated_counts)

# Extract only the varaibles we want to annotate in the heatmap
lung_annotation <- 
  colData(lung_annotated_counts) %>% 
  magrittr::extract(, c(5,7)) %>% 
  as.data.frame()

# Create matrix with top 5000 variance genes
top_genes <- 
  lung_count_matrix %>%
  rowVars() %>% 
  order(decreasing = T) %>% 
  magrittr::extract(seq_len(5000))

# Make heatmap of top 5000 variance genes
library(RColorBrewer)
colors <- 
  colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

pheatmap::pheatmap(
  lung_count_matrix[top_genes,] %>% t(),
  annotation_row = lung_annotation,
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
  col = colors
)





# MATERNAL LIVER 
# Subset matrix with only maternal liver samples
liver_annotated_counts <- 
  annotated_counts[, annotated_counts$tissue == "liver"]

maternal_liver_annotated_counts <-
liver_annotated_counts[, liver_annotated_counts$maternal_fetal == "maternal"]

maternal_liver_count_matrix <- 
  assay(maternal_liver_annotated_counts)

maternal_liver_annotation <- 
  colData(liver_annotated_counts) %>% 
  magrittr::extract(, c(5,7)) %>% 
  as.data.frame()

# Create matrix with top 5000 variance genes
top_genes <- 
  maternal_liver_count_matrix %>%
  rowVars() %>% 
  order(decreasing = T) %>% 
  magrittr::extract(seq_len(5000))

# Make heatmap of top 5000 variance genes
library(RColorBrewer)
colors <- 
  colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

pheatmap::pheatmap(
  maternal_liver_count_matrix[top_genes,] %>% t(),
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
  annotation_row = maternal_liver_annotation,
  col = colors
)

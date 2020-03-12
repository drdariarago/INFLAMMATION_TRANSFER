library(here)
library(magrittr)
library(SummarizedExperiment)
library(DESeq2)

# Loading vsd object and integrating summarized experiment, to obtain all variables
variance_stabilized_counts <-
  readRDS(here('results/tximeta/vsd.Rdata')) %>%
  DESeqDataSet(se = ., design = ~ tissue + maternal_fetal + exposure + timepoint) %>%
  vst(object = ., blind = T) %>%
  assay() 

# Create matrix with top 1000 variance genes
top_1000_genes <- 
  rowVars(variance_stabilized_counts) %>% 
  order(decreasing = T) %>% 
  magrittr::extract(seq_len(1000))

#Making heatmap of top 1000 variance genes
library(RColorBrewer)
colors <- 
  colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

top_1000_heatmap <-
  pheatmap::pheatmap(
    variance_stabilized_counts[top_1000_genes,] %>% 
    t(),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "average",
    cluster_rows = T,
    cluster_cols = T,
    display_numbers = F,
    cex= 0.5,
    col = colors
  )

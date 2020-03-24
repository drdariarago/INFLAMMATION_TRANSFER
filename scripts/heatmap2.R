library(here)
library(pheatmap)
library(magrittr)
library(rlist)
library(DESeq2)
library(tidyverse)

limma_coefs <- 
  readRDS(here("results/limma/limma_coefs.Rdata"))

placenta_fetal_tbl <- 
  limma_coefs[["placenta_fetal"]] %>%
  filter(treatment == "exposureLPS") %>%
  as_tibble()

timepoint_annotation <-
  placenta_fetal_tbl["timepoint"] %>%
  as.data.frame() 

placenta_fetal_wider <-
placenta_fetal_tbl %>%
  group_by(gene_id) %>%
    pivot_wider(
      names_from = timepoint,
      values_from = LogFC) %>%
column_to_rownames("gene_id") %>%
  select(-c(treatment, AveExpr))
  
placenta_fetal_matrix <-
  as.matrix(placenta_fetal_wider)
 
library(RColorBrewer)
colors <- 
  colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

pheatmap::pheatmap(
  placenta_fetal_matrix %>% t(),
  #annotation_row = timepoint_annotation,
  #clustering_distance_rows = "euclidean",
  #clustering_distance_cols = "euclidean",
  #clustering_method = "average",
  scale = "column",
  cluster_rows = F,
  cluster_cols = T,
  display_numbers = F,
  cex= 0.9,
  show_colnames = F
)






library(here)
library(pheatmap)
library(magrittr)
library(rlist)
library(DESeq2)
library(tidyverse)
library(stringr)

# Read limma coefficients and choose LPS exposure
limma_coefs <- 
  readRDS(here("results/limma/limma_coefs.Rdata"))

placenta_fetal_tbl <- 
  limma_coefs[["placenta_fetal"]] %>%
  filter(treatment == "exposureLPS") %>%
  as.data.frame()

#timepoint_annotation <-
 # placenta_fetal_tbl["timepoint"] %>%
  #as.data.frame()

# Transform into wider format and delete gene_id suffix
placenta_fetal_wider <-
  placenta_fetal_tbl %>%
  group_by(gene_id) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = LogFC) %>%
  select(-c(treatment, AveExpr))

placenta_fetal_wider$gene_id <-
  gsub('\\..*', '', 
      placenta_fetal_wider$gene_id)

# Read in ENSEMBL annotated results and extract gene_id
annotated_results <- 
        readRDS(here("results/limma_results/annotated_results.Rdata"))

placenta_fetal_ensembl_gene_id <-
  annotated_results[["placenta_fetal"]] %>%
  magrittr::extract(c("ensembl_gene_id"))

# Right join the two data frames and make gene_id rownames
placenta_fetal_join <-
  right_join(x= placenta_fetal_wider, 
            y = placenta_fetal_ensembl_gene_id, 
            by = c("gene_id" = "ensembl_gene_id"))
  
placenta_fetal_joined <-
  placenta_fetal_join %>%
  column_to_rownames("gene_id") 

# Make matrix
placenta_fetal_matrix <-
  as.matrix(placenta_fetal_joined)

# Heatmap!
colors <- 
  colorRampPalette(
    RColorBrewer::brewer.pal(9, "PiYG")
  )(255)

pheatmap::pheatmap(
  placenta_fetal_matrix %>% t(),
  color = colors,
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  scale = "column",
  cluster_rows = F,
  cluster_cols = T,
  display_numbers = F,
  cex= 0.9,
  show_colnames = F,
  main = "placenta fetal p = 0.05, 1008 genes"
)



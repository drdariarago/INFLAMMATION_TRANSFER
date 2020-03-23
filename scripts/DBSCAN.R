## Summarize limma results via DBSCAN projected on 2D MDS
library(dbscan)
library(tidyverse)

# Import data
fold_change_list <-
  read_rds(path = snakemake@input[["coefs"]])

# Convert to tissue-specific list matrices, scaled by 

fold_change_matrices <-
  fold_change_list %>% 
  map(
    filter,
    treatment == 'exposureLPS'
    ) %>% 
  map(
    group_by,
    gene_id
  ) %>% 
  map(
    mutate,
    LogFC = scale(x = LogFC, scale = TRUE, center = TRUE)
  ) %>% 
  map(
    pivot_wider,
    id_cols = gene_id,
    names_from = timepoint,
    values_from = LogFC
  ) %>% 
  map(
    column_to_rownames,
    var = 'gene_id'
  ) %>% 
  map(
    as.matrix
  )

# Calculate euclidean distance between genes
fold_distances <- 
  fold_change_matrices %>% 
  map(
    dist,
    method = 'euclidean'
  )

# Summarize via MDS
projected_genes <- 
  fold_distances %>% 
  map(
    cmdscale,
    k = 2
  )

# DBSCAN
scan_results <-
  fold_distances %>% 
  map(
    hdbscan, 
    minPts = 10
  )

write_rds(x = scan_results, path = snakemake@output[["clusters"]])

# Create hull plots
pdf(file = snakemake@output[["plots"]])

pmap(
  .l = list(
    projection = projected_genes,
    clusters = scan_results,
    plot_names = names(projected_genes)
  ),
  .f = function(projection, clusters, plot_names){
    hullplot(
      x = projection,
      cl = clusters,
      main = plot_names
    )
  }
)

dev.off()
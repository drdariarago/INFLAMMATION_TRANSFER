library(tidyverse)
library(here)
library(purrr)
library(pheatmap)

limma_coefs <- 
  readRDS(here("results/limma/limma_coefs.Rdata"))

# transform results 
non_annotated_results <-
limma_coefs %>%
  map(.x = ., .f = ~ filter(.x, treatment == "exposureLPS")) %>%
  map(.x = ., .f = ~ select(.x, -c(treatment))) %>% 
  map(.x = ., .f = ~ mutate(.data = .x, gene_id = gsub('\\..*', '', gene_id)))

# Read annotated results
annotated_results <- 
  readRDS(here::here("results/limma_results/annotated_results.Rdata")) %>% 
  map(.x = ., 
      .f = ~ pivot_longer( 
        data = .x, 
        cols = c("timepoint2", "timepoint5", "timepoint12", "timepoint24"), 
        names_to = "timepoint", values_to = "significant"
      )
  ) %>% 
  map(.x = ., 
      .f = ~ dplyr::select(.data = .x, -3)
  )

 #Join the two lists and transform into matrix
all_results <-
  map2(.x= annotated_results, .y = non_annotated_results,
       .f = ~ full_join(
         .x, .y,
         by = c("ensembl_gene_id" = "gene_id", "timepoint" = "timepoint")
       )
) %>%
  map(.x = ., .f = ~ select(.x, -c(mgi_symbol, significant, AveExpr))) %>%
  map(.x = ., .f = ~ group_by(.x, timepoint)) %>%
  map(.x = ., .f = ~ pivot_wider(
    data = .x,
    names_from = timepoint,
    values_from = LogFC)) %>%
  map(.x = ., .f = ~ column_to_rownames(.x, var = 'ensembl_gene_id')) %>%
  map(.x = ., .f = as.matrix) #%>%
  #map(.x = ., .f = t)
 
#heatmap!
colors <- 
  colorRampPalette(
    RColorBrewer::brewer.pal(9, "PiYG")
  )(255)

# heatmap_all <-
map(
  .x = all_results,
  .f = ~ pheatmap(
  mat = .x,
  color = colors,
  clustering_distance_rows = "correlation",
  clustering_method = "average",
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  display_numbers = F,
  cex= 0.9,
  show_colnames = T
  ))

# Saving heatmaps - not functioning before heatmaps are completely fine
imap(
  .x = heatmap_all,
  .f = ~ ggsave(
    plot = .x,
    path = here::here("results/heatmaps_all_tissues.png"),
    filename = paste0(.y, "_heatmaps_all_tissues.png"),
    width = 6.67, height = 7.5, units = 'in',
    device = "png"
  )
)

#Old 'standard' heatmap
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




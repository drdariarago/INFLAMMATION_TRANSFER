## Plot frequency of intersecting responses for timepoint by tissue

library(UpSetR)
library(tidyverse)
library(magrittr)
library(janitor)

ALPHA = snakemake@params[["ALPHA"]]
MIN_LOGFC = snakemake@params[["MIN_LOGFC"]]


filtered_results <- 
  snakemake@input[[1]] %>%
  # here::here("results/limma_maternal_lung/fold_change_summary.rds") %>% 
  readRDS %>% 
  filter(
    q_value < ALPHA, 
    abs(logFC) > MIN_LOGFC,
    exposure == "response"
  ) %>% 
  # Remove contrasts between chorion and decidua in placenta limma models
  { if ("maternal" %in% colnames(.) ) {
    filter(., maternal == "shared")
  } else . 
  } 

counts_data_frame <-
  filtered_results %>% 
  select(
    ensembl_gene_id, timepoint
  ) %>%
  janitor::tabyl( ensembl_gene_id, timepoint ) %>% 
  as.data.frame

gene_annotation <-
  filtered_results %>% 
    arrange(ensembl_gene_id) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(
      mgi_symbol = first(mgi_symbol),
      direction = case_when(
        all(logFC > 0) ~ "upregulated",
        all(logFC < 0) ~ "downregulated",
        TRUE ~ "mixed"
      ) %>% as.factor
    )
  
annotated_results <-
  left_join(
    gene_annotation, 
    counts_data_frame
  ) %>% 
  filter(
    direction != "mixed"
  ) %>% 
  as.data.frame()

#### Plot and save the upsetr plots ####

A4 = c(8.27, 11.69)
pdf(file = snakemake@output[[1]], width = A4[2]/2, height = A4[1]/2, onefile = FALSE)

upset(
  annotated_results,
  nintersects = NA,
  sets = c("timepoint24", "timepoint12", "timepoint5", "timepoint2"),
  keep.order = T,
  order.by = "freq",
  empty.intersections = T,
  queries = list(
    list(
      query = elements, params = list("direction", "downregulated"), active = T, color = 'lightgrey'
    )
  )
)

dev.off()

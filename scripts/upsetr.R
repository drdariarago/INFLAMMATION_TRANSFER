## Plot frequency of intersecting responses for timepoint by tissue

library(UpSetR)
library(tidyverse)
library(magrittr)
library(janitor)

ALPHA = snakemake@params[["ALPHA"]]
MIN_LOGFC = snakemake@params[["MIN_LOGFC"]]

annotated_results <- 
  snakemake@input[[1]] %>%
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
  } %>% 
  select(
    ensembl_gene_id, timepoint
  ) %>%
  janitor::tabyl( ensembl_gene_id, timepoint ) %>% 
  as.data.frame

#### Plot and save the upsetr plots ####

A4 = c(8.27, 11.69)
pdf(file = snakemake@output[[1]], width = A4[2]/2, height = A4[1]/2, onefile = FALSE)

upset(
  annotated_results,
  nintersects = NA,
  sets = c("timepoint24", "timepoint12", "timepoint5", "timepoint2"),
  keep.order = T,
  order.by = "freq",
  empty.intersections = T
)

dev.off()

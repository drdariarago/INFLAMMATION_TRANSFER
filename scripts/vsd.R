library(DESeq2)  
library(magrittr)

expression_results <- readRDS(snakemake@input[["expression_results"]])
expression_results$timepoint <- as.factor(expression_results$timepoint)

DESeq_results <- 
  DESeqDataSet(
    expression_results, 
    design = ~ tissue + maternal_fetal + exposure + timepoint
  ) 

filtered_genes <-
  DESeq_results %>% 
  counts() %>% 
  rowSums() %>% 
  is_greater_than(snakemake@params[['min_counts']])

DESeq_filtered <- DESeq_results[filtered_genes,]

variance_stabilized_counts <- vst(DESeq_filtered, blind = FALSE)

readr::write_rds(variance_stabilized_counts, snakemake@output[["variance_stabilized_counts"]])
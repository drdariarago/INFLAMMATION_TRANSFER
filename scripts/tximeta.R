### Import data from salmon raw output into R
# Note: this script uses tximport as installed from the github repository:
# the bioconductor version for R3.5 has bugs

# Load packagess
library(magrittr)
library(tximeta)
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)

# Import metadata and add filenames and paths

quant_file_data <- tibble(
  files = snakemake@input[["salmon_dirs"]] %>% paste0(., "/quant.sf"),
  names = str_extract(
    string = files, pattern = "[:upper:][:alpha:]{1,2}i?[:digit:]{1,2}X?") %>% 
    str_replace(pattern = "lL?i", replacement = "Li")
)

column_data <- 
  readRDS(file = snakemake@input[["metadata"]]) %>% 
  left_join(quant_file_data, ., by = c("names" = "sample_id"))

# Retrieve feature metadata (transcript names)
# NOTE: We lose all features that lack metadata at this step
setTximetaBFC(dir = "data/tximetaBFC")

transcript_data <- tximeta(
  coldata = column_data,
  type = "salmon",
  txOut = TRUE
)

# Convert to gene level
gene_data <- summarizeToGene(
  object = transcript_data,
  varReduce = TRUE,
  ignoreTxVersion = FALSE,
  ignoreAfterBar = FALSE,
  countsFromAbundance = 'lengthScaledTPM'
)

## Save output

# As summarizedExperiment
saveRDS(object = gene_data, file = snakemake@output[[1]])

# As matrix of variance stabilized counts
vsd = DESeqDataSet(se = gene_data, design = ~ tissue + maternal_fetal + exposure) %>%
  vst(object = ., blind = T) %>%
  assay() %>%
  write.csv(file = snakemake@output[[2]])
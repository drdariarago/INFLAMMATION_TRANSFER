# Check initial DE expression of genes over time

library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)

# Import data
experiment_data <-
  here::here("results/tximeta/expression_results.Rdata") %>%
  readRDS() %>%
  {DGEList(
    counts  = assays(.) %>% extract2("counts"),
    samples = colData(.),
    genes   = rowData(.),
    group   =
      colData(.) %>%
      as_tibble %>%
      transmute(group = paste(tissue, timepoint, maternal_fetal, exposure, sep = ":")) %>%
      pull(group)
  )}

# Remove TiO2 samples
experiment_data <- 
  experiment_data[,experiment_data$samples$exposure != "TiO2"] 

experiment_data$samples <- droplevels(experiment_data$samples)

# Perform library normalization
experiment_data <- 
  experiment_data %>% 
  calcNormFactors(method = "TMM")

# Create list of datasets (one per tissue)
experiment_data$samples$tissue2 <-
  paste(experiment_data$samples$tissue, 
        experiment_data$samples$maternal_fetal, 
        sep = "_") %>% 
  as.factor()

tissue_types = unique(experiment_data$samples$tissue2)

experiment_data$samples$timepoint <-
  as.factor(experiment_data$samples$timepoint)

experiment_list <-
  lapply(
    tissue_types, 
    function(tissue){
      experiment_data[,experiment_data$samples$tissue2 == tissue]
    }
  ) %>% 
  set_names(tissue_types)

# Specify design matrix
design_list <-
  lapply(
    experiment_list,
    function(experiment){
      model.matrix(
        object = formula( ~ timepoint + timepoint:exposure) ,
        contrasts.arg = list(
          timepoint = "contr.treatment"),
        data = experiment$samples
      )}
  )

# Filter low expression genes and apply voom transformation

filtered_data <-
  lapply(
    tissue_types, 
    function(tissue){
      filterByExpr(
        y = experiment_list[[tissue]],
        design = design_list[[tissue]],
        min.count = 5,
        min.total.count = 1
      ) %>% 
        experiment_list[[tissue]][., , keep.lib.sizes = T] %>% 
        voom(
          design = design_list[[tissue]],
          plot = T,
          lib.size = .$lib.size
        )
    }
  ) %>% 
  set_names(tissue_types)

# Check initial DE expression of genes over time
library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)

# Import data
experiment_data <-
  snakemake@input[[1]] %>%
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
  )} %>% 
  .[,.$samples$exposure != "TiO2"] %>% 
  calcNormFactors(method = "TMM")

# Create convenience variables
# exposure removes TiO2 factor level
# tissue2 is necessary to split the dataset
# timepoint needs to be converted to factor to be used in the contrast design

experiment_data$samples <-
  experiment_data$samples %>% 
  mutate(
    exposure = droplevels(exposure),
    tissue2 = paste(tissue, maternal_fetal, sep = "_") %>% as.factor(),
    timepoint = as.factor(experiment_data$samples$timepoint)
  )

# Create list of datasets (one per tissue)
tissue_types = unique(experiment_data$samples$tissue2)

experiment_list <-
  map(
    .x = tissue_types,
    .f = ~ experiment_data[, experiment_data$samples$tissue2 == .x]
  ) %>% 
  set_names(tissue_types)

# Specify design matrix
design_list <-
  map(
    .x = experiment_list, 
    .f = 
      ~ model.matrix(
        object = formula( ~ timepoint + timepoint:exposure),
        contrasts.arg = list( timepoint = "contr.treatment"),
        data = .x$samples
    )
  )

# Filter low expression genes and apply voom transformation
pdf(file = snakemake@output[["plots"]])

filtered_data <-
  map2(
    .x = experiment_list,
    .y = design_list,
    .f = filterByExpr,
    min.count = 5,
    min.total.count = 1
  ) %>% 
  map2(
    .x = .,
    .y = experiment_list,
    .f = ~ .y[.x, , keep.lib.sizes = T]
  ) %>% 
  map2(
    .x = .,
    .y = design_list,
    .f = ~ voom(
      counts = .x,
      design = .y,
      lib.size = .x$lib.size,
      plot = TRUE
    )
  )

dev.off()

# Fit linear models
linear_model_lists <-
  map2(
    .x = filtered_data,
    .y = design_list,
    .f = lmFit
  ) %>%
  map(., eBayes)

write_rds(linear_model_lists, snakemake@output[["fitted_models"]])

# Summarize results
limma_coefs <- map(
  .x = linear_model_lists, 
  .f = topTable, 
  number = Inf
) %>% 
  map(
    .x = .,
    .f = rename_at,
    .vars = vars(matches("timepoint[0-9]{1,2}$")),
    .funs = list(~ paste0(., ".exposureCTRL"))
  ) %>% 
  map(
    .x = ., 
    .f = select, 
    ends_with(c("CTRL", "LPS", "Expr"))
  ) %>% 
  map(
    .x = .,
    .f = pivot_longer,
    cols = ends_with(c("CTRL", "LPS")),
    names_to = c("timepoint", "treatment"),
    names_pattern = "([[:alnum:]]*)\\.([[:alnum:]]*)",
    values_to = "LogFC"
  ) 

write_rds(limma_coefs, snakemake@output[["limma_coefs"]])
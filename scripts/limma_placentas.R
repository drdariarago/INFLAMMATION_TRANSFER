## Compare expression between maternal and foetal placentas

library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)

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
  )} %>% 
  .[,.$samples$exposure != "TiO2"] %>% 
  calcNormFactors(method = "TMM")


# Subset only the placentas
placentas_data <- 
  experiment_data[, experiment_data$samples$tissue == "placenta"]

placentas_data$samples <-
  placentas_data$samples %>% 
  mutate(
    exposure = droplevels(exposure),
    tissue = tissue %>% as.factor() %>% droplevels(),
    maternal = maternal_fetal %>% as.factor(),
    timepoint = factor(x = placentas_data$samples$timepoint, levels = c(2,5,12,24))
  )

# Design matrix
# Set intercept on maternal placenta, calculate coefs for maternal and foetal separately at each stage
# Calculate response for maternal and compare with foetal

design <- 
  model.matrix(
    object = formula( ~ maternal * timepoint + timepoint:exposure + maternal:timepoint:exposure ),
    contrasts.arg = list( timepoint = "contr.treatment"),
    data = placentas_data$samples
  )


# Filter low expression genes and apply voom transformation

filtered_data <-
  filterByExpr(
    y = placentas_data, 
    design = design, 
    min.count = 15,
    min.total.count = 1
  ) %>% 
  placentas_data[., , keep.lib.sizes = T] %>% 
  voom(
    counts = .,
    design = design,
    lib.size = .$lib.size,
    plot = TRUE
  )

linear_model <-
  lmFit(
    object = filtered_data, 
    design = filtered_data$design
    ) %>% 
  eBayes()

write_rds(linear_model, here::here("results/limma_placentas/linear_model.Rdata"))

# Summarize results
limma_coefs <-
  linear_model %>% 
  topTable(number = Inf) %>% 
  dplyr::select(contains("LPS")) %>% 
  rownames_to_column("gene_id") %>% 
  rename_at(
    .vars = vars(matches("^timepoint")), 
    .funs =  ~ paste0("maternal.", .)
  ) %>% 
  pivot_longer(
    cols = contains("LPS"),
    names_to = c("maternal", "timepoint"),
    names_pattern = "([a-z]*)[.]timepoint([0-9]{1,2})",
    values_to = "expression"
  ) %>% 
  mutate(
    timepoint = factor(x = timepoint, levels = c(2,5,12,24)),
    maternal = factor(x = maternal, levels = c("maternal", "maternalfetal"), labels = c("maternal", "fetal"))
  )

limma_qvals <-
  linear_model$p.value %>% 
  as_tibble(rownames = "gene_id") %>% 
  dplyr::select(matches(".*(LPS|gene).*")) %>% 
  rename_at(
    .vars = vars(matches("^timepoint")), 
    .funs =  ~ paste0("maternal.", .)
  ) %>% 
  pivot_longer(
    cols = contains("LPS"),
    names_to = c("maternal", "timepoint"),
    names_pattern = "([a-z]*)[:.]timepoint([0-9]{1,2})",
    values_to = "pval"
  ) %>% 
  mutate(
    timepoint = factor(x = timepoint, levels = c(2,5,12,24)),
    maternal = factor(x = maternal, levels = c("maternal", "maternalfetal"), labels = c("maternal", "fetal")),
    qval = fdrtool::fdrtool(x = pval, statistic = "pvalue") %$% lfdr
  )

limma_results <- full_join(limma_coefs, limma_qvals)
write_rds(x = limma_results, path = here::here("results/limma_placentas/limma_placenta_results.Rdata"))

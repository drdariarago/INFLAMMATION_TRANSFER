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
    object = formula( ~ 0 + timepoint * exposure * maternal),
    contrasts.arg = list( 
      timepoint = contr.treatment(n = 4, base = 2),
      maternal = "contr.sum", 
      exposure = "contr.sum"
    ),
    data = placentas_data$samples
  )

design_check <-
  design %>% 
  set_rownames(placentas_data$samples[ ,1])

colnames(design_check)

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

# Check distribution of q-values across different types of factors

qvalues <- 
  linear_model %$% 
  p.value %>% 
  as_tibble(x = ., rownames = 'gene_id') %>% 
  pivot_longer(data = ., cols = -contains('gene'), names_to = 'contrast', values_to = 'p_value') %>% 
  mutate(
    q_value = fdrtool::fdrtool(x = p_value, statistic = 'pvalue') %$% qval,
    contrast = factor(
      x = contrast, 
      levels = c(
        "maternal1", 
        "exposure1",
        "exposure1:maternal1",
        paste("timepoint", c(2,5,12,24), sep = ""),
        paste("timepoint", 1:4, ":exposure1", sep = ""),
        paste("timepoint", 1:4, ":maternal1", sep = ""),
        paste("timepoint", 1:4, ":exposure1", ":maternal1", sep = "")
      )),
    maternal = factor(x = grepl(pattern = "maternal", x = contrast), labels = c('other', 'maternal')),
    exposure = factor(x = grepl(pattern = "exposure", x = contrast), labels = c('other', 'exposure'))
  )

# Plot q-value distribution

q_distribution <-
  qvalues %>% 
  ggplot(data = ., mapping = aes(x = q_value, fill = contrast)) +
  geom_histogram(position = 'dodge', binwidth = 0.08) +
  facet_grid(exposure ~ maternal, switch = 'y') +
  scale_fill_manual(values = viridis::viridis(16)) 

q_distribution

ggsave(filename = here::here('results/limma_placentas/q_value_distribution.pdf'))

# Zoom to focus on interaction terms

q_distribution +
  coord_cartesian(ylim = c(0,6000)) 

ggsave(filename = here::here('results/limma_placentas/q_value_distribution_zoomed.pdf'))

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

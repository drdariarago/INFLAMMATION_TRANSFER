library(tidyverse)
library(magrittr)

#Using expression results because VSD gives less confidences for genes with low counts

experiment_data <-
  readRDS(here::here("results/tximeta/expression_results.Rdata"))

#https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#subsetting

# Subset just counts

#make list for each tissue
tissue_list <-
  list(
    maternal_lung = experiment_data[, experiment_data$tissue == 'lung'],
    all_placenta = experiment_data[, experiment_data$tissue == 'placenta'],
    fetal_liver = experiment_data[, experiment_data$tissue == 'liver' & 
                                    experiment_data$maternal_fetal == 'fetal'],
    maternal_liver = experiment_data[, experiment_data$tissue == 'liver' & 
                                       experiment_data$maternal_fetal == 'maternal']
  )

tissue_matrix <-
  tissue_list %>% 
  map(.x = ., .f = ~ assay(.x)) %>%
  map(.x = ., .f = ~ as.matrix(.x)) %>% 
  map(.f = t)

# make PCAs with map for each tissue (placentas together)
# PCs as columns, samples as rows

# doing the pca..
pca <-
  tissue_matrix %>%
  map( ~ prcomp(.))
  
# extracting x       
pca_data <-
  pca %>%
       map(~ .x$x[,1:10] %>%
         as_tibble(rownames = "sample_id")) %>%
  bind_rows(.id = 'contrast')

# extracting proportion of variance
pca_variance <-
  pca %>%
  map(~ summary(.x)$importance[2,1:10] * 100) %>% 
  modify(ceiling) %>% 
  modify(~ purrr::set_names(x = paste(names(.x), ":", .x, "%", sep = ""), nm = names(.x))) %>% 
  purrr::transpose() %>% 
  modify(unlist) 

# Extract metadata: sample_id, exposure and timepoint
meta <- 
  read_csv(here::here("results/metadata/sample_metadata.csv")) %>%
  extract(c("sample_id", "exposure", "timepoint", "tissue", "maternal_fetal"))

# Join PCA data with metadata, exclude TiO2 samples 

pca_joined <-
  pca_data %>%
  full_join(x = ., y = meta,
            by = c("sample_id")) %>%
  .[1:414,] %>%
  .[.$exposure %in% c("LPS", "ctr"),]

# PCAs
pca_labels <-
  purrr::set_names(
    x = paste(
      names(pca_variance[[1]]),
      pca_variance$PC1,
      pca_variance$PC2,
      sep = "\n"
    ), 
    nm = names(pca_variance[[1]])
  )

ggplot(pca_joined) +
  geom_point(aes(
    x = PC1, y = PC2, col = maternal_fetal, shape = exposure)) +
  facet_grid(
    timepoint ~ contrast,
    labeller = labeller(contrast = pca_labels))

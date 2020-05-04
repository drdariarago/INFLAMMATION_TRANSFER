library(tidyverse)
library(magrittr)
library(DESeq2)

#https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#subsetting

#Using expression results because VSD gives less confidences for genes with low counts

experiment_data <-
  readRDS(here::here("results/tximeta/vsd.Rdata"))

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
  left_join(x = ., y = meta,
            by = c("sample_id")) %>%
  filter(
    exposure %in% c("LPS", "ctr")
  )

# PCAs
pca_labels <-
  purrr::set_names(
    x = paste(
      names(pca_variance[[1]]),
      pca_variance$PC3,
      pca_variance$PC4,
      sep = "\n"
    ), 
    nm = names(pca_variance[[1]])
  )

# Separating by timepoint 
ggplot(pca_joined) +
  geom_point(
    aes(
      x = PC3, y = PC4, shape = maternal_fetal, color = exposure),
    alpha = 0.8, size = 3) +
  facet_grid(
    timepoint ~ contrast,
    labeller = labeller(contrast = pca_labels)) + 
  scale_colour_brewer(type = 'qual', palette = 7)


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

# Separating by exposure, timepoint discrete
ggplot(pca_joined) +
  geom_point(
    aes(
      x = PC1, y = PC2, shape = maternal_fetal, color = timepoint),
    alpha = 0.8, size = 3) +
  facet_grid(
    exposure ~ contrast,
    labeller = labeller(contrast = pca_labels)) + 
  scale_colour_viridis_c()

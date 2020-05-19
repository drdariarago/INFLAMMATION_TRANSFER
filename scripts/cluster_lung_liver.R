## Cluster lung and maternal liver response via difference in the area under normalized curves

library(tidyverse)
library(magrittr)
library(cluster)

## Import log FC data
model_results <- 
  here::here("results/limma_maternal_liver/fold_change_summary.rds") %>% 
  read_rds() %>% 
  filter(
    exposure == 'response'
  ) 

## Filter only genes that are significant across all timepoints
## Then center (subtract mean), scale (divide by sd) and convert to positive (add constant)
normalized_fold_change_matrix <-
  model_results %>% 
  group_by(ensembl_gene_id) %>% 
  filter(any(q_value < 0.05 & logFC > 0.5)) %>% 
  mutate(
    logFC = ifelse(q_value < 0.01, logFC, 0.001)
  ) %>%
  pivot_wider(id_cols = timepoint, names_from = ensembl_gene_id, values_from = logFC)  %>% 
  column_to_rownames(var = "timepoint") %>% 
  as.matrix() %>% 
  exp() %>% 
  scale(x = ., center = TRUE, scale = FALSE) %>% 
  add(10) %>% 
  t() %>% 
  as_tibble(x = ., rownames = 'ensembl_gene_id')


normalized_fold_change_matrix %>% 
  pivot_longer(data = ., cols = -ensembl_gene_id, names_to = "timepoint", values_to = 'normalized_fc') %>% 
  ggplot(data = ., mapping = aes(x = timepoint, y = normalized_fc, group = ensembl_gene_id)) +
  geom_line(alpha = 0.05) 

## Functionalize following section, just use a vector 1:4 for each base calculation
normalized_distance_matrices <-
  normalized_fold_change_matrix %>%  
  map(
    .x = ., 
    .f = ~ .x
  ) %>% 
  magrittr::extract(
    grepl(pattern = "timepoint[0-9]{1,2}", x = names(.))
  ) %>%
  map(
    .x = .,
    .f = ~ outer(X = .x, Y = .x, FUN = `-`) 
  ) %>% 
  map(
    .f = ~ matrix(data = .x, nrow = nrow(normalized_fold_change_matrix),
                  dimnames = list(
                    normalized_fold_change_matrix$ensembl_gene_id,
                    normalized_fold_change_matrix$ensembl_gene_id)
    )
  ) 

# trapezoid area formula for each time step, proportional to difference in hours
# caveat: if the 2 sides have different sign we had an inversion: divide the area by another 2 to get the area of 2 triangles instead
delta_hours = c(1,1,1) #c(3,7,12)
timepoints = names(normalized_distance_matrices)

areas <-
  normalized_distance_matrices %>% 
  map2(
    .x = .[1:3],
    .y = .[-1],
    .f = ~ .x + .y
  ) %>% 
  set_names(
    x = ., nm = c("interval1", "interval2", "interval3")
  ) %>% 
  map2(
    .x = ., 
    .y = delta_hours,
    .f = ~ (.x/2) * .y
  ) %>% 
  map(
    .x = .,
    .f = ~ ifelse(.x <= 0, abs(.x/2), .x)
  )

# Create final distance
delta_areas <- reduce(.x = areas, .f = `+`)

# Create groups (test with basic hclust)
tree <-
  delta_areas %>% 
  as.dist() %>% 
  hclust(d = ., method = "single") 
plot(tree)

hclust_clusters <- 
  tree %>%
  cutree(k = 20) 

# Use PAM
pam_clusters <-
  delta_areas %>% 
  pam(x = ., diss = TRUE, k = 8, cluster.only = TRUE)


# Merge clusters with previous data and plot mean vs spread of expression for each
pam_clusters %>% 
  as_tibble(rownames = "ensembl_gene_id") %>% 
  full_join(x = ., y = model_results, by = "ensembl_gene_id") %>% 
  filter(!is.na(value)) %>%  
  group_by(ensembl_gene_id) %>% 
  mutate(
    logFC = scale(logFC, center = FALSE, scale = FALSE)
  ) %>% 
  ggplot(data = ., 
         mapping = aes(x = timepoint, y = logFC, group = ensembl_gene_id)) +
  geom_line(alpha = 0.1) + facet_wrap(~ value) +
  geom_line(data = . %>% group_by(value, timepoint) %>% summarise(logFC = median(logFC)),
            aes(group = value), col = 'hotpink'
  ) + 
  coord_cartesian(ylim = c(-1,2))


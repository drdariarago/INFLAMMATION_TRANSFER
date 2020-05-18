## Cluster lung and maternal liver response via difference in the area under normalized curves

library(tidyverse)
library(magrittr)
library(pracma)

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
  filter(all(q_value < 0.05)) %>% 
  pivot_wider(id_cols = timepoint, names_from = ensembl_gene_id, values_from = logFC)  %>% 
  column_to_rownames(var = "timepoint") %>% 
  as.matrix() %>% 
  scale(x = ., center = TRUE, scale = FALSE) %>% 
  add(10) %>% 
  t()

normalized_fold_change_matrix %>% 
  as_tibble(x = ., rownames = 'ensembl_gene_id') %>% 
  pivot_longer(data = ., cols = -ensembl_gene_id, names_to = "timepoint", values_to = 'logFC') %>% 
  ggplot(data = ., mapping = aes(x = timepoint, y = logFC, group = ensembl_gene_id)) +
  geom_line(alpha = 0.1) 

## Functionalize following section, just use a vector 1:4 for each base calculation
# Calculate base1 
time1 <- outer(normalized_fold_change_matrix[,1], normalized_fold_change_matrix[,1], FUN = '-') 
# Calculate base2 
time2 <- outer(normalized_fold_change_matrix[,2], normalized_fold_change_matrix[,2], FUN = '-') 
# Calculate base3 
time3 <- outer(normalized_fold_change_matrix[,3], normalized_fold_change_matrix[,3], FUN = '-') 
# Calculate base 4
time4 <- outer(normalized_fold_change_matrix[,4], normalized_fold_change_matrix[,4], FUN = '-') 

# trapezoid area formula for each time step, proportional to difference in hours
# caveat: if the 2 sides have different sign we had an inversion: divide the area by another 2 to get the area of 2 triangles instead
delta_hours = c(3,7,12)

areas <-
  list(
    area1 = time1 + time2,
    area2 = time2 + time3,
    area3 = time3 + time4
  ) %>% 
  map(.x = ., ~ .x/2) %>% 
  map2(
    .x = ., 
    .y = delta_hours,
    .f = ~ .x * .y
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
  cutree(h = 2.5) 

# Use PAM
pam_clusters <-
  delta_areas %>% 
  pam(x = ., diss = TRUE, k = 9, cluster.only = TRUE)


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

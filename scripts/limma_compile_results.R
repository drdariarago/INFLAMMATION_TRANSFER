### Compile all limma results in a single table

library(tidyverse)
library(magrittr)

fold_change_summary_maternal <-
  snakemake@input[["limma_results"]]  %>%
  setNames(object = ., nm = snakemake@params[["models"]]) %>% 
  map_dfr(.f = read_rds, .id = "tissue") %>% 
  mutate(
    tissue = as.factor(tissue)
  )

write_csv(
  x = fold_change_summary_maternal,
  path = snakemake@output[["maternal"]]
)

fold_change_summary <-
  fold_change_summary_maternal %>% # Remove chorion/decidua contrasts
  filter(is.na(maternal) | maternal == "shared") %>% 
  select(-maternal)

write_csv(
  x = fold_change_summary,
  path = snakemake@output[["summary"]]
)
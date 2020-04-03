## SAMPLE METADATA
# Create metadata tables for all samples

library(tidyverse)
library(here)

sample_metadata_list <- 
  c("ml", "mli", "fli", "fp", "mp") %>% 
  paste0(., "_metadata.csv") %>% 
  paste0("data/00_metadata/", .) %>% 
  here(.) %>% 
  lapply(X = ., FUN = read_delim, col_names = TRUE, delim = ';')

#Using Rbind to combine the different metadata data-frame/vectors by columns.

sample_metadata_table <- 
  plyr::ldply(sample_metadata_list)

# Make new columns
sample_metadata_tibble <- 
  sample_metadata_table %>% 
  as_tibble(.) %>% 
  transmute(
    animal_no = 
      str_extract(string = Animal_no, pattern = "[0-9]{1,}X?$"),
    sample_id = 
      str_replace(string = Animal_no, pattern = "lLi", replacement = "Li"),
    tissue = 
      str_extract(string = sample_id, pattern = "^[M,F][L,P]i?") %>%
      factor(x=., 
             levels = c("ML","MP","FP","MLi","FLi"), 
             labels = c("lung","placenta","placenta","liver", "liver")),
    maternal_fetal = 
      str_extract(string = sample_id, pattern = "^[M,F]") %>% 
      factor(x = ., 
             levels = c("M","F"), 
             labels = c("maternal","fetal")),
    exposure = 
      Exposure %>% 
      factor(x = ., 
             levels = c("nanopurevand","nanopure","Control","ctr","Ctr","LPS","TiO2"), 
             labels = c("ctr","ctr","ctr","ctr","ctr","LPS","TiO2")),
    rna_date = 
      rna_date,
    timepoint =
      timepoint
  )

# Saving the r script to RDS

saveRDS(object = sample_metadata_tibble, file = snakemake@output[["rdata"]])

# Saving it as csv file

write.csv(x = sample_metadata_tibble, file = snakemake@output[["csv"]], row.names = FALSE)


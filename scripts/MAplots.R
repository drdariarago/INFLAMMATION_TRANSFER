library(tidyverse)
library(ggrepel)

## Test MA - just for placenta fetal (look at MA_responses_all_tissues.R 
## for a version for all tissues)

limma_coefs <- 
  readRDS(here::here("results/limma/limma_coefs.Rdata"))

placenta_fetal_tbl <- 
  limma_coefs[["placenta_fetal"]] %>%
  filter(treatment == "exposureLPS") %>%
  as_tibble() %>%
  select(-c(treatment))

placenta_fetal_tbl$gene_id <-
  gsub('\\..*', '', 
       placenta_fetal_tbl$gene_id)

annotated_results <- 
  readRDS(here::here("results/limma_results/annotated_results.Rdata")) 

placenta_fetal_annotated <-
  annotated_results[["placenta_fetal"]] 

placenta_fetal_longer <-
  placenta_fetal_annotated %>%
  pivot_longer( 
    cols = c("timepoint2", "timepoint5", "timepoint12", "timepoint24"), 
    names_to = "timepoint")

placenta_fetal_all <-
  full_join(placenta_fetal_longer, placenta_fetal_tbl, 
            by = c("ensembl_gene_id" = "gene_id", "timepoint" = "timepoint"))

done <-  
  placenta_fetal_all %>%
  mutate(
    value = 
      ifelse( (is.na(value) | abs(LogFC) < 0.6), 0, value) %>% 
      factor(x = ., 
             levels = c(0, 1), 
             labels = c("not_significant", "significant")),
    timepoint = 
      factor(timepoint, 
             levels = paste0("timepoint", c(2, 5, 12, 24)), 
             labels = paste0(c(2, 5, 12, 24), "hours after exposure" ))
  )

placenta_fetal_MA <-
  ggplot(done, 
         aes(x = AveExpr, y = LogFC)) +
  geom_point(alpha = 0.3) +
  geom_point(data = subset(done, value == "significant"), aes(x = AveExpr, y = LogFC), col = 'hotpink') +
  geom_label_repel(data = subset(done, value == "significant" & AveExpr > 5), aes(x = AveExpr, y = LogFC, label = mgi_symbol), box.padding = 0.5) +
  facet_grid(timepoint ~ .) + 
  ggtitle("Placenta fetal all genes") +
  scale_y_continuous(breaks = c(seq(from = -3, to = 3, by = 1)))

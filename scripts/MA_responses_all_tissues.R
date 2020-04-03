library(tidyverse)
library(ggrepel)

limma_coefs <- 
  readRDS(here::here("results/limma/limma_coefs.Rdata"))

placenta_fetal_tbl <- 
  limma_coefs %>% 
  map(.x = ., .f = ~ filter(.x, treatment == "exposureLPS")) %>% 
  map(.x = ., .f = ~ select(.x, -c(treatment))) %>% 
  map(.x = ., .f = ~ mutate(.data = .x, gene_id = gsub('\\..*', '', gene_id)))

annotated_results <- 
  readRDS(here::here("results/limma_results/annotated_results.Rdata")) %>% 
  map(.x = ., 
      .f = ~ pivot_longer( 
        data = .x, 
        cols = c("timepoint2", "timepoint5", "timepoint12", "timepoint24"), 
        names_to = "timepoint", values_to = "significant"
        )
  ) %>% 
  map(.x = ., 
      .f = ~ dplyr::select(.data = .x, -3)
  )

placenta_fetal_all <-
  map2(.x = annotated_results, .y = placenta_fetal_tbl, 
       .f = ~ full_join(
         .x, .y,
         by = c("ensembl_gene_id" = "gene_id", "timepoint" = "timepoint")
       )
  ) %>% 
  map(.x = ., 
      .f =  ~ mutate(
        .data  = .x, 
        significant = 
          ifelse( (is.na(significant) | abs(LogFC) < 1.0 | AveExpr < 1), 0, significant) %>% 
          factor(x = ., 
                 levels = c(0, 1), 
                 labels = c("not_significant", "significant")),
        timepoint = 
          factor(timepoint, 
                 levels = paste0("timepoint", c(2, 5, 12, 24)), 
                 labels = paste0(c(2, 5, 12, 24), " hours after exposure" ))
      )
  ) 

MA_plots <- 
  imap(
    .x = placenta_fetal_all, 
    .f = ~ ggplot(
      data = .x, 
      aes(x = AveExpr, y = LogFC)) +
      geom_point(alpha = 0.3) +
      geom_point(
        data = subset(.x, significant == "significant"), 
        aes(x = AveExpr, y = LogFC), col = 'hotpink') +
      geom_label_repel(
        data = subset(.x, significant == "significant" & AveExpr > median(AveExpr) & abs(LogFC) > median(abs(LogFC))), 
        aes(x = AveExpr, y = LogFC, label = mgi_symbol), box.padding = 0.5) +
      facet_grid(timepoint ~ .) + 
      xlim(c(0,10)) +
      ggtitle(label = .y) +
      scale_y_continuous(breaks = c(seq(from = -3, to = 3, by = 1)))
  )

imap(
  .x = MA_plots,
  .f = ~ ggsave(
    plot = .x,
    path = here::here("results/MA_responses_all_tissues"),
    filename = paste0(.y, "_MA_plot.png"),
    width = 6.67, height = 7.5, units = 'in',
    device = "png"
  )
)

## Plot frequency of intersecting responses for timepoint by tissue

library(UpSetR)
library(tidyverse)
library(magrittr)

annotated_results <- 
  readRDS(here::here("results/limma_results/annotated_results.Rdata")) %>% 
  map(.x = ., 
      .f = ~ dplyr::select(.data = .x, -7)
  )  


set_plots <- imap(
  .x = annotated_results,
  .f = UpSetR::upset, 
  nintersects = NA,
  sets = c("timepoint24", "timepoint12", "timepoint5", "timepoint2"),
  keep.order = T,
  order.by = "freq",
  empty.intersections = T,
  mainbar.y.max = 1500
)

pdf(file = here::here("results/upsetr/maternal_lung.pdf"), width = 6.67, height = 7.5)
set_plots[c("lung_maternal")]
dev.off()

pdf(file = here::here("results/upsetr/maternal_liver.pdf"), width = 6.67, height = 7.5)
set_plots[c("liver_maternal")]
dev.off()

pdf(file = here::here("results/upsetr/maternal_placenta.pdf"), width = 6.67, height = 7.5)
set_plots[c("placenta_maternal")]
dev.off()

pdf(file = here::here("results/upsetr/fetal_placenta.pdf"), width = 6.67, height = 7.5)
set_plots[c("placenta_fetal")]
dev.off()

pdf(file = here::here("results/upsetr/fetal_liver.pdf"), width = 13.33, height = 7.5)
set_plots[c("liver_fetal")]
dev.off()

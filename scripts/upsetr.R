## Plot frequency of intersecting responses for timepoint by tissue

library(UpSetR)
library(tidyverse)
library(magrittr)
library(janitor)

ALPHA = 0.05
MIN_LOGFC = 0.5
tissue <- c("fetal_liver", "maternal_liver")

annotated_results <- 
  glue::glue("results/limma_{tissue}/fold_change_summary.rds") %>% 
  here::here() %>% 
  map( readRDS ) %>% 
  set_names( tissue ) %>% 
  map(
    .f = filter,
    q_value < ALPHA, 
    abs(logFC) > MIN_LOGFC,
    exposure == "response",
    mgi_symbol != ""
  ) %>% 
  map(
    .f = select,
    mgi_symbol, timepoint
  ) %>%
  map(
    .f = ~ janitor::tabyl( .x, mgi_symbol, timepoint)
  ) %>% 
  map( as.data.frame)


set_plots <- 
  imap(
    .x = annotated_results,
    .f = upset,
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

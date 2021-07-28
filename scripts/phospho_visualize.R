### Visualize phospho data results

library(tidyverse)
library(magrittr)

site_results <-
  here::here("results/phospho_import/sitewise_results.rds") %>% 
  read_rds()

## MA plots for sites

site_results %>% 
  ggplot(
    aes(
      x = (m_ctrl + m_LPS) / 2,
      y = logFC
    )
  ) +
  geom_hex() +
  viridis::scale_fill_viridis(trans = "log", breaks = c(1,5,25,125), direction = -1) +
  scale_x_log10() +
  # scale_y_continuous(limits = c(-2,2)) +
  facet_grid(tissue ~ timepoint) +
  geom_point(
    data = site_results %>% filter(FDR < 0.05, abs(logFC) > 1 ),
    alpha = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = site_results %>% filter(FDR < 0.05, logFC > 1 ),
    aes(label = mgi_symbol)
  ) +
  geom_rug(
    data = site_results %>% filter(FDR < 0.05, abs(logFC) > 1 ),
    alpha = 0.1
  )

ggsave(here::here("results/phospho_visualize/site_ma_plot.pdf"))

## Volcano plot for sites

site_results %>% 
  ggplot(
    aes(
      x = logFC,
      y = FDR
    )
  ) +
  geom_hex() + 
  viridis::scale_fill_viridis(trans = "log", breaks = c(4^(1:7)), direction = -1) +
  facet_grid(tissue ~ timepoint) +
  scale_y_continuous(trans = 'log10', limits = c(1e-06, 1), oob = scales::squish) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::squish) +
  geom_hline(yintercept = 0.05, col = 'hotpink') +
  ggrepel::geom_text_repel(
    data = site_results %>% filter(FDR < 0.001, abs(logFC) > 0.75 ),
    aes(label = mgi_symbol)
  ) +
  ggtitle(label = "Volcano plot for site-wise phosphorilation")

ggsave(filename = here::here("results/phospho_visualize/site_volcano_plot.pdf"), 
       width = 15, height = 10)

## Same but on averaged gene results

gene_results <-
  here::here("results/phospho_import/genewise_results.rds") %>%
  read_rds()

## Volcano plot for genes

gene_results %>%
  ggplot(
    aes(
      x = min_log_fc,
      y = q_value_avg
    )
  ) +
  geom_hex() +
  viridis::scale_fill_viridis(trans = "log", breaks = c(3^(1:6)), direction = -1) +
  facet_grid(tissue ~ timepoint) +
  geom_hline(yintercept = 0.05) +
  scale_y_continuous(trans = 'log10', limits = c(1e-14, 1), oob = scales::squish) +
  scale_x_continuous(limits = c(-2,2), oob = scales::squish) +
  geom_hline(yintercept = 0.05, col = 'hotpink') +
  ggrepel::geom_text_repel(
    data = gene_results %>% filter(q_value_avg < 1e-3, abs(min_log_fc) > 1 ),
    aes(label = mgi_symbol)
  ) +
  ggtitle(label = "Volcano plot for gene-wise phosphorilation")

ggsave(filename = here::here("results/phospho_visualize/gene_volcano_plot.pdf"), 
       width = 15, height = 10)

gene_results %>%
  ggplot(
    aes(
      x = min_log_fc,
      y = q_value_min
    )
  ) +
  geom_hex() +
  viridis::scale_fill_viridis(trans = "log", breaks = c(3^(1:6)), direction = -1) +
  facet_grid(tissue ~ timepoint) +
  geom_hline(yintercept = 0.05) +
  scale_y_continuous(trans = 'log10', limits = c(1e-7, 1), oob = scales::squish) +
  scale_x_continuous(limits = c(-2,2), oob = scales::squish) +
  geom_hline(yintercept = 0.05, col = 'hotpink') +
  ggrepel::geom_text_repel(
    data = gene_results %>% filter(q_value_avg < 1e-04, abs(min_log_fc) > 1 ),
    aes(label = mgi_symbol)
  ) +
  ggtitle(label = "Volcano plot for gene-wise phosphorilation")


# Import raw results and plot PCA of samples

data_path <- here::here("data/proteomics/20210202_second_run/MaxQuant_table_all4exp_Jan_2021.xlsx")

phospho_data <-
  data_path %>% 
  readxl::excel_sheets() %>% 
  setNames(object = ., nm = .) %>% 
  map(
    .id = "contrast_id",
    .x = .,
    .f = ~ readxl::read_xlsx(
      data_path,
      sheet = .x
    ) %>% 
      select(
        "LP_position", starts_with("Reporter")
      ) %>% 
      magrittr::set_names(
        janitor::make_clean_names(
          string = names(.) %>% 
            gsub(pattern = "Reporter Intensity_", replacement = ""))
      ) %>% 
      column_to_rownames("lp_position")
  ) 

phospho_data <-
  here::here("results/phospho_normalize/normalized_results_list.rds") %>% 
  readRDS()

shared_sites <-
  phospho_data %>% 
  map(.f = rownames) %>% 
  reduce(.f = intersect)

shared_pca_data <-
  phospho_data %>% 
  set_names( 
    names(.) %>% 
      janitor::make_clean_names(parsing_option = 3) %>% 
      str_extract(pattern = "(liver|placenta)")
    ) %>% 
  map(
    .f = ~ .[rownames(.x) %in% shared_sites,]
  ) %>%
  imap(
    .f = ~ colnames(.x) %>% 
      janitor::make_clean_names() %>% 
      str_extract(string = ., pattern = "(lps|control)_[5,2]h_[:digit:]") %>% 
      paste(.y, ., sep = "_") %>% 
      set_colnames(.x, .)
  ) %>%
  reduce(.f = cbind)
# %>% 
#   filter(
#     across( everything(), .fns = ~ var(.x) > 0)
#   )
  
phospho_pca_data <-
  shared_pca_data %>% 
  t %>% 
  prcomp(scale = TRUE)

summary(phospho_pca_data)

pca_plot_data <-
  phospho_pca_data %$% 
  x %>% 
  as_tibble(rownames = "sample_id") %>% 
  mutate(
    sample_id = tolower(sample_id),
    timepoint = str_extract(string = sample_id, pattern = "[2,5]h"),
    tissue = str_extract(string = sample_id, pattern = "(placenta|liver)"),
    treatment = str_extract(string = sample_id, pattern = "control|lps")
  ) 

pca_plot_data %>% 
  ggplot(
    aes(x = PC2, y = PC3, col = treatment, fill = treatment, group = tissue)
  ) +
  ggforce::geom_voronoi_tile(alpha = 0.25) +
  geom_point(alpha = 1) +
  facet_grid( tissue ~ timepoint ) +
  scale_color_brewer(type = 'qual', name = "treatment") + 
  scale_fill_brewer(type = 'qual', name = "treatment") +
  theme_minimal() +
  scale_x_continuous(name = "PC2\n[8% variance]") +
  scale_y_continuous(name = "PC3\n[3% variance]") 

ggsave(filename = here::here('results/phospho_visualize/pca_plot.pdf'))

## Create dendrogram of samples
dendrogram <-
  shared_pca_data %>% 
  set_colnames(
    paste(
      ifelse(grepl("Liver", x = colnames(.) ), "liv", "pla"),
      ifelse(grepl("lps", x = colnames(.) ), "lps", "ctrl"),
      ifelse(grepl("5h", x = colnames(.) ), "5h", "2h"),
      sep = "_"
    )
  ) %>% 
  t %>% 
  scale %>% 
  dist %>% 
  hclust(method = "average") %>% 
  as.dendrogram()

dendextend::labels_colors(dendrogram) <- 
  grepl("2h", labels(dendrogram)) + ( 1 + grepl("lps", labels(dendrogram)) ) * 2

pdf(file = here::here('results/phospho_visualise/dendrogram.pdf'))
plot(dendrogram)
dev.off()

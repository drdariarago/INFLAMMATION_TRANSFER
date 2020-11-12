library(tidyverse)
library(magrittr)
library(gprofiler2)

# Import parameters from snakemake
min_log_fc = snakemake@params[['min_log_fc']]
min_log_base_counts = 1
max_fdr = snakemake@params[["max_fdr"]]
direction = snakemake@params[['up_or_down']]
model = snakemake@params[['model']]

# Import compiled limma results and subset tissue of interest
limma_results <-
  # here::here("results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>%
  snakemake@input[['limma_results']] %>%
  read_csv() %>% 
  filter(tissue == model)

# Plot thresholds vs values for QC
limma_results %>% 
  ggplot(aes(x = logFC)) + 
  geom_density() + 
  facet_wrap(exposure ~ timepoint, scales = "free") + 
  geom_vline(
    data = data.frame(x = min_log_base_counts, exposure = "baseline"),
    aes(xintercept = x), col = 'red') +
  geom_vline(
    data = data.frame(x = min_log_fc, exposure = "response"),
    aes(xintercept = x), col = 'blue') +
  ggtitle(
    glue::glue(
      "Distribution and thresholds for base and fold change values in {model}"
      )
  )

ggsave(
  filename = snakemake@output[['threshold_plot']], 
  width = 297, height = 210, units = 'mm'
)

# Set background gene expression as all genes with baseline log counts > 1
background_genes <-
  limma_results %>% 
  filter(
    exposure == "baseline",
    logFC > min_log_base_counts
  ) %>% 
  pull(ensembl_gene_id) %>% 
  unique()

# Filter log fold change above threshold in the right direction
filtered_results <-
  if (direction == "upregulated") {
    limma_results %>% 
      filter(exposure == "response") %>% 
      filter( logFC > min_log_fc ) %>% 
      arrange( -logFC )
  } else if (direction == "downregulated") {
    limma_results %>% 
      filter(exposure == "response") %>% 
      filter( logFC < -min_log_fc ) %>% 
      arrange( logFC )
  } else {
    print("Direction needs to be either 'downregulated' or 'upregulated'")
  }

# Convert to list of ranked gene IDs
ranked_gene_list<-
  filtered_results$timepoint %>% 
  unique %>% 
  purrr::set_names() %>% 
  map(
    .f = ~ filtered_results[filtered_results$timepoint == .x, "ensembl_gene_id"]
  ) %>% 
  map(.f = pull)

# Run GO enrichment
gost_result_list <-
  map(
    .x = ranked_gene_list,
    .f = ~ gprofiler2::gost(
      query = .x,
      ordered_query = TRUE, 
      multi_query = FALSE,
      organism = "mmusculus",
      sources = c("GO"), 
      domain_scope = "custom", 
      custom_bg = background_genes,
      correction_method = "gSCS",
      user_threshold = max_fdr,
      significant = FALSE,
      measure_underrepresentation = FALSE
    )
  )

write_rds(x = gost_result_list, file = snakemake@output[["raw_results"]])
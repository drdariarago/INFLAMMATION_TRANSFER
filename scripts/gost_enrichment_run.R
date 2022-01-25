library(tidyverse)
library(magrittr)
library(gprofiler2)

# Import parameters from snakemake
# min_log_base_counts = 1
direction = snakemake@params[['up_or_down']]
model = snakemake@params[['model']]

# Import compiled limma results and subset tissue of interest
limma_results <-
  # here::here("results/limma_compile_results/limma_results_no_maternal_contrasts.csv") %>%
  snakemake@input[['limma_results']] %>%
  read_csv() %>% 
  filter(tissue == model) 

# Set background gene expression as all expressed genes across all timepoints
background_genes <-
  limma_results %>% 
  filter(
    exposure == "baseline",
    logFC > 0
  ) %>% 
  dplyr::pull(ensembl_gene_id) %>% 
  unique()

# Filter log fold change above threshold in the right direction
filtered_results <-
  limma_results %>% 
  filter(
    exposure == "response",
    # Ensure that we select only genes present in background
    ensembl_gene_id %in% background_genes
  ) %>% 
  arrange( log10(q_value) )

directional_results <-
  if (direction == "upregulated") {
    filtered_results[filtered_results$logFC > 0,]
  } else if (direction == "downregulated") {
    filtered_results[filtered_results$logFC < 0,]
  } else {
    print("Direction needs to be either 'downregulated' or 'upregulated'")
  }

# Convert to list of ranked gene IDs
ranked_gene_list<-
  directional_results$timepoint %>% 
  unique %>% 
  purrr::set_names() %>% 
  map(
    .f = ~ directional_results[directional_results$timepoint == .x, "ensembl_gene_id"]
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
      user_threshold = 0.05,
      significant = FALSE,
      measure_underrepresentation = FALSE
    )
  )

write_rds(x = gost_result_list, file = snakemake@output[["raw_results"]])
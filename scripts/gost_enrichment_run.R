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
  filter(tissue == model) %>% 
  filter(exposure == "response")

# Set background gene expression as all responder genes across all timepoints
background_genes <-
  limma_results %>% 
  pull(ensembl_gene_id) %>% 
  unique()

# Filter log fold change above threshold in the right direction
filtered_results <-
  if (direction == "upregulated") {
    limma_results %>% 
      arrange( -log10(q_value) * -sign(logFC) )
  } else if (direction == "downregulated") {
    limma_results %>% 
      arrange( -log10(q_value) * sign(logFC) )
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
      user_threshold = 0.05,
      significant = FALSE,
      measure_underrepresentation = FALSE
    )
  )

write_rds(x = gost_result_list, file = snakemake@output[["raw_results"]])
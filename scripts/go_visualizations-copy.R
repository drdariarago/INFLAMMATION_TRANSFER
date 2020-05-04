library(tidyverse)
library(magrittr)

## DOWNREFULATED PLACENTA
down_genes_csv <-
  read_csv(here::here("results/gost_enrichment_format/downregulated_response.csv")) 

down_tissue_list <-
  down_genes_csv %>% 
  # Order timepoint according to our correct scheme and rename in human readable terms
  mutate(
    timepoint = factor(
      x = timepoint, 
      levels = c("2_hrs", "5_hrs", "12_hrs", "24_hrs"),
      labels = paste( c(2,5,12,25), "hrs post\nexposure", sep = " "))
  ) %>% 
  list(
    lung_mat = filter(.data = ., contrast == "lung_maternal"),
    liver_mat = filter(.data = ., contrast == "liver_maternal"),
    liver_fet = filter(.data = ., contrast == "liver_fetal"),
    placenta_all_down = filter(.data = ., contrast %in% c("placenta_maternal", "placenta_fetal"))
  ) %>% 
  # This section arranges the GO IDs by term size
  # I do this by specifying term_name as a factor with order set by the term_size
  map(
    .x = .,
    .f = ~ mutate(
      .data = .x,
      term_name = factor(
        x = term_name, 
        levels = arrange(.x, term_size) %>% pull(term_name) %>% unique()
      )
    )
  )

# ggplot
down_tissue_list %>% 
  # Filter terms of interest
  # Two examples: term type (BP, MF or CC) and minimum term_size
  map(.x = ., 
      .f = ~ filter(
        .data = .x, 
        source == "GO:BP",
        term_size > 10
      ) 
  ) %>% 
  map(.x = .,
      # this line picks the top N GOs, ordered by term size. 
      # Use increments of 5 since each GO is 5 lines (one per timepoint)
      .f = ~.x[1:100,]) %>% 
  imap(
  .x = ., 
  .f = ~ ggplot(data = .x, mapping = aes(x = timepoint, y = term_name)) +
    geom_point(aes(
      size = term_size,
      color = significant)) +
    theme_bw(base_size = 14) +
    scale_colour_brewer(type = 'qual', palette = 2) +
    ggtitle(label = paste("Downregulated GO terms in", .y)) +
    facet_grid(~ contrast)
)

##UPREGULATED PLACENTA

## Repeat as above, with a different dataset :)

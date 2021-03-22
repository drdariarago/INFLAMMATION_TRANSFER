### Filter expressed receptors in target tissues and rank DE ligands in messenger tissues

#### Load libraries and data ####

library(tidyverse)
library(magrittr)
library(ggplot2)
library(patchwork)
library(ggrepel)

a4 = c(297, 210)

# Check for presence of cytokines

cytokines <-
  read_tsv(file = "results/match_orthologs/human_mouse_ligands_receptors.txt") %>% 
  filter(
    grepl(pattern = "^(I|Cxc)l.*", secreted_gene_symbol_mouse)
  ) %>% 
  dplyr::pull(secreted_ensembl_gene_mouse) %>% 
  unique()

receptor_ligand_pairs <-
  read_tsv(file = "results/match_orthologs/human_mouse_ligands_receptors.txt") %>% 
  select(matches("ensembl.*mouse")) %>% 
  set_names(c("ligand", "receptor")) %>% 
  filter(!is.na(receptor), !is.na(ligand)) %>% 
  distinct()

fold_change_summary <-
  c("fetal_liver", "maternal_liver", "maternal_lung", "placentas") %>% 
  setNames(., .) %>% 
  map(~ paste("results/limma_", .x, "/fold_change_summary.rds", sep = "")) %>% 
  map_dfr(read_rds, .id = "tissue") %>% 
  mutate(
    tissue = as.factor(tissue)
  ) %>% # Remove chorion/decidua contrasts
  filter(is.na(maternal) | maternal == "shared") %>% 
  select(-maternal)

# Define which tissues we scan for ligands (messenger tissues) and receptors (target tissues)
# Also define other global criteria

messenger_tissues <- c("maternal_lung", "maternal_liver")
target_tissues <- c("placentas", "fetal_liver")
timepoints <- paste("timepoint", c(2,5,12,24), sep = "")
log_fc_threshold <- 0.5
q_value_threshold <-0.05

# Subset only receptors that are expressed (more than 8 counts baseline in at least one timepoint of interest) 
# in the target tissues and the ligands that target them

expressed_receptors <-
  fold_change_summary %>% 
  filter(tissue %in% target_tissues) %>%
  filter(timepoint %in% timepoints) %>% 
  filter(exposure == "baseline", logFC > 3) %>% 
  select(expressed_receptor = ensembl_gene_id) %>% 
  distinct() 

# Subset ligands differentially expressed in messenger tissues at 2 and 5 hours (logFC > 0.5, q-val < 0.05)
# Remove ligands whose receptors are not expressed in target tissues

expressed_ligands <-
  fold_change_summary %>% 
  filter(tissue %in% messenger_tissues) %>% 
  filter(timepoint %in% timepoints ) %>% 
  filter(exposure == "response") %>%   
  filter(q_value < q_value_threshold, abs(logFC) > log_fc_threshold ) %>% 
  select(expressed_ligand = ensembl_gene_id ) %>% 
  unique() 

# Merge expressed receptor and differentially expressed ligands, discard interactions without a partner

expressed_receptor_ligand_pairs <-
  receptor_ligand_pairs %>% 
  filter(
    ligand %in% expressed_ligands$expressed_ligand,
    receptor %in% expressed_receptors$expressed_receptor
  )
  
# Extract annotation for all genes in the final ligand-receptor pairs

ligand_receptor_annotation <-
  fold_change_summary %>% 
  filter(
    timepoint %in% timepoints
  ) %>% 
  select( ensembl_gene_id, tissue, timepoint, exposure, log_fc = logFC, q_value ) %>% 
  pivot_wider( names_from = exposure, values_from = c(log_fc, q_value) )

# Annotate receptor-ligand pairs in long format
annotated_receptor_ligand_pairs <-
  expressed_receptor_ligand_pairs %>% 
  mutate(
    pair_id = paste(ligand, receptor, sep = "_")
  ) %>% 
  pivot_longer(
    cols = c("ligand", "receptor"), 
    names_to = "gene_class", 
    values_to = "ensembl_gene_id"
  ) %>% 
  left_join(
    x = ., y = ligand_receptor_annotation,
    by = "ensembl_gene_id"
  ) %>% # Filter only receptor data from target tissues and ligand data from messenger tissues
  filter(gene_class == "receptor" & tissue %in% target_tissues | gene_class == "ligand" & tissue %in% messenger_tissues) %>% 
  filter( # Filter out receptors with no expression (baseline q > 0.05) and ligands with no response (baseline q > 0.05 or response fc < 0.5)
    gene_class == "receptor" & q_value_baseline < q_value_threshold | gene_class == "ligand" & q_value_response < q_value_threshold & log_fc_response < log_fc_threshold
  )

# Annotate if the ligands have an expressed interacting receptor and vice versa or if the ligand/receptor is orphan
annotated_receptor_ligand_pairs_orphaned <-
  annotated_receptor_ligand_pairs %>% 
  select(pair_id, gene_class, timepoint) %>% 
  distinct() %>% 
  group_by( pair_id, timepoint ) %>% 
  summarize(
    paired = ifelse( n() > 1, "paired", "orphan") %>% as.factor()
  ) %>% 
  left_join(
    x = annotated_receptor_ligand_pairs,
    y = ., 
    by = c('pair_id', 'timepoint')
  )

# Import gene names
gene_names <- 
  read_rds(here::here("results/download_gene_data/gene_names.Rdata")) %>% 
  select(-description)

named_receptor_ligand_pairs <- 
  right_join(
    x = gene_names,
    y = annotated_receptor_ligand_pairs_orphaned,
    by = "ensembl_gene_id"
  ) %>% 
  as_tibble()

# # Graph interactions
# # Create 4 plots: lung to placenta, lung to liver, liver to placenta, liver to liver
# # use purrr map wiht 2 lists (messengers, targets) to subset and iterate ggplot
# 
# interaction_plots <-
#   crossing(
#   messenger_tissues, target_tissues
# ) %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   as.list() %>% 
#   map(.x = ., .f = ~ as.vector(.x)) %>% 
#   map(
#   .x = .,
#   .f = ~ filter( named_receptor_ligand_pairs, tissue %in% .x, paired == "paired")
# ) %>% 
#   map(
#     .x = .,
#     .f = ~ ggplot(
#       data = .x,
#       aes(
#         x = gene_class, y = log_fc_response, col = tissue, alpha = - log10(q_value_response)
#       )
#     ) + 
#       geom_point() +
#       geom_line(aes(group = pair_id)) +
#       facet_wrap(~ timepoint) +
#       scale_y_continuous(name = "log Fold Change") +
#       scale_x_discrete(name = NULL)  +
#       scale_alpha_continuous(name = "-log10 qvalue", breaks = c(2,5,10,15))
#   )
# 
# interaction_plots[[1]] + interaction_plots[[2]] + interaction_plots[[3]] + interaction_plots[[4]]
# 
# ggsave(filename = here::here("results/receptor_ligand_filter/receptor_plot.pdf"), 
#        width = 297, height = 210, units = "mm")

# Add labels to top ligands/receptors (gene symbols)


# Graph baseline expression of receptors on x and fold change of ligand on y

wide_receptor_ligand_data <-
  receptor_ligand_pairs %>% 
  # add ligand data
  left_join(
    x = ., y = ligand_receptor_annotation,
    by = c("ligand" = "ensembl_gene_id")
  ) %>% 
  rename(
    log_fc_baseline_ligand = log_fc_baseline,
    log_fc_response_ligand = log_fc_response,
    q_value_baseline_ligand = q_value_baseline,
    q_value_response_ligand = q_value_response,
    tissue_ligand = tissue
  ) %>% 
  # add receptor data
  left_join(
    x = ., y = ligand_receptor_annotation,
    by = c("receptor" = "ensembl_gene_id", "timepoint" = "timepoint")
  ) %>% 
  rename(
    log_fc_baseline_receptor = log_fc_baseline,
    log_fc_response_receptor = log_fc_response,
    q_value_baseline_receptor = q_value_baseline,
    q_value_response_receptor = q_value_response,
    tissue_receptor = tissue
  ) %>% 
  # filter only receptor data from target tissues and ligand data from messenger tissues
  filter(
    tissue_ligand %in% c(messenger_tissues, "placentas")
  ) %>% 
  filter(
    tissue_receptor %in% target_tissues
  ) %>% 
  # add gene names to ligands
  right_join(
    x = gene_names, y = .,
    by = c( "ensembl_gene_id" = "receptor")
  ) %>% 
  # add gene names to receptors
  right_join(
    x = gene_names, y = .,
    by = c( "ensembl_gene_id" = "ligand"), 
    suffix = c("_ligand", "_receptor")
  ) %>% 
  as_tibble() %>% 
  rename(
    "ensembl_gene_id_ligand" = "ensembl_gene_id"
  ) %>% 
  mutate(
    pair_id = paste(mgi_symbol_ligand, mgi_symbol_receptor, sep = ":")
  )
  
## Save formatted data

write_csv(x = wide_receptor_ligand_data, path = here::here("results/receptor_ligand_filter/receptor_ligand_wide_table.csv"))

wide_receptor_ligand_data %>% 
  ggplot(
    aes(x = log_fc_baseline_receptor, y = log_fc_response_ligand, col = q_value_response_ligand < 0.05)
  ) +
  geom_point(alpha = 0.5) +
  facet_grid( timepoint ~ tissue_ligand + tissue_receptor) +
  ggrepel::geom_text_repel(
    data = filter(
      wide_receptor_ligand_data, 
      abs(log_fc_response_ligand) > 2, 
      log_fc_baseline_receptor > 2
    ),
    aes(label = pair_id)
  ) +
  scale_color_brewer(name = "Significant Ligand DE", type = 'qual', palette = 2) +
  scale_x_continuous(name = "log normalized receptor counts in target tissue") +
  scale_y_continuous(name = "log fold change of ligands in messenger tissue") +
  theme_bw() +
  ggtitle(label = "Absolute expression of receptors compared with\nfold change in expression of their ligands")

ggsave(filename = here::here("results/receptor_ligand_filter/receptor_base_ligand_fold_plot_full.pdf"), 
       units = "mm", width = a4[1]*2, height = a4[2]*2)

wide_receptor_ligand_data %>% 
  filter(tissue_ligand == "maternal_lung") %>% 
  ggplot(
    aes(x = log_fc_baseline_receptor, y = log_fc_response_ligand, col = q_value_response_ligand < 0.05)
  ) +
  geom_point(alpha = 0.5) +
  facet_grid( timepoint ~ tissue_receptor) +
  ggrepel::geom_text_repel(
    data = 
      filter(
        wide_receptor_ligand_data, 
        tissue_ligand == "maternal_lung", 
        abs(log_fc_response_ligand) > 2, 
        log_fc_baseline_receptor > 2
      ),
    aes(label = pair_id)
  ) +
  scale_color_brewer(name = "Significant Ligand DE", type = 'qual', palette = 2) +
  scale_x_continuous(name = "log normalized receptor counts in target tissue") +
  scale_y_continuous(name = "log fold change of ligands in messenger tissue") +
  coord_cartesian(xlim = c(0,10)) +
  theme_bw() +
  ggtitle(label = "Absolute expression of receptors compared with fold change in expression of their ligands in maternal lung")

ggsave(filename = here::here("results/receptor_ligand_filter/receptor_base_ligand_fold_plot_lung.pdf"), 
       units = "mm", width = a4[1], height = a4[2])

## Same plot, but track IL6:OSm through all timepoints

wide_receptor_ligand_data %>% 
  filter(
    tissue_ligand == "maternal_lung", 
    grepl("Il6", mgi_symbol_ligand) | grepl("Osm", mgi_symbol_receptor) |  grepl("Lifr", mgi_symbol_receptor) | grepl("Il6", mgi_symbol_receptor)
  ) %>% 
  ggplot(
    aes(x = log_fc_baseline_receptor, y = log_fc_response_ligand, col = q_value_response_ligand < 0.05)
  ) +
  geom_point(alpha = 0.5) +
  facet_grid( tissue_receptor ~ timepoint ) +
  ggrepel::geom_text_repel(
    aes(label = pair_id)
  ) +
  scale_color_brewer(name = "Significant Ligand DE", type = 'qual', palette = 2) +
  scale_x_continuous(name = "log normalized receptor counts in target tissue") +
  scale_y_continuous(name = "log fold change of ligands in messenger tissue") +
  # coord_cartesian(xlim = c(0,10)) +
  theme_bw() +
  ggtitle(label = "Absolute expression of receptors compared with fold change in expression of their ligands in maternal lung")

ggsave(filename = here::here("results/receptor_ligand_filter/IL6_osm_plot_lung.pdf"), 
       units = "mm", width = a4[1], height = a4[2])

## Check if receptors have a pattern over time when exposed to LPS
wide_receptor_ligand_data %>% 
  select(
    timepoint, tissue_receptor, 
    ensembl_gene_id_receptor, mgi_symbol_receptor,
    log_fc_baseline_receptor, log_fc_response_receptor, 
    q_value_response_receptor
  ) %>% 
  mutate(
    lps_expression = log_fc_baseline_receptor + log_fc_response_receptor
  ) %>% 
  unique() %>% 
  filter(
    grepl("Osm", mgi_symbol_receptor) |  grepl("Lifr", mgi_symbol_receptor) | grepl("Il6", mgi_symbol_receptor)
  ) %>% 
  ggplot(
    aes( x = timepoint, y = lps_expression, 
         group = mgi_symbol_receptor, 
         label = mgi_symbol_receptor, 
         col = q_value_response_receptor < 0.1)
  ) +
  geom_line() +
  facet_wrap(~ tissue_receptor) +
  ggrepel::geom_text_repel() +
  scale_color_brewer(name = "Significant (q<0.1)\nreceptor\nresponse to LPS", type = 'qual', palette = 2) +
  theme_bw()
  
ggsave(filename = here::here("results/receptor_ligand_filter/IL6_Osm_receptor_responses.pdf"), 
       units = "mm", width = a4[1], height = a4[2])

## Check if receptors have a pattern over time when exposed to LPS
wide_receptor_ligand_data %>% 
  select(
    timepoint, tissue_receptor, 
    ensembl_gene_id_receptor, mgi_symbol_receptor,
    log_fc_baseline_receptor, log_fc_response_receptor, 
    q_value_response_receptor
  ) %>% 
  mutate(
    lps_expression = log_fc_baseline_receptor + log_fc_response_receptor
  ) %>% 
  unique() %>% 
  filter(
    grepl("Osm", mgi_symbol_receptor) |  grepl("Lifr", mgi_symbol_receptor) | grepl("Il6", mgi_symbol_receptor)
  ) %>% 
  ggplot(
    aes( x = timepoint, y = log_fc_response_receptor, 
         group = mgi_symbol_receptor, 
         label = mgi_symbol_receptor, 
         col = q_value_response_receptor < 0.1)
  ) +
  geom_line() +
  facet_wrap(~ tissue_receptor) +
  ggrepel::geom_text_repel() +
  scale_color_brewer(name = "Significant (q<0.1)\nreceptor\nresponse to LPS", type = 'qual', palette = 2) +
  theme_bw()

# Create table for manual selection of interesting ligands
interesting_ligand_table <-
  wide_receptor_ligand_data %>% 
  filter(
    log_fc_baseline_receptor > 3,
    log_fc_baseline_ligand > 0.5,
    log_fc_response_ligand > 0.5,
    q_value_response_ligand < q_value_threshold,
    timepoint %in% paste( "timepoint", c(2,5), sep = "")
  ) %>% 
  dplyr::select(
    starts_with("mgi"), tissue_ligand, tissue_receptor, timepoint, log_fc_baseline_ligand, log_fc_response_ligand, starts_with("ensembl")
  ) %>% 
  group_by(mgi_symbol_ligand, tissue_ligand) %>% 
  summarise(
    mgi_symbol_receptors = unique(mgi_symbol_receptor) %>%  paste(., collapse = ":"),
    timepoint = timepoint,
    log_fc_baseline_ligand = mean(log_fc_baseline_ligand),
    log_fc_response_ligand = log_fc_response_ligand,
    ensembl_gene_id_ligand = ensembl_gene_id_ligand
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = timepoint, 
    values_from = log_fc_response_ligand
  ) %>% 
  select(
    mgi_symbol = mgi_symbol_ligand, tissue_ligand, 
    log_baseline = log_fc_baseline_ligand, log_fc_2_hrs = timepoint2, log_fc_5_hrs = timepoint5,
    mgi_symbol_receptors, ensembl_gene_id_ligand
  )

write_csv(
  x = interesting_ligand_table, 
  file = here::here("results/receptor_ligand_filter/interesting_ligands.csv")
)

# Final format: one row per ligand. Must contain expression in mat lung/liver at 2&5 hrs, plus base counts
# Include only ligands with at least 1 receptor with log base counts > 3
# Include only ligands with significant expression in at least one timepoint and logFC > 1, and base counts > 3


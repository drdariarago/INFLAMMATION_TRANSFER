### Filter expressed receptors in target tissues and rank DE ligands in messenger tissues

#### Load libraries and data ####

library(tidyverse)
library(magrittr)
library(ggplot2)

receptor_ligand_pairs <-
  read_tsv(file = "results/match_orthologs/human_mouse_ligands_receptors.txt") %>% 
  select(matches("ensembl.*mouse")) %>% 
  set_names(c("ligand", "receptor"))

fold_change_summary <-
  c("fetal_liver", "maternal_liver", "maternal_lung", "placentas") %>% 
  setNames(., .) %>% 
  map(~ paste("results/limma_", .x, "/fold_change_summary.rds", sep = "")) %>% 
  map_dfr(read_rds, .id = "tissue") %>% 
  mutate(
    tissue = as.factor(tissue)
  )

# Define which tissues we scan for ligands (messenger tissues) and receptors (target tissues)
# Also define other global criteria

messenger_tissues <- c("maternal_lung", "maternal_liver")
target_tissues <- c("placentas", "fetal_liver")
timepoints <- paste("timepoint", c(2,5), sep = "")
log_fc_threshold <- 0.5
q_value_threshold <-0.05

# Subset only receptors that are expressed in the target tissues and the ligands that target them

expressed_receptors <-
  fold_change_summary %>% 
  filter(tissue %in% target_tissues) %>% 
  select(expressed_receptor = ensembl_gene_id) %>% 
  distinct() %>% 
  inner_join(x = ., y = receptor_ligand_pairs, by = c("expressed_receptor" = "receptor"))

# Subset ligands differentially expressed in messenger tissues at 2 and 5 hours (logFC > 0.5, q-val < 0.05)
# Remove ligands whose receptors are not expressed in target tissues

expressed_receptor_ligand_pairs <-
  fold_change_summary %>% 
  filter( tissue %in% messenger_tissues) %>% 
  filter( timepoint %in% timepoints ) %>% 
  filter( exposure == "response") %>%   
  filter( q_value < q_value_threshold, abs(logFC) > log_fc_threshold ) %>% 
  select( differentially_expressed_ligand = ensembl_gene_id ) %>% 
  inner_join( x = ., y = expressed_receptors, by = c("differentially_expressed_ligand" = "ligand")) %>% 
  unique()

# Annotated subset of DE ligand - expressed receptor pairs with ligand expression data

annotated_curated_ligands <-
  fold_change_summary %>% 
  filter(
    tissue %in% messenger_tissues, 
    timepoint %in% timepoints
  ) %>% 
  select( ensembl_gene_id, tissue, timepoint, exposure, logFC ) %>% 
  left_join( 
    x = expressed_receptor_ligand_pairs, y = ., 
    by = c("differentially_expressed_ligand" = "ensembl_gene_id") 
  ) %>% 
  select( - expressed_receptor) %>%
  distinct() %>% 
  pivot_wider( names_from = exposure, values_from = logFC) %>% 
  filter(abs(response) > log_fc_threshold) %>% 
  group_by( tissue, timepoint) %>% 
  mutate(
    ligand_rank = rank(- response),
    ligand_fc = response
  )

# Rank ligands at time 2 by logFC
# Rank receptors by mean logFC of their ligands
# Annotate average expression of receptors and ligands
# Group ligands by receptor
# plot

ranked_receptor_ligand_pairs <-
  annotated_curated_ligands %>% 
  inner_join(x = expressed_receptor_ligand_pairs, y = ., by = "differentially_expressed_ligand") %>% 
  as_tibble() %>% 
  select( - baseline, - response) %>% 
  group_by(tissue, timepoint, expressed_receptor) %>% 
  mutate(
    receptor_rank = mean(ligand_rank),
    receptor_fc = mean(ligand_fc),
    n_ligands = n_distinct(differentially_expressed_ligand)
    ) %>% 
  arrange(timepoint, tissue, - receptor_fc) %>% 
  ungroup()

# Import gene names
gene_names <- 
  read_rds(here::here("results/download_gene_data/gene_names.Rdata")) %>% 
  select(-description)

named_receptor_ligand_pairs <- 
  ranked_receptor_ligand_pairs %>% 
  right_join(
    x = gene_names, y = ., by = c("ensembl_gene_id" = "differentially_expressed_ligand")
  ) %>% 
  rename(ligand_ensembl_id = ensembl_gene_id, ligand_name = mgi_symbol) %>% 
  right_join(
    x = gene_names, y = ., by = c("ensembl_gene_id" = "expressed_receptor")
  ) %>% 
  rename(receptor_ensembl_id = ensembl_gene_id, receptor_name = mgi_symbol) %>% 
  select(
    timepoint, tissue, 
    name_receptor = receptor_name, rank_receptor = receptor_rank, fc_receptor = receptor_fc, 
    n_ligands, name_ligand = ligand_name, rank_ligand = ligand_rank, fc_ligand = ligand_fc
  )
  
# Create human readable output table
named_receptor_ligand_pairs %>% 
  group_by(name_receptor) %>% 
  summarise(
    avg_receptor_rank = mean(rank_receptor)
  ) %>% 
  ungroup() %>% 
  transmute(
    name_receptor = name_receptor,
    receptor_avg_rank = rank(avg_receptor_rank)
  ) %>% 
  left_join(x = ., y = named_receptor_ligand_pairs, by = "name_receptor") %>% 
  arrange(receptor_avg_rank, name_ligand, timepoint, tissue) %>% 
  write_csv(here::here("results/receptor_ligand_filter/receptor_ligand_fc_table.csv"))

# Graph interactions

long_paired_data <-
  named_receptor_ligand_pairs %>% 
  select(- n_ligands) %>% 
  mutate(
    pairs = paste(name_ligand, name_receptor, sep = "_")
  ) %>% 
  pivot_longer(
    cols = contains("_"),
    names_to = c(".value", "ligand"),
    names_sep = "_",
    values_drop_na = FALSE
  ) %>% 
  mutate(
    ligand = factor(x = ligand, levels = c("ligand", "receptor"))
  )
  
long_paired_data %>% 
  ggplot(
    aes(x = ligand, y = fc, group = pairs)
  ) + 
  geom_point(alpha = 0.1) +
  geom_line(alpha = 0.1) +
  facet_grid(tissue ~ timepoint) +
  scale_y_continuous(name = "log Fold Change of Ligands") +
  scale_x_discrete(name = NULL) +
  ggtitle(label = "log Fold Change of ligands in messenger tissues and average of the ligands\nthat bind to each receptor in target tissues")

ggsave(filename = here::here("results/receptor_ligand_filter/receptor_plot.pdf"))

# Add labels to top ligands/receptors (gene symbols)
# Annotate if receptors are placenta or liver

# MA plot of curated ligands by tissue and time
annotated_curated_ligands %>% 
  ggplot( aes( x = baseline, y = response) )+
  geom_point() +
  facet_grid( tissue ~ timepoint )

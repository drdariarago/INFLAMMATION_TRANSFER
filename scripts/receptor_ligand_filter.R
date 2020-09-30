### Filter expressed receptors in target tissues and rank DE ligands in messenger tissues

#### Load libraries and data ####

library(tidyverse)
library(magrittr)
library(ggplot2)
library(patchwork)

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
timepoints <- paste("timepoint", c(2,5), sep = "")
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
  )


# Graph interactions
# Create 4 plots: lung to placenta, lung to liver, liver to placenta, liver to liver
# use purrr map wiht 2 lists (messengers, targets) to subset and iterate ggplot

interaction_plots <-
  crossing(
  messenger_tissues, target_tissues
) %>% 
  t() %>% 
  as.data.frame() %>% 
  as.list() %>% 
  map(.x = ., .f = ~ as.vector(.x)) %>% 
  map(
  .x = .,
  .f = ~ filter( named_receptor_ligand_pairs, tissue %in% .x, paired == "paired")
) %>% 
  map(
    .x = .,
    .f = ~ ggplot(
      data = .x,
      aes(
        x = gene_class, y = log_fc_response, col = tissue, alpha = - log10(q_value_response)
      )
    ) + 
      geom_point() +
      geom_line(aes(group = pair_id)) +
      facet_wrap(~ timepoint) +
      scale_y_continuous(name = "log Fold Change") +
      scale_x_discrete(name = NULL)  +
      scale_alpha_continuous(name = "-log10 qvalue", breaks = c(2,5,10,15))
  )

interaction_plots[[1]] + interaction_plots[[2]] + interaction_plots[[3]] + interaction_plots[[4]]

ggsave(filename = here::here("results/receptor_ligand_filter/receptor_plot.pdf"), 
       width = 297, height = 210, units = "mm")

# Add labels to top ligands/receptors (gene symbols)

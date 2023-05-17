#### Libraries ####
library(tidyverse)
library(broom)
library(multidplyr)
cluster <- new_cluster(parallel::detectCores() - 1) 
cluster_library(cluster, c('ape', 'tidytree', 'tibble', 'dplyr', 'stringr'))


#### Data ####
cafe_trees <- read_lines('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/errorModel/Base_asr.tre') %>%
  tibble(tree = .) %>%
  mutate(FamilyID = str_extract(tree, 'OG[0-9]+'),
         tree = str_remove(tree, str_c(' *TREE ', FamilyID, ' = '))) %>%
  filter(!is.na(tree)) %>%
  rowwise(FamilyID) %>%
  partition(cluster) %>%
  summarise(read.tree(text = tree) %>%
              as_tibble() %>%
              mutate(cafe_nodeID = str_extract(label, '(([a-z_])+)?\\<[0-9]+\\>'),
                     label = str_replace(label, '\\<[0-9]+\\>', str_c('<', node, '>'))) %>%
              select(label, cafe_nodeID)) %>%
  collect %>%
  ungroup %>%
  mutate(is_significant = str_detect(label, '\\*'),
         ortho_nodeID = str_extract(label, '\\<[0-9]+\\>') %>% str_remove_all(c('\\<|\\>')) %>% as.integer,
         n_gene = str_extract(label, '[0-9]+$') %>% as.integer()) %>%
  rename(Orthogroup = FamilyID) %>%
  mutate(Orthogroup = str_c('N0.H', Orthogroup))

expansion_contractions <- read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/errorModel/Base_change.tab', 
                                     delim = '\t', show_col_types = FALSE) %>%
  rename(Orthogroup = FamilyID) %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta')

kegg_paths <- read_csv('../intermediate_files/kegg_orthogroup_pathways.csv.gz',
                       show_col_types = FALSE) %>%
  filter(major_category != 'Human Diseases') %>%
  select(-gene_id, -species) %>%
  distinct %>%
  nest(kegg_paths = c(kegg_path_id, name:rel_pathway)) %>%
  select(-kegg_type)



#### Join data together ####
full_join(cafe_trees,
          expansion_contractions,
          by = c('Orthogroup', 'cafe_nodeID' = 'node')) %>%
  left_join(kegg_paths,
            by = 'Orthogroup') %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction')) %>%
  filter(is_significant) %>%
  select(Orthogroup, cafe_nodeID, ortho_nodeID, n_gene, delta, change, kegg_gene, kegg_paths) %>%
  unnest(kegg_paths, keep_empty = TRUE) %>%
  group_by(cafe_nodeID, change, major_category, minor_category, name, description) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            sum_delta = sum(delta),
            .groups = 'drop') %>%
  
  filter(str_detect(cafe_nodeID, 'acer')) %>% 
  select(-n_ortho) %>%
  pivot_wider(names_from = change,
              values_from = sum_delta,
              values_fill = 0) %>%
  
  filter(str_detect(minor_category, 'Immune'))



full_join(cafe_trees,
          expansion_contractions,
          by = c('Orthogroup', 'cafe_nodeID' = 'node')) %>%
  left_join(kegg_paths,
            by = 'Orthogroup') %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction')) %>%
  filter(is_significant) %>%
  select(Orthogroup, cafe_nodeID, ortho_nodeID, n_gene, delta, change, kegg_gene, kegg_paths) %>%
  filter(str_detect(cafe_nodeID, 'acer')) %>%
  filter(!is.na(kegg_gene)) %>%
  filter(map_lgl(kegg_paths, ~any(str_detect(.x$minor_category, 'Immune')))) %>% 
  arrange(kegg_gene) %>% View
  
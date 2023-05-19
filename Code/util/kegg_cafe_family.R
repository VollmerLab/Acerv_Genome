library(tidyverse)
maxDiff <- 200

kegg_paths <- read_csv('../../intermediate_files/kegg_orthogroup_pathways.csv.gz',
                       show_col_types = FALSE) %>%
  filter(major_category != 'Human Diseases') %>%
  select(-gene_id, -species) %>%
  distinct %>%
  nest(kegg_paths = c(kegg_path_id, name:rel_pathway)) %>%
  select(-kegg_type)

#Unique KEGGS per pathway
for_cafe <- read_delim('../../../Bioinformatics/Phylogenomics/Phylogenetic_Hierarchical_Orthogroups/N0.tsv', show_col_types = FALSE) %>%
  rename(Orthogroup = HOG) %>%
  select(-OG, -`Gene Tree Parent Clade`) %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'species',
               values_to = 'gene_list') %>%
  filter(!is.na(gene_list)) %>%
  rowwise(Orthogroup, species) %>%
  summarise(gene_id = str_split(gene_list, ', ') %>%
              unlist, 
            .groups = 'drop') %>%
  filter(!is.na(gene_id)) %>%
  group_by(Orthogroup, species) %>%
  summarise(n = n(), .groups = 'drop_last')  %>%
  ungroup %>%
  left_join(kegg_paths, 
            by = 'Orthogroup') %>%
  
  
  unnest(kegg_paths) %>%
  select(-kegg_gene, -name:-rel_pathway) %>%
  group_by(kegg_path_id, kegg_orthology, species) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            n_gene = sum(n),
            .groups = 'drop') %>%
  group_by(kegg_path_id, species) %>%
  summarise(n = n_distinct(kegg_orthology),
            .groups = 'drop_last') %>%

# %>%
#   select(-kegg_gene:-prop_top_hits) %>%
#   unnest(kegg_paths) %>%
#   select(kegg_path_id, species, n) %>%
#   group_by(kegg_path_id, species) %>%
#   summarise(n = sum(n),
#             .groups = 'drop_last') %>%
  mutate(max_diff = max(n) - min(n),
         n_species = n_distinct(species)) %>%
  pivot_wider(names_from = species,
              values_from = n,
              values_fill = 0L) %>%
  mutate(Desc = '(null)', .before = everything()) %>%
  mutate(kegg_path_id = str_remove(kegg_path_id, 'path:')) %>%
  rename('Family ID' = kegg_path_id) %>%
  ungroup %>%
  # arrange(-max_diff)
  
  filter(n_species > 1,
         max_diff <= maxDiff) %>%
  select(-max_diff, -n_species)

write_delim(for_cafe, delim = '\t', file = '../../intermediate_files/kegg_families.txt')


#### Gene Copies per KEGG pathway
for_cafe <- read_delim('../../../Bioinformatics/Phylogenomics/Phylogenetic_Hierarchical_Orthogroups/N0.tsv', show_col_types = FALSE) %>%
  rename(Orthogroup = HOG) %>%
  select(-OG, -`Gene Tree Parent Clade`) %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'species',
               values_to = 'gene_list') %>%
  filter(!is.na(gene_list)) %>%
  rowwise(Orthogroup, species) %>%
  summarise(gene_id = str_split(gene_list, ', ') %>%
              unlist, 
            .groups = 'drop') %>%
  filter(!is.na(gene_id)) %>%
  group_by(Orthogroup, species) %>%
  summarise(n = n(), .groups = 'drop_last')  %>%
  ungroup %>%
  left_join(kegg_paths, 
            by = 'Orthogroup') %>%
  
  
  unnest(kegg_paths) %>%
  select(-kegg_gene, -name:-rel_pathway) %>%
  group_by(kegg_path_id, kegg_orthology, species) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            n_gene = sum(n),
            .groups = 'drop') %>%
  
  group_by(kegg_path_id, species) %>%
  summarise(n = sum(n_gene),
            .groups = 'drop_last') %>%
  mutate(max_diff = max(n) - min(n),
         n_species = n_distinct(species)) %>%
  pivot_wider(names_from = species,
              values_from = n,
              values_fill = 0L) %>%
  mutate(Desc = '(null)', .before = everything()) %>%
  mutate(kegg_path_id = str_remove(kegg_path_id, 'path:')) %>%
  rename('Family ID' = kegg_path_id) %>%
  ungroup %>%
  # arrange(-max_diff)
  
  filter(n_species > 1,
         max_diff <= maxDiff) %>%
  select(-max_diff, -n_species)

write_delim(for_cafe, delim = '\t', file = '../../intermediate_files/kegg_families_gene_level.txt')


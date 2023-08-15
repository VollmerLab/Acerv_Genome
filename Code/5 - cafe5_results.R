#### Libraries ####
library(tidyverse)
library(magrittr)
library(broom)

#### Functions ####
row_fisher <- function(n_sig, n_genome, total_sig, total_genome, direction = 'two.sided'){
  #http://yulab-smu.top/clusterProfiler-book/chapter2.html#over-representation-analysis
  data.frame(significant = c(n_sig, total_sig - n_sig), not_significant = c(n_genome - n_sig, total_genome - total_sig - (n_genome - n_sig))) %>%
    fisher.test(alternative = direction) %>%
    tidy %>%
    select(-alternative, -method)
}

#### Data ####
cafe5_results <- full_join(read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/errorModel/Base_count.tab', delim = '\t',
                     show_col_types = FALSE) %>%
            pivot_longer(cols = -FamilyID,
                         names_to = 'tree_node',
                         values_to = 'n_gene'),
          
          read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/errorModel/Base_change.tab', delim = '\t',
                     show_col_types = FALSE) %>%
            pivot_longer(cols = -FamilyID,
                         names_to = 'tree_node',
                         values_to = 'gene_change'),
          by = c('FamilyID', 'tree_node')) %>%
  
  full_join(read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/errorModel/Base_branch_probabilities.tab', delim = '\t',
                       show_col_types = FALSE) %>%
              select(-...49, -`<46>`) %>%
              rename(FamilyID = `#FamilyID`) %>%
              pivot_longer(cols = -FamilyID,
                           names_to = 'tree_node',
                           values_to = 'node_probability'),
            by = c('FamilyID', 'tree_node')) %>%
  left_join(read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/errorModel/Base_family_results.txt', delim = '\t',
                       show_col_types = FALSE) %>%
              rename(FamilyID = `#FamilyID`,
                     family_probability = pvalue) %>%
              select(-`Significant at 0.05`),
            by = 'FamilyID') %>%
  
  rename(Orthogroup = FamilyID) %>%
  mutate(species = str_extract(tree_node, '[a-z_]+')) %>%
  filter(!is.na(species)) %>%
  select(-tree_node) %>%
  select(Orthogroup, species, n_gene, gene_change, family_probability, node_probability)


orthogroup_kegg <- read_csv('../intermediate_files/kegg_gene_pathways.csv.gz', 
                            show_col_types = FALSE, guess_max = 1e5)

#### What expanded/contracted in Acer ####
acer_changes <- full_join(cafe5_results, 
          orthogroup_kegg,
          by = c('Orthogroup')) %>%
  filter(species == 'acer',
         !is.na(kegg_gene),
         family_probability < 0.05,
         node_probability < 0.05) %>%
  janitor::remove_empty(which = 'cols') %T>%
  write_csv('../Results/acer_orthogroupChange.csv')


all_acer_paths <- full_join(cafe5_results, 
          orthogroup_kegg,
          by = c('Orthogroup')) %>%
  filter(species == 'acer',
         !is.na(kegg_gene)) %>%
  pivot_longer(cols = -c(Orthogroup:prop_top_hits),
               names_to = c('major', 'minor'),
               names_pattern = '(.*) - (.*)',
               values_to = 'pathway', 
               values_drop_na = TRUE) %>%
  rowwise %>%
  mutate(pathway = str_split(pathway, '; ')) %>%
  unnest(pathway) %>%
  group_by(major, minor) %>%
  summarise(total_paths = n_distinct(pathway),
            .groups = 'drop')
  

acer_expand_contract <- acer_changes %>%
  pivot_longer(cols = -c(Orthogroup:prop_top_hits),
               names_to = c('major', 'minor'),
               names_pattern = '(.*) - (.*)',
               values_to = 'pathway', 
               values_drop_na = TRUE) %>%
  rowwise %>%
  mutate(pathway = str_split(pathway, '; ')) %>%
  unnest(pathway) %>%
  group_by(major, minor) %>%
  summarise(expansion = n_distinct(pathway[gene_change > 0]),
            contraction = n_distinct(pathway[gene_change < 0]),
            .groups = 'drop') 

pathway_overrep <- full_join(acer_expand_contract, 
          all_acer_paths,
          by = c('major', 'minor')) %>%
  filter(total_paths > 5) %>%
  
  mutate(across(c(expansion, contraction), ~replace_na(., 0L))) %>%
  mutate(total_expansion = sum(expansion),
         total_contraction = sum(contraction),
         total_total = sum(total_paths)) %>%
  rowwise(major, minor) %>%
  mutate(expansion_test = row_fisher(expansion, total_paths, total_expansion, total_total, 'greater'),
            contraction_test = row_fisher(contraction, total_paths, total_contraction, total_total, 'greater')) %>%
  select(-total_expansion, -total_contraction, -total_total) %>%
  ungroup

pathway_overrep %>%
  filter(expansion_test$p.value < 0.05 | contraction_test$p.value < 0.05) %T>%
  print %>%
  select(-expansion_test) 
  

gene_changes <- acer_changes %>%
  pivot_longer(cols = -c(Orthogroup:prop_top_hits),
               names_to = c('major', 'minor'),
               names_pattern = '(.*) - (.*)',
               values_to = 'pathway', 
               values_drop_na = TRUE) %>%
  rowwise %>%
  mutate(pathway = str_split(pathway, '; ')) %>%
  unnest(pathway) %>%
  inner_join(pathway_overrep %>%
               filter(expansion_test$p.value < 0.05 | contraction_test$p.value < 0.05) %>%
               select(minor),
             by = 'minor') %>%
  mutate(delta = if_else(gene_change > 0, 'expand', 'contract')) %>%
  select(major, minor, pathway, kegg_gene, kegg_orthology, delta) %>%
  distinct %>%
  mutate(major_minor = str_c(major, minor, sep = ' - '),
         .keep = 'unused') %>%
  pivot_wider(names_from = major_minor, values_from = pathway, values_fn = ~str_c(., collapse = ';;')) %>% 
  arrange(kegg_gene) 
write_csv('../Results/expansion_contractions_ora.csv')
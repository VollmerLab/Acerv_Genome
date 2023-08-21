##Rerun cafe based on orthofinder results 
#Rows are KEGG orthologs with numbers filled in being the number of HOGs
# - rework orthofinder step (8 in markdown)
# - make spawn cafe5 using HOGS and KEGGS
# - need to rerun KEGG/Entap

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
  
  rename(kegg_orthology = FamilyID) %>%
  mutate(species = str_extract(tree_node, '[a-z_]+')) %>%
  filter(!is.na(species)) %>%
  select(-tree_node) %>%
  select(kegg_orthology, species, n_gene, gene_change, family_probability, node_probability)

orthogroup_kegg <- read_csv('../../Bioinformatics/Phylogenomics/Annotations/species_kegg_orthogroup.csv.gz', 
                         show_col_types = FALSE, guess_max = 1e5) %>%
  select(kegg_gene, kegg_orthology, contains(' - ')) %>%
  distinct %>%
  janitor::remove_empty(which = c('rows', 'cols'))

#### What expanded/contracted in Acer ####
acer_kegg_ortho <- cafe5_results %>%
  filter(species == 'acer') %>%
  mutate(node_probability = p.adjust(node_probability, 'fdr')) %>%
  full_join(orthogroup_kegg,
            by = c('kegg_orthology')) %>%
  filter(species == 'acer',
         !is.na(kegg_gene))

# acer_kegg_ortho <- full_join(cafe5_results, 
#           orthogroup_kegg,
#           by = c('kegg_orthology')) %>%
#   filter(species == 'acer',
#          !is.na(kegg_gene)) 

acer_kegg_ortho %>%
  filter(kegg_orthology == 'ko:K00002')

acer_kegg_ortho %>%
  filter(node_probability < 0.05) %>%
  janitor::remove_empty(which = 'cols') %T>%
  write_csv('../Results/acer_orthogroupChange.csv')

acer_kegg_ortho %>%
  filter(node_probability < 0.05) %>%
  janitor::remove_empty(which = 'cols') %>%
  filter(if_all(contains(' - '), is.na))


acer_ortho <- acer_kegg_ortho %>%
  pivot_longer(cols = -c(kegg_orthology:kegg_gene),
             names_to = c('major_minor'),
             values_to = 'pathway', 
             values_drop_na = TRUE) %>%
  group_by(across(kegg_orthology:major_minor)) %>%
  reframe(pathway = str_split(pathway, ';') %>% unlist %>% str_trim) %>%
  select(kegg_gene, kegg_orthology, major_minor, pathway, n_gene, gene_change, family_probability, node_probability) %>%
  separate(major_minor, 
         sep = ' - ', 
         into = c('major', 'minor'), 
         extra = 'merge') %>%
  filter(!minor %in% c('Cellular community - prokaryotes',
                       'Information processing in viruses'))


pathway_overrep <- acer_ortho %>%
  # filter(minor == 'Nervous system') %>%
  
  group_by(major, minor, pathway) %>%
  summarise(n_expand = n_distinct(kegg_orthology[node_probability < 0.05 & gene_change > 0]),
            n_contract = n_distinct(kegg_orthology[node_probability < 0.05 & gene_change < 0]),
            path_ortho = n_distinct(kegg_orthology),
            .groups = 'rowwise') %>%
  bind_cols(summarise(acer_ortho,
                      total_expansion = n_distinct(kegg_orthology[node_probability < 0.05 & gene_change > 0]),
                      total_contraction = n_distinct(kegg_orthology[node_probability < 0.05 & gene_change < 0]),
                      total_ortho = n_distinct(kegg_orthology))) %>%
  # left_join(group_by(acer_ortho, major, minor) %>%
  #             summarise(total_expansion = n_distinct(Orthogroup[node_probability < 0.05 & gene_change > 0]),
  #                       total_contraction = n_distinct(Orthogroup[node_probability < 0.05 & gene_change < 0]),
  #                       total_ortho = n_distinct(Orthogroup),
  #                       .groups = 'drop'),
  #           by = c('major', 'minor')) %>%
  mutate(expansion_test = row_fisher(n_expand, path_ortho, total_expansion, total_ortho, 'greater') %>%
           rename_with(~str_c('expansion', ., sep = '_')),
         contraction_test = row_fisher(n_contract, path_ortho, total_contraction, total_ortho, 'greater') %>%
           rename_with(~str_c('contraction', ., sep = '_'))) %>%
  ungroup %>%
  unnest(c(expansion_test, contraction_test)) %>%
  # group_by(major, minor) %>%
  mutate(expansion_fdr = p.adjust(expansion_p.value, method = 'fdr'),
         contraction_fdr = p.adjust(contraction_p.value, method = 'fdr')) %>%
  ungroup


#https://en.wikipedia.org/wiki/Fisher%27s_method
minor_cat_significance <- pathway_overrep %>%
  # filter(minor == 'Nervous system') %>%
  group_by(major, minor) %>%
  summarise(df = 2 * n(),
            expansion_x2 = -2 * sum(log(expansion_p.value)),
            contraction_x2 = -2 * sum(log(contraction_p.value)),
            .groups = 'drop') %>%
  mutate(expansion_p.value = pchisq(expansion_x2, df, lower.tail = FALSE),
         contraction_p.value = pchisq(contraction_x2, df, lower.tail = FALSE)) %>%
  mutate(expansion_fdr = p.adjust(expansion_p.value, method = 'fdr'),
         contraction_fdr = p.adjust(contraction_p.value, method = 'fdr')) %>%
  filter(expansion_fdr < 0.05 | contraction_fdr < 0.05)
minor_cat_significance
write_csv(minor_cat_significance, '../Results/expand_contract_minor_pathways.csv')

minor_cat_significance %>%
  pivot_longer(cols = c(starts_with('expansion'), starts_with('contraction')),
               names_to = c('direction', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  filter(fdr < 0.05) %>%
  arrange(direction, major) %>%
  mutate(across(c(p.value, fdr), ~scales::pvalue(.))) %>%
  filter(direction != 'expansion') %>%
  select(-major) 


overrep_pathways <- pathway_overrep %>%
  inner_join(minor_cat_significance %>%
               select(major, minor, contains('fdr')) %>%
               mutate(across(contains('fdr'), ~. < 0.01)) %>%
               pivot_longer(cols = contains('fdr'),
                            names_to = 'overrep_direction',
                            names_pattern = '(.*)_fdr',
                            values_to = 'significance') %>%
               filter(significance) %>%
               select(-significance),
             by = c('major', 'minor'),
             relationship = "many-to-many") %>%
  # filter(minor == 'Nervous system')
  filter((contraction_p.value < 0.05 & overrep_direction == 'contraction') | 
           (expansion_p.value < 0.05 & overrep_direction == 'expansion'))
write_csv(overrep_pathways, '../Results/expand_contract_pathways.csv')


overrep_pathways %>%
  select(-total_expansion:-total_ortho) %>%
  pivot_longer(cols = c(starts_with('expansion'), starts_with('contraction')),
               names_to = c('direction', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  arrange(direction, major) %>%
  filter(direction == overrep_direction) %>%
  # filter(direction != 'expansion') %>%
  select(-major, -overrep_direction, -direction, -starts_with('conf')) %>%
  mutate(across(c(p.value, fdr), ~scales::pvalue(.))) %>%
  # slice(-1:-10) %>%
  # pull(pathway)
  select(-minor)
  # select(minor, pathway)
  identity()

kegg_change <- overrep_pathways %>%
  select(major, minor, pathway, contains('p.value')) %>%
  mutate(across(contains('p.value'), ~. < 0.05)) %>%
  pivot_longer(cols = contains('p.value'),
               names_to = 'overrep_direction',
               names_pattern = '(.*)_p.value',
               values_to = 'significance') %>%
  filter(significance) %>%
  # filter(overrep_direction == 'expansion')
  
  select(-significance) %>%
  distinct %>%
  left_join(acer_ortho,
            by = c('major', 'minor', 'pathway'),
            relationship = "many-to-many") %>%
  # filter((overrep_direction == 'contraction' & gene_change < 0) |
  #          (overrep_direction == 'expansion' & gene_change > 0)) %>%
  filter(node_probability < 0.05) %>%
  select(-overrep_direction) %>% 
  distinct %>%
  pivot_wider(names_from = c('major', 'minor'),
              values_from = pathway,
              names_sep = ' - ',
              values_fn = ~str_c(., collapse = ';;'))
write_csv(kegg_change, '../Results/expand_contract_keggs.csv')

kegg_change %>%
  mutate(kegg_gene = str_c(kegg_gene, ' (', kegg_orthology, ')'),
         .keep = 'unused') %>%
  select(-family_probability, -node_probability) %>%
  mutate(kegg_gene = str_c(kegg_gene, ': n ortho = ', n_gene, ' (delta = ', gene_change, ')'),
         .keep = 'unused') %>% 
  filter(!is.na(`Organismal Systems - Immune system`)) %>%
  filter(!str_detect(`Organismal Systems - Immune system`, 'NOD')) %>%
  pull(kegg_gene)


acer_ortho %>%
  pivot_wider(names_from = c('major', 'minor'),
              values_from = pathway,
              names_sep = ' - ',
              values_fn = ~str_c(., collapse = ';;')) %>%
  mutate(kegg_gene = str_c(kegg_gene, ' (', kegg_orthology, ')'),
         .keep = 'unused') %>%
  # mutate(node_probability = p.adjust(node_probability, 'fdr')) %>%
  filter(node_probability < 0.05) %>%
  janitor::remove_empty(which = 'cols') 


#### Over-represented in changers (up or down) ####
pathway_overrep <- acer_ortho %>%
  # filter(minor == 'Nervous system') %>%
  
  group_by(major, minor, pathway) %>%
  summarise(expand =  n_distinct(kegg_orthology[node_probability < 0.05 & gene_change > 0]),
            contract = n_distinct(kegg_orthology[node_probability < 0.05 & gene_change < 0]),
            n_change = n_distinct(kegg_orthology[node_probability < 0.05]),
            path_ortho = n_distinct(kegg_orthology),
            .groups = 'rowwise') %>%
  bind_cols(summarise(acer_ortho,
                      total_change = n_distinct(kegg_orthology[node_probability < 0.05]),
                      total_ortho = n_distinct(kegg_orthology))) %>%
  mutate(change_test = row_fisher(n_change, path_ortho, total_change, total_ortho, 'greater')) %>%
  ungroup %>%
  unnest(c(change_test)) %>%
  # group_by(major, minor) %>%
  mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
  ungroup

pathway_overrep %>%
  filter(fdr < 0.05) %>% 
  mutate(across(c(p.value, fdr), ~scales::pvalue(.))) %>%
  # pull(pathway)
  # select(major, minor) 
  select(pathway, expand, contract, path_ortho, estimate, p.value, fdr) 


kegg_change <- pathway_overrep %>%
  select(major, minor, pathway, contains('fdr')) %>%
  filter(fdr < 0.05) %>%
  select(-fdr) %>%
  distinct %>%
  left_join(acer_ortho,
            by = c('major', 'minor', 'pathway'),
            relationship = "many-to-many") %>%
  # filter((overrep_direction == 'contraction' & gene_change < 0) |
  #          (overrep_direction == 'expansion' & gene_change > 0)) %>%
  filter(node_probability < 0.05) %>%
  distinct %>%
  pivot_wider(names_from = c('major', 'minor'),
              values_from = pathway,
              names_sep = ' - ',
              values_fn = ~str_c(., collapse = ';;'))


kegg_change %>%
  mutate(kegg_gene = str_c(kegg_gene, ' (', kegg_orthology, ')'),
         .keep = 'unused') %>%
  select(-family_probability, -node_probability) %>%
  mutate(kegg_gene = str_c(kegg_gene, ': n ortho = ', n_gene, ' (delta = ', gene_change, ')'),
         .keep = 'unused') %>% View
  
  # filter(!is.na(`Cellular Processes - Cell growth and death`))
  filter(!is.na(`Environmental Information Processing - Signal transduction`))

  
  
  filter(!is.na(`Organismal Systems - Immune system`)) %>%
  filter(!str_detect(`Organismal Systems - Immune system`, 'NOD')) %>%
  pull(kegg_gene)
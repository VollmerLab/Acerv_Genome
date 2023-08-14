#### Libraries ####
library(tidyverse)

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
full_join(cafe5_results, 
          orthogroup_kegg,
          by = c('Orthogroup')) %>%
  filter(species == 'acer',
         !is.na(kegg_gene)) %>%
  write_csv('../intermediate_files/acer_orthogroupChange.csv')



species_kegg %>%
  filter(species == 'amic_shinzato')

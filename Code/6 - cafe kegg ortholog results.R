##TODO - wtf is the increase/decrease number and why doesnt it match up with the significant counts at each node...

#### Libraries ####
library(tidyverse)
library(magrittr)
library(treedataverse)
library(deeptime)
library(gt)

library(multidplyr)
cluster <- new_cluster(parallel::detectCores() - 1) 
cluster_library(cluster, c('ape', 'tidytree', 'tibble', 'dplyr', 'stringr'))

#### KEGG Pathway data ####
paths_particular_interest <- c('map04620', 'map04624', 'map00260', 'map00270', 'map01100', 'map01110', 'map01230', 'map01120')

kegg_paths <- read_csv('../intermediate_files/kegg_orthogroup_pathways.csv.gz',
                       show_col_types = FALSE) %>%
  filter(major_category != 'Human Diseases') %>%
  select(kegg_orthology, kegg_gene, kegg_path_id, name, description, minor_category, major_category, pathway_map, rel_pathway) %>%
  distinct %>%
  mutate(kegg_path_id = str_remove(kegg_path_id, 'path:'),
         kegg_orthology = str_remove(kegg_orthology, 'ko:')) %>%
  nest(path_data = c(kegg_path_id:rel_pathway))

#### Tree ####
the_tree <- read_rds('../intermediate_files/updated_tree.rds')

simple_tree <- mutate(the_tree, label = species) %>%
  ggtree()

#### Significance of KEGG orthologs ####
pathway_changes <- read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_kegg_koCopy/cafeOut/errorModel/Base_family_results.txt',
                              delim = '\t', show_col_types = FALSE) %>%
  rename(FamilyID = `#FamilyID`) %>%
  full_join(kegg_paths,
            by = c('FamilyID' = 'kegg_orthology'))

pathway_changes %>%
  unnest(path_data) %>%
  # filter(minor_category == 'Development and regeneration') %>%
  group_by(major_category, minor_category) %>%
  summarise(sig_paths = n_distinct(FamilyID[pvalue < 0.05 & !is.na(pvalue)]),
            tested_paths = n_distinct(FamilyID[!is.na(pvalue)]),
            total_paths = n_distinct(FamilyID),
            untested_paths = n_distinct(FamilyID[is.na(pvalue)]),
            .groups = 'drop') %>%
  select(-untested_paths) %>%
  mutate(value = sig_paths / tested_paths) %>%
  # pivot_longer(cols = ends_with('paths')) %>%
  # mutate(name = str_remove(name, '_paths'),
  #        name = case_when(name == 'sig' ~ 'Significant',
  #                         name == 'tested' ~ 'Tested',
  #                         name == 'total' ~ 'Total'),
  #        name = factor(name, levels = rev(c('Significant', 'Tested', 'Total')))) %>%
  # filter(tested_paths < sig_paths)
  
  ggplot(aes(x = value, y = minor_category)) +
  geom_linerange(aes(xmin = 0, xmax = value),
                 position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~major_category, scales = 'free_y') +
  scale_x_continuous(labels = scales::percent_format()) +
  guides(colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  labs(x = 'KEGG Orthologs (%)',
       y = NULL,
       colour = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        axis.text = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 14),
        legend.position = 'bottom')
ggsave('../Results/significant_change_ko.png', height = 7, width = 15, scale = 1)

filter(pathway_changes, FamilyID %in% paths_particular_interest)

major_categories <- unnest(kegg_paths, path_data) %>%
  pull(major_category) %>%
  unique %>%
  sort

pathway_changes %>%
  filter(pvalue < 0.05) %>%
  select(-`Significant at 0.05`) %>%
  unnest(path_data) %>%
  select(-pathway_map, -rel_pathway, -description) %>%
  mutate(kegg_pathway = str_c(name, ' (', kegg_path_id, ')'), .keep = 'unused') %>%
  mutate(category = str_c(major_category, ' - ', minor_category), .keep = 'unused') %>%
  mutate(gene_name = str_c(kegg_gene, ' (', FamilyID, ')'), .keep = 'unused') %>%
  pivot_wider(names_from = category, values_from = kegg_pathway,
              values_fill = NA_character_,
              values_fn = ~str_c(., collapse = '; ')) %>%
  select(gene_name, pvalue, order(colnames(.))) %>%
  arrange(gene_name) %>%
  gt(rowname_col = 'gene_name') %>%
  tab_header(title = 'KEGG Orthologs with Significant Changes in number of gene copies across Anthozoan Tree') %>%
  tab_stubhead('KEGG Ortholog') %>%
  tab_spanner(label = major_categories[1], columns = starts_with(major_categories[1])) %>%
  tab_spanner(label = major_categories[2], columns = starts_with(major_categories[2])) %>%
  tab_spanner(label = major_categories[3], columns = starts_with(major_categories[3])) %>%
  tab_spanner(label = major_categories[4], columns = starts_with(major_categories[4])) %>%
  tab_spanner(label = major_categories[5], columns = starts_with(major_categories[5])) %>%
  cols_label_with(columns = where(is.character), ~str_remove_all(., str_c(major_categories, collapse = '|'))) %>%
  cols_label_with(columns = where(is.character), ~str_remove_all(., '^ - ')) %>%
  fmt(columns = 'pvalue', fns = scales::pvalue_format()) %>%
  sub_missing(columns = where(is.character), missing_text = '-') %T>%
  gtsave('../Results/significantly_evolving_ko.html') 

#### Read in Cafe Results - KEGG Orthologs ####
cafe_trees <- read_lines('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_kegg_koCopy/cafeOut/errorModel/Base_asr.tre') %>%
  str_subset('^ *TREE') %>%
  tibble(tree = .) %>%
  mutate(FamilyID = str_extract(tree, 'K[0-9]+'),
         tree = str_remove(tree, str_c(' *TREE ', FamilyID, ' = '))) %>%
  # slice(1:10) %>%
  rowwise(FamilyID) %>%
  partition(cluster) %>%
  summarise(read.tree(text = tree) %>%
              as_tibble() %>%
              mutate(cafe_nodeID = str_extract(label, '(([a-z_])+)?\\<[0-9]+\\>'),
                     label = str_replace(label, '\\<[0-9]+\\>', str_c('<', node, '>')),
                     n_keggs = str_extract(label, '[0-9]+$') %>% as.integer) %>%
              select(label, cafe_nodeID, n_keggs)) %>%
  collect %>%
  ungroup %>%
  mutate(is_significant = str_detect(label, '\\*'),
         ortho_nodeID = str_extract(label, '\\<[0-9]+\\>') %>% str_remove_all(c('\\<|\\>')) %>% as.integer,
         n_gene = str_extract(label, '[0-9]+$') %>% as.integer()) 


gene_family_change <- read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_kegg_koCopy/cafeOut/errorModel/Base_clade_results.txt',
                                 delim = '\t', show_col_types = FALSE) %>%
  rename(cafe_nodeID = `#Taxon_ID`) %>%
  left_join(cafe_trees %>%
              select(ortho_nodeID, cafe_nodeID) %>%
              distinct,
            by = 'cafe_nodeID')


complete_family_species_cafe <- full_join(
  read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_kegg_koCopy/cafeOut/errorModel/Base_change.tab', 
             delim = '\t', show_col_types = FALSE) %>%
    pivot_longer(cols = -FamilyID,
                 names_to = 'cafe_nodeID',
                 values_to = 'gene_change'),
  
  read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_kegg_koCopy/cafeOut/errorModel/Base_branch_probabilities.tab', 
             na = 'N/A', delim = '\t', show_col_types = FALSE) %>%
    select(-...49) %>%
    rename(FamilyID = `#FamilyID`) %>%
    pivot_longer(cols = -FamilyID,
                 names_to = 'cafe_nodeID',
                 values_to = 'p'),
  
  by = c("FamilyID", "cafe_nodeID")
) %>%
  rename(node_p = p) %>%
  # filter(p < 0.05) %>%
  # select(-p) %>%
  
  # mutate(pos_neg = c('N', 'P')[(gene_change > 0) + 1]) %>%
  # pivot_wider(names_from = 'pos_neg',
  #             values_from = 'gene_change', values_fill = 0L) %>%
  
  left_join(read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_kegg_koCopy/cafeOut/errorModel/Base_family_results.txt',
                       delim = '\t', show_col_types = FALSE),
            by = c('FamilyID' = '#FamilyID')) %>%
  select(-`Significant at 0.05`) %>%
  rename(family_p = pvalue) %>%
  left_join(cafe_trees %>%
              select(ortho_nodeID, cafe_nodeID) %>%
              distinct,
            by = 'cafe_nodeID')

#### Annotate tree with Cafe results ####
bipartition_colours <- function(x) c("red", "yellow", "green", "purple", "blue", 'black')
cafe_tree_plot <- the_tree %>%
  left_join(select(gene_family_change, -cafe_nodeID),
            by = c('node' = 'ortho_nodeID')) %>%
  ggtree(layout = 'rectangular') +
  
  geom_tiplab(aes(label = str_c(species, ' (', Increase, '/', Decrease, ')')), hjust = 0,
              fontface = "italic") +
  geom_range('CI_date', colour = 'red', size = 3, alpha = 0.3) +
  # geom_nodelab(aes(label = boot_support), vjust = -1, hjust = 1.3) +
  geom_nodepoint(data = . %>% filter(!node %in% c(25)),
                 aes(fill = bipartition_support), size = 3,
                 shape = 21) +
  
  geom_nodelab(aes(label = scales::comma(Increase)), 
               colour = 'darkgreen', vjust = 1, hjust = 1.3) +
  geom_nodelab(aes(label = scales::comma(Decrease)), 
               colour = 'darkred', vjust = 1, hjust = -0.05) +
  
  coord_geo(xlim = c(-425, 0), ylim = c(0, Ntip(the_tree)), 
            neg = TRUE, abbrv = TRUE,
            dat = 'periods', clip = "off") +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn", 
               palette = bipartition_colours,
               breaks = seq(0, 1, by = 0.2),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps") +
  scale_x_continuous(limits = c(-425, 10), labels = abs) +
  guides(fill = guide_coloursteps(title.position = 'top', barwidth = 15, title.hjust = 0.5)) +
  labs(x = 'MYA',
       fill = 'Bipartition Support') +
  theme_tree2() +
  theme(plot.margin = margin(10, 150, 10, 10),
        panel.background = element_blank(),
        plot.background = element_rect(colour = 'black', linewidth = 1),
        legend.position = c(0.3, 0.8),
        legend.direction = "horizontal",
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 16),
        axis.text.x = element_text(colour = 'black', size = 12),
        axis.title.x = element_text(colour = 'black', size = 16))
cafe_tree_plot <- revts(cafe_tree_plot)
cafe_tree_plot
ggsave('../Results/cafe_time_tree_ko.png', plot = cafe_tree_plot, height = 7, width = 10)


#### What KOs expanded/contracted at nodes of interest ####
nodes_of_interest <- tibble(node = c(28, 30, 32, 14),
                            node_description = c('anemone_coral', 'robust_complex', 
                                                 'acroporid_split',
                                                 'acropora cervicornis')) %>%
  mutate(node_description = fct_inorder(node_description))

cafe_trees %>%
  filter(is_significant) %>%
  left_join(select(the_tree, node, orthofinder_code, species),
            by = c('ortho_nodeID' = 'node')) %>%
  rename(node = ortho_nodeID) %>%
  inner_join(nodes_of_interest,
             by = 'node') %>%
  left_join(pathway_changes,
            by = 'FamilyID') %>%
  select(FamilyID, node_description, kegg_gene, path_data) %>%
  unnest(path_data) %>%
  select(-pathway_map, -rel_pathway, -description) %>%
  mutate(kegg_pathway = str_c(name, ' (', kegg_path_id, ')'), .keep = 'unused') %>%
  mutate(category = str_c(major_category, ' - ', minor_category), .keep = 'unused') %>%
  mutate(gene_name = str_c(kegg_gene, ' (', FamilyID, ')'), .keep = 'unused') %>%
  pivot_wider(names_from = category, values_from = kegg_pathway,
              values_fill = NA_character_,
              values_fn = ~str_c(., collapse = '; ')) %>%
  select(gene_name, node_description, order(colnames(.))) %>%
  arrange(node_description, gene_name) %>%
  group_by(node_description) %>%
  mutate(n = n()) %>%
  ungroup %>%
  mutate(node_description = str_c(node_description, ' (', n, ')'), .keep = 'unused') %>%
  gt(rowname_col = 'gene_name', groupname_col = 'node_description') %>%
  tab_header(title = 'KEGG Orthologs with Significant Changes at Important Nodes of Tree') %>%
  tab_stubhead('KEGG Ortholog') %>%
  tab_spanner(label = major_categories[1], columns = starts_with(major_categories[1])) %>%
  tab_spanner(label = major_categories[2], columns = starts_with(major_categories[2])) %>%
  tab_spanner(label = major_categories[3], columns = starts_with(major_categories[3])) %>%
  tab_spanner(label = major_categories[4], columns = starts_with(major_categories[4])) %>%
  tab_spanner(label = major_categories[5], columns = starts_with(major_categories[5])) %>%
  cols_label_with(columns = where(is.character), ~str_remove_all(., str_c(major_categories, collapse = '|'))) %>%
  cols_label_with(columns = where(is.character), ~str_remove_all(., '^ - ')) %>%
  # fmt(columns = 'pvalue', fns = scales::pvalue_format()) %>%
  sub_missing(columns = where(is.character), missing_text = '-') %T>%
  gtsave('../Results/ko_changes_nodes.html') 

cafe_trees %>%
  filter(is_significant) %>%
  left_join(select(the_tree, node, orthofinder_code, species),
            by = c('ortho_nodeID' = 'node')) %>%
  rename(node = ortho_nodeID) %>%
  inner_join(nodes_of_interest,
             by = 'node') %>%
  left_join(pathway_changes,
            by = 'FamilyID') %>%
  select(FamilyID, node_description, kegg_gene, path_data) %>%
  unnest(path_data) %>%
  select(-pathway_map, -rel_pathway, -description) %>%
  mutate(kegg_pathway = str_c(name, ' (', kegg_path_id, ')'), .keep = 'unused') %>%
  mutate(category = str_c(major_category, ' - ', minor_category), .keep = 'unused') %>%
  mutate(gene_name = str_c(kegg_gene, ' (', FamilyID, ')'), .keep = 'unused') %>%
  group_by(node_description) %>%
  summarise(n = n_distinct(gene_name))

#### Heatmap of Increase : Decrease at tree tips ####
species_order <- get_taxa_name(simple_tree)
alpha <- 0; beta <- 0

species_change_prop <- cafe_trees %>%
  left_join(select(the_tree, node, orthofinder_code, species),
            by = c('ortho_nodeID' = 'node')) %>%
  full_join(select(complete_family_species_cafe, FamilyID, gene_change, node_p, family_p, ortho_nodeID),
            by = c('FamilyID', 'ortho_nodeID')) %>%
  filter(!is.na(species)) %>%
  filter(is_significant) %>%
  left_join(pathway_changes,
            by = 'FamilyID') %>%
  select(species, orthofinder_code, gene_change, FamilyID, kegg_gene, family_p, node_p, path_data) %>%
  unnest(path_data) %>%
  select(-pathway_map, -rel_pathway, -description) %>%
  mutate(change_type = case_when(gene_change < 0 ~ 'Decrease',
                                 gene_change > 0 ~ 'Increase',
                                 TRUE ~ 'No Change')) %>%
  group_by(species, orthofinder_code, major_category, minor_category, change_type) %>%
  summarise(n_gene = sum(gene_change),
            n_ko = n_distinct(FamilyID),
            n_path = n_distinct(kegg_path_id),
            .groups = 'drop') %>%
  select(-n_gene, -n_path) %>%
  pivot_wider(names_from = change_type,
              values_from = n_ko, 
              values_fill = 0L) %>%
  mutate(species = factor(species, levels = rev(species_order)),
         total = Decrease + Increase,
         prob_increase = (Increase + alpha) / (total + alpha + beta)) 


species_change_prop %>%
  
  ggplot(aes(x = minor_category, y = species, fill = prob_increase)) +
  geom_raster() +
  facet_grid(~major_category, scales = 'free_x', shrink = TRUE,
             labeller = as_labeller(~str_replace_all(., ' ', '\n')), 
             switch = 'x') +
  scale_x_discrete(position = 'top') +
  scale_fill_gradient2(midpoint = 0.5, labels = scales::percent_format(),
                       breaks = seq(0, 1, by = 0.25), limits = c(0, 1), 
                       na.value = 'black') +
  guides(fill = guide_colorbar(ticks.colour = 'black')) +
  labs(x = NULL,
       y = NULL,
       fill = '% Significant\nKO Expanding') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', face = 'italic'),
        strip.background = element_blank(),
        strip.placement = 'outside',
        axis.text.x = element_text(colour = 'black', angle = 90, hjust = 0),
        panel.background = element_rect(colour = 'black'))
  
species_change_prop %>%
  filter(str_detect(minor_category, 'Immune')) %>% View
  


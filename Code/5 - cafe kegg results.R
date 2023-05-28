#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9266397/ - maybe useful model
##TODO - add robust vs complex clade markings to tree
##TODO - add clade markings for within acropora(?)
##TODO - Make axis title (x MYA) 
##TODO - add cafe +/- to tip labels

#Cafe results. If "family p" is non-significant then no test for significance of individual nodes


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
  select(kegg_path_id, name, description, minor_category, major_category, pathway_map, rel_pathway) %>%
  distinct %>%
  mutate(kegg_path_id = str_remove(kegg_path_id, 'path:'))

kegg_paths

filter(kegg_paths, kegg_path_id %in% paths_particular_interest)

#range for isthmus of panama - Bacon et al 2015 and O'Dea et al 2016
panama_range <- -c(20, 3)

#range for tethys sea
tethys_range <- -c(50) #need to find good citation for range

bipartition_colours <- function(x) c("red", "yellow", "green", "purple", "blue", 'black')

#### Make Preliminary time calibrated tree ####
species_code_rename <- read_csv('../../Bioinformatics/Phylogenomics/tree_data_sources.csv', show_col_types = FALSE) %>%
  select(species, orthofinder_code)


the_tree <- read.beast('../../Bioinformatics/Phylogenomics/Time Calibration/tree.orthofinder.result.date.nexus') %>% 
  mutate(orthofinder_code = label,
         label = str_remove(label , '_.*'),
         label = if_else(str_detect(label, '[0-9]+'), as.character(round(as.numeric(label), 2)), label),
         date = round(-1 * date),
         bipartition_support = if_else(!isTip, as.numeric(label), NA_real_)) %>%
  left_join(species_code_rename,
            by = 'orthofinder_code') 

the_tree %>%
  ggtree() +
  geom_nodelab(aes(label = node)) +
  geom_tiplab()

write_rds(the_tree, '../intermediate_files/updated_tree.rds')

#### Get Divergence Dates ####
simple_tree_plot <- the_tree %>%
  ggtree() +
  geom_nodelab(aes(label = node)) +
  # geom_tiplab(aes(label = node))
  geom_tiplab()

viewClade(simple_tree_plot, MRCA(simple_tree_plot, "ayon", "amil"))

as_tibble(the_tree) %>%
  #25 = anemones vs corals
  #28 = robust vs complex (acropora = complex)
  # ^ https://reefs.com/new-classification-stony-corals/
  # ^ https://link.springer.com/article/10.1186/s13059-018-1552-8
  #30 = montipora vs acropora 
  #32 = initial diversification of acropora
  #38 = acerv split
  filter(node %in% c(25, 28, 30, 32, 38)) %>% 
  unnest(c(CI_height, CI_date)) %>%
  select(node, date, CI_height) %>%
  group_by(node, date) %>%
  mutate(version = if_else(CI_height == min(CI_height), 'lower', 'upper')) %>%
  ungroup %>%
  pivot_wider(names_from = version,
              values_from = CI_height)

#### Time Calibrated Plot ####
tree_plot <- the_tree %>%
  ggtree(layout = 'rectangular') +
  
  # geom_rect(ymin = -Inf, ymax = Inf, 
  #           xmin = panama_range[1], 
  #           xmax = panama_range[2]) +
  # geom_rect(ymin = -Inf, ymax = Inf, 
  #           xmin = tethys_range[1], 
  #           xmax = tethys_range[2]) +
  
  geom_tiplab(aes(label = species), hjust = 0) +
  geom_range('CI_date', colour = 'red', size = 3, alpha = 0.3) +
  # geom_nodelab(aes(label = boot_support), hjust = 0.1) +
  geom_nodepoint(data = . %>% filter(!node %in% c(25)), #, 35
                 aes(fill = bipartition_support), size = 3,
                 shape = 21) +
  
  # geom_cladelab(node = 46, label="test label", hjust = 0, offset = 0) +
  
  coord_geo(xlim = c(-425, 0), ylim = c(0, Ntip(the_tree)), 
            neg = TRUE, abbrv = TRUE,
            dat = 'periods', clip = "off") +
  scale_x_continuous(limits = c(-425, 10), labels = abs) +
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn", 
               palette = bipartition_colours,
               breaks = seq(0, 1, by = 0.2),
               limits = c(0, 1),
               show.limits = TRUE, 
               guide = "colorsteps") +
  guides(fill = guide_coloursteps(title.position = 'top')) +
  labs(x = 'MYA',
       fill = 'Bipartition Support') +
  theme_tree2() +
  theme(plot.margin = margin(10, 150, 10, 10),
        panel.background = element_blank(),
        plot.background = element_rect(colour = 'black', linewidth = 1),
        legend.position = c(0.2, 0.8),
        legend.direction = "horizontal",
        axis.text.x = element_text(colour = 'black', size = 14),
        axis.title.x = element_text(colour = 'black', size = 18))
time_tree_plot <- revts(tree_plot)
time_tree_plot
ggsave('../Results/time_tree.png', plot = time_tree_plot, height = 7, width = 10)
viewClade(time_tree_plot, MRCA(time_tree_plot, "ayon", "amil"))

#### Significance of KEGG paths ####
pathway_changes <- read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut/errorModel/Base_family_results.txt',
           delim = '\t', show_col_types = FALSE) %>%
  rename(FamilyID = `#FamilyID`) %>%
  full_join(kegg_paths,
            by = c('FamilyID' = 'kegg_path_id'))

pathway_changes %>%
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
  labs(x = 'KEGG Pathways (%)',
       y = NULL,
       colour = NULL) +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        axis.text = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 14),
        legend.position = 'bottom')
ggsave('../Results/significant_change_pathways.png', height = 7, width = 15, scale = 1)

filter(pathway_changes, FamilyID %in% paths_particular_interest)

pathway_changes %>%
  filter(pvalue < 0.05) %>%
  select(major_category, minor_category, name, FamilyID, pvalue) %>%
  arrange(major_category, minor_category, name) %>%
  gt(rowname_col = 'major_category') %>%
  tab_header(title = 'KEGG Pathways with Significant Changes in number of unique KEGG Orthologs across Anthozoan Tree') %>%
  tab_stubhead('Major KEGG Category') %>%
  cols_label(minor_category = 'Minor KEGG Category',
             name = 'KEGG Pathway',
             pvalue = 'p-value') %>%
  fmt(columns = 'pvalue', fns = scales::pvalue_format()) %>%
  cols_merge(columns = c('name', 'FamilyID'), pattern = "{1} ({2})") %T>%
  gtsave('../Results/significantly_evolving_pathways.html')

#### Read in Cafe Results - KEGG Paths ####
cafe_trees <- read_lines('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut/errorModel/Base_asr.tre') %>%
  str_subset('^ *TREE') %>%
  tibble(tree = .) %>%
  mutate(FamilyID = str_extract(tree, 'map[0-9]+'),
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


gene_family_change <- read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut/errorModel/Base_clade_results.txt',
                                 delim = '\t', show_col_types = FALSE) %>%
  rename(cafe_nodeID = `#Taxon_ID`) %>%
  left_join(cafe_trees %>%
              select(ortho_nodeID, cafe_nodeID) %>%
              distinct,
            by = 'cafe_nodeID')


complete_family_species_cafe <- full_join(
  read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut/errorModel/Base_change.tab', 
             delim = '\t', show_col_types = FALSE) %>%
    pivot_longer(cols = -FamilyID,
                 names_to = 'cafe_nodeID',
                 values_to = 'gene_change'),
  
  read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut/errorModel/Base_branch_probabilities.tab', 
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
  
  left_join(read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut/errorModel/Base_family_results.txt',
                       delim = '\t', show_col_types = FALSE),
            by = c('FamilyID' = '#FamilyID')) %>%
  select(-`Significant at 0.05`) %>%
  rename(family_p = pvalue) %>%
  left_join(cafe_trees %>%
              select(ortho_nodeID, cafe_nodeID) %>%
              distinct,
            by = 'cafe_nodeID')

#### Annotate tree with Cafe results ####
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
ggsave('../Results/cafe_time_tree.png', plot = cafe_tree_plot, height = 7, width = 10)

#### Trees for Individual Pathways ####
# tree <- individual_pathway_trees$tree[[1]]
make_tree_plot <- function(tree){
  tree %>%
    mutate(species = str_c(species, ' (', n_gene, ')')) %>%
    ggtree(layout = 'rectangular') +
    
    geom_tiplab(aes(label = species), hjust = 0,
                fontface = "italic") +
    geom_tippoint(aes(colour = is_significant)) +
    geom_nodelab(aes(label = scales::comma(n_gene)),
                 hjust = 1.5, vjust = -0.5) +
    geom_nodepoint(aes(colour = is_significant)) +
    scale_x_continuous(limits = c(-10, 425)) +
    scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
    guides(colour = 'none') +
    theme_tree2()
}

individual_pathway_trees <- cafe_trees %>%
  nest(data = -c(FamilyID)) %>%
  left_join(pathway_changes,
            by = 'FamilyID') %>%
  rowwise %>%
  mutate(tree = list(left_join(the_tree, 
                               data,
                               by = c('node' = 'ortho_nodeID')))) %>%
  ungroup %>%
  # filter(FamilyID %in% paths_particular_interest) %>%
  rowwise %>%
  mutate(tree_plot = list(make_tree_plot(tree) + labs(title = name, subtitle = scales::pvalue(pvalue)))) %>%
  ungroup

individual_pathway_trees %>%
  filter(FamilyID == 'map00430') %>%
  pull(tree_plot) %>%
  pluck(1)

#### What pathways expanded/contracted at nodes of interest ####
simple_tree_plot
nodes_of_interest <- tibble(node = c(28, 30, 32, 14),
                           node_description = c('anemone_coral', 'robust_complex', 
                                           'acroporid_split',
                                           'acropora cervicornis')) %>%
  mutate(node_description = fct_inorder(node_description))

node_changes <- cafe_trees %>%
  filter(is_significant) %>%
  left_join(select(the_tree, node, orthofinder_code, species),
            by = c('ortho_nodeID' = 'node')) %>%
  rename(node = ortho_nodeID) %>%
  inner_join(nodes_of_interest,
             by = 'node') %>%
  left_join(pathway_changes,
            by = 'FamilyID') %>%
  select(node_description, name, minor_category, major_category, FamilyID) %>%
  arrange(node_description)
  # count(node_description)
  
count(node_changes, node_description)
count(node_changes, FamilyID)
  
node_changes %>%
  gt()

individual_pathway_trees %>%
  filter(FamilyID == 'map00983') %>%
  pull(tree_plot) %>%
  pluck(1)


cafe_trees %>%
  filter(is_significant) %>%
  left_join(select(the_tree, node, orthofinder_code, species),
            by = c('ortho_nodeID' = 'node')) %>%
  rename(node = ortho_nodeID) %>%
  inner_join(nodes_of_interest,
             by = 'node') %>%
  left_join(pathway_changes,
            by = 'FamilyID') %>%
  select(-pathway_map, -rel_pathway, -description) %>%
  mutate(kegg_pathway = str_c(name, ' (', FamilyID, ')'), .keep = 'unused') %>%
  select(node_description, kegg_pathway, major_category, minor_category) %>%
  pivot_wider(names_from = node_description, values_from = kegg_pathway,
              values_fill = NA_character_,
              values_fn = ~str_c(., collapse = '; ')) %>%
  select(major_category, minor_category, anemone_coral, robust_complex, acroporid_split, `acropora cervicornis`) %>%
  arrange(major_category, minor_category) %>%
  gt(rowname_col = 'minor_category', groupname_col = 'major_category') %>%
  tab_header(title = 'KEGG Pathways with Significant Changes at Important Nodes of Tree') %>%
  tab_stubhead('KEGG Pathway') %>%
  cols_label(anemone_coral = 'Anemone - Scleractinian',
             robust_complex = 'Robust - Complex',
             acroporid_split = 'Acroporid Diversification',
             `acropora cervicornis` = md('*Acropora cervicornis*')) %>%
  # fmt(columns = 'pvalue', fns = scales::pvalue_format()) %>%
  sub_missing(columns = where(is.character), missing_text = '-') %T>%
  gtsave('../Results/pathway_changes_nodes.html') 

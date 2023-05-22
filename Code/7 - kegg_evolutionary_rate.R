##TODO - why are numbers of families expand/contract so huge with gamma kegg...
##TODO - What is the difference between posterior probability in 'Gamma_family_likelihoods.txt' and pvalue in Gamma_family_results.txt and which to use...
##        - in family likelihoods for each path they sum to 1 across all rate groups...so it seesm like that posterior must be testing "fit" of lambda not expansion of catgegory...
##TODO - change analysis such that first there are significant families/nodes. and there are also some families which significantly differ from the base evolutionary rate - breakdown of what those families are
##TODO - REDO/INCORPORATE # of total gene copies along with the unique keggs

#### Libraries ####
library(tidyverse)
library(magrittr)
library(treedataverse)
library(deeptime)
library(brms)
library(tidybayes)
library(patchwork)
library(gt)

library(multidplyr)
cluster <- new_cluster(parallel::detectCores() - 1) 
cluster_library(cluster, c('ape', 'tidytree', 'tibble', 'dplyr', 'stringr'))

paths_particular_interest <- c('map04620', 'map04624', 'map00260', 'map00270', 'map01100', 'map01110', 'map01230', 'map01120', 'map04361')
cafe_folder <- '../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggCopies/cafeOut'

#### KEGG Pathway data ####
kegg_paths <- read_csv('../intermediate_files/kegg_orthogroup_pathways.csv.gz',
                       show_col_types = FALSE) %>%
  filter(major_category != 'Human Diseases') %>%
  select(kegg_path_id, name, description, minor_category, major_category, pathway_map, rel_pathway) %>%
  distinct %>%
  mutate(kegg_path_id = str_remove(kegg_path_id, 'path:'))

the_tree <- read_rds('../intermediate_files/updated_tree.rds')

species_pathway_keggs <- read_delim('../intermediate_files/kegg_families.txt', show_col_types = FALSE) %>%
  select(-Desc) %>% 
  rename(FamilyID = 'Family ID') %>%
  pivot_longer(cols = -FamilyID,
               names_to = 'species',
               values_to = 'n_kegg')

#### Redo Cafe but with various evolutionary rate groups ####
evolution_rate <- list.files(cafe_folder, 
                             pattern = 'Gamma_results.txt|Base_results.txt', recursive = TRUE,
                             full.names = TRUE) %>%
  tibble(file = .) %>%
  filter(str_detect(file, 'K[0-9]+')) %>%
  mutate(k_group = str_extract(file, 'K[0-9]+') %>% str_remove('K') %>% as.integer) %>%
  rowwise(k_group) %>%
  summarise(line = read_lines(file, n_max = 2),
            .groups = 'keep') %>%
  mutate(value = if_else(row_number() == 1, 'likelihood', 'lambda')) %>%
  ungroup %>%
  mutate(line = str_extract(line, '[0-9\\.]+') %>% as.numeric()) %>%
  pivot_wider(names_from = 'value', values_from = 'line') 

evolution_rate %>%
  filter(k_group != 7) %>%
  ggplot(aes(x = k_group, y = likelihood)) +
  geom_line() +
  geom_point(aes(colour = k_group == 4), show.legend = FALSE) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
  labs(x = 'Number of Evolutionary Rate Groups',
        y = '-log(likelihood)') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'))
ggsave('../Results/number_gamma_categories.png', height = 5, width = 5)

evolution_rate %>%
  ggplot(aes(x = k_group, y = lambda)) +
  geom_line() +
  geom_point()

evolution_rate %>%
  ggplot(aes(x = likelihood, y = lambda)) +
  geom_line() +
  geom_point()

# k_categories <- evolution_rate$k_group[which.max(evolution_rate$likelihood)]
k_categories <- 1
global_lambda <- evolution_rate$lambda[evolution_rate$k_group == k_categories]

#### Families with Significant Gene Evolution ####
pathway_changes <- read_delim(str_c(cafe_folder, '/K', k_categories, 
                                    if_else(k_categories == 1, '/Base_family_results.txt', '/Gamma_family_results.txt')),
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
  gtsave('../Results/significantly_evolving_pathways.pdf')  %T>%
  gtsave('../Results/significantly_evolving_pathways.png')


#### What Pathways are evolving at a signficantly different rate than others ####
othrogroup_evolutionary_rates <- list.files(cafe_folder, 
                                            pattern = 'Gamma_family_likelihoods.txt|Base_family_likelihoods.txt', 
                                            recursive = TRUE,
                                            full.names = TRUE) %>%
  tibble(file = .) %>%
  filter(str_detect(file, 'K[0-9]+')) %>%
  mutate(k_group = str_extract(file, 'K[0-9]+') %>% str_remove('K') %>% as.integer) %>%
  
  #Pick number of Ks based on minimizing -loglikelihood
  filter(k_group == k_categories) %>%
  rowwise(k_group) %>%
  summarise(data = list(read_delim(file, delim = '\t', show_col_types = FALSE) %>% janitor::clean_names()),
            .groups = 'drop') %>%
  select(data) %>%
  unnest(data) %>%
  filter(significant != 'N/S') %>%
  rename(FamilyID = number_family_id) %>%
  full_join(kegg_paths,
            by = c('FamilyID' = 'kegg_path_id')) %>%
  mutate(pvalue = 1 - posterior_probability)

othrogroup_evolutionary_rates %>%
  filter(!is.na(gamma_cat_mean)) %>%
  mutate(rate_group = if_else(gamma_cat_mean == min(gamma_cat_mean), 'Slower', 'Faster')) %>%
  select(major_category, minor_category, name, FamilyID, rate_group) %>%
  
  arrange(major_category, minor_category, name) %>%
  gt(rowname_col = 'major_category', groupname_col = 'rate_group') %>%
  tab_header(title = 'KEGG Pathways Significantly Faster or Slower Evolutionary Rates than background in Anthozoans') %>%
  tab_stubhead('Major KEGG Category') %>%
  cols_label(minor_category = 'Minor KEGG Category',
             name = 'KEGG Pathway') %>%
  cols_merge(columns = c('name', 'FamilyID'), pattern = "{1} ({2})") %T>%
  gtsave('../Results/differential_evolutionary_rate.pdf')  %T>%
  gtsave('../Results/differential_evolutionary_rate.png')

  
#### Redo Cafe evolutionary Rates ####
cafe_trees <- read_lines(str_c(cafe_folder, '/K', k_categories, 
                               if_else(k_categories == 1, '/Base_asr.tre', '/Gamma_asr.tre'))) %>%
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


gene_family_change <- read_delim(str_c(cafe_folder, '/K', k_categories,
                                       if_else(k_categories == 1, '/Base_clade_results.txt', '/Gamma_clade_results.txt')),
                                 delim = '\t', show_col_types = FALSE) %>%
  rename(cafe_nodeID = `#Taxon_ID`) %>%
  left_join(cafe_trees %>%
              select(ortho_nodeID, cafe_nodeID) %>%
              distinct,
            by = 'cafe_nodeID')

complete_family_species_cafe <- full_join(
  read_delim(str_c(cafe_folder, '/K', k_categories, 
                   if_else(k_categories == 1, '/Base_change.tab', '/Gamma_change.tab')), 
             delim = '\t', show_col_types = FALSE) %>%
    pivot_longer(cols = -FamilyID,
                 names_to = 'cafe_nodeID',
                 values_to = 'gene_change'),
  
  read_delim(str_c(cafe_folder, '/K', k_categories, 
                   if_else(k_categories == 1, '/Base_branch_probabilities.tab', '/Gamma_branch_probabilities.tab')), 
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
  
  left_join(read_delim(str_c(cafe_folder, '/K', k_categories, 
                             if_else(k_categories == 1, '/Base_family_results.txt', '/Gamma_family_results.txt')),
                       delim = '\t', show_col_types = FALSE),
            by = c('FamilyID' = '#FamilyID')) %>%
  select(-`Significant at 0.05`) %>%
  rename(family_p = pvalue) %>%
  left_join(cafe_trees %>%
              select(ortho_nodeID, cafe_nodeID) %>%
              distinct,
            by = 'cafe_nodeID')

#### Make Tree with Evolutionary Categories ####
bipartition_colours <- function(x) c("red", "yellow", "green", "purple", "blue", 'black')

cafe_tree_plot <- the_tree %>%
  left_join(select(gene_family_change, -cafe_nodeID),
            by = c('node' = 'ortho_nodeID')) %>% 
  ggtree(layout = 'rectangular') +
  
  geom_tiplab(aes(label = scales::comma(Increase)), hjust = 0,
              colour = 'darkgreen') +
  geom_tiplab(aes(label = scales::comma(Decrease)), hjust = -1,
              colour = 'darkred') +
  
  geom_tiplab(aes(label = species), hjust = -0.5,
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

#### Characterize Node Changes ####
gene_family_change

complete_family_species_cafe %>%
  filter(str_detect(cafe_nodeID, '<44>')) %>%
  left_join(kegg_paths, by = c('FamilyID' = 'kegg_path_id')) %>%
  filter(str_detect(minor_category, 'Immune'))
  # filter(!is.na(node_p)) %>%
  arrange(node_p)
  filter(node_p < 0.05)


#### Write output to use in future ####
as_tibble(the_tree) %>%
  select(parent, node, label, orthofinder_code, species) %>%
  rename(ortho_parentID = parent,
         ortho_nodeID = node) %>%
  full_join(complete_family_species_cafe, 
            by = 'ortho_nodeID') %>%
  write_csv('../intermediate_files/geneChange_tree_out.csv')


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
  filter(FamilyID %in% paths_particular_interest | pvalue < 0.05) %>%
  rowwise %>%
  mutate(tree_plot = list(make_tree_plot(tree) + 
                            labs(title = name, subtitle = scales::pvalue(pvalue, add_p = TRUE)))) %>%
  ungroup


individual_pathway_trees$tree_plot[[3]]

individual_pathway_trees %>%
  filter(FamilyID %in% paths_particular_interest) %>%
  pull(tree_plot) %>%
  wrap_plots()


individual_pathway_trees %>%
  filter(str_detect(minor_category, 'Immune')) %>%
  pull(tree_plot) %>%
  wrap_plots()

individual_pathway_trees %>%
  count(minor_category)

individual_pathway_trees %>%
  # filter(FamilyID == 'map04361')
  filter(str_detect(minor_category, 'Develop')) %>%
  pull(tree_plot) %>%
  wrap_plots()


individual_pathway_trees %>%
  group_by(major_category, minor_category) %>%
  summarise(plot = list(wrap_plots(tree_plot)),
            n = n()) %>%
  filter(minor_category == 'Amino acid metabolism') %>%
  pull(plot) %>%
  pluck(1)

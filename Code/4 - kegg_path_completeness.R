#### Libraries ####
library(treedataverse)
library(tidyverse)
library(broom)
library(patchwork)

#### Data ####
orthofinder_tree <- read.newick('../../Bioinformatics/Phylogenomics/Species_Tree/SpeciesTree_rooted.txt') %>% 
  as_tibble %>%
  as.phylo() %>%
  ggtree(ladderize = TRUE, branch.length = 'none', right = TRUE) +
  geom_tiplab(as_ylab = TRUE, color = 'black')
species_order <- get_taxa_name(orthofinder_tree)

species_names <- read_csv('../../Bioinformatics/Phylogenomics/tree_data_sources.csv', 
         show_col_types = FALSE) %>%
  select(species, orthofinder_code) %>%
  filter(orthofinder_code %in% species_order) %>%
  mutate(orthofinder_code = factor(orthofinder_code, levels = species_order)) %>%
  arrange(orthofinder_code) %>%
  pull(species) %>%
  set_names(species_order)

full_pathway_data <- read_csv('../intermediate_files/pathdata_for_discovery.csv', show_col_types = FALSE)

kegg_paths <- read_csv('../intermediate_files/kegg_orthogroup_pathways.csv.gz', show_col_types = FALSE) %>%
  select(kegg_path_id, name, description, major_category, minor_category) %>%
  distinct

filter(kegg_paths, kegg_path_id %in% str_c('path:', 'map00270'))

#### Read in Model & Results ####
path_complete_model <- read_rds('../intermediate_files/kegg_pathway_completeness_fixed.rds')

major_cat_emmeans <- read_rds('../intermediate_files/kegg_pathway_completeness_fixed_major.rds.gz')
minor_cat_emmeans <- read_rds('../intermediate_files/kegg_pathway_completeness_fixed_minor.rds.gz')

#### Investigate Model ####
plot(path_complete_model)
pp_check(path_complete_model)

summary(path_complete_model)

#### Major Category Comparisons ####
summary(major_cat_emmeans, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  mutate(species = factor(species, levels = rev(species_order))) %>%
  ggplot(aes(y = species, x = prob, xmin = lower.HPD, xmax = upper.HPD)) +
  geom_pointrange() +
  scale_x_continuous(labels = scales::percent_format(suffix = '')) +
  scale_y_discrete(labels = species_names) +
  labs(y = NULL, x = 'Pathway Completeness (%)') +
  facet_wrap(~ major_category, nrow = 1, 
             labeller = as_labeller(function(x) str_replace_all(x, ' ', '\n'))) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(face = 'italic'),
        axis.title = element_text(colour = 'black', size = 16),
        strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'))

#### Minor Categories ####
if(file.exists('../intermediate_files/minor_category_path_completeness.csv')){
  minor_cats <- read_csv('../intermediate_files/minor_category_path_completeness.csv', show_col_types = FALSE)
} else {
  minor_cats <- summary(minor_cat_emmeans, type = 'response') %>%
    broom::tidy(conf.int = TRUE)
  write_csv(minor_cats, '../intermediate_files/minor_category_path_completeness.csv')
}


minor_cat_plots <- mutate(minor_cats, species = factor(species, levels = rev(species_order))) %>%
  nest_by(major_category) %>%
  mutate(plot = list(ggplot(data, aes(y = species, x = prob, xmin = lower.HPD, xmax = upper.HPD)) +
                       geom_pointrange() +
                       scale_x_continuous(labels = scales::percent_format(suffix = '')) +
                       scale_y_discrete(labels = species_names) +
                       labs(y = NULL, x = 'Pathway Completeness (%)',
                            title = major_category) +
                       facet_wrap(~ minor_category, 
                                  labeller = as_labeller(function(x) str_replace_all(x, ' ', '\n'))) +
                       theme_classic() +
                       theme(axis.text = element_text(colour = 'black', size = 12),
                             axis.text.y = element_text(face = 'italic'),
                             axis.title = element_text(colour = 'black', size = 16),
                             strip.background = element_blank(),
                             panel.background = element_rect(colour = 'black')))) 
minor_cat_plots$plot[[4]]
minor_cat_plots$plot[[5]]


#### Individual Paths ####
library(tidybayes)

immune_paths <- full_pathway_data %>%
  filter(str_detect(minor_category, 'Immune')) %>%
  # select(kegg_path_id, species) %>%
  # distinct %>%
  mutate(total_keggs = 1) %>%
  add_epred_draws(path_complete_model) %>%
  point_interval() %>%
  left_join(kegg_paths,
            by = 'kegg_path_id')

  
immune_paths$name %>% unique

immune_paths %>%
  filter(str_detect(name, 'Toll')) %>%
  ggplot(aes(y = species, x = .epred, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  scale_x_continuous(labels = scales::percent_format(suffix = '')) +
  # scale_y_discrete(labels = species_names) +
  labs(y = NULL, x = 'Pathway Completeness (%)') +
  facet_wrap(~ name, nrow = 1) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(face = 'italic'),
        axis.title = element_text(colour = 'black', size = 16),
        strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'))


full_pathway_data %>%
  filter(kegg_path_id %in% str_c('path:', 'map00270')) %>%
  # select(kegg_path_id, species) %>%
  # distinct %>%
  mutate(total_keggs = 1) %>%
  add_epred_draws(path_complete_model) %>%
  point_interval() %>%
  left_join(kegg_paths,
            by = 'kegg_path_id') %>%
  ggplot(aes(y = species, x = .epred, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_point(aes(x = prop_complete), colour = 'red') +
  scale_x_continuous(labels = scales::percent_format(suffix = '')) +
  # scale_y_discrete(labels = species_names) +
  labs(y = NULL, x = 'Pathway Completeness (%)') +
  facet_wrap(~ name, nrow = 1) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(face = 'italic'),
        axis.title = element_text(colour = 'black', size = 16),
        strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'))

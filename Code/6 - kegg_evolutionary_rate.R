##TODO - why are numbers of families expand/contract so huge with gamma kegg...

#### Libraries ####
library(tidyverse)
library(magrittr)
library(treedataverse)
library(deeptime)
library(brms)
library(tidybayes)
library(patchwork)

library(multidplyr)
cluster <- new_cluster(parallel::detectCores() - 1) 
cluster_library(cluster, c('ape', 'tidytree', 'tibble', 'dplyr', 'stringr'))

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
evolution_rate <- list.files('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut', 
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
k_categories <- 4
global_lambda <- evolution_rate$lambda[evolution_rate$k_group == k_categories]

#### Evolutionary Rates of each KEGG pathway & significance ####
othrogroup_evolutionary_rates <- list.files('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut', 
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
  group_by(number_family_id) %>%
  
  #For each family choose maximum likelihood gamma mean value
  filter(likelihood_of_category == max(likelihood_of_category)) %>%
  ungroup %>%
  rename(FamilyID = number_family_id) %>%
  full_join(kegg_paths,
            by = c('FamilyID' = 'kegg_path_id')) %>%
  mutate(pvalue = 1 - posterior_probability) %>%
  nest(data = -c(gamma_cat_mean)) %>%
  arrange(gamma_cat_mean) %>%
  mutate(gamma_category = LETTERS[row_number()],
         gamma_category = if_else(is.na(gamma_cat_mean), NA_character_, gamma_category)) %>%
  unnest(data)

othrogroup_evolutionary_rates %>%
  # filter(minor_category == 'Development and regeneration') %>%
  group_by(major_category, minor_category) %>%
  summarise(sig_paths = n_distinct(FamilyID[pvalue < 0.05 & !is.na(pvalue)]),
            tested_paths = n_distinct(FamilyID[!is.na(pvalue)]),
            total_paths = n_distinct(FamilyID),
            untested_paths = n_distinct(FamilyID[is.na(pvalue)]),
            .groups = 'drop') %>%
  select(-untested_paths) %>%
  pivot_longer(cols = ends_with('paths')) %>%
  mutate(name = str_remove(name, '_paths'),
         name = case_when(name == 'sig' ~ 'Significant',
                          name == 'tested' ~ 'Tested',
                          name == 'total' ~ 'Total'),
         name = factor(name, levels = rev(c('Significant', 'Tested', 'Total')))) %>%
  # filter(tested_paths < sig_paths)
  
  ggplot(aes(x = value, y = minor_category, colour = name)) +
  geom_linerange(aes(xmin = 0, xmax = value),
                 position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~major_category, scales = 'free_y') +
  guides(colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  labs(x = 'Number of KEGG Pathways',
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

#### Composition of Evolutionary Rate Groups ####
# data <- tmp; catagories <- unique(kegg_paths$major_category); group_var <- 'major_category'; n_var <- 'n_keggs'
format_brms_data <- function(data, catagories, group_var = 'major_category', n_var = 'n_keggs'){
  pivot_wider(data, 
              names_from = all_of(group_var), 
              values_from = all_of(n_var),
              values_fill = 0L) %>%
    mutate(Y = cbind(!!!syms(catagories)),
           total = rowSums(Y)) %>%
    select(-all_of(catagories))
}

brms_count_model_fitting <- function(data){
  brm(total ~ gamma_category + (1 | species),
      family = poisson(),
      prior = prior(normal(0, 10), 
                    class = 'b'), 
      data = data,
      sample_prior = 'yes',
      chains = 4,
      cores = 4,
      iter = 2000,
      warmup = 1000,
      backend = 'cmdstanr')
}

brms_composition_model_fitting <- function(data){
  brm(Y | trials(total) ~ gamma_category + (1 | species),
      family = multinomial(),
      prior = prior(normal(0, 10), 
                    class = 'b'), 
      data = data,
      sample_prior = 'yes',
      chains = 4,
      cores = 4,
      iter = 2000,
      warmup = 1000,
      backend = 'cmdstanr')
}

count_plot <- function(data, model){
  data %>%
    select(gamma_category) %>%
    distinct %>%
    mutate(total = 1) %>%
    add_epred_draws(model, re_formula = NA) %>%
    point_interval() %>%
    ungroup %>%
    
    ggplot(aes(x = gamma_category, y = .epred, ymin = .lower, ymax = .upper)) +
    geom_pointrange() +
    # facet_wrap(~ .category) +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(y = 'Orthogroups (#)',
         x = 'Evolutionary Rate (λ)') +
    theme_classic() +
    theme(axis.title = element_text(colour = 'black', size = 16),
          axis.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 16),
          legend.text = element_text(colour = 'black', size = 12),
          panel.background = element_rect(colour = 'black'))
}

composition_plot <- function(data, model){
  data %>%
    select(gamma_category) %>%
    distinct %>%
    mutate(total = 1) %>%
    add_epred_draws(model, re_formula = NA) %>%
    point_interval() %>%
    ungroup %>%
    
    ggplot(aes(x = gamma_category, y = .epred, ymin = .lower, ymax = .upper, 
               colour = .category)) +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = 'Orthogroups (%)',
         x = 'Evolutionary Rate (λ)') +
    theme_classic() +
    theme(axis.title = element_text(colour = 'black', size = 16),
          axis.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 16),
          legend.text = element_text(colour = 'black', size = 12),
          panel.background = element_rect(colour = 'black'),
          legend.position = 'bottom')
}

major_cat_data <- select(othrogroup_evolutionary_rates, FamilyID, gamma_cat_mean, gamma_category, name, 
       minor_category, major_category) %>%
  filter(!is.na(gamma_category)) %>%
  mutate(gamma_category = factor(gamma_category, ordered = TRUE)) %>%
  left_join(species_pathway_keggs,
            by = 'FamilyID') %>%
  
  group_by(major_category, species, gamma_category) %>%
  summarise(n_paths = n_distinct(FamilyID),
            n_keggs = sum(n_kegg),
            .groups = 'drop') %>%
  select(-n_paths) %>%
  format_brms_data(unique(kegg_paths$major_category))

major_count_model <- brms_count_model_fitting(major_cat_data)
major_comp_model <- brms_composition_model_fitting(major_cat_data)


count_plot(major_cat_data, major_count_model) / 
  (composition_plot(major_cat_data, major_comp_model) + 
     labs(colour = NULL) + guides(colour = guide_legend(ncol = 2))) &
  plot_annotation(title = 'Major KEGG Categories')
ggsave('../Results/evolutionary_rate_major_keggs.png', height = 12, width = 7)

minor_category_models <- select(othrogroup_evolutionary_rates, FamilyID, gamma_cat_mean, gamma_category, name, 
       minor_category, major_category) %>%
  filter(!is.na(gamma_category)) %>%
  mutate(gamma_category = factor(gamma_category, ordered = TRUE)) %>%
  left_join(species_pathway_keggs,
            by = 'FamilyID') %>%
  
  group_by(major_category, minor_category, species, gamma_category) %>%
  summarise(n_paths = n_distinct(FamilyID),
            n_keggs = sum(n_kegg),
            .groups = 'drop') %>%
  select(-n_paths) %>%
  nest(data = -c(major_category)) %>%
  left_join(kegg_paths %>%
              group_by(major_category) %>%
              summarise(internal_categories = list(unique(minor_category))),
            by = 'major_category') %>%
  rowwise %>%
  mutate(data = list(format_brms_data(data, internal_categories, group_var = 'minor_category', n_var = 'n_keggs'))) %>%
  mutate(count_model = list(brms_count_model_fitting(data)),
         composition_model = list(brms_composition_model_fitting(data)))

minor_cat_plots <- minor_category_models %>%
  mutate(count_plot = list(count_plot(data, count_model)),
         composition_plot = list(composition_plot(data, composition_model) + 
                                   labs(colour = 'Minor KEGG Category'))) %>%
  rowwise() %>%
  mutate(combined_plot = list(wrap_plots(count_plot + labs(title = str_replace_all(major_category, ' ', '\n')), 
                                         composition_plot + 
                                           labs(colour = NULL) + 
                                           guides(colour = guide_legend(ncol = 1)), 
                                         ncol = 1))) %>%
  ungroup %>%
  summarise(plot = list(wrap_plots(combined_plot, ncol = nrow(.)))) %>%
  pull(plot) %>%
  pluck(1)

ggsave('../Results/evolutionary_rate_minor_keggs.png', plot = minor_cat_plots, width = 18, height = 12)

#### Redo Cafe evolutionary Rates ####
cafe_trees <- read_lines(str_c('../../Bioinformatics/Phylogenomics/Time Calibration/cafe_keggPaths/cafeOut/K', k_categories, '/Gamma_asr.tre')) %>%
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
                     label = str_replace(label, '\\<[0-9]+\\>', str_c('<', node, '>'))) %>%
              select(label, cafe_nodeID)) %>%
  collect %>%
  ungroup %>%
  mutate(is_significant = str_detect(label, '\\*'),
         ortho_nodeID = str_extract(label, '\\<[0-9]+\\>') %>% str_remove_all(c('\\<|\\>')) %>% as.integer,
         n_gene = str_extract(label, '[0-9]+$') %>% as.integer()) 


gene_family_change <- read_delim(str_c('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/K', k_categories, '/Gamma_clade_results.txt'),
                                 delim = '\t', show_col_types = FALSE) %>%
  rename(cafe_nodeID = `#Taxon_ID`) %>%
  left_join(cafe_trees %>%
              select(ortho_nodeID, cafe_nodeID) %>%
              distinct,
            by = 'cafe_nodeID')


complete_family_species_cafe <- full_join(
  read_delim(str_c('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/K', k_categories, '/Gamma_change.tab'), 
             delim = '\t', show_col_types = FALSE) %>%
    pivot_longer(cols = -FamilyID,
                 names_to = 'cafe_nodeID',
                 values_to = 'gene_change'),
  
  read_delim(str_c('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/K', k_categories, '/Gamma_branch_probabilities.tab'), 
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
  
  left_join(read_delim(str_c('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/K', k_categories, '/Gamma_family_results.txt'),
                       delim = '\t', show_col_types = FALSE),
            by = c('FamilyID' = '#FamilyID')) %>%
  select(-`Significant at 0.05`) %>%
  rename(family_p = pvalue) %>%
  left_join(cafe_trees %>%
              select(ortho_nodeID, cafe_nodeID) %>%
              distinct,
            by = 'cafe_nodeID')

#### Make Tree with Evolutionary Categories ####
cafe_tree_plot <- the_tree %>%
  left_join(select(gene_family_change, -cafe_nodeID),
            by = c('node' = 'ortho_nodeID')) %>%
  ggtree(layout = 'rectangular') +
  
  geom_tiplab(aes(label = species), hjust = 0,
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


#### Write output to use in future ####

as_tibble(the_tree) %>%
  select(parent, node, label, orthofinder_code, species) %>%
  rename(ortho_parentID = parent,
         ortho_nodeID = node) %>%
  full_join(complete_family_species_cafe, 
            by = 'ortho_nodeID') %>%
  write_csv('../intermediate_files/geneChange_tree_out.csv')


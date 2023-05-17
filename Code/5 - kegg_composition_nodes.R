#### Libraries ####
library(tidyverse)
library(broom)

#### Functions ####
fisher_posthoc <- function(data, group_col, addition = 1L){
  data %>%
    mutate(not_contract = sum(contraction) - contraction,
           not_expand = sum(expansion) - expansion,
           not_none = sum(none) - none) %>%
    mutate(across(where(is.integer), ~. + addition)) %>%
    rename(n_contract = contraction,
           n_expand = expansion,
           n_none = none) %>%
    
    pivot_longer(cols = -all_of(group_col), 
                 names_pattern = '(.*)_(.*)',
                 names_to = c('.value', 'type')) %>%
    nest(data = -c(all_of(group_col))) %>%
    rowwise(group_col) %>%
    summarise(filter(data, type != 'expand') %>%
                column_to_rownames('type') %>%
                fisher.test() %>% tidy %>%
                select(-contains('conf'), -ends_with('method'),
                       -ends_with('alternative')) %>%
                rename_with(~str_c('contract', ., sep = '_')),
              
              filter(data, type != 'contract') %>%
                column_to_rownames('type') %>%
                fisher.test() %>% tidy %>%
                select(-contains('conf'), -ends_with('method'),
                       -ends_with('alternative')) %>%
                rename_with(~str_c('expand', ., sep = '_')),
              .groups = 'drop') #the estimate is the odds of a category being expanded/contracted vs the base case of no change
  
}

#### Data ####
expansion_contractions <- read_delim('../../Bioinformatics/Phylogenomics/Time Calibration/cafeOut/errorModel/Base_change.tab', 
                                     delim = '\t', show_col_types = FALSE) %>%
  rename(Orthogroup = FamilyID)

kegg_paths <- read_csv('../intermediate_files/kegg_orthogroup_pathways.csv.gz',
                       show_col_types = FALSE) %>%
  filter(major_category != 'Human Diseases')

#### Test for differences in composition of each node between expansions/contractions compared to base case of none in terms of composition 
major_expansion_contraction <- expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction')) %>%
  filter(!is.na(major_category)) %>%
  group_by(node, change, major_category) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  pivot_wider(names_from = 'change',
              values_from = 'n_ortho',
              values_fill = 0L) %>%
  nest_by(node) %>%
  mutate(chisq = list(column_to_rownames(data, 'major_category') %>%
                        chisq.test()),
         tidy(chisq)) %>%
  ungroup %>%
  mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
  rowwise %>%
  mutate(post_hoc = list(fisher_posthoc(data, 'major_category'))) %>%
  ungroup


major_expansion_contraction %>%
  filter(fdr < 0.05) %>%
  select(-data:-fdr) %>%
  unnest(post_hoc) %>%
  filter(str_detect(node, 'acer'))


#### Test within minor categories ####
minor_expansion_contraction <- expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction')) %>%
  filter(!is.na(major_category)) %>%
  group_by(node, change, major_category, minor_category) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  pivot_wider(names_from = 'change',
              values_from = 'n_ortho',
              values_fill = 0L) %>%
  nest_by(node, major_category) %>%
  mutate(chisq = list(column_to_rownames(data, 'minor_category') %>%
                        chisq.test()),
         tidy(chisq)) %>%
  ungroup %>%
  group_by(major_category) %>%
  mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
  ungroup %>%
  rowwise %>%
  mutate(post_hoc = list(fisher_posthoc(data, 'minor_category'))) %>%
  ungroup

minor_expansion_contraction %>%
  filter(str_detect(node, '^[a-zA-Z]')) %>%
  filter(fdr < 0.05) %>%
  slice(1) %>%
  unnest(post_hoc)


minor_expansion_contraction %>%
  filter(str_detect(node, 'acer')) %>%
  select(-data, -chisq) %>%
  unnest(post_hoc) %>%
  filter(fdr < 0.05) %>%
  select(-statistic:-fdr)

expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction')) %>%
  filter(str_detect(node, 'acer'),
         minor_category == 'Immune system') %>% 
  select(change, major_category, minor_category, kegg_path_id, name, description) %>%
  distinct %>%
  filter(change == 'contraction')

#### Heatmap ####
minor_expansion_contraction %>%
  filter(str_detect(node, '^[a-z]')) %>%
  select(node, major_category, post_hoc) %>%
  unnest(post_hoc) %>%
  mutate(across(ends_with('p.value'), p.adjust, method = 'fdr'),
         node = str_remove(node, '_.*'),
         node = str_remove(node, '<[0-9]+>')) %>%
  pivot_longer(cols = c(starts_with('contract'), starts_with('expand')),
               names_to = c('direction', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  ggplot(aes(x = log(estimate))) +
  geom_histogram() +
  geom_vline(xintercept = log(1))

minor_expansion_contraction %>%
  # filter(str_detect(node, '^[a-z]')) %>%
  select(node, major_category, post_hoc) %>%
  unnest(post_hoc) %>%
  # mutate(node = str_remove(node, '_.*'),
  #        node = str_remove(node, '<[0-9]+>')) %>%
  group_by(major_category) %>%
  mutate(across(ends_with('p.value'), p.adjust, method = 'fdr')) %>%
  ungroup %>%
  pivot_longer(cols = c(starts_with('contract'), starts_with('expand')),
               names_to = c('direction', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  # filter(major_category == 'Cellular Processes') %>%
  ggplot(aes(y = node, x = interaction(direction, minor_category), 
             fill = log(estimate))) +
  geom_raster() +
  geom_point(data = . %>% filter(p.value < 0.1)) +
  # scale_alpha_discrete(values = c('TRUE' = 1, 'FALSE' = 0.5)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1, angle = 90)) +
  scale_fill_gradient2(midpoint = log(1)) +
  facet_grid(~ major_category, scales = 'free_x')
  

minor_expansion_contraction %>%
  # filter(str_detect(node, '^[a-z]')) %>%
  select(node, major_category, post_hoc) %>%
  unnest(post_hoc) %>%
  # mutate(node = str_remove(node, '_.*'),
  #        node = str_remove(node, '<[0-9]+>')) %>%
  group_by(major_category) %>%
  mutate(across(ends_with('p.value'), p.adjust, method = 'fdr')) %>%
  ungroup %>%
  pivot_longer(cols = c(starts_with('contract'), starts_with('expand')),
               names_to = c('direction', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  filter(p.value < 0.05) %>%
  
  left_join(expansion_contractions %>%
              pivot_longer(cols = -Orthogroup,
                           names_to = 'node',
                           values_to = 'delta') %>%
              left_join(kegg_paths %>%
                          select(-gene_id, -species) %>%
                          distinct(),
                        by = c('Orthogroup')) %>%
              mutate(change = case_when(delta == 0 ~ 'none',
                                        delta > 0 ~ 'expand',
                                        delta < 0 ~ 'contract')) %>%
              select(Orthogroup, node, change, name, major_category, minor_category) %>%
              distinct %>%
              group_by(major_category, minor_category, node, change) %>%
              summarise(pathway = str_c(str_replace_na(name) %>% unique, collapse = ';; '),
                        .groups = 'drop') ,
            by = c('major_category', 'minor_category', 'node', 'direction' = 'change')) %>% View
   
  
expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expand',
                            delta < 0 ~ 'contract')) %>%
  select(Orthogroup, node, change, name, major_category, minor_category) %>%
  distinct %>%
  group_by(major_category, minor_category, node, change) %>%
  summarise(pathway = str_c(str_replace_na(name), collapse = ';; '),
            .groups = 'drop') %>%
  
  filter(node == '<14>',
         major_category == 'Environmental Information Processing',
         minor_category == 'Membrane transport')



  left_join(select(kegg_paths, major_category, minor_category, name) %>%
              distinct %>%
              group_by(major_category, minor_category) %>%
              summarise(pathways = str_c(name, collapse = ';; '),
                        .groups = 'drop'),
            by = c('major_category', 'minor_category')) %>% View




expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction')) %>%
  select(Orthogroup, node, change, name, major_category, minor_category) %>%
  distinct %>%
  group_by(major_category, minor_category, node, change) %>%
  summarise(pathway = str_c(name, collapse = ';; '),
            .groups = 'drop')


#### Make Composition Plots ####
expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction'),
         change = factor(change, levels = c('contraction', 'none', 'expansion'))) %>%
  filter(!is.na(major_category)) %>%
  group_by(node, change, major_category) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  # pivot_wider(names_from = 'change',
  #             values_from = 'n_ortho',
  #             values_fill = 0L) %>%
  filter(str_detect(node, '^[a-zA-Z]')) %>%
  
  ggplot(aes(x = change, y = n_ortho, fill = major_category)) +
  geom_bar(position = "fill", stat = "identity", colour = 'black',
           show.legend = FALSE) +
  facet_grid(node ~ .)


expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction'),
         change = factor(change, levels = c('contraction', 'none', 'expansion'))) %>%
  filter(!is.na(major_category)) %>%
  group_by(node, change, major_category, minor_category) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  # pivot_wider(names_from = 'change',
  #             values_from = 'n_ortho',
  #             values_fill = 0L) %>%
  filter(str_detect(node, '^[a-zA-Z]')) %>%
  
  ggplot(aes(x = change, y = n_ortho, fill = minor_category)) +
  geom_bar(position = "fill", stat = "identity", colour = 'black',
           show.legend = TRUE) +
  facet_grid(node ~ major_category, labeller = as_labeller(function(x) str_replace_all(x, ' ', '\n') %>% str_remove('_shinzato|_fuller') %>% str_remove('<[0-9]+>'))) +
  scale_y_continuous(labels = scales::percent_format(suffix = ''), expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  scale_x_discrete(labels = ~str_to_title(.), expand = c(0, 0)) +
  guides(fill = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  labs(x = NULL, 
       y = 'Orthogroups (%)',
       fill = 'Minor KEGG Category') +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'),
        # panel.background = element_rect(colour = 'black'),
        # plot.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        legend.position = 'bottom')


#### Bayes Model ####
tmp_data <- expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction'),
         change = factor(change, levels = c('contraction', 'none', 'expansion'))) %>%
  filter(!is.na(major_category)) %>%
  group_by(node, change, major_category) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  pivot_wider(names_from = major_category,
              values_from = n_ortho,
              values_fill = 0L) %>%
  mutate(Y = cbind(!!!syms(unique(kegg_paths$major_category))),
         total = rowSums(Y)) %>%
  select(-all_of(unique(kegg_paths$major_category)))


library(brms)

test_model <- brm(Y | trials(total) ~ change + (1 | node),
                  family = multinomial(),
                  prior = prior(normal(0, 10), 
                                class = 'b'), 
                  data = tmp_data,
                  sample_prior = 'yes',
                  chains = 4,
                  cores = 4,
                  iter = 2000,
                  warmup = 1000,
                  backend = 'cmdstanr')

test_model

tmp_data %>%
  select(change) %>%
  distinct %>%
  mutate(total = 1) %>%
  add_epred_draws(test_model, re_formula = NA) %>%
  point_interval() %>%
  ungroup %>%
  
  ggplot(aes(x = change, y = .epred, ymin = .lower, ymax = .upper, 
             colour = .category)) +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(colour = 'Minor Category',
       y = 'Orthogroups (%)',
       x = 'Node Change') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 16),
        axis.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 16),
        legend.text = element_text(colour = 'black', size = 12),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'bottom')


#### Bayes Minor Categories ####
format_brms_data <- function(data, catagories){
  pivot_wider(data, 
              names_from = minor_category, 
              values_from = n,
              values_fill = 0L) %>%
    mutate(Y = cbind(!!!syms(catagories)),
           total = rowSums(Y)) %>%
    select(-all_of(catagories))
}

brms_count_model_fitting <- function(data){
  brm(total ~ change + (1 | node),
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
  brm(Y | trials(total) ~ change + (1 | node),
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
    select(change) %>%
    distinct %>%
    mutate(total = 1) %>%
    add_epred_draws(model, re_formula = NA) %>%
    point_interval() %>%
    ungroup %>%
    
    ggplot(aes(x = change, y = .epred, ymin = .lower, ymax = .upper)) +
    geom_pointrange() +
    # facet_wrap(~ .category) +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(colour = 'Major Category',
         y = 'Orthogroups (#)',
         x = 'Node Change') +
    theme_classic() +
    theme(axis.title = element_text(colour = 'black', size = 16),
          axis.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 16),
          legend.text = element_text(colour = 'black', size = 12),
          panel.background = element_rect(colour = 'black'))
}

composition_plot <- function(data, model){
  data %>%
    select(change) %>%
    distinct %>%
    mutate(total = 1) %>%
    add_epred_draws(model, re_formula = NA) %>%
    point_interval() %>%
    ungroup %>%
    
    ggplot(aes(x = change, y = .epred, ymin = .lower, ymax = .upper, 
               colour = .category)) +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(colour = 'Minor Category',
         y = 'Orthogroups (%)',
         x = 'Node Change') +
    theme_classic() +
    theme(axis.title = element_text(colour = 'black', size = 16),
          axis.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 16),
          legend.text = element_text(colour = 'black', size = 12),
          panel.background = element_rect(colour = 'black'),
          legend.position = 'bottom')
}

minor_changes_data <- expansion_contractions %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'node',
               values_to = 'delta') %>%
  left_join(kegg_paths %>%
              select(-gene_id, -species) %>%
              distinct(),
            by = c('Orthogroup')) %>%
  mutate(change = case_when(delta == 0 ~ 'none',
                            delta > 0 ~ 'expansion',
                            delta < 0 ~ 'contraction'),
         change = factor(change, levels = c('contraction', 'none', 'expansion'))) %>%
  filter(!is.na(major_category)) %>%
  group_by(node, change, major_category, minor_category) %>%
  summarise(n = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  nest(data = -c(major_category)) %>%
  rowwise %>%
  mutate(minor_categories = list(unique(data$minor_category))) %>%
  mutate(data = list(format_brms_data(data, minor_categories)))


minor_changes_data$data[[1]]


test_models <- minor_changes_data %>%
  mutate(count_model = list(brms_count_model_fitting(data)),
         composition_model = list(brms_composition_model_fitting(data)))


test_comp <- test_models %>%
  mutate(count_plot = list(count_plot(data, count_model) + labs(title = major_category)),
         comp_plot = list(composition_plot(data, composition_model) + labs(title = major_category)))

(wrap_plots(test_comp$count_plot, nrow = 1)) / (wrap_plots(test_comp$comp_plot, nrow = 1))

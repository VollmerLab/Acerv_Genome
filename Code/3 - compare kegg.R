##TODO - clean up code to make voting not on a gene by gene level but an orthogroup by orthogroup level
##TODO - Assign each orthogroup a KEGG ID
##TODO - analysis of completeness of pathways by species/major/minor/pathway

#### Libraries ####
library(multcomp)
library(multcompView)
library(KEGGREST)
library(tidyverse)
library(magrittr)
library(broom)
library(emmeans)
library(cmdstanr)
library(brms)
library(tidybayes)
library(patchwork)

#### Read in Data ####
missing_genes <- c('K01697', 'K10150')

gff_file <- microseq::readGFF('../../Bioinformatics/genome_annotation/k2_structuralAnnotations.gff3') %>%
  filter(Type %in% c('gene', 'mRNA')) %>%
  mutate(gene_id = if_else(str_detect(Attributes, 'Parent'),
                           str_extract(Attributes, 'Parent=Acer_[0-9]+') %>% str_remove('Parent='),
                           str_extract(Attributes, 'ID=Acer_[0-9]+') %>% str_remove('ID='))) %>%
  group_by(gene_id) %>%
  filter(n() == 2) %>%
  ungroup %>%
  select(-Seqid, -Source, -Start, -End, -Score, -Strand, -Phase) %>%
  pivot_wider(names_from = Type, values_from = Attributes) %>%
  mutate(gene = str_extract(gene, 'Name=.*') %>% str_remove('Name='),
         mRNA = str_remove(mRNA, str_c('ID=.*', gene_id)) %>% str_remove('^;'))


species_stats <- read_delim('../../Bioinformatics/Phylogenomics/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv', 
                            delim = '\t', n_max = 10) %>%
  rename(param = 1) %>%
  pivot_longer(cols = -param,
               names_to = 'species') %>%
  pivot_wider(names_from = param, values_from = value)

species_stats %>%
  ggplot(aes(y = species, x = `Percentage of unassigned genes`)) +
  geom_col()

species_stats %>%
  filter(species == 'acer') %>%
  select(-species) %>%
  pivot_longer(cols = everything())

kegg_annotations <- read_csv('../../Bioinformatics/Phylogenomics/Annotations/kegg_annotations.csv.gz', 
                             show_col_types = FALSE)

all_orthogroups <- read_delim('../../Bioinformatics/Phylogenomics/Phylogenetic_Hierarchical_Orthogroups/N0.tsv', show_col_types = FALSE) %>%
  dplyr::select(-OG, -`Gene Tree Parent Clade`) %>%
  dplyr::rename(Orthogroup = HOG) %>%
  pivot_longer(cols = -Orthogroup,
               names_to = 'species',
               values_to = 'gene_list') %>%
  filter(!is.na(gene_list)) %>%
  rowwise(Orthogroup, species) %>%
  summarise(gene_id = str_split(gene_list, ', ') %>%
              unlist, 
            .groups = 'drop') %>%
  filter(!is.na(gene_id))

kegg_annotations %>%
  filter(kegg_orthology %in% str_c('ko:', missing_genes)) %>%
  count(species)

duplications <- read_delim('../../Bioinformatics/Phylogenomics/Gene_Duplication_Events/Duplications.tsv', 
                           delim = '\t', show_col_types = FALSE)

kegg_db <- map(c('pathway', 'brite', 'module'), 
               ~keggLink(.x, "ko") %>%
                 enframe(name = 'ko', value = .x)) %>%
  reduce(full_join, by = "ko") %>%
  pivot_longer(cols = -ko,
               names_to = 'kegg_class',
               values_to = 'path_id') %>%
  filter(!is.na(path_id)) %>%
  distinct %>%
  filter((kegg_class == 'pathway' & str_detect(path_id, 'map')) |
           kegg_class != 'pathway')

#### Identify top hit for each gene ####
gene_kegg_vote <- function(data, group_vars){
  #For each geneID from up to top 10 protein hits with KEGG orthology
  # pick modal kegg orthology term
  # If there is a tie then: pick minimum evalue
  # if still a tie then: pick max bitscore
  # if still a tie (e.g. 'Acer_00000041-RA') pick at random
  group_by(data, across(all_of(c(group_vars, 'kegg_orthology')))) %>%
    mutate(n = n()) %>%
    ungroup %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(prop_top_hits = n / n()) %>%
    filter(n == max(n)) %>%
    filter(evalue == min(evalue)) %>%
    filter(bitscore == max(bitscore)) %>%
    sample_n(1) %>%
    ungroup %>%
    select(-n)
}

ortho_kegg_vote <- all_orthogroups %>%
  left_join(kegg_annotations,
            by = c('species', 'gene_id' = 'qseqid')) %>%
  filter(!is.na(kegg_orthology)) %>%
  gene_kegg_vote('Orthogroup') %>%
  select(-species:-uniprot_id)

ortho_kegg_vote %>%
  filter(kegg_orthology %in% str_c('ko:', missing_genes))

ortho_kegg_vote %>%
  ggplot(aes(x = prop_top_hits)) +
  geom_histogram(bins = 30) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = 'Blast hit KO Agreement - must have a KO to include (%)')

species_ortho_vote <- all_orthogroups %>%
  left_join(kegg_annotations,
            by = c('species', 'gene_id' = 'qseqid')) %>%
  select(Orthogroup, species, gene_id) %>%
  distinct %>%
  inner_join(ortho_kegg_vote,
             by = 'Orthogroup') 

species_ortho_vote %>%
  filter(kegg_orthology %in% str_c('ko:', missing_genes)) %>%
  count(species)

species_ortho_vote %>%
  count(species) %>%
  left_join(select(species_stats, species, 'Number of genes in orthogroups'),
            by = 'species') %>%
  rename(total_genes = `Number of genes in orthogroups`,
         kegg_genes = n) %>%
  mutate(prop_kegg = kegg_genes / total_genes) %>%
  
  ggplot(aes(x = prop_kegg, y = species)) +
  geom_col() +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = 'Genes with concensus KEGG orthology (%)')


prop_kegg_model <- species_ortho_vote %>%
  count(species) %>%
  left_join(select(species_stats, species, 'Number of genes in orthogroups'),
            by = 'species') %>%
  rename(total_genes = `Number of genes in orthogroups`,
         kegg_genes = n) %>%
  mutate(prop_kegg = kegg_genes / total_genes) %>%
  glm(prop_kegg ~ species, data = ., family = binomial(link = 'logit'), weights = total_genes)

emmeans(prop_kegg_model, ~species, type = 'response') %>%
  cld(Letters = LETTERS, reversed = TRUE) %>%
  as_tibble() %>%
  mutate(.group = str_trim(.group)) %>%
  mutate(species = fct_reorder(species, prob)) %>%
  ggplot(aes(x = prob, xmin = asymp.LCL, xmax = asymp.UCL,
             y = species)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(xend = 0, yend = species), 
               linetype = 'dashed', colour = 'gray50') +
  geom_pointrange() +
  geom_text(aes(x = Inf, label = .group), hjust = 1) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = 'Genes with concensus KEGG orthology (%)',
       y = NULL) +
  theme_classic()


#### KEGG Pathway Identification ####
if(file.exists('../intermediate_files/kegg_orthogroup_pathways.csv.gz')){
  kegg_paths <- read_csv('../intermediate_files/kegg_orthogroup_pathways.csv.gz', show_col_types = FALSE)
} else {
  kegg_paths <- species_ortho_vote %>%
    filter(!is.na(kegg_orthology)) %>%
    select(Orthogroup, gene_id, species, starts_with('kegg'), prop_top_hits) %>%
    pivot_longer(cols = c(starts_with('kegg'), -kegg_orthology, -kegg_gene),
                 names_prefix = 'kegg_',
                 names_to = 'kegg_type',
                 values_to = 'kegg_path_id') %>%
    filter(!is.na(kegg_path_id)) %>%
    mutate(kegg_path_id = str_split(kegg_path_id, ';;')) %>%
    unnest(kegg_path_id) %>%
    filter(kegg_type == 'pathway') %>%
    
    nest_by(kegg_type, kegg_path_id) %>%
    mutate(kegg_api = possibly(keggGet, otherwise = list(NULL), quiet = FALSE)(kegg_path_id)) %>%
    
    unnest_wider(kegg_api) %>%
    select(-contains(c('ENTRY', 'REFERENCE', 'DBLINKS', 'DISEASE', 'MODULE', 'DRUG', 'KO_PATHWAY'))) %>%
    mutate(name_length = lengths(NAME)) %>%
    rowwise %>%
    mutate(DESCRIPTION = str_c(DESCRIPTION, collapse = ';;;;'),
           DESCRIPTION = if_else(DESCRIPTION == "", NA_character_, DESCRIPTION),
           
           #May not want to do this - related pathways - network??
           REL_PATHWAY = str_c(REL_PATHWAY, collapse = ';;;;'),
           REL_PATHWAY = if_else(REL_PATHWAY == "", NA_character_, REL_PATHWAY)) %>%
    mutate(NAME = case_when(name_length == 2 ~ list(str_c(unlist(NAME), collapse = ': ')), 
                            name_length == 0 ~ list(NA_character_),
                            TRUE ~ list(NAME))) %>%
    ungroup %>%
    select(-name_length) %>%
    unnest(NAME) %>%
    janitor::clean_names() %>%
    separate(class, sep = '; ', into = c('major_category', 'minor_category')) %>%
    unnest(data)
  
  write_csv(kegg_paths, '../intermediate_files/kegg_orthogroup_pathways.csv.gz')
}

#
#### KEGG Compositional Description of Acerv ####
summarized_keggs <- kegg_paths %>%
  filter(major_category != 'Human Diseases') %>%
  group_by(major_category, minor_category, species) %>%
  summarise(n_kegg = n_distinct(kegg_orthology),
            n_gene = n_distinct(gene_id),
            n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') 

summarized_keggs %>%
  filter(species == 'acer') %>%
  pivot_longer(cols = starts_with('n_'),
               names_to = 'what',
               values_to = 'value',
               names_prefix = 'n_') %>%
  ggplot(aes(x = major_category, y = value, fill = minor_category)) +
  # geom_col(position = 'fill') +
  geom_col() +
  facet_wrap(~ what)

#### Get orthogroups unique/duplicated in Acerv ####
unique_acer_kegg <- all_orthogroups %>%
  group_by(Orthogroup) %>%
  filter(n_distinct(species) == 1) %>%
  ungroup %>%
  filter(species == 'acer') %>%
  select(-species) %>%
  distinct %>%
  left_join(kegg_paths, by = c('Orthogroup', 'gene_id')) %>%
  filter(!is.na(kegg_orthology))
#The two genes in one orthgroup with KEGG paths are involved in cocaine addiction...so ignoring

duplicate_acer_kegg <- duplications %>%
  filter(`Species Tree Node` == 'acer') %>%
  select(Orthogroup, starts_with('Genes')) %>%
  pivot_longer(cols = starts_with('Genes')) %>%
  select(-name) %>%
  rowwise(Orthogroup) %>%
  summarise(gene_id = str_split(value, ', ') %>% 
              unlist %>% str_remove('^acer_'),
            .groups = 'drop') %>%
  select(-Orthogroup) %>%
  distinct %>%
  left_join(kegg_paths, by = c('gene_id')) %>%
  filter(!is.na(kegg_orthology)) %>%
  distinct()

all_acer_kegg <- all_orthogroups %>%
  filter(species == 'acer') %>%
  select(-species) %>%
  
  left_join(kegg_paths, by = c('Orthogroup', 'gene_id')) %>%
  filter(!is.na(kegg_orthology)) %>%
  distinct 

#### Compositional Differences Comparing Duplicated orthogroups within Acerv to the background ####
bind_rows(Duplicate = duplicate_acer_kegg,
          Overall = all_acer_kegg,
          .id = 'class') %>%
  mutate(class = factor(class, levels = c('Duplicate', 'Overall'))) %>%
  group_by(class, major_category, minor_category) %>%
  summarise(n_gene = n_distinct(gene_id),
            n_ortho = n_distinct(Orthogroup),
            n_kegg = n_distinct(kegg_orthology),
            .groups = 'drop') %>%
  filter(major_category != 'Human Diseases') %>%
  mutate(major_category = str_replace_all(major_category, ' ', '\n')) %>%
  
  ggplot(aes(x = class, y = n_gene, fill = minor_category)) +
  geom_bar(position = "fill", stat = "identity", colour = 'black',
           show.legend = TRUE) +
  # geom_bar(stat = "identity", colour = 'black') +
  facet_grid(~major_category) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  guides(fill = guide_legend(title.position = 'top', title.hjust = 0.5)) +
  labs(x = NULL, 
       y = 'Genes (%)',
       fill = 'Minor KEGG Category') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        legend.position = 'bottom')
ggsave('../Results/acerv_kegg_composition.png', height = 7, width = 15)

#### Test if there is an overrepresentation of duplicated genes in the major/minor kegg categories compared the composition of Acerv ####
model_program <- cmdstan_model('fisher_test.stan')
run_model <- function(x, n, sample_prior_only){
  stan_data <- list(
    beta_a = 1,
    beta_b = 1,
    n_total_category = n[1],
    n_total_background = n[2],
    n_genes_category = x[1],
    n_genes_background = x[2],
    sample_prior_only = sample_prior_only
  )
  
  model <- model_program$sample(data = stan_data, refresh = 0)
  model$summary() %>%
    as_tibble() %>%
    filter(variable == 'OR')
}

major_duplicate_test <- bind_rows(Duplicate = duplicate_acer_kegg,
                                  Overall = all_acer_kegg,
                                  .id = 'class') %>%
  mutate(class = factor(class, levels = c('Duplicate', 'Overall'))) %>%
  group_by(class, major_category) %>%
  summarise(n_gene = n_distinct(gene_id),
            .groups = 'drop') %>%
  filter(major_category != 'Human Diseases') %>%
  pivot_wider(names_from = 'class',
              values_from = 'n_gene', 
              values_fill = 0L) %>%
  # mutate(across(where(is.integer), ~. + 1L)) %>%
  column_to_rownames('major_category') %T>%
  print %>%
  chisq.test() %T>%
  print

minor_duplicate_test <- bind_rows(Duplicate = duplicate_acer_kegg,
                                  Overall = all_acer_kegg,
                                  .id = 'class') %>%
  mutate(class = factor(class, levels = c('Duplicate', 'Overall'))) %>%
  filter(major_category != 'Human Diseases') %>%
  group_by(class, major_category, minor_category) %>%
  summarise(n_gene = n_distinct(gene_id),
            .groups = 'drop') %>%
  pivot_wider(names_from = 'class',
              values_from = 'n_gene', 
              values_fill = 0L) %>%
  # mutate(across(where(is.integer), ~. + 1L)) %>%
  nest_by(major_category) %>%
  mutate(chisq = list(column_to_rownames(data, 'minor_category') %>%
                        chisq.test()),
         tidy(chisq),
         prop_postHoc = list(mutate(data,
                                    total_dupe = sum(Duplicate),
                                    total_overall = sum(Overall)) %>%
                               rowwise(minor_category) %>%
                               summarise(tidy(prop.test(x = c(Duplicate, Overall), 
                                                        n = c(total_dupe, total_overall))),
                                         .groups = 'drop')),
         
         fisher_postHoc = list(mutate(data, not_dupe = sum(Duplicate) - Duplicate,
                                      not_overall = sum(Overall) - Overall) %>%
                                 mutate(across(where(is.integer), ~. + 1L)) %>%
                                 rename(n_dupe = Duplicate,
                                        n_overall = Overall) %>%
                                 pivot_longer(cols = -minor_category, 
                                              names_pattern = '(.*)_(.*)',
                                              names_to = c('.value', 'type')) %>%
                                 nest(data = -c(minor_category)) %>%
                                 rowwise(minor_category) %>%
                                 summarise(column_to_rownames(data, 'type') %>%
                                             fisher.test() %>% tidy, 
                                           .groups = 'drop') %>%
                                 ungroup),
         
         bayes_postHoc = list(mutate(data,
                                     total_dupe = sum(Duplicate),
                                     total_overall = sum(Overall)) %>%
                                rowwise(minor_category) %>%
                                summarise(run_model(x = c(Duplicate, Overall), 
                                                    n = c(total_dupe, total_overall),
                                                    sample_prior_only = FALSE),
                                          .groups = 'drop'))) %>%
  ungroup

major_duplicate_test #there is a significant difference in the gene composition of major categories of KEGG function in the duplicate genes compared with acerv overall 

minor_duplicate_test #in particular those differences are observed within the environmental information processing, genetic information processing, metabolism, and organismal systems categories. Cellular processes are equivilantly represented between the background and duplicated orthogroups

minor_duplicate_test$fisher_postHoc[[2]] #excess of signalling molecules and dearth of membrane transport proteins
minor_duplicate_test$fisher_postHoc[[3]] #prop/fisher agree: there is an excess of folding/sorting/degradation, and Replication/Repair. Dearth of Transcription & translation
minor_duplicate_test$fisher_postHoc[[4]] #prop/fisher agree: excess of Glycan biosynthesis and metabolism
minor_duplicate_test$fisher_postHoc[[5]] #prop/fisher agree: dearth of circulartory, and excess of development/regeneration and immune

paths_over_or_under_represented_duplicates <- minor_duplicate_test %>%
  filter(p.value < 0.05) %>%
  select(major_category, fisher_postHoc) %>%
  unnest(fisher_postHoc) %>%
  filter(p.value < 0.05) %>%
  select(major_category, minor_category, estimate) %>%
  left_join(duplicate_acer_kegg,
            by = c('major_category',
                   'minor_category')) %>% 
  select(major_category, minor_category, estimate, kegg_path_id, name, description) %>%
  distinct 
  
#### Compositional Differences among species ####
format_brms_data <- function(data, category_name, values_name, catagories){
  pivot_wider(data, 
              names_from = all_of(category_name), 
              values_from = all_of(values_name),
              values_fill = 0L) %>%
    mutate(Y = cbind(!!!syms(catagories))) %>%
    select(-all_of(catagories)) %>%
    mutate(total = rowSums(Y))
}

brms_count_model_fitting <- function(data){
  brm(total ~ (1 | species),
      family = poisson(), 
      data = data,
      sample_prior = 'yes',
      chains = 4,
      cores = 4,
      iter = 2000,
      warmup = 1000,
      backend = 'cmdstanr')
}

brms_composition_model_fitting <- function(data){
  brm(Y | trials(total) ~ (1 | species),
      family = multinomial(),
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
    select(species) %>%
    distinct %>%
    add_epred_draws(model, re_formula = ~(1 | species)) %>%
    point_interval() %>%
    ungroup %>%
    mutate(species = fct_reorder(species, .epred)) %>%
    
    ggplot(aes(y = species, x = .epred, xmin = .lower, xmax = .upper)) +
    
    geom_pointrange() +
    # facet_wrap(~ .category) +
    scale_x_continuous(labels = scales::comma_format()) +
    labs(colour = 'Major Category',
         x = 'Orthogroups (#)',
         y = NULL) +
    theme_classic() +
    theme(axis.title = element_text(colour = 'black', size = 16),
          axis.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 16),
          legend.text = element_text(colour = 'black', size = 12),
          panel.background = element_rect(colour = 'black'))
}

composition_plot <- function(data, model){
  data %>%
    select(species) %>%
    distinct %>%
    mutate(total = 1) %>%
    add_epred_draws(model, re_formula = ~(1 | species)) %>%
    point_interval() %>%
    ungroup %>%
    
    ggplot(aes(y = species, x = .epred, xmin = .lower, xmax = .upper, 
               colour = .category)) +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(colour = 'Minor Category',
         y = NULL,
         x = 'Orthogroups (%)') +
    guides(colour = guide_legend(title.position = 'top', title.hjust = 0.5)) +
    theme_classic() +
    theme(axis.title = element_text(colour = 'black', size = 16),
          axis.text = element_text(colour = 'black', size = 12),
          legend.title = element_text(colour = 'black', size = 16),
          legend.text = element_text(colour = 'black', size = 12),
          panel.background = element_rect(colour = 'black'),
          legend.position = 'bottom')
}

majorCat_composition_data <- kegg_paths %>%
  filter(major_category != 'Human Diseases') %>%
  group_by(major_category, species) %>%
  summarise(n_kegg = n_distinct(kegg_orthology),
            n_gene = n_distinct(gene_id),
            n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  select(-n_kegg, -n_gene) %>%
  format_brms_data('major_category', 'n_ortho', unique(summarized_keggs$major_category))
  
majorCat_count_model <- brms_count_model_fitting(majorCat_composition_data)
count_plot(majorCat_composition_data, majorCat_count_model)

majorCat_composition_model <- brms_composition_model_fitting(majorCat_composition_data)
composition_plot(majorCat_composition_data, majorCat_composition_model) + labs(colour = 'Major Category')

job::job({
  minorCat_models <- summarized_keggs %>%
    nest_by(major_category) %>%
    mutate(minor_categories = list(unique(data$minor_category))) %>%
    mutate(data = list(select(data, -n_kegg, -n_gene) %>%
                         format_brms_data('minor_category', 'n_ortho', minor_categories))) %>%
    mutate(count_model = list(brms_count_model_fitting(data)),
           composition_model = list(brms_composition_model_fitting(data))) %>%
    mutate(count_plots = list(count_plot(data, count_model) + labs(x = str_c(major_category, ': Orthogroups (#)'))),
           comp_plots = list(composition_plot(data, composition_model) + labs(x = str_c(major_category, ': Orthogroups (%)')))) %>%
    ungroup()
})

minorCat_models$count_plots[[4]]

#### KEGG Pathway completeness ####


kegg_paths %>%
  group_by(kegg_type, kegg_path_id, major_category, minor_category,
           species) %>%
  summarise(n_keggs = n_distinct(kegg_orthology),
            .groups = 'drop')


path_complete_data <- full_join(kegg_paths %>%
            group_by(kegg_type, kegg_path_id, major_category, minor_category,
                     species) %>%
            summarise(n_keggs = n_distinct(kegg_orthology),
                      .groups = 'drop'),
          
          kegg_paths %>%
            group_by(kegg_type, kegg_path_id, major_category, minor_category) %>%
            summarise(total_keggs = n_distinct(kegg_orthology),
                      .groups = 'drop'),
          by = c("kegg_type", "kegg_path_id", "major_category", "minor_category")) %>%
  mutate(prop_complete = n_keggs / total_keggs) %>%
  filter(major_category != 'Human Diseases')

full_pathway_data <- expand_grid(kegg_path_id = unique(path_complete_data$kegg_path_id),
                                 species = unique(path_complete_data$species)) %>%
  left_join(select(path_complete_data, kegg_path_id, 
                   major_category, minor_category, total_keggs) %>%
              distinct,
            by = 'kegg_path_id') %>%
  left_join(select(path_complete_data, -kegg_type, 
                   -major_category, -minor_category,
                   -total_keggs),
            by = c('kegg_path_id', 'species')) %>%
  mutate(n_keggs = if_else(is.na(n_keggs), 0L, n_keggs),
         prop_complete = if_else(is.na(prop_complete), 0, prop_complete))

full_pathway_data %>%
  group_by(major_category, species) %>%
  summarise(mean_complete = mean(prop_complete),
            se_complete = sd(prop_complete) / sqrt(n()),
            n_kegg = sum(n_keggs),
            total_kegg = sum(total_keggs),
            .groups = 'drop') %>% #count(species)
  filter(species %in% c('acer', 'adig_shinzato', 'spis', 'nvec')) %>%
  
  ggplot(aes(y = species, x = mean_complete, xmin = mean_complete - se_complete, xmax = mean_complete + se_complete)) +
  geom_pointrange() + 
  # geom_point(aes(x = n_kegg / total_kegg), colour = 'red') +
  facet_wrap(~ major_category)

ggplot(data = tibble(x = -4:4), aes(x = x)) +
  geom_function(fun = dnorm, args = list(mean = 0, sd = 1), size = 1) +
  theme_classic()

ggplot(data = tibble(x = 0:4), aes(x = x)) +
  geom_function(fun = dgamma, args = list(shape = 2, rate = 1), size = 1, colour = 'black') +
  geom_function(fun = dgamma, args = list(shape = 4, rate = 2), size = 1, colour = 'red') +
  theme_classic()

#STOP HERE - MOVE TO 4

if(!file.exists('../intermediate_files/kegg_pathway_completeness.rds')){
  write_csv(full_pathway_data, '../intermediate_files/pathdata_for_discovery.csv')
  #Then run util/fit_path_completeness_model_discovery.R on Discovery to fit the model 
} else {
  path_complete_model <- read_rds('../intermediate_files/kegg_pathway_completeness_fixed.rds')
  
  major_cat_emmeans <- read_rds('../intermediate_files/kegg_pathway_completeness_fixed_major.rds.gz')
  minor_cat_emmeans <- read_rds('../intermediate_files/kegg_pathway_completeness_fixed_minor.rds.gz')
}

#Need to shift initial post-processing onto Discovery

plot(path_complete_model)
pp_check(path_complete_model)

summary(path_complete_model) #

## Plot Major Categories
summary(major_cat_emmeans, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  ggplot(aes(y = species, x = prob, xmin = lower.HPD, xmax = upper.HPD)) +
  geom_pointrange() +
  
  geom_point(data = full_pathway_data %>%
               group_by(species, major_category) %>%
               summarise(prop_complete = mean(prop_complete),
                         .groups = 'drop'),
             inherit.aes = FALSE, colour = 'red',
             aes(x = prop_complete, y = species)) +
  geom_point(data = full_pathway_data %>%
               group_by(species, major_category) %>%
               summarise(total_keggs = sum(total_keggs),
                         n_keggs = sum(n_keggs),
                         .groups = 'drop') %>%
               mutate(prop_complete = n_keggs / total_keggs),
             inherit.aes = FALSE, colour = 'blue',
             aes(x = prop_complete, y = species)) +
  
  facet_wrap(~ major_category)



## Plot random set of minor categories 
minor_cats_plot <- sample(unique(full_pathway_data$minor_category), 12)

summary(minor_cat_emmeans, type = 'response') %>%
  broom::tidy(conf.int = TRUE) %>%
  filter(minor_category %in% minor_cats_plot) %>%
  ggplot(aes(y = species, x = prob, xmin = lower.HPD, xmax = upper.HPD)) +
  geom_pointrange(position = position_dodge(0.5)) +
  
  geom_point(data = full_pathway_data %>%
               filter(minor_category %in% minor_cats_plot) %>%
               group_by(species, major_category, minor_category) %>%
               summarise(prop_complete = mean(prop_complete),
                         .groups = 'drop'),
             inherit.aes = FALSE, colour = 'red',
             aes(x = prop_complete, y = species)) +
  geom_point(data = full_pathway_data %>%
               filter(minor_category %in% minor_cats_plot) %>%
               group_by(species, major_category, minor_category) %>%
               summarise(total_keggs = sum(total_keggs),
                         n_keggs = sum(n_keggs),
                         .groups = 'drop') %>%
               mutate(prop_complete = n_keggs / total_keggs),
             inherit.aes = FALSE, colour = 'blue',
             aes(x = prop_complete, y = species)) +
  
  facet_wrap(major_category ~ minor_category)




#Plot 12 random paths
paths_plot <- sample(unique(full_pathway_data$kegg_path_id), 12)
full_pathway_data %>%
  select(species, major_category, minor_category, kegg_path_id) %>%
  distinct %>%
  filter(kegg_path_id %in% paths_plot) %>%
  mutate(total_keggs = 1) %>%
  add_epred_draws(path_complete_model, 
                  re_formula = ~(1 + species | major_category / minor_category / kegg_path_id)) %>%
  point_interval() %>%
  ungroup %>%
  # mutate(species = fct_reorder(species, .epred)) %>%
  ggplot(aes(y = species, x = .epred, xmin = .lower, xmax = .upper)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_point(data = full_pathway_data %>%
               filter(kegg_path_id %in% paths_plot),
             inherit.aes = FALSE, colour = 'red',
             aes(x = prop_complete, y = species)) +
  facet_wrap(major_category + minor_category ~ kegg_path_id)




#### Compare Species for particular pathways ####
cysteine_metabolism_related_paths <- c('map00260', 'map00270', 'map01100', 'map01110', 'map01230', 'map00920', 'map01120') 

full_pathway_data %>%
  select(species, major_category, minor_category, kegg_path_id) %>%
  distinct %>%
  filter(kegg_path_id %in% str_c('path:', cysteine_metabolism_related_paths)) %>%
  mutate(total_keggs = 1) %>%
  add_epred_draws(path_complete_model, 
                  re_formula = ~(1 + species | major_category / minor_category / kegg_path_id)) %>%
  point_interval() %>%
  ungroup %>%
  left_join(select(kegg_paths, kegg_path_id, name, description) %>%
              distinct,
            by = 'kegg_path_id') %>%
  # mutate(species = fct_reorder(species, .epred)) %>%
  ggplot(aes(y = species, x = .epred, xmin = .lower, xmax = .upper)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_point(data = full_pathway_data %>%
               filter(kegg_path_id %in% str_c('path:', cysteine_metabolism_related_paths)) %>%
               left_join(select(kegg_paths, kegg_path_id, name, description) %>%
                           distinct,
                         by = 'kegg_path_id'),
             inherit.aes = FALSE, colour = 'red',
             aes(x = prop_complete, y = species)) +
  facet_wrap( ~ name)

full_pathway_data %>%
  select(species, major_category, minor_category, kegg_path_id) %>%
  distinct %>%
  filter(str_detect(minor_category, 'Immune')) %>%
  mutate(total_keggs = 1) %>%
  add_epred_draws(path_complete_model, 
                  re_formula = ~(1 + species | major_category / minor_category / kegg_path_id)) %>%
  point_interval() %>%
  ungroup %>%
  left_join(select(kegg_paths, kegg_path_id, name, description) %>%
              distinct,
            by = 'kegg_path_id') %>%
  # mutate(species = fct_reorder(species, .epred)) %>%
  ggplot(aes(y = species, x = .epred, xmin = .lower, xmax = .upper)) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_point(data = full_pathway_data %>%
               filter(str_detect(minor_category, 'Immune')) %>%
               left_join(select(kegg_paths, kegg_path_id, name, description) %>%
                           distinct,
                         by = 'kegg_path_id'),
             inherit.aes = FALSE, colour = 'red',
             aes(x = prop_complete, y = species)) +
  facet_wrap( ~ name)

#















































#### Compare Number of KEGGs in Orthogroups across Taxa ####
# vote what each orthogroup is based on # of sp

all_orthogroups %>%
  # filter(Orthogroup == 'N0.HOG0000020') %>%
  filter(Orthogroup == 'N0.HOG0000050') %>%
  left_join(kegg_annotations,
            by = c('species', 'gene_id' = 'qseqid')) %>%
  filter(!is.na(kegg_orthology))

full_join(all_orthogroups,
          voted_kegg,
          by = c('species', 'gene_id' = 'qseqid')) %>%
  filter(!is.na(kegg_orthology)) %>%
  # filter(Orthogroup == 'N0.HOG0000020')
  
  filter(Orthogroup == 'N0.HOG0000050')



all_orthogroups %>%
  group_by(Orthogroup) %>%
  filter(n_distinct(species) >= 23,
         any(species == 'acer')) %>%
  
  full_join(voted_kegg,
            by = c('species', 'gene_id' = 'qseqid')) %>%
  filter(!is.na(kegg_orthology)) %>%
  filter(any(species == 'acer')) %>%
  ungroup %>%
  
  filter(Orthogroup == 'N0.HOG0000785')
  
  

# all_orthogroups %>%

# 1. just ahve blast hit - is there more info in blast hit not just top hit
#     ID - vote based on modal kegg 

# 2. don't collapse within a species - use orthogroup and for anything with orthogroup do voting - majority rule

# 3. within species vote and then correct with 


# compare option 1 v 2 for Acerv "gap"


kegg_paths_described %>%
  filter(str_detect(minor_category, 'Immune')) %>%
  select(species, name, n_keggs, total_keggs) %>%
  pivot_wider(names_from = 'species', values_from = 'n_keggs') %>%
  filter(str_detect(name, 'Toll'))

kegg_paths_described %>%
  filter(str_detect(minor_category, 'Immune')) %>%
  
  ggplot(aes(x = prop_complete, y = species)) +
  geom_boxplot()


kegg_paths_described %>%
  filter(str_detect(minor_category, 'Immune')) %>% 
  filter(prop_complete < 0.5)



kegg_paths_described %>%
  filter(str_detect(minor_category, 'Amino acid metabolism')) %>%
  select(species, name, n_keggs, total_keggs) %>%
  pivot_wider(names_from = 'species', values_from = 'n_keggs')




kegg_paths_described %>%
  filter(str_detect(minor_category, 'Amino acid metabolism')) %>%
  ggplot(aes(x = prop_complete, y = species)) +
  geom_boxplot()



















all_orthogroups %>%
  filter(species == 'acer') %>%
  mutate(gene_id = str_remove(gene_id, '-RA')) %>%
  full_join(gff_file,
            by = 'gene_id') %>%
  filter(is.na(Orthogroup)) %>%
  select(gene_id) %>%
  left_join(kegg_annotations,
            by = c('gene_id' = 'qseqid')) %>% 
  filter(!is.na(swissprot_name))

  
filter(kegg_annotations, is.na(kegg_orthology))
2617226/3976874


anti_join(kegg_annotations, species_ortho_vote,
          by = c('species', 'qseqid' = 'gene_id'))

species_ortho_vote %>%
  filter(gene_id == 'rna-XM_001628103.2')

all_orthogroups %>%
  filter(gene_id == 'rna-XM_001628103.2')


gene_swissprot_vote <- function(data, group_vars){
  #For each geneID from up to top 10 protein hits with KEGG orthology
  # pick modal kegg orthology term
  # If there is a tie then: pick minimum evalue
  # if still a tie then: pick max bitscore
  # if still a tie (e.g. 'Acer_00000041-RA') pick at random
  group_by(data, across(all_of(c(group_vars, 'swissprot_id')))) %>%
    mutate(n = n()) %>%
    ungroup %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(prop_top_hits = n / n()) %>%
    filter(n == max(n)) %>%
    filter(evalue == min(evalue)) %>%
    filter(bitscore == max(bitscore)) %>%
    sample_n(1) %>%
    ungroup %>%
    select(-n)
}

group_vars <- c('Orthogroup')

tmp <- all_orthogroups %>%
  left_join(kegg_annotations,
            by = c('species', 'gene_id' = 'qseqid')) %>%
  gene_swissprot_vote('Orthogroup') %>%
  select(-species, -gene_id)




species_ortho_vote

swiss_species_ortho_vote <- all_orthogroups %>%
  left_join(kegg_annotations,
            by = c('species', 'gene_id' = 'qseqid')) %>%
  select(Orthogroup, species, gene_id) %>%
  distinct %>%
  inner_join(tmp,
             by = 'Orthogroup') 

swiss_species_ortho_vote %>%
  count(species) %>%
  left_join(select(species_stats, species, 'Number of genes in orthogroups'),
            by = 'species') %>%
  rename(total_genes = `Number of genes in orthogroups`,
         kegg_genes = n) %>%
  mutate(prop_kegg = kegg_genes / total_genes) %>%
  
  ggplot(aes(x = prop_kegg, y = species)) +
  geom_col() +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = 'Genes with concensus KEGG orthology (%)')

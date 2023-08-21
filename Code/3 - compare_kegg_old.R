##Goals
# 1 - use the orthogroup assignments to get KEGG genes for all ACERV & other genes (where possible)
# 2 - characterize annotations of ACERV compare with other species
#       - amount with kegg annotation
#       - breakdown of kegg pathways



#### Libraries ####
library(multcomp)
library(KEGGREST)
library(tidyverse)
library(broom)
library(emmeans)
library(patchwork)
library(ggmosaic)
library(chisq.posthoc.test)

#### Functions ####
get_kegg_path <- function(pathways){
  tibble(kegg_api = keggGet(unlist(str_split(pathways, ';;')))) %>%
    unnest_wider(kegg_api) %>%
    select(contains(c('ENTRY', 'NAME', 'CLASS'))) %>%
    mutate(name_length = lengths(NAME)) %>%
    rowwise %>%
    mutate(NAME = case_when(name_length == 2 ~ list(str_c(unlist(NAME), collapse = ': ')), 
                            name_length == 0 ~ list(NA_character_),
                            TRUE ~ list(NAME))) %>%
    ungroup %>%
    select(-name_length) %>%
    unnest(NAME) %>%
    janitor::clean_names() %>%
    separate(class, sep = '; ', into = c('major_category', 'minor_category')) 
}

process_posthoc <- function(post){
  post %>%
    as_tibble() %>%
    pivot_longer(cols = -c(Dimension, Value),
                 names_to = 'species',
                 values_to = 'value',
                 values_transform = ~str_remove(., '\\*') %>% as.numeric) %>%
    pivot_wider(names_from = Value,
                values_from = value) %>%
    rename(fdr = 'p values')
}

#### Data ####
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
  reframe(gene_id = str_split(gene_list, ', ') %>%
              unlist) %>%
  filter(!is.na(gene_id))

#### Identify each Orthogroup ####
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
  ggplot(aes(x = prop_top_hits)) +
  geom_histogram(bins = 30) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = 'Blast hit KO Agreement - must have a KO to include (%)')

#### Get Kegg Pathway data ####
kegg_pathways <- ortho_kegg_vote %>%
  select(kegg_orthology, kegg_pathway) %>%
  filter(!is.na(kegg_pathway)) %>%
  distinct %>%
  rowwise(kegg_orthology) %>%
  reframe(kegg_path = str_split(kegg_pathway, ';;') %>% unlist) %>%
  distinct %>%
  nest(data = c(kegg_orthology)) %>%
  mutate(grouping = row_number() %/% 10) %>%
  group_by(grouping) %>%
  reframe(kegg_paths = str_c(kegg_path, collapse = ';;'),
          data = list(data)) %>%
  rowwise %>%
  mutate(kegg_pathways = list(possibly(get_kegg_path, otherwise = list(NULL), quiet = FALSE)(kegg_paths))) %>%
  ungroup %>%
  select(data, kegg_pathways) %>%
  unnest(c(data, kegg_pathways)) %>%
  unnest(data) %>%
  distinct %>%
  filter(!str_detect(major_category, 'Human')) %>%
  mutate(name = str_c(name, ' (', entry, ')')) %>%
  select(-entry) %>%
  
  rename(path = name,
         major = major_category,
         minor = minor_category) %>%
  group_by(kegg_orthology, minor, major) %>%
  summarise(all_paths = str_c(path, collapse = '; '),
            .groups = 'drop') %>%
  mutate(minor = str_c(major, minor, sep = ' - '), .keep = 'unused') %>%
  arrange(minor) %>%
  
  pivot_wider(names_from = minor,
              values_from = all_paths,
              names_vary = 'slowest')


select(all_orthogroups, Orthogroup) %>%
  distinct %>%
  full_join(ortho_kegg_vote,
            by = 'Orthogroup') %>%
  full_join(kegg_pathways, 
          by = 'kegg_orthology') %>%
  select(-kegg_pathway, -kegg_brite, -kegg_module) %>%
  arrange(str_extract(Orthogroup, '[0-9]+') %>% as.integer) %>%
  write_csv('../intermediate_files/kegg_gene_pathways.csv.gz')

#### Join KEGGS to species ####
species_kegg <- all_orthogroups %>%
  left_join(ortho_kegg_vote,
            by = 'Orthogroup') %>%
  left_join(kegg_pathways, 
            by = 'kegg_orthology')
  
#### Just Acer KEGGS ####
species_kegg %>%
  group_by(Orthogroup) %>%
  filter(n_distinct(species) == 1) %>%
  filter(species == 'acer') %>%
  ungroup %>%
  filter(!is.na(kegg_gene)) %>%
  select(-Orthogroup:-gene_id, -prop_top_hits) %>%
  distinct %>%
  janitor::remove_empty(which = 'cols')


#### Check Cystine biosynthesis ####
shouldnt_have <- c('K01697', 'K10150')
may_have <- c('K00641', 'K00651', 'K10764')

species_kegg %>%
  filter(str_detect(kegg_orthology, str_c(c(shouldnt_have, may_have), collapse = '|'))) %>%
  select(species, kegg_gene, kegg_orthology) %>%
  distinct %>%
  mutate(exists = TRUE) %>%
  pivot_wider(names_from = 'species',
              values_from = 'exists',
              values_fill = FALSE) 


summary_stats_kegg <- species_kegg %>%
  group_by(species) %>%
  summarise(n_gene = n_distinct(gene_id),
            n_gene_w_kegg = n_distinct(gene_id[!is.na(kegg_orthology)]),
            prop_gene = n_gene_w_kegg / n_gene, 
            
            unique_kegg = n_distinct(kegg_orthology),
            
            n_ortho = n_distinct(Orthogroup),
            n_ortho_w_kegg = n_distinct(Orthogroup[!is.na(kegg_orthology)]),
            prop_ortho = n_ortho_w_kegg / n_ortho,
            
            across(all_of(colnames(kegg_pathways)[-1]), ~sum(!is.na(.)))) %>%
  mutate(is_acroporid = (str_detect(species, '^a') & species != 'aten'))

#### Prop Genes with KEGG ####
gene_prop_kegg <- glm(cbind(n_gene_w_kegg, n_gene - n_gene_w_kegg) ~ species, 
                       data = summary_stats_kegg,
                       family = binomial(link = 'logit')) 
car::Anova(gene_prop_kegg)

emmeans(gene_prop_kegg, ~species, type = 'response') %>% 
  cld(Letters = LETTERS, alpha = 0.05) %>%
  tidy %>%
  mutate(is_acroporid = (str_detect(species, '^a') & species != 'aten')) %>%
  mutate(.group = str_trim(.group),
         species = fct_reorder(species, prob)) %>%
  
  ggplot(aes(y = species, x = prob, colour = is_acroporid)) +
  geom_segment(xend = 0, aes(yend = species), linetype = 'dashed') +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  geom_text(aes(x = Inf, label = .group), hjust = 1, size = 2.5) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
  scale_y_discrete(labels = ~str_remove(., '_(.*)')) +
  guides(colour = 'none') +
  labs(x = 'KEGG Annotated Genes (%)',
       y = NULL) +
  theme_classic()


#### Prop Orthogroups with KEGG ####
ortho_prop_kegg <- glm(cbind(n_ortho_w_kegg, n_ortho - n_ortho_w_kegg) ~ species, 
    data = summary_stats_kegg,
    family = binomial(link = 'logit')) 
car::Anova(ortho_prop_kegg)

emmeans(ortho_prop_kegg, ~species, type = 'response') %>% 
  cld(Letters = LETTERS, alpha = 0.05) %>%
  tidy %>%
  mutate(is_acroporid = (str_detect(species, '^a') & species != 'aten')) %>%
  mutate(.group = str_trim(.group),
         species = fct_reorder(species, prob)) %>%
  
  ggplot(aes(y = species, x = prob, colour = is_acroporid)) +
  geom_segment(xend = 0, aes(yend = species), linetype = 'dashed') +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  geom_text(aes(x = Inf, label = .group), hjust = 1, size = 2) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
  scale_y_discrete(labels = ~str_remove(., '_(.*)')) +
  guides(colour = 'none') +
  labs(x = 'KEGG Annotated Orthogroups (%)',
       y = NULL) +
  theme_classic()
  
summary_stats_kegg %>%
  mutate(species = fct_reorder(species, unique_kegg)) %>%
  
  ggplot(aes(y = species, x = unique_kegg)) +
  geom_point()

#### Composition of KEGG pathways ####
species_paths <- species_kegg %>%
  pivot_longer(cols = all_of(colnames(kegg_pathways)[-1]),
               names_to = 'major_minor',
               values_to = 'pathway',
               values_drop_na = TRUE) %>%
  # select(gene_id, pathway) %>%
  group_by(across(Orthogroup:major_minor)) %>%
  reframe(pathway = str_split(pathway, ';') %>% unlist %>% str_trim) %>%
  select(species, Orthogroup, gene_id, kegg_gene, kegg_orthology, major_minor, pathway) %>%
  separate(major_minor, 
           sep = ' - ', 
           into = c('major', 'minor'), 
           extra = 'merge') %>%
  filter(!minor %in% c('Cellular community - prokaryotes',
                       'Information processing in viruses'))

path_matrix <- species_paths %>%
  filter(str_detect(species, '^a'),
         species != 'aten') %>%
  group_by(species, major, minor, pathway) %>%
  summarise(n_gene = n_distinct(gene_id),
            .groups = 'drop') %>%
  pivot_wider(names_from = species,
              values_from = n_gene,
              values_fill = 0) %>%
  mutate(pathway = str_c(major, minor, pathway, sep = ';;'), .keep = 'unused', .before = everything()) %>%
  column_to_rownames('pathway')

chisq.test(path_matrix, simulate.p.value = FALSE)
path_posthoc <- chisq.posthoc.test(path_matrix, method = 'fdr') %>%
  process_posthoc %>%
  filter(fdr < 0.05) %>%
  separate(Dimension, into = c('major', 'minor', 'pathway'), sep = ';;') 
# count(major, minor, pathway) %>%
# arrange(-n)


#### Do same as above but with different groupings ####
summarize_pathways <- function(data, count_var, group_var){
  group_by(data, !!!syms(unique(c('species', group_var)))) %>%
    summarise(n = n_distinct(!!sym(count_var)),
              .groups = 'drop') %>%
    pivot_wider(names_from = species,
                values_from = n,
                values_fill = 0) %>%
    mutate(pathway = str_c(!!!syms(str_subset(group_var, 'species', negate = TRUE)), sep = ';;'), .keep = 'unused', .before = everything()) %>%
    column_to_rownames('pathway')
}

tst <- expand_grid(counting = c('gene_id', 'Orthogroup', 'kegg_orthology'),
                   all_levels = list('major', c('major', 'minor'), c('major', 'minor', 'pathway')),
                   species_inclusion = c('all', 'acropora')) %>%
  rowwise(counting, all_levels, species_inclusion) %>%
  mutate(level = all_levels[length(all_levels)],
         data = if_else(species_inclusion == 'acropora', 
                        list(filter(species_paths, str_detect(species, '^a'),
                                    species != 'aten')), 
                        list(species_paths))) %>%
  mutate(data = list(summarize_pathways(data, counting, all_levels))) %>%
  mutate(chisq = list(chisq.test(data)),
         tidy(chisq),
         posthoc = list(chisq.posthoc.test(data, method = 'fdr') %>%
                          process_posthoc %>%
                          filter(fdr < 0.05) %>%
                          separate(Dimension, into = c('major', 'minor', 'pathway'), sep = ';;')))



#### All KEGG Orthologs ####
kegg_ortholog_matrix <- species_paths %>%
  filter(str_detect(species, '^a'),
         species != 'aten') %>%
  group_by(kegg_gene, kegg_orthology, species) %>%
  summarise(n_gene = n_distinct(gene_id),
            .groups = 'drop_last') %>%
  filter(min(n_gene) > 5) %>%
  ungroup %>%
  pivot_wider(names_from = species,
              values_from = n_gene,
              values_fill = 0) %>%
  # filter(!if_any(where(is.numeric), ~. == 0)) %>%
  mutate(kegg_gene = str_c(kegg_gene, kegg_orthology, sep = ';;'), .keep = 'unused', .before = everything()) %>%
  column_to_rownames('kegg_gene')

dim(kegg_ortholog_matrix)

chisq.test(kegg_ortholog_matrix, simulate.p.value = FALSE)
path_posthoc <- chisq.posthoc.test(kegg_ortholog_matrix, method = 'fdr') %>%
  process_posthoc %>%
  filter(fdr < 0.05)

species_paths %>%
  inner_join(filter(path_posthoc, species == 'acer'),
             by = c('species', 'major', 'minor', 'pathway')) %>%
  select(-Residuals, -fdr) %>%
  pivot_wider(names_from = c(major, minor),
              values_from = pathway,
              names_sep = ';;',
              values_fn = ~str_c(., sep = ';;')) %>%
  group_by(across(c(everything(), -gene_id, -Orthogroup))) %>%
  summarise(n_gene = n_distinct(gene_id),
            n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  filter(!is.na(`Organismal Systems;;Immune system`)) %>% 
  janitor::remove_empty(which = 'cols') %>%
  View
  pull(kegg_gene)

#### Pathway Completeness ####
complete_pathways <- map(c('pathway'), 
                         ~keggLink(.x, "ko") %>%
                           enframe(name = 'ko', value = .x)) %>%
  reduce(full_join, by = "ko") %>%
  distinct %>%
  filter(str_detect(pathway, 'map')) %>%
  group_by(pathway) %>%
  summarise(total_kegg = n_distinct(ko)) %>%
  mutate(pathway = str_remove(pathway, '^path:'))

complete_cnidarian_pathways <- species_paths %>%
  # select(major_minor, pathway, kegg_gene, kegg_orthology) %>%
  # distinct %>%
  group_by(major, minor, pathway) %>%
  summarise(total_species = n_distinct(species),
            total_cnidarian = n_distinct(kegg_gene),
            .groups = 'drop')

complete_cnidarian_pathways %>%
  arrange(total_cnidarian) %>%
  filter(total_species != 24)

species_pathway_completeness <- species_paths %>%
  group_by(species, major, minor, pathway) %>%
  summarise(n_kegg = n_distinct(kegg_gene),
            n_ortho = n_distinct(Orthogroup),
            n_gene = n_distinct(gene_id),
            .groups = 'drop') %>%
  left_join(complete_cnidarian_pathways,
            by = c('major', 'minor', 'pathway')) %>%
  mutate(path_id = str_extract(pathway, 'map[0-9]+')) %>%
  full_join(expand_grid(complete_pathways, species = unique(species_paths$species)),
            by = c('path_id' = 'pathway', 'species')) %>%
  select(-path_id) %>%
  
  mutate(percent_species = n_kegg / total_cnidarian,
         percent_cnidarian = total_cnidarian / total_kegg) %>%
  filter(total_species == 24 & total_cnidarian > 5) %>%
  
  select(species, major, minor, pathway, 
         total_cnidarian, n_kegg) %>%
  mutate(major_minor_path = str_c(major, minor, pathway, sep = ';;'),
         .keep = 'unused', .after = 'species')
write_csv(species_pathway_completeness, '../intermediate_files/pathway_completeness.csv.gz')

#Run this model on Discovery #
if(file.exists('../intermediate_files/species_major.csv.gz')){
  major_response <- read_csv('../intermediate_files/species_major.csv.gz', show_col_types = FALSE)
  minor_response <- read_csv('../intermediate_files/species_minor.csv.gz', show_col_types = FALSE)
  path_response <- read_csv('../intermediate_files/species_path.csv.gz', show_col_types = FALSE)
} else {
  library(tidyverse)
  library(brglm2)
  library(emmeans)
  library(broom)
  
  species_pathway_completeness <- read_csv('/scratch/j.selwyn/pathway_completeness.csv.gz') %>%
    separate(major_minor_path, into = c('major', 'minor', 'path'), sep = ';;', remove = FALSE) 
  
  if(file.exists('/scratch/j.selwyn/pathway_completeness_model.rds.gz')){
    completeness_model <- read_rds('/scratch/j.selwyn/pathway_completeness_model.rds.gz')
  } else {
    completeness_model <- glm(cbind(n_kegg, total_cnidarian - n_kegg) ~ species * major_minor_path, 
                              data = species_pathway_completeness,
                              family = binomial(link = 'logit'),
                              method = "brglm_fit")
    write_rds(completeness_model, '/scratch/j.selwyn/pathway_completeness_model.rds.gz')
  }
  
  car::Anova(completeness_model)
  
  minor_cats <- species_pathway_completeness %>%
    select(major_minor_path, major, minor, path) %>%
    distinct
  
  major_cats <- minor_cats %>%
    select(major, minor) %>%
    distinct
  
  completeness_grid <- ref_grid(completeness_model, ~species * major_minor_pathway) %>%
    add_grouping('minor', 'major_minor_path', 
                 factor(minor_cats$minor)) %>%
    add_grouping('major', 'minor',
                 factor(major_cats$major))
  
  major_response <- emmeans(completeness_grid, ~species * major, type = 'response')
  write_csv(tidy(major_response, conf.int = TRUE), '/scratch/j.selwyn/species_major.csv.gz')
  
  minor_response <- emmeans(completeness_grid, ~species * minor, type = 'response')
  write_csv(tidy(minor_response, conf.int = TRUE), '/scratch/j.selwyn/species_minor.csv.gz')
  
  path_response <- emmeans(completeness_grid, ~species * major_minor_path, type = 'response')
  write_csv(tidy(minor_response, conf.int = TRUE), '/scratch/j.selwyn/species_path.csv.gz')
}


library(tidytext)
major_response %>%
  mutate(species = reorder_within(species, prob, major)) %>%
  ggplot(aes(y = species, x = prob, xmin = conf.low, xmax = conf.high)) +
  geom_segment(aes(yend = species), xend = 0, linetype = 'dotted') +
  geom_pointrange() +
  facet_wrap(~major, scales = 'free_y') +
  scale_y_reordered() +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1))


major_response %>%
  mutate(species = fct_reorder(species, prob)) %>%
  ggplot(aes(y = species, x = prob, xmin = conf.low, xmax = conf.high,
             colour = major)) +
  geom_segment(aes(yend = species), xend = 0, linetype = 'dotted',
               position = position_dodge(0.75)) +
  geom_pointrange(position = position_dodge(0.75)) +
  scale_x_continuous(labels = scales::percent_format(), 
                     limits = c(0, 1)) + 
  scale_y_discrete(labels = ~str_remove(., '_.*'))


minor_response %>%
  mutate(species = reorder_within(species, prob, major)) %>%
  ggplot(aes(y = species, x = prob, xmin = conf.low, xmax = conf.high,
             colour = minor)) +
  geom_segment(aes(yend = species), xend = 0, linetype = 'dotted',
               position = position_dodge(0.75)) +
  geom_pointrange(position = position_dodge(0.75)) +
  scale_x_continuous(labels = scales::percent_format(), 
                     limits = c(0, 1)) + 
  scale_y_reordered() +
  facet_wrap(~major, scales = 'free_y')

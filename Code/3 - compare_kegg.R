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

#### Composition of Major KEGG categories ####
kegg_gene_composition <- species_kegg %>%
  select(species, gene_id, all_of(colnames(kegg_pathways)[-1])) %>%
  distinct %>%
  pivot_longer(cols = -c(gene_id, species), 
               values_drop_na = TRUE) %>%
  separate(name, into = c('major', 'minor', 'extra'), sep = ' - ', fill = 'right') %>%
  group_by(species, major) %>%
  summarise(value = n(),
            .groups = 'drop') %>%
  
  pivot_wider(names_from = 'major',
              values_from = 'value')

major_chisq <- kegg_gene_composition %>%
  column_to_rownames('species') %>%
  chisq.test()

major_chisq

as_tibble(major_chisq$observed - major_chisq$expected, rownames = 'species') %>%
  pivot_longer(cols = -species) %>%
  ggplot(aes(x = name, y = species, fill = value)) +
  geom_tile()

species_kegg %>%
  select(species, gene_id, all_of(colnames(kegg_pathways)[-1])) %>%
  distinct %>%
  pivot_longer(cols = -c(gene_id, species), 
               values_drop_na = TRUE) %>%
  separate(name, into = c('major', 'minor', 'extra'), sep = ' - ', fill = 'right') %>%
  mutate(major = str_replace(major, ' ', '\n')) %>%
  
  ggplot() +
  geom_mosaic(aes(x = product(major), fill = major, conds = product(species))) +
  guides(fill = 'none') +
  scale_x_productlist(labels = ~str_remove(., '_(.*)'), position = 'bottom') +
  scale_y_productlist(labels = ~str_replace(., ' ', '\n'), position = 'right') +
  coord_flip() +
  theme_classic() +
  theme(axis.title = element_blank())

#### Just look at acroporids major pathway breakdown ####
major_acro_chisq <- kegg_gene_composition %>%
  filter(str_detect(species, '^a'), species != 'aten') %>%
  column_to_rownames('species') %>%
  chisq.test()

major_acro_chisq

kegg_gene_composition %>%
  filter(str_detect(species, '^a'), species != 'aten', species != 'amil_fuller') %>%
  column_to_rownames('species') %>%
  chisq.test()

kegg_gene_composition %>%
  select(-`Environmental Information Processing`, -`Organismal Systems`) %>%
  filter(str_detect(species, '^a'), species != 'aten') %>%
  column_to_rownames('species') %>%
  chisq.test()

as_tibble(major_acro_chisq$observed - major_acro_chisq$expected, rownames = 'species') %>%
  pivot_longer(cols = -species) %>%
  ggplot(aes(x = name, y = species, fill = value)) +
  geom_tile()

species_kegg %>%
  filter(str_detect(species, '^a'), species != 'aten', species != 'amil_fuller') %>%
  select(species, gene_id, all_of(colnames(kegg_pathways)[-1])) %>%
  distinct %>%
  pivot_longer(cols = -c(gene_id, species), 
               values_drop_na = TRUE) %>%
  separate(name, into = c('major', 'minor', 'extra'), sep = ' - ', fill = 'right') %>%
  
  ggplot() +
  geom_mosaic(aes(x = product(major), fill = major, conds = product(species))) +
  guides(fill = 'none') +
  scale_x_productlist(labels = ~str_remove(., '_(.*)'), position = 'bottom') +
  scale_y_productlist(labels = ~str_replace(., ' ', '\n'), position = 'right') +
  coord_flip() +
  theme_classic() +
  theme(axis.title = element_blank())


kegg_gene_composition %>%
  mutate(Y = cbind(!!!syms(colnames(kegg_gene_composition)[-1]))) %>%
  select(species, Y) %>%
  mutate(total = rowSums(Y))

#### Minor Categories withing Major ####
minor_kegg_cats <- species_kegg %>%
  select(species, gene_id, all_of(colnames(kegg_pathways)[-1])) %>%
  # filter(str_detect(species, '^a'), species != 'aten') %>%
  distinct %>%
  pivot_longer(cols = -c(gene_id, species), 
               values_drop_na = TRUE) %>%
  separate(name, into = c('major', 'minor', 'extra'), sep = ' - ', fill = 'right') %>%
  nest(minor_breakdown = -c(major)) %>%
  rowwise 

minor_kegg_breakdown <- minor_kegg_cats %>%
  mutate(species_breakdown = list(group_by(minor_breakdown, species, minor) %>% 
                                    summarise(n_gene = n_distinct(gene_id),
                                              .groups = 'drop') %>%
                                    pivot_wider(names_from = minor, 
                                                values_from = n_gene)),
         
         column_to_rownames(species_breakdown, 'species') %>%
           chisq.test() %>%
           tidy,
         
         minor_category_tests = list(pivot_longer(species_breakdown, cols = -species,
                        names_to = 'minor_category') %>%
           nest_by(minor_category) %>%
           summarise(column_to_rownames(data, 'species') %>%
                       chisq.test() %>%
                       tidy, 
                     .groups = 'drop')))
  
minor_kegg_breakdown %>%
  ungroup %>%
  slice(4) %>%
  select(major, minor_category_tests) %>%
  unnest(minor_category_tests) %>%
  select(minor_category, parameter, statistic, p.value)

minor_kegg_breakdown %>%
  select(major, minor_category_tests) %>%
  unnest(minor_category_tests) %>%
  select(minor_category, parameter, statistic, p.value) %>%
  mutate(fdr = p.adjust(p.value, 'fdr')) %>% 
  filter(fdr < 0.05)


minor_kegg_breakdown_acropora <- minor_kegg_cats %>%
  mutate(species_breakdown = list(group_by(minor_breakdown, species, minor) %>% 
                                    summarise(n_gene = n_distinct(gene_id),
                                              .groups = 'drop') %>%
                                    filter(str_detect(species, '^a'), species != 'aten', species != 'amil_fuller') %>%
                                    pivot_wider(names_from = minor, 
                                                values_from = n_gene)),
         
         column_to_rownames(species_breakdown, 'species') %>%
           chisq.test() %>%
           tidy,
         
         minor_category_tests = list(pivot_longer(species_breakdown, cols = -species,
                                                  names_to = 'minor_category') %>%
                                       nest_by(minor_category) %>%
                                       summarise(column_to_rownames(data, 'species') %>%
                                                   chisq.test() %>%
                                                   tidy, 
                                                 .groups = 'drop')))

minor_kegg_breakdown_acropora %>%
  ungroup %>%
  slice(5) %>%
  select(major, minor_category_tests) %>%
  unnest(minor_category_tests) %>%
  select(minor_category, parameter, statistic, p.value)

minor_kegg_breakdown_acropora %>%
  select(major, minor_category_tests) %>%
  unnest(minor_category_tests) %>%
  select(minor_category, parameter, statistic, p.value) %>%
  mutate(fdr = p.adjust(p.value, 'fdr')) %>% 
  filter(fdr < 0.05)

minor_kegg_breakdown_acropora %>%
  mutate(mosaic_plot = list(ggplot(data = mutate(minor_breakdown, species = str_remove(species, '_.*$'))) +
                              geom_mosaic(aes(x = product(minor, species), fill = minor)) +
                              theme(legend.position = 'none',
                                    axis.title = element_blank())),
         
         mosaic_plot = list(mosaic_plot + 
                              labs(title = str_replace(major, ' ', '\n')))) %>%
  pull(mosaic_plot) %>%
  wrap_plots()
  

minor_kegg_breakdown_acropora %>%
  filter(major == 'Organismal Systems') %>%
  select(minor_breakdown) %>%
  unnest(minor_breakdown) %>%
  # filter(minor == 'Immune system') %>%
  group_by(species, minor) %>%
  summarise(n = n_distinct(gene_id),
            .groups = 'drop') %>%
  mutate(minor = if_else(minor == 'Immune system', 'Immune', 'Other')) %>%
  group_by(species, minor) %>%
  summarise(n = sum(n),
            .groups = 'drop') %>%
  pivot_wider(names_from = 'minor',
              values_from = n) %>%
  mutate(prop_immune = Immune / (Immune + Other)) %>%
  mutate(species = fct_reorder(species, prop_immune)) %>%
  ggplot(aes(y = species, x = prop_immune)) +
  geom_point()

#### Individual Pathway ORA ####
species_paths <- species_kegg %>%
  pivot_longer(cols = all_of(colnames(kegg_pathways)[-1]),
               names_to = 'major_minor',
               values_to = 'pathway',
               values_drop_na = TRUE) %>%
  # select(gene_id, pathway) %>%
  group_by(across(Orthogroup:major_minor)) %>%
  reframe(pathway = str_split(pathway, ';') %>% unlist %>% str_trim) %>%
  select(species, Orthogroup, gene_id, kegg_gene, kegg_orthology, major_minor, pathway)

path_matrix <- species_paths %>%
  group_by(species, major_minor, pathway) %>%
  summarise(n_gene = n_distinct(gene_id),
            .groups = 'drop') %>%
  pivot_wider(names_from = species,
              values_from = n_gene,
              values_fill = 0) %>%
  separate(major_minor, 
           sep = ' - ', 
           into = c('major', 'minor'), 
           extra = 'merge') %>%
  filter(if_all(where(is.integer), ~. > 5))

path_matrix %>%
  select(where(is.integer)) %>%
  chisq.test()

path_matrix %>%
  nest_by(major, minor, pathway) %>%
  summarise(select(data, starts_with('a'), -aten) %>%
              chisq.test() %>%
              tidy,
            .groups = 'drop') %>%
  mutate(fdr = p.adjust(p.value, 'fdr')) %>%
  filter(fdr < 0.05) %>% View


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
  group_by(major_minor, pathway) %>%
  summarise(total_species = n_distinct(species),
            total_cnidarian = n_distinct(kegg_gene),
            .groups = 'drop')

complete_pathways %>%
  arrange(total_kegg) %>%
  filter(total_species != 24)

species_pathway_completeness <- species_paths %>%
  group_by(species, major_minor, pathway) %>%
  summarise(n_kegg = n_distinct(kegg_gene),
            n_ortho = n_distinct(Orthogroup),
            n_gene = n_distinct(gene_id),
            .groups = 'drop') %>%
  left_join(complete_cnidarian_pathways,
            by = c('major_minor', 'pathway')) %>%
  separate(major_minor, 
           sep = ' - ', 
           into = c('major', 'minor'), 
           extra = 'merge') %>%
  mutate(path_id = str_extract(pathway, 'map[0-9]+')) %>%
  full_join(expand_grid(complete_pathways, species = unique(species_paths$species)),
            by = c('path_id' = 'pathway', 'species')) %>%
  select(-path_id) %>%
  
  
  mutate(percent_species = n_kegg / total_kegg,
         percent_cnidarian = total_cnidarian / total_kegg) %>%
  filter(total_species == 24 & total_kegg > 5) %>%
  filter(!minor %in% c('Cellular community - prokaryotes',
                       'Information processing in viruses'))


## For any group of pathways is there 
species_pathway_completeness %>%
  mutate(diff = percent_species - percent_cnidarian) %>%
  select(species:pathway, diff) %>%
  pivot_wider(names_from = species,
              values_from = diff)


species_pathway_completeness %>%
  select(-species, -n_kegg:-n_gene, -percent_complete) %>%
  distinct %>%
  mutate(percent_cnidarian = total_cnidarian / total_kegg) %>%
  group_by(major, minor) %>%
  summarise(mean_complete = mean(percent_cnidarian),
            se_complete = sd(percent_cnidarian) / sqrt(n()),
            .groups = 'drop') %>%
  
  ggplot(aes(x = major, y = mean_complete, fill = major)) +
  geom_boxplot()
  
  


species_pathway_completeness %>%
  group_by(species, major, minor) %>%
  summarise(mean_complete = mean(percent_complete),
            se_complete = sd(percent_complete) / sqrt(n()),
            .groups = 'drop') %>%
  filter(str_detect(species, 'acer|adig|spis|epal|nvec')) %>% 
  
  ggplot(aes(x = species, y = mean_complete, fill = species)) +
  geom_boxplot() +
  facet_wrap(~major)

library(tidytext)
species_pathway_completeness %>%
  group_by(species, major, minor) %>%
  summarise(mean_complete = mean(percent_complete),
            se_complete = sd(percent_complete) / sqrt(n()),
            .groups = 'drop') %>%
  mutate(species = reorder_within(species, mean_complete, major)) %>%
  
  ggplot(aes(y = species, x = mean_complete)) +
  geom_boxplot() +
  scale_y_reordered(labels = ~str_remove(., '_(.*)$')) +
  facet_wrap(~major, scales = 'free_y')


species_pathway_completeness %>%
  # filter(minor == 'Amino acid metabolism') %>%
  select(-n_kegg:-total_kegg) %>%
  pivot_wider(names_from = 'species',
              values_from = 'percent_complete') %>%
  
  filter(if_any(where(is.numeric), ~is.na(.))) %>% 
  pull(pathway)


tst_model <- lm(percent_complete ~ species * minor,
   data = species_pathway_completeness)

summary(tst_model)
anova(tst_model)

species_pathway_completeness %>%
  # filter(minor == 'Amino acid metabolism') %>%
  select(-n_kegg:-total_kegg) %>%
  pivot_wider(names_from = 'species',
              values_from = 'percent_complete')

tst_emmeans <- emmeans(tst_model, ~species * major / minor)

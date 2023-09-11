#### Libraries ####
library(multcomp)
library(tidyverse)
library(broom)
library(emmeans)
library(patchwork)
library(ggmosaic)
library(chisq.posthoc.test)

#### Functions ####
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
species_kegg <- read_csv('../../Bioinformatics/Phylogenomics/Annotations/species_kegg_orthogroup.csv.gz', 
                         show_col_types = FALSE, guess_max = 1e5)

functional_annotations <- read_rds('../../Bioinformatics/genome_annotation/k2_functionalAnnotations.rds')

#### Count Genes & percent with a kegg ####
species_kegg %>%
  group_by(species) %>%
  summarise(n_gene = n_distinct(gene_id),
            n_kegg = n_distinct(gene_id[!is.na(kegg_gene)]),
            
            n_ortho = n_distinct(Orthogroup),
            n_kegg_ortho = n_distinct(Orthogroup[!is.na(kegg_gene)])) %>%
  mutate(pct_gene = n_kegg / n_gene,
         pct_ortho = n_kegg_ortho / n_ortho) %>%
  select(species, pct_ortho)

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

species_kegg %>%
  filter(kegg_orthology %in% c('ko:K03654', 'ko:K07497', 'ko:K08378')) %>%
  count(Orthogroup, species) %>%
  pivot_wider(names_from = species, 
              values_from = n) %>% View

species_kegg %>%
  group_by(Orthogroup) %>%
  filter(n_distinct(species) == 1) %>%
  filter(species == 'acer') %>%
  ungroup 

unique_acer <- species_kegg %>%
  group_by(Orthogroup) %>%
  filter(n_distinct(species) == 1) %>%
  filter(species == 'acer') %>%
  ungroup %>%
  select(Orthogroup, gene_id) %>%
  distinct %>%
  mutate(gene_id = str_remove(gene_id, '-RA$')) %>%
  left_join(functional_annotations,
            by = c('gene_id' = 'qseqid')) %>% 
  filter(!is.na(has_swissprot))

swissprot_unique <- unique_acer %>%
  filter(!is.na(swissprot_name)) %>%
  select(Orthogroup, swissprot_id, swissprot_name) %>%
  distinct %>%
  mutate(swissprot = str_c(swissprot_name, ' (', swissprot_id, ')'),
         .keep = 'unused') %>%
  group_by(Orthogroup) %>%
  mutate(id_number = str_c('swissprot_', row_number())) %>%
  pivot_wider(names_from = id_number, values_from = swissprot) %>%
  ungroup

eggnog_unique <- unique_acer %>%
  select(Orthogroup, gene_id, entap_results) %>%
  unnest(entap_results) %>% 
  select(Orthogroup, gene_id, `EggNOG Predicted Gene`, `EggNOG Description`) %>%
  filter(!is.na(`EggNOG Predicted Gene`)) %>%
  select(-gene_id) %>%
  mutate(eggnog_pred = str_c(`EggNOG Description`, ' (', `EggNOG Predicted Gene`, ')'),
         .keep = 'unused') %>%
  distinct %>%
  group_by(Orthogroup) %>%
  mutate(id_number = str_c('eggnog_', row_number())) %>%
  pivot_wider(names_from = id_number, values_from = eggnog_pred) %>%
  ungroup

pfam_unique <- unique_acer %>%
  select(Orthogroup, gene_id, interproscan_results) %>%
  unnest(interproscan_results) %>%
  filter(analysis == 'Pfam') %>%
  select(Orthogroup, signature_accession, signature_description) %>%
  distinct() %>%
  filter(!signature_accession %in% c('PF18701', 'PF20209')) %>%
  mutate(pfam = str_c(signature_description, ' (', signature_accession, ')'),
         .keep = 'unused') %>%
  distinct %>%
  group_by(Orthogroup) %>%
  mutate(id_number = str_c('pfam_', row_number())) %>%
  pivot_wider(names_from = id_number, values_from = pfam) %>%
  ungroup

kegg_unique <- species_kegg %>%
  group_by(Orthogroup) %>%
  filter(n_distinct(species) == 1) %>%
  filter(species == 'acer') %>%
  ungroup %>%
  filter(!is.na(kegg_gene)) %>%
  select(Orthogroup, kegg_gene, kegg_orthology) %>%
  distinct %>%
  mutate(kegg = str_c(kegg_gene, ' (', kegg_orthology, ')'),
         .keep = 'unused') %>%
  group_by(Orthogroup) %>%
  mutate(id_number = str_c('kegg_', row_number())) %>%
  pivot_wider(names_from = id_number, values_from = kegg) %>%
  ungroup

full_join(swissprot_unique,
          eggnog_unique,
          by = join_by(Orthogroup)) %>%
  full_join(pfam_unique,
            by = join_by(Orthogroup)) %>%
  full_join(kegg_unique,
            by = join_by(Orthogroup)) %>%
  write_csv('../Results/unique_kegg_acer.csv')

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
              values_fill = FALSE) %>% View


#### Check KEGG Assignments ####
summary_stats_kegg <- species_kegg %>%
  group_by(species) %>%
  summarise(n_gene = n_distinct(gene_id),
            n_gene_w_kegg = n_distinct(gene_id[!is.na(kegg_orthology)]),
            prop_gene = n_gene_w_kegg / n_gene, 
            
            unique_kegg = n_distinct(kegg_orthology),
            
            n_ortho = n_distinct(Orthogroup),
            n_ortho_w_kegg = n_distinct(Orthogroup[!is.na(kegg_orthology)]),
            prop_ortho = n_ortho_w_kegg / n_ortho) %>%
  mutate(is_acroporid = (str_detect(species, '^a') & species != 'aten'))

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

acer_vs_acropora <- c(-1/16, -1/16, 1, -1/16, -1/16, -1/16, -1/16, -1/16, 
                      -1/16, -1/16, -1/16, -1/16, -1/16, -1/16, -1/16, 
                      0, -1/16, -1/16, 0, 0, 0, 0, 0, 0)
emmeans(gene_prop_kegg, ~species, type = 'link') %>%
  contrast(method = list(acer_vs_acropora = acer_vs_acropora))

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

acer_vs_acropora <- c(-1/16, -1/16, 1, -1/16, -1/16, -1/16, -1/16, -1/16, 
                      -1/16, -1/16, -1/16, -1/16, -1/16, -1/16, -1/16, 
                      0, -1/16, -1/16, 0, 0, 0, 0, 0, 0)

acropora_vs_out <- c(1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 
                     1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 
                     -1/5, 1/17, 1/17, -1/5, 0, 0, -1/5, -1/5, -1/5)

acropora_vs_montipora <- c(1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 
                           1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 1/17, 
                           0, 1/17, 1/17, 0, -1/2, -1/2, 0, 0, 0)

emmeans(ortho_prop_kegg, ~species, type = 'link') %>%
  contrast(method = list(acer_vs_acropora = acer_vs_acropora,
                         acropora_vs_out = acropora_vs_out,
                         acropora_vs_montipora = acropora_vs_montipora))


emmeans(ortho_prop_kegg, ~species, type = 'response') %>%
  tidy(conf.int = TRUE) %>%
  filter(species == 'acer')

summary_stats_kegg %>%
  mutate(species = fct_reorder(species, unique_kegg)) %>%
  
  ggplot(aes(y = species, x = unique_kegg)) +
  geom_point()

#### Composition of KEGG pathways ####
species_paths <- species_kegg %>%
  pivot_longer(cols = all_of(contains(' - ')),
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

#### Distribution of Acer paths ####
acer_pathways <- species_paths %>%
  filter(species == 'acer') %>%
  group_by(major, minor, pathway) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  filter(n_ortho > 0) %>%
  mutate(pathway = str_c(major, minor, pathway, sep = ';;'), .keep = 'unused', .before = everything()) %>%
  column_to_rownames('pathway') 

chisq.test(acer_pathways)

species_paths %>%
  filter(species == 'acer') %>%
  group_by(major) %>%
  summarise(n_ortho = n_distinct(Orthogroup),
            .groups = 'drop') %>%
  mutate(prop = n_ortho / sum(n_ortho))

#### Across all species ####
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
fisher.test(path_matrix, simulate.p.value = TRUE)

path_posthoc <- chisq.posthoc.test(path_matrix, method = 'fdr') %>%
  process_posthoc %>%
  filter(fdr < 0.05) %>%
  separate(Dimension, into = c('major', 'minor', 'pathway'), sep = ';;') 
# count(major, minor, pathway) %>%
# arrange(-n)

path_posthoc %>%
  filter(species == 'acer')

path_posthoc %>%
  count(species) %>%
  arrange(-n)

path_posthoc %>%
  count(major, minor, pathway) %>%
  arrange(-n)

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

all_kegg_comps <- expand_grid(counting = c('gene_id', 'Orthogroup', 'kegg_orthology'),
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

all_kegg_comps %>%
  ungroup %>%
  filter(level == 'pathway',
         counting != 'kegg_orthology') %>%
  select(counting, species_inclusion,
         parameter, statistic, p.value)


all_kegg_comps %>%
  ungroup %>%
  filter(counting == 'gene_id',
         species_inclusion == 'all',
         level == 'pathway') %>%
  select(posthoc) %>%
  unnest(posthoc) %>%
  # filter(species == 'acer') %>% pull(pathway)
  # count(species) %>% arrange(-n)
  # filter(species == 'amil_fuller')
  
  # filter(str_detect(species, '^a'),
  #        species != 'aten', species != 'amil_fuller') %>%
  # # count(species) %>% arrange(-n)
  # count(major, minor, pathway) %>%
  # filter(n == 15) %>% pull(pathway)
  
  count(major, minor, pathway) %>%
  arrange(-n)

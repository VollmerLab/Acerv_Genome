if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  blast_kegg_annotations <- args[1]
  othrogroup_file <- args[2]
  outdir <- args[3]
} else {
  blast_kegg_annotations <- '/scratch/j.selwyn/phylogenetics/annotations/kegg_annotations.csv.gz'
  othrogroup_file <- '/scratch/j.selwyn/phylogenetics/orthofinder/Results_Aug18/Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
  outdir <- '/scratch/j.selwyn/phylogenetics'
}

#### Libraries ####
suppressMessages(suppressWarnings(library(KEGGREST)))
suppressMessages(suppressWarnings(library(tidyverse)))

#### Data ####
kegg_annotations <- read_csv(blast_kegg_annotations, show_col_types = FALSE)

all_orthogroups <- read_delim(othrogroup_file, show_col_types = FALSE) %>%
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

#### Identify Orthogroups in each KEGG ####
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
  left_join(select(kegg_annotations, -species),
            by = c('gene_id' = 'qseqid'),
            multiple = "all") %>%
  filter(!is.na(kegg_orthology)) %>%
  gene_kegg_vote('Orthogroup') %>%
  select(-species:-uniprot_id)

kegg_vote_agreement <- ortho_kegg_vote %>%
  ggplot(aes(x = prop_top_hits)) +
  geom_histogram(bins = 30) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = 'Blast hit KO Agreement - must have a KO to include (%)') +
  theme_classic()
ggsave(str_c(outdir, 'kegg_vote_agreement.png', sep = '/'), height = 7, width = 7)

#### Get Kegg Pathway data ####
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

kegg_pathways <- ortho_kegg_vote %>%
  select(kegg_orthology, kegg_pathway) %>%
  filter(!is.na(kegg_pathway)) %>%
  distinct %>%
  rowwise(kegg_orthology) %>%
  summarise(kegg_path = str_split(kegg_pathway, ';;') %>% unlist,
          .groups = 'drop') %>%
  distinct %>%
  nest(data = c(kegg_orthology)) %>%
  mutate(grouping = row_number() %/% 10) %>%
  group_by(grouping) %>%
  summarise(kegg_paths = str_c(kegg_path, collapse = ';;'),
          data = list(data),
          .groups = 'drop') %>%
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

#### Join Kegg and Orthogroups to species ####
species_kegg <- all_orthogroups %>%
  left_join(ortho_kegg_vote,
            by = 'Orthogroup') %>%
  left_join(kegg_pathways, 
            by = 'kegg_orthology')

write_csv(species_kegg, str_c(outdir, 'species_kegg_orthogroup.csv.gz', sep = '/'))

#### Convert to Orthofinder like output for easy CAFE5 conversion ####
orthofinder_annotations <- select(species_kegg, Orthogroup, species, kegg_orthology) %>%
  filter(!is.na(kegg_orthology)) %>%
  distinct %>%
  pivot_wider(names_from = species, 
              values_from = 'Orthogroup',
              values_fn = ~str_c(., collapse = ', ')) %>%
  rename(HOG = kegg_orthology) %>%
  mutate(OG = NA_character_,
         "Gene Tree Parent Clade" = NA_character_,
         .after = 'HOG')
write_csv(orthofinder_annotations, str_c(outdir, 'orthofinder_annotations.csv.gz', sep = '/'))

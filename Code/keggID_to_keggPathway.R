library(KEGGREST)
library(tidyverse)

kegg_annotations <- read_csv('C:/Users/jdsel/Documents/Google Drive/Research/Vollmer Lab PostDoc/Gene_Expression/intermediate_files/kegg_annotations.csv.gz', 
                             show_col_types = FALSE)



#### KEGG Pathway Identification ####
if(file.exists('../intermediate_files/kegg_orthogroup_pathways.csv.gz')){
  kegg_paths <- read_csv('../intermediate_files/kegg_orthogroup_pathways.csv.gz', show_col_types = FALSE)
} else {
  kegg_paths <- kegg_annotations %>%
    filter(!is.na(kegg_orthology)) %>%
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

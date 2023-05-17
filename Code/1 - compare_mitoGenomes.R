library(tidyverse)
library(magrittr)
library(rentrez)
library(read.gb)

all_mito_ids <- read_csv('../Data/Acropora mtGenomes.csv')


process_ncbi <- function(ID){
  if(file.exists(ID)){
    raw_out <- read.gb(ID) 
  } else {
    tmp_file <- tempfile()
    write_lines(entrez_fetch(db="nuccore", id = ID, rettype = "gbwithparts"), tmp_file)
    raw_out <- read.gb(tmp_file)
  }
  
  processed_out <- raw_out[[1]]$FEATURES[names(raw_out[[1]]$FEATURES) %in% c('tRNA', 'rRNA', 'gene')] %>%
    map_dfr(~slice(.x, 1:2) %>%
              mutate(Location = c('location', 'product')) %>%
              # rename(type = Location) %>%
              pivot_wider(names_from = Location,
                          values_from = Qualifier), 
            .id = 'type')
  file.remove(tmp_file)
  processed_out
}


all_features <- all_mito_ids %>%
  mutate(Title = case_when(str_detect(Title, 'aff') ~ str_extract(Title, 'aff. [a-z]+'),
                           TRUE ~ str_extract(Title, 'Acropora [a-z]+')),
         Title = str_remove_all(Title, 'aff. |Acropora '),
         Title = if_else(Title == 'tenuis' & EntrezUID == '2102348512', 'tenuis2', Title)) %>%
  select(Title, EntrezUID) %>%
  mutate(EntrezUID = as.character(EntrezUID)) %>%
  add_row(Title = 'Acerv', EntrezUID = '../../Bioinformatics/mitochondrial_genome/k2_mitogenome.gbf', .before = 0) %>%
  rowwise(Title) %>%
  summarise(process_ncbi(EntrezUID), 
            .groups = 'drop')


all_features %>%
  select(-type) %>%
  distinct %>%
  group_by(Title) %>%
  mutate(order = row_number()) %>%
  ungroup %>%
  select(-location) %>%

  # filter(Title == 'tenuis', product == 'COX1')
  
  pivot_wider(names_from = 'Title',
              values_from = 'product') %T>%
  write_csv('../intermediate_files/mitochondrial_order.csv')

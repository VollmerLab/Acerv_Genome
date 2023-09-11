library(Biostrings)
library(tidyverse)

mtgenomes <- readDNAStringSet('../../intermediate_files/acropora_mtgenomes.fasta')
lengths(mtgenomes)

renamed_names <- tibble(name = names(mtgenomes)) %>%
  rowwise %>%
  mutate(new_name = if_else(str_locate(name, '_')[1] == 5, 
                            name, NA_character_)) %>%
  ungroup %>%
  mutate(species_code = str_extract(name, '(Acropora|Montipora) [a-z]+') %>%
           str_to_lower() %>%
           str_remove('(cropora|ontipora) ') %>%
           str_sub(start = 1, end = 4),
         accession = str_extract(name, '^.*\\.') %>% str_remove('\\..*')) %>%
  mutate(update_name = str_c(species_code, accession, sep = '_'),
         .keep = 'unused') %>%
  mutate(new_name = coalesce(new_name, update_name), 
         .keep = 'unused') 
write_csv(renamed_names, '../../intermediate_files/mitochondria_rename.csv')

names(mtgenomes) <- renamed_names$new_name


mtgenomes %>%
  # replaceAmbiguities(new = "N") %>%
  unique() %>%
  magrittr::extract(str_detect(names(.), '^asp', negate = TRUE)) %>%
  magrittr::extract(str_detect(names(.), '^acer_M', negate = TRUE)) %>%
  magrittr::extract(str_detect(names(.), '^apro', negate = TRUE)) %>%
  magrittr::extract(str_detect(names(.), '^apal', negate = TRUE)) %>%
  writeXStringSet(filepath = '../../intermediate_files/reduced_acropora_mtgenomes.fa')




#### Plot ####

library(treedataverse)
library(tidyverse)


rename_sequences <- read_csv('../../intermediate_files/mitochondria_rename.csv')

read.newick('../../intermediate_files/species.treefile') %>% 
  root('meff_NC_040137') %>% 
  # drop.tip('meff_NC_040137') %>%
  # drop.tip('acer_MK574927') %>%
  as.polytomy(feature='node.label', fun=function(x) as.numeric(x) < 75) %>%
  as_tibble %>%
  left_join(rename_sequences, 
            by = c('label' = 'new_name')) %>%
  mutate(label = case_when(str_detect(label, '^[0-9]+$') ~ label,
                           label == '' ~ label,
                           str_detect(label, 'thisStudy') ~ 'Acropora cervicornis (this study)',
                           label == name ~ str_sub(name, 1, 4),
                           TRUE ~ str_extract(name, '(Acropora|Montipora) [a-z]+'))) %>% 
  # filter(is.na(label2)) %>% View
  # mutate(label = str_remove(label, '_.*')) %>%
  as.treedata() %>%
  ggtree(layout = 'rectangular', ladderize = TRUE,
         branch.length = "branch.length") +
  geom_tippoint(aes(colour = label), show.legend = FALSE) +
  geom_tiplab() +
  geom_nodelab()



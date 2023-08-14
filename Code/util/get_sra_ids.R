library(tidyverse)
library(rentrez)
library(xml2)

shinzato_bioproj <- 'PRJDB8519'

#### Functions ####
get_bioproj_sra <- function(proj_id){
  bio_ids <- entrez_search(db = 'bioproject', term = proj_id)
  
  link_ids <- entrez_link(dbfrom = 'bioproject', id = bio_ids$ids, db = 'sra')
  
  entrez_fetch(db = 'sra', id = link_ids$links$bioproject_sra, rettype = 'XML') %>%
    process_sra_out()
}

process_sra_out <- function(sra_out){
 tibble(sequence_id = read_xml(sra_out) %>%
          xml_find_all('EXPERIMENT_PACKAGE') %>%
          xml_find_all('EXPERIMENT') %>%
          xml_find_all('IDENTIFIERS') %>%
          xml_find_all('PRIMARY_ID') %>%
          xml_text(),
        
        species = read_xml(sra_out) %>%
          xml_find_all('EXPERIMENT_PACKAGE') %>%
          xml_find_all('SAMPLE') %>%
          xml_find_all('SAMPLE_NAME') %>%
          xml_find_all('SCIENTIFIC_NAME') %>%
          xml_text(),
        
        taxon_id = read_xml(sra_out) %>%
          xml_find_all('EXPERIMENT_PACKAGE') %>%
          xml_find_all('SAMPLE') %>%
          xml_find_all('SAMPLE_NAME') %>%
          xml_find_all('TAXON_ID') %>%
          xml_text(),
        
        sample_id = read_xml(sra_out) %>%
          xml_find_all('EXPERIMENT_PACKAGE') %>%
          xml_find_all('SAMPLE') %>%
          xml_find_all('IDENTIFIERS') %>%
          xml_find_all('EXTERNAL_ID') %>%
          xml_text(),
        
        sra_id = read_xml(sra_out) %>%
          xml_find_all('EXPERIMENT_PACKAGE') %>%
          xml_find_all('RUN_SET') %>%
          xml_find_all('RUN') %>%
          xml_find_all('IDENTIFIERS') %>%
          xml_find_all('PRIMARY_ID') %>%
          xml_text()) 
}




#### Data ####
bioproj_data <- tribble(
  ~author, ~bioproject,
  'shinzato', 'PRJDB8519'
)


#### Download SRA metadata ####
sequence_ids <- bioproj_data %>%
  rowwise %>%
  reframe(get_bioproj_sra(bioproject)) %>%
  filter(str_detect(species, 'Acropora'))

write_csv(sequence_ids, '../../intermediate_files/sra_ids.csv')
  

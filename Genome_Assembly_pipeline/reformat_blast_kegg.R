#coral OMA2
# GetKEGG.Rmd
# KEGG will lock you out after too long
# uniprot2keggKO.Rdata (has already downloaded maps)
# load blast hits
# get protein ID for all blast hits & cross reference what is already retrevied
# look at all organism specific keggs
# can do top 5 blast hits & then pick best with an ortholog
# 

if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  blast_dir <- args[1]
  protein_dir <- args[2]
} else {
  blast_dir <- '/scratch/j.selwyn/phylogenetics/annotations/blast_array'
  protein_dir <- '/scratch/j.selwyn/phylogenetics/proteins'
  # blast_dir <- '../../Bioinformatics/tmp/small.tsv'
}

#### Libraries ####
library(rlang, lib.loc = '/home/j.selwyn/R/x86_64-pc-linux-gnu-library/4.2') #need this version for multidplyr
suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(KEGGREST)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(multidplyr)))

#### Make Parallel Cluster ####
cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c('KEGGREST', 'dplyr', 'tibble', 'purrr'))

#### Functions ####
# dbentries <- full_kegg$data[[10]]$kegg_orgID
# dbentries[5] <- 'blahblah'
keggGet_jds <- function (dbentries, option = c("aaseq", "ntseq", "mol", "kcf", 
                                               "image", "kgml")) 
{
  if (length(dbentries) > 10) 
    warning(paste("More than 10 inputs supplied, only the first", 
                  "10 results will be returned."))
  db_names <- dbentries
  dbentries <- paste(dbentries, collapse = "+")
  url <- sprintf("%s/get/%s", KEGGREST:::.getRootUrl(), dbentries)
  if (!missing(option)) {
    url <- sprintf("%s/%s", url, option)
    if (option == "image") 
      return(content(GET(url), type = "image/png"))
    if (option %in% c("aaseq", "ntseq")) {
      t <- tempfile()
      cat(KEGGREST:::.getUrl(url, KEGGREST:::.textParser), file = t)
      if (option == "aaseq") 
        return(readAAStringSet(t))
      else if (option == "ntseq") 
        return(readDNAStringSet(t))
    }
    if (option %in% c("mol", "kcf", "kgml")) 
      return(KEGGREST:::.getUrl(url, KEGGREST:::.textParser))
  }
  if (grepl("^br:", dbentries[1])) 
    return(KEGGREST:::.getUrl(url, KEGGREST:::.textParser))
  out <- KEGGREST:::.getUrl(url, KEGGREST:::.flatFileParser)
  names(out) <- purrr::map_chr(out, ~stringr::str_c(names(.x$ORGANISM), .x$ENTRY, sep = ':'))
  out
}

kegg_ortho_get <- function(data){
  keggGet_jds(data$kegg_orgID) %>%
    enframe(name = 'kegg_orgID', value = 'kegg_info') %>%
    left_join(data, ., by = 'kegg_orgID')
}

# data <- kegg_orthologies; uniprot_column <- 'uniprot_id'
get_kegg_ortho <- function(data, uniprot_column){
  # initial_keggs <- pull(data, !!sym(uniprot_column)) %>%
  #   keggConv("genes", ., querySize = 100) %>%
  #   enframe(name = uniprot_column, value = 'kegg_orgID')
  
  cluster_copy(cluster, c('uniprot_column', 'keggGet_jds', 'kegg_ortho_get'))
  
  initial_keggs <- data %>%
    select(!!sym(uniprot_column)) %>%
    mutate(groupings = (row_number() - 1) %/% 100) %>%
    nest_by(groupings) %>%
    mutate(data = list(pull(data, !!sym(uniprot_column)))) %>%
    partition(cluster) %>%
    mutate(kegg_out = list(keggConv("genes", data, querySize = 100) %>%
                             enframe(name = uniprot_column, 
                                     value = 'kegg_orgID'))) %>%
    collect() %>%
    ungroup %>%
    select(kegg_out) %>%
    unnest(kegg_out)
  message('Finished getting KEGG orgIDs: ', Sys.time())
  
  kegg_with_info <- initial_keggs %>%
    filter(!is.na(kegg_orgID)) %>%
    #need to split into groups of 10 to submit to KEGG
    
    mutate(groupings = (row_number() - 1) %/% 10) %>%
    nest_by(groupings) %>% 
    
    #In each of the kegg groups of 10 get the kegg info
    partition(cluster) %>%
    # mutate(data = list(data %>%
    #                      mutate(kegg_info = keggGet(kegg_orgID)))) %>% 
    mutate(data = list(possibly(kegg_ortho_get, otherwise = NULL, 
                                quiet = FALSE)(data))) %>%
    collect() %>%
    ungroup %>%
    unnest(data) 
  write_rds(kegg_with_info, str_c(blast_dir, '/../', 'kegg_annotations.rds'))
  
  full_kegg <- kegg_with_info %>%
    select(-groupings, -kegg_orgID) %>%
    hoist(kegg_info, 'ORTHOLOGY') %>%
    select(-kegg_info) %>%
    unnest(ORTHOLOGY) %>%
    mutate(kegg_orthology = possibly(names, otherwise = NA_character_)(ORTHOLOGY)) %>%
    rename(kegg_gene = ORTHOLOGY) %>%
    distinct %>%
    bind_rows(filter(initial_keggs, is.na(kegg_orgID)) %>% 
                select(-kegg_orgID))
  message('Finished getting KEGG Orthologs: ', Sys.time())
  
  data %>%
    left_join(filter(full_kegg, !is.na(kegg_orthology)),
              by = uniprot_column)
}

#### Set up KEGG Database ####
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
           kegg_class != 'pathway') %>%
  pivot_wider(names_from = kegg_class,
              values_from = path_id,
              names_prefix = 'kegg_',
              values_fn = ~str_c(., collapse = ';;')) 
message('Finished make KEGG DB: ', Sys.time())

#### Read In Annotations ####
swissprot <- list.files(path = blast_dir, pattern = 'blastp$', recursive = TRUE,
           full.names = TRUE) %>%
  tibble(file = .) %>%
  rowwise %>%
  summarise(read_delim(file, delim = '\\t',
                       col_names = c('qseqid', 'swissprot_id', 'swissprot_name',
                                     'pident', 'length', 
                                     'mismatch', 'gapopen', 'qstart', 'qend', 
                                     'sstart', 'send', 'evalue', 'bitscore'),
                       col_types = cols(.default = col_number(),
                                        qseqid = col_character(),
                                        swissprot_id = col_character(),
                                        swissprot_name = col_character())),
            .groups = 'drop') %>%
  mutate(swissprot_name = str_remove_all(swissprot_name, swissprot_id) %>% str_remove('^\\|\\|') %>% str_trim)
message('Finished reading in Swissprot data: ', Sys.time())

#### Get KEGG Orthologies for all blast hits ####
kegg_orthologies <- swissprot %>%
  nest(coral_data = -c(swissprot_id, swissprot_name)) %>%
  mutate(uniprot_id = str_extract(swissprot_id, '\\|.*\\|') %>% str_remove_all('\\|') %>% str_c('up:', .)) %>%
  # sample_n(100) #%>% #for testing
  
  get_kegg_ortho('uniprot_id') %>%
  unnest(coral_data)
message('Finished getting KEGG orthology: ', Sys.time())

#### Get KEGG pathways associated with orthologies ####
kegg_paths <- kegg_orthologies %>%
  mutate(kegg_orthology = str_c('ko:', kegg_orthology)) %>%
  left_join(kegg_db,
            by = c('kegg_orthology' = 'ko')) 
message('Finished merging pathways: ', Sys.time())

#### Output File ####
protein_names <- list.files(protein_dir, pattern = 'fasta$', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(species = str_extract(file, '[a-z_]+.fasta') %>% str_remove('.fasta')) %>% 
  rowwise(species) %>%
  summarise(gene_id = readAAStringSet(file) %>%
              names(),
            .groups = 'drop') %>%
  mutate(gene_id = str_extract(gene_id, '[a-zA-Z0-9_\\.\\-]+'))

kegg_paths %>%
  left_join(protein_names,
            by = c('qseqid' = 'gene_id')) %>%
  select(qseqid, species, starts_with('swissprot'), everything()) %>% 
  write_csv(str_c(blast_dir, '/../', 'kegg_annotations.csv.gz'))

if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  save_dir <- args[1]
} else {
  save_dir <- 'C:/Users/jdsel/Documents/Google Drive/Research/Vollmer Lab PostDoc/Bioinformatics/Phylogenomics/acropora_proteomes'
}

library(tidyverse)
library(rvest)

base_website <- 'https://marinegenomics.oist.jp/gallery/gallery/index'


#### Functions ####
download_species <- function(session_id, gff_or_fasta, species, outdir){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  
  if(str_detect(species, 'Montipora')){
    gff_or_fasta <- if_else(gff_or_fasta == 'gff', 'gtf', 'scaffold.fa')
  }
  
  website_location <- session_follow_link(session_id, species) %>%
    session_follow_link('Downloads') %>%
    session_follow_link(str_c(gff_or_fasta, '.gz'))
  
  # base_name <- str_extract(website_location$url, '[a-z]+.(g[ft]f|fasta).gz')
  base_name <- str_extract(website_location$url, '[a-z]+[\\.a-z]+.gz') %>%
    str_replace_all(c('scaffold.fa' = 'fasta', 'gtf' = 'gff')) %>%
    str_replace_all(c('.fasta.gz' = '_shinzato.fasta.gz', '.gff.gz' = '_shinzato.gff.gz'))
  outfile <- str_c(outdir, base_name, sep = '/')
  
  download.file(website_location$url, outfile)
  outfile
}


#find all "Acropora" and then navigate to downloads and download the predicted proteins
available_acroporids <- read_html(base_website) %>%
  html_nodes('a') %>% 
  html_text2() %>%
  str_subset('Acropora|Montipora') %>%
  unique %>%
  tibble(species_names = .)

website_session <- session(base_website)

downloaded_acroporids <- available_acroporids %>%
  # filter(str_detect(species_names, 'Montipora')) %>%
  expand_grid(file_choice = c('gff', 'fasta')) %>%
  rowwise %>%
  mutate(downloaded_file = possibly(download_species, otherwise = 'not available', quiet = FALSE)(website_session, file_choice, species_names, save_dir)) %>%
  pivot_wider(names_from = file_choice, values_from = downloaded_file)




##TODO - if caribbean is closest then confirm if it is east of the panama line. then keep caribbean otherwise flip to pacific
##TODO - genomes not easily found on ncbi - e.g. Acropora loripes?
##TODO - too many acropora genomes not enough other genera/class genomes?

library(tidyverse)
library(magrittr)
library(rentrez)
library(taxize)
library(rgbif)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

world <- ne_countries(scale = "small", returnclass = "sf")
caribbean_coord <- st_sfc(st_point(c(-78.6569, 21.4691)), crs = 'WGS84')
indopacific_coord <- st_sfc(st_point(c(122.0, 2.8)), crs = 'WGS84')
panama_line <- st_sfc(st_linestring(matrix(c(c(-76.5, 6.3), c(-105.3, 20.2)), nrow = 2, byrow = TRUE)), crs = 'WGS84')

genomes_used <- read_csv('../../Bioinformatics/Phylogenomics/tree_data_sources.csv', show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  filter(used_orthofinder_final == 'y') %>%
  select(species, orthofinder_code)

# Coral Species Distributions 
#http://www.coralsoftheworld.org
# https://www.frontiersin.org/articles/10.3389/fmars.2014.00081/full - paper
caribbean_total <- 73
indopacific_total <- 758

# Coral Genomes from NCBI 
coral_genomes <- entrez_search(db = "genome", term = "scleractinia", retmax = 50) %$%
  entrez_link(dbfrom = 'genome', id = ids, db = 'taxonomy') %$%
  
  classification(links$genome_taxonomy, db = 'ncbi') %>%
  map(as_tibble) %>%
  bind_rows(.id = 'genome_id') %>%
  filter(rank == 'species') %>%
  select(-rank, -id) %>%
  rename(species = name) %>%
  filter(!str_detect(species, 'sp\\.')) %>%
  filter(species != 'Desmophyllum pertusum') %>%
  
  full_join(genomes_used,
          by = 'species') %>% 
  filter(!orthofinder_code %in% c('acer', 'nvec', 'aten', 'edia'))

coral_genomes %>%
  arrange(species) %>% View

# Get distribution data for each genome species
coral_locations <- coral_genomes %>%
  rowwise(species) %>%
  mutate(sampling_data = list(occ_data(scientificName = species))) %>%
  ungroup %>%
  rowwise(species) %>%
  summarise(select(sampling_data$data,starts_with('decimal'), continent, 
                   stateProvince, countryCode, country, waterBody) %>%
              filter_at(vars(everything()), any_vars(!is.na(.))),
            .groups = 'keep')

coral_locations %>%
  ungroup %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>% 
  group_by(species) %>%
  sample_frac(0.1) %>%
  ungroup %>%
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'),
           crs = 'WGS84') %>%
  ggplot(data = .) +
  geom_sf(data = world, fill= "gray", colour = NA) +
  geom_sf(aes(colour = species), alpha = 1, show.legend = FALSE) +
  facet_wrap(~species)


coral_classification <- coral_locations %>%
  summarise(across(c(countryCode, country), ~list(unique(.))),
            mid_lat = median(decimalLatitude, na.rm = TRUE),
            mid_lon = median(decimalLongitude, na.rm = TRUE)) %>%
  st_as_sf(coords = c('mid_lon', 'mid_lat')) %>%
  st_set_crs('WGS84') %>%
  mutate(dist_carib = st_distance(geometry, caribbean_coord)[,1],
         dist_indopacific = st_distance(geometry, indopacific_coord)[,1]) %>%
  as_tibble() %>%
  # sample_frac(0.1) %>%
  rowwise() %>%
  mutate(closer = c('carib', 'indo')[which.min(c_across(c(dist_carib, dist_indopacific)))]) %>%
  ungroup 

#Fix manually a couple errors
coral_classification <- coral_classification %>%
  mutate(closer = if_else(species %in% c('Montipora capitata', 'Pocillopora damicornis'),
                          'indo', closer))

coral_classification %>%
  st_as_sf %>%
  ggplot(data = .) +
  geom_sf(data = world, fill= "gray", colour = NA) +
  geom_sf(aes(colour = species), alpha = 1, show.legend = TRUE) +
  facet_wrap(~closer)

#### Are Caribbean Corals underrepresented in genomes ####
coral_classification %>%
  # filter(species != 'Acropora palmata') %>%
  group_by(closer) %>%
  summarise(n = n_distinct(species)) %>%
  pivot_wider(names_from = closer,
              values_from = n) %>%
  mutate(total = carib + indo) %$%
  prop.test(x = carib, n = total, p = caribbean_total / (caribbean_total + indopacific_total))

prop.test(x = c(2, 73), n = c(34, 758))

tibble(category = c('Genome', 'Species'),
       caribbean = c(2, 73),
       total = c(34, 758)) %>%
  rowwise(category) %>%
  summarise(broom::tidy(prop.test(caribbean, total)),
            .groups = 'drop') %>%
  ggplot(aes(x = category, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange() +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = NULL, 
       y = 'Caribbean (%)') +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 12))


tibble(category = c('Genome', 'Species'),
       caribbean = c(2, 73),
       total = c(34, 758)) %>%
  mutate(prop = caribbean / total) %>%
  glm(prop ~ category, weights = total, family = binomial(), data = .) %>%
  car::Anova()

binom.test(x = 2, n = 34, p = 73 / 758)

coral_classification %>%
  group_by(closer) %>%
  summarise(n = n_distinct(species)) %>%
  mutate(no_genome = c(caribbean_total, indopacific_total),
         no_genome = no_genome - n) %>%
  column_to_rownames('closer') %>%
  fisher.test()
  
#### Get Genome Papers ####
genome_papers <- entrez_search(db = "genome", term = "scleractinia", retmax = 50) %$%
  entrez_link(dbfrom = 'genome', id = ids, db = 'pubmed') %$%
  str_c('https://pubmed.ncbi.nlm.nih.gov/', links$genome_pubmed)
write_lines(genome_papers, '../intermediate_files/genome_paper_urls.txt')
#Select all urls and then drag them onto firefox to open: https://www.reddit.com/r/firefox/comments/pa10ga/open_a_list_of_urls_with_firefox/

#### Bayes for fun ####
library(cmdstanr)
model_program <- cmdstan_model('fisher_test.stan')

estBetaParams <- function(mu, var) {
  #https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams(0.2, 0.1)

ggplot(data = data.frame(x = c(0, 1)), aes(x)) +
  stat_function(fun = dbeta, n = 101, args = list(shape1 = 1, shape2 = 1), colour = 'black') +
  stat_function(fun = dbeta, n = 101, args = list(shape1 = 10, shape2 = 2), colour = 'red') +
  stat_function(fun = dbeta, n = 101, args = list(shape1 = 2, shape2 = 10), colour = 'blue') +
  stat_function(fun = dbeta, n = 101, args = list(shape1 = 0.12, shape2 = 0.48), colour = 'purple')


model <- model_program$sample(data = list(beta_a = 0.12,
                                          beta_b = 0.48,
                                          n_total_category = 34,
                                          n_total_background = 758,
                                          n_genes_category = 2,
                                          n_genes_background = 73,
                                          sample_prior_only = FALSE), 
                              refresh = 0)
model$summary() %>%
  as_tibble()

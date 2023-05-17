# srun -t 24:00:00 --nodes=1 --cpus-per-task=24 --mem=250G --pty /bin/bash
# module load singularity/3.5.3
# INLA="/work/vollmer/software/inla.sif"
# singularity shell --bind /work,/scratch,/tmp ${INLA}
# R

library(tidyverse)
library(brms)

model_choice <- 'fixed'

data_file <- '/scratch/j.selwyn/comparative_genomes/pathdata_for_discovery.csv'
out_model <- str_c('/scratch/j.selwyn/comparative_genomes/kegg_pathway_completeness', 
                   model_choice, sep = '_')
the_data <- read_csv(data_file, show_col_types = FALSE)

#
if(model_choice == 'random'){
  the_form <- bf(n_keggs | trials(total_keggs) ~ species + 
                   (1 + species | major_category / 
                      minor_category / kegg_path_id))
  
  the_prior <- prior(normal(0, 1), class = 'b') +
    prior(gamma(4, 2), class = 'sd') +
    prior(lkj(2), class = 'cor')
} else {
  the_form <- bf(n_keggs | trials(total_keggs) ~ 
                   species * (kegg_path_id))
  #major_category / minor_category / #can't include in model directly because of issue with long vectors
  #maybe try including as just + major_category + minor_category - therefore the interaction with species is only with pathID??
  the_prior <- prior(normal(0, 1), class = 'b')
}


path_complete_model <- brm(the_form,
                           prior = the_prior,
                           family = beta_binomial(link = 'logit',
                                                  link_phi = 'log'),
                           
                           data = the_data,
                           sample_prior = 'yes',
                           
                           chains = 4,
                           iter = 2000,
                           warmup = 1000,
                           
                           control = list(adapt_delta = if_else(model_choice == 'random', 0.95, 0.8), 
                                          max_treedepth = 10),
                           
                           cores = 4,
                           threads = if_else(model_choice == 'random', 6L, 6L),
                           backend = 'cmdstanr',
                           
                           file = out_model,
                           file_refit = 'on_change')

if(model_choice == 'fixed'){
  # path_complete_model <- read_rds(str_c(out_model, '.rds'))
  library(emmeans)
  path_ref_grid <- ref_grid(path_complete_model)
  # path_ref_grid@levels$kegg_path_id
  minor_category_fct <- select(the_data, minor_category, kegg_path_id) %>%
    distinct %>%
    arrange(kegg_path_id) %>% #they are in alphabetical order
    pull(minor_category) %>%
    fct_inorder()
    
  major_category_fct <- select(the_data, major_category, minor_category, kegg_path_id) %>%
    distinct %>%
    arrange(kegg_path_id) %>%
    mutate(major_category = fct_inorder(major_category)) %>%
    select(major_category, minor_category) %>%
    distinct %>%
    pull(major_category)
  
  path_ref_grid_new <- add_grouping(path_ref_grid, 'minor_category', 'kegg_path_id', 
                                    minor_category_fct) %>%
    add_grouping('major_category', 'minor_category', major_category_fct)
  
  major_cat_emmeans <- emmeans(path_ref_grid_new, ~species * major_category)
  write_rds(major_cat_emmeans, str_c(out_model, '_major.rds.gz'))
  
  minor_cat_emmeans <- emmeans(path_ref_grid_new, ~species * minor_category)
  write_rds(minor_cat_emmeans, str_c(out_model, '_minor.rds.gz'))
  
  # minor_cat_emmeans <- read_rds(str_c(out_model, '_minor.rds.gz'))
  minor_cats <- summary(minor_cat_emmeans, type = 'response') %>%
    broom::tidy(conf.int = TRUE)
  
  for(CAT in levels(minor_category_fct)){
    message('\n')
    message('Starting to get emmeans for minor category: ', CAT, '\n')
    message(Sys.time())
    kegg_path_emmeans <- emmeans(path_ref_grid_new, ~species * kegg_path_id,
                                 at = list(minor_category = CAT))
    write_rds(kegg_path_emmeans, str_c(out_model, '_keggPath', 
                                       str_replace_all(CAT, ' ', '.'), 
                                       '.rds.gz'))
    
    message('Finished getting emmeans for minor category: ', CAT, '\n')
    message(Sys.time())
    message('\n')
  }
  
}

library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)


busco_data <- read_csv('../Data/summarized_buscos.csv', show_col_types = FALSE) %>%
  select(-complete) %>%
  mutate(Y = cbind(`Complete Single Copy`, `Complete Duplicated`, Fragmented, Missing)) %>%
  select(species, n, Y)


multi_model_busco <- brm(Y | trials(n) ~ (1 | species),
           family = multinomial(),
           data = busco_data,
           chains = 4,
           cores = 4,
           iter = 2000,
           warmup = 1000,
           backend = 'cmdstanr')

busco_data %>%
  select(species) %>%
  mutate(n = 1) %>%
  add_epred_draws(multi_model_busco) %>%
  filter(.category != 'Missing') %>%
  group_by(species, .draw) %>%
  mutate(cum_pred = cumsum(.epred)) %>%
  ungroup %>%
  select(species, .category, .epred, cum_pred) %>%
  group_by(species, .category) %>%
  point_interval() %>%
  ungroup %>%
  
  ggplot(aes(x = species, y = cum_pred, ymin = cum_pred.lower, ymax = cum_pred.upper)) +
  geom_col(data = busco_data %>%
             mutate(Y = Y / n) %>%
             select(-n) %>%
             rowwise(species) %>%
             summarise(as_tibble(Y),
                       .groups = 'drop') %>%
             pivot_longer(cols = -species) %>%
             mutate(name = fct_inorder(name) %>% fct_rev()), 
           inherit.aes = FALSE,
           aes(x = species, y = value, fill = name)) +
  geom_pointrange() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = ~str_remove_all(.x, '_.*')) +
  labs(x = NULL,
       y = 'BUSCOs (%)',
       fill = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.text = element_text(colour = 'black'))




dirich_model_busco <- brm(Y / n | weights(n) ~ (1 | species),
                         family = dirichlet(),
                         data = busco_data,
                         chains = 4,
                         cores = 4,
                         iter = 2000,
                         warmup = 1000,
                         backend = 'cmdstanr')


busco_data %>%
  select(species) %>%
  mutate(n = 1) %>%
  add_epred_draws(multi_model_busco) %>%
  filter(.category != 'Missing') %>%
  group_by(species, .draw) %>%
  mutate(cum_pred = cumsum(.epred)) %>%
  ungroup %>%
  select(species, .category, .epred, cum_pred) %>%
  group_by(species, .category) %>%
  point_interval() %>%
  ungroup %>%
  
  ggplot(aes(x = species, y = cum_pred, ymin = cum_pred.lower, ymax = cum_pred.upper)) +
  geom_col(data = busco_data %>%
             mutate(Y = Y / n) %>%
             select(-n) %>%
             rowwise(species) %>%
             summarise(as_tibble(Y),
                       .groups = 'drop') %>%
             pivot_longer(cols = -species) %>%
             mutate(name = fct_inorder(name) %>% fct_rev()), 
           inherit.aes = FALSE,
           aes(x = species, y = value, fill = name)) +
  geom_pointrange() +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = ~str_remove_all(.x, '_.*')) +
  labs(x = NULL,
       y = 'BUSCOs (%)',
       fill = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom')


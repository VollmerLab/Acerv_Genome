data {
  real<lower=0> beta_a; 
  real<lower=0> beta_b; 
  int<lower=0> n_total_category; 
  int<lower=0> n_total_background; //the total genes in major category in the background
  int<lower=0> n_genes_category; //the number of genes in minor category in the category (unique or duplicated)
  int<lower=0> n_genes_background; //the number of genes in minor category in the background
  int<lower=0, upper=1> sample_prior_only;
}
parameters { 
  real<lower=0, upper=1> theta_category;
  real<lower=0, upper=1> theta_background;
}
model {
  theta_category ~ beta(beta_a, beta_b);
  theta_background ~ beta(beta_a, beta_b);
  
  if (sample_prior_only != 1) {
    n_genes_category ~ binomial(n_total_category, theta_category);
    n_genes_background ~ binomial(n_total_background, theta_background);
  }
}
generated quantities {
  real diff;
  real OR;
  diff = theta_category - theta_background;
  OR = (theta_category / (1 - theta_category)) / (theta_background / (1 - theta_background));
}

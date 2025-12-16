data {
  int<lower = 0> K_pos;
  array[K_pos] int<lower = 0> N_pos;
  array[K_pos] int<lower = 0> n_pos;

  int<lower = 0> K_neg;
  array[K_neg] int<lower = 0> N_neg;
  array[K_neg] int<lower = 0> n_neg;

  int<lower = 0> K_unk;
  array[K_unk] int<lower = 0> N_unk;
  array[K_unk] int<lower = 0> n_unk;

  real<lower = 0> sigma_sigma_logit_sens;
  real<lower = 0> sigma_sigma_logit_spec;
}
parameters {
  real mu_logit_sens;
  real mu_logit_spec;
  real<lower = 0> sigma_logit_sens;
  real<lower = 0> sigma_logit_spec;

  vector<offset = mu_logit_sens,
         multiplier = sigma_logit_sens>[K_pos] logit_sens;
  vector<offset = mu_logit_spec,
         multiplier = sigma_logit_spec>[K_neg] logit_spec;
  vector<offset = mu_logit_sens,
         multiplier = sigma_logit_sens>[K_unk] logit_sens_unk;
  vector<offset = mu_logit_spec,
         multiplier = sigma_logit_spec>[K_unk] logit_spec_unk;

  real<lower = 0, upper = 1> pi;
}
transformed parameters {
  vector[K_pos] sens = inv_logit(logit_sens);
  vector[K_neg] spec = inv_logit(logit_spec);
  vector[K_unk] sens_unk = inv_logit(logit_sens_unk);
  vector[K_unk] spec_unk = inv_logit(logit_spec_unk);
}
model {
  mu_logit_sens ~ normal(4, 2);
  sigma_logit_sens ~ normal(0, sigma_sigma_logit_sens);
  logit_sens ~ normal(mu_logit_sens, sigma_logit_sens);
  n_pos ~ binomial_logit(N_pos, logit_sens);

  mu_logit_spec ~ normal(4, 2);
  sigma_logit_spec ~ normal(0, sigma_sigma_logit_spec);
  logit_spec ~ normal(mu_logit_spec, sigma_logit_spec);
  n_neg ~ binomial_logit(N_neg, logit_spec);

  pi ~ uniform(0, 1);
  logit_sens_unk ~ normal(mu_logit_sens, sigma_logit_sens);
  logit_spec_unk ~ normal(mu_logit_spec, sigma_logit_spec);
  n_unk ~ binomial(N_unk, pi * sens_unk + (1 - pi) * (1 - spec_unk));
}

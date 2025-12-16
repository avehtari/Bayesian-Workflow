data {
  int y_sample;
  int n_sample;
  int J_spec;
  int J_sens;
  array[J_spec] int y_spec;
  array[J_spec] int n_spec;
  array[J_sens] int y_sens;
  array[J_sens] int n_sens;
  real logit_spec_prior_scale;
  real logit_sens_prior_scale;
}
parameters {
  real<lower=0, upper=1> p;
  real mu_logit_spec;
  real mu_logit_sens;
  real<lower=0> sigma_logit_spec;
  real<lower=0> sigma_logit_sens;
  vector<offset=mu_logit_spec, multiplier=sigma_logit_spec>[J_spec] logit_spec;
  vector<offset=mu_logit_sens, multiplier=sigma_logit_sens>[J_sens] logit_sens;
}
transformed parameters {
  vector[J_spec] spec = inv_logit(logit_spec);
  vector[J_sens] sens = inv_logit(logit_sens);
  real p_sample = p * sens[1] + (1 - p) * (1 - spec[1]);
}
model {
  y_sample ~ binomial(n_sample, p_sample);
  y_spec ~ binomial(n_spec, spec);
  y_sens ~ binomial(n_sens, sens);
  logit_spec ~ normal(mu_logit_spec, sigma_logit_spec);
  logit_sens ~ normal(mu_logit_sens, sigma_logit_sens);
  sigma_logit_spec ~ normal(0, logit_spec_prior_scale);
  sigma_logit_sens ~ normal(0, logit_sens_prior_scale);
  mu_logit_spec ~ normal(4, 2); // weak prior on mean of dist of specificity
  mu_logit_sens ~ normal(4, 2); // weak prior on mean of dist of sensitivity
}
generated quantities {
  real log_lik = binomial_lpmf(y_sample | n_sample, p_sample);
  real lprior = normal_lpdf(mu_logit_spec | 4, 2) +
                normal_lpdf(mu_logit_sens | 4, 2); 
}

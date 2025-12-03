data {
  int<lower = 0> y_sample;
  int<lower = 0> n_sample;
  int<lower = 0> J_specificity;
  int<lower = 0> J_sensitivity;
  array[J_specificity] int<lower = 0> y_specificity;
  array[J_specificity] int<lower = 0> n_specificity;
  array[J_sensitivity] int<lower = 0> y_sensitivity;
  array[J_sensitivity] int<lower = 0> n_sensitivity;
  real<lower = 0> logit_specificity_prior_scale;
  real<lower = 0> logit_sensitivity_prior_scale;
}
parameters {
  real<lower = 0, upper = 1> prevalence;
  real mu_logit_specificity;
  real mu_logit_sensitivity;
  real<lower = 0> sigma_logit_specificity;
  real<lower = 0> sigma_logit_sensitivity;
  vector<offset = mu_logit_specificity, multiplier = sigma_logit_specificity>[J_specificity] logit_specificity;
  vector<offset = mu_logit_sensitivity, multiplier = sigma_logit_sensitivity>[J_sensitivity] logit_sensitivity;
}
transformed parameters {
  vector[J_specificity] specificity = inv_logit(logit_specificity);
  vector[J_sensitivity] sensitivity = inv_logit(logit_sensitivity);
  real prevalence_sample = prevalence * sensitivity[1] + (1 - prevalence) * (1 - specificity[1]);
}
model {
  y_sample ~ binomial(n_sample, prevalence_sample);
  y_specificity ~ binomial(n_specificity, specificity);
  y_sensitivity ~ binomial(n_sensitivity, sensitivity);
  logit_specificity ~ normal(mu_logit_specificity, sigma_logit_specificity);
  logit_sensitivity ~ normal(mu_logit_sensitivity, sigma_logit_sensitivity);
  sigma_logit_specificity ~ normal(0, logit_specificity_prior_scale);
  sigma_logit_sensitivity ~ normal(0, logit_sensitivity_prior_scale);
  mu_logit_specificity ~ normal(4, 2); // weak prior on mean of distribution of specificity
  mu_logit_sensitivity ~ normal(4, 2); // weak prior on mean of distribution of sensitivity
}
generated quantities {
  real log_lik = binomial_lpmf(y_sample | n_sample, prevalence_sample);
  real lprior = normal_lpdf(mu_logit_specificity | 4, 2) +
                normal_lpdf(mu_logit_sensitivity | 4, 2); 
}

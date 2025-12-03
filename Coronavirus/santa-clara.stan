data {
  int<lower = 0> y_sample;
  int<lower = 0> n_sample;
  int<lower = 0> y_specificity;
  int<lower = 0> n_specificity;
  int<lower = 0> y_sensitivity;
  int<lower = 0> n_sensitivity;
}
parameters {
  real<lower=0, upper = 1> prevalence;
  real<lower=0, upper = 1> specificity;
  real<lower=0, upper = 1> sensitivity;
}
transformed parameters {
  real prevalence_sample = prevalence * sensitivity +
                           (1 - prevalence) * (1 - specificity);
}
model {
  y_sample ~ binomial(n_sample, prevalence_sample);
  y_specificity ~ binomial(n_specificity, specificity);
  y_sensitivity ~ binomial(n_sensitivity, sensitivity);
}
generated quantities {
  real log_lik = binomial_lpmf(y_sample | n_sample, prevalence_sample);
  real lprior = binomial_lpmf(y_specificity | n_specificity, specificity) +
       	        binomial_lpmf(y_sensitivity | n_sensitivity, sensitivity);   
}

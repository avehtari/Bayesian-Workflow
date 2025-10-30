data {
  int N;   // number of observations
  int J;   // number of students
  int K;   // number of items on exam
  array[N] int<lower=0, upper=J> student;
  array[N] int<lower=0, upper=K> item;
  array[N] int<lower=0, upper=1> y;
  real mu_mu_beta;
  real<lower=0> sigma_mu_beta, mu_sigma_alpha, mu_sigma_beta, mu_sigma_gamma;
}
parameters {
  real mu_beta;
  real<lower=0> sigma_alpha, sigma_beta, sigma_gamma;
  vector<offset=0, multiplier=sigma_alpha>[J] alpha;
  vector<offset=mu_beta, multiplier=sigma_beta>[K] beta;
  vector<offset=1, multiplier=sigma_gamma>[K] gamma;
}
model {
  alpha ~ normal(0, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta);
  gamma ~ normal(1, sigma_gamma);
  mu_beta ~ normal(mu_mu_beta, sigma_mu_beta);
  sigma_alpha ~ exponential(1/mu_sigma_alpha);
  sigma_beta ~ exponential(1/mu_sigma_beta);
  sigma_gamma ~ exponential(1/mu_sigma_gamma);
  y ~ bernoulli(0.25 + 0.75*inv_logit(gamma[item] .* (alpha[student] - beta[item])));
}

data {
  int N;   // number of observations
  int J;   // number of students
  int K;   // number of items on exam
  array[N] int<lower=0, upper=J> student;
  array[N] int<lower=0, upper=K> item;
  array[N] int<lower=0, upper=1> y;
  vector[J] x;
  real mu_mu_a, mu_mu_b;
  real<lower=0> sigma_mu_a, sigma_mu_b, mu_sigma_a, mu_sigma_b;
}
transformed data {
  vector[J] x_adj = (x - mean(x))/sd(x);
}
parameters {
  real mu_a, mu_b;
  real<lower=0> sigma_a, sigma_b;
  vector<offset=mu_a, multiplier=sigma_a>[K] a;
  vector<offset=mu_b, multiplier=sigma_b>[K] b;
}
model {
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  mu_a ~ normal(mu_mu_a, sigma_mu_a);
  mu_b ~ normal(mu_mu_b, sigma_mu_b);
  sigma_a ~ exponential(1/mu_sigma_a);
  sigma_b ~ exponential(1/mu_sigma_b);
  y ~ bernoulli(0.25 + 0.75*inv_logit(a[item] + b[item] .* x_adj[student]));
}

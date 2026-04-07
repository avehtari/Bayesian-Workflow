data {
  int J;
  array[J] int<lower=0, upper=1> y;
  vector[J] x;
  real mu_a, mu_b;
  real<lower=0> sigma_a, sigma_b;
}
parameters {
  real a, b;
}
model {
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ bernoulli(0.25 + 0.75*inv_logit(a + b*x));
}

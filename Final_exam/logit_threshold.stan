data {
  int J;
  array[J] int<lower=0, upper=1> y;
  vector[J] x;
  real mu_a, mu_b;
  real<lower=0> sigma_a, sigma_b;
}
transformed data {
  vector[J] x_adj = (x - mean(x))/sd(x);
}
parameters {
  real a, b;
  real<lower=0, upper=1> p0;
}
model {
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ bernoulli(p0 + (1-p0)*inv_logit(a + b*x_adj));
}

data {
  int N;
  vector[N] x;
  vector[N] y;
  real mu_a, mu_b;
  real<lower=0> sigma_a, sigma_b;
}
parameters {
  real a, b;
  real<lower=0> sigma;
}
model {
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ normal(a + b*x, sigma);
}

data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real a, b;
  real<lower=0> sigma;
}
model {
  y ~ normal(a + b*x, sigma);
}

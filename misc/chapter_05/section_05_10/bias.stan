data {
  int N;
  vector[N] y;
}
parameters {
  real theta, bias;
  real<lower=0> sigma;
}
model {
  theta ~ normal(0, 1);
  y ~ normal(theta + bias, sigma);
  bias ~ cauchy(0, 0.1);
}

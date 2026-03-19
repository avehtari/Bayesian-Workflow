data {
  int N;
  vector[N] y;
}
parameters {
  real theta;
  real<lower=0> sigma;
}
model {
  theta ~ cauchy(0, 1);
  y ~ normal(theta, sigma);
}

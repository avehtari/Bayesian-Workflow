// logistic regression
data {
  int<lower=0> N;
  int<lower=0> M;
  vector[N] y;
  matrix[N,M] x;
}
parameters {
  real alpha;
  vector[M] beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  y ~ normal_id_glm(x, alpha, beta, sigma);
}

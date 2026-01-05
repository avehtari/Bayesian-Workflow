// logistic regression
data {
  int<lower=0> N;
  int<lower=0> M;
  array[N] int<lower=0,upper=1> y;
  matrix[N,M] x;
}
parameters {
  real alpha;
  vector[M] beta;
  real gamma;
}
model {
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  y ~ bernoulli_logit_glm(x, alpha, beta);
}

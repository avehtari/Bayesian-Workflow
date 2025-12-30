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
}
model {
  alpha ~ cauchy(0, 10);
  beta ~ cauchy(0, 10);
  y ~ bernoulli_logit_glm(x, alpha, beta);
}

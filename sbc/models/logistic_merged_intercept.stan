data {
  int<lower=0> N_obs;
  array[N_obs] int<lower=0, upper=1> y;
  int<lower=1> N_predictors;
  matrix[N_obs, N_predictors] X;
}
parameters {
  vector[N_predictors] beta;
}
model {
  target += bernoulli_logit_lpmf(y | X * beta);
  target += normal_lpdf(beta[1] | 0, 2);
  target += normal_lpdf(beta | 0, 1);
}

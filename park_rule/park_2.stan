data {
  int<lower=0> N, J, K, L;
  array[N] int<lower=0, upper=1> y;
  array[N] int<lower=1, upper=J> respondent;
  array[N] int<lower=1, upper=K> item;
  matrix[N,L] X;
}
parameters {
  vector[L] b;
  real<lower=0> sigma_respondent, sigma_item;
  vector[J] a_respondent;
  vector[K] a_item;
}
model {
  a_respondent ~ normal(0, sigma_respondent);
  a_item ~ normal(0, sigma_item);
  y ~ bernoulli_logit(X*b + a_respondent[respondent] + a_item[item]);
}

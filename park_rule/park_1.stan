data {
  int<lower=0> N, J, K, L;
  array[N] int<lower=0, upper=1> y;
  array[N] int<lower=1, upper=J> respondent;
  array[N] int<lower=1, upper=K> item;
  matrix[N, L] X;
}
parameters {
  real a;
  vector[L] b;
  real<lower=0> sigma_respondent, sigma_item;
  vector<multiplier=sigma_respondent>[J] a_respondent;
  vector<multiplier=sigma_item>[K] a_item;
}
model {
  a_respondent ~ normal(0, sigma_respondent);
  a_item ~ normal(0, sigma_item);
  b ~ normal(0, 1);
  {sigma_respondent, sigma_item} ~ normal(0, 3);
  y ~ bernoulli_logit_glm(X, a + a_respondent[respondent] + a_item[item], b);
}

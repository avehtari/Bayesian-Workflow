data {
  int<lower=0> N, J, K, L;
  array[N] int<lower=0, upper=1> y;
  array[N] int<lower=1, upper=J> respondent;
  array[N] int<lower=1, upper=K> item;
  matrix[N, L] X;
}
transformed data {
  matrix[N, L] X_c;
  for (l in 1:L) {
    X_c[, l] = X[, l] - mean(X[, l]);
  }
}
parameters {
  real a;
  vector[L] b;
  real<lower=0> sigma_respondent, sigma_item;
  sum_to_zero_vector[J] a_respondent;
  sum_to_zero_vector[K] a_item;
}
model {
  a_respondent ~ normal(0, sigma_respondent);
  a_item ~ normal(0, sigma_item);
  b ~ normal(0, 1);
  {sigma_respondent, sigma_item} ~ normal(0, 3);
  y ~ bernoulli_logit_glm(X_c, a + a_respondent[respondent] + a_item[item], b);
}

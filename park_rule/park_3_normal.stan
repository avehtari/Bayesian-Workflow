data {
  int<lower=0> N, J, K, L;
  vector[N] y;
  array[N] int<lower=1, upper=J> respondent;
  array[N] int<lower=1, upper=K> item;
  matrix[N,L] X;
}
parameters {
  vector[L] b;
  real<lower=0> sigma_respondent, sigma_item, sigma_y;
  sum_to_zero_vector[J] a_respondent;
  sum_to_zero_vector[K] a_item;
}
model {
  a_respondent ~ normal(0, sigma_respondent);
  a_item ~ normal(0, sigma_item);
  y ~ normal(X*b + a_respondent[respondent] + a_item[item], sigma_y);
}

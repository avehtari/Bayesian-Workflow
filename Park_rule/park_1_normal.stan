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
  vector<offset=0, multiplier=sigma_respondent>[J] a_respondent;
  vector<offset=0, multiplier=sigma_item>[K] a_item;
}
model {
  a_respondent ~ normal(0, sigma_respondent);
  a_item ~ normal(0, sigma_item);
  y ~ normal(X*b + a_respondent[respondent] + a_item[item], sigma_y);
}

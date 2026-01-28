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
  sum_to_zero_vector[J] z_respondent;
  sum_to_zero_vector[K] z_item;
}
transformed parameters {
  vector[J] a_respondent = z_respondent * sigma_respondent;
  vector[K] a_item = z_item * sigma_item;
}
model {
  z_respondent ~ std_normal();
  z_item ~ std_normal();
  b ~ normal(0, 1);
  {sigma_respondent, sigma_item} ~ normal(0, 3);
  y ~ bernoulli_logit_glm(X, a + a_respondent[respondent] + a_item[item], b);
}

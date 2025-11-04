data {
  int N;
  vector[N] y;
  int J;
  int K;
  array[N] int<lower=1, upper=J> movie;
  array[N] int<lower=1, upper=K> rater;
}
parameters {
  vector[J] alpha;
  vector[K] beta;
  real mu;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[J] a = mu + sigma_a * alpha;
}
model {
  y ~ normal(a[movie] - sigma_b * beta[rater], sigma_y);
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  mu ~ normal(3, 5);
  sigma_a ~ normal(0, 5);
  sigma_b ~ normal(0, 5);
  sigma_y ~ normal(0, 5);
}

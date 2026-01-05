data {
  int N;
  vector[N] y;
  int J;
  array[N] int<lower=1, upper=J> movie;
}
parameters {
  vector<lower=0, upper=5>[J] theta;
  real<lower=0> sigma_y;
}
model {
  theta ~ normal(3, 1);
  sigma_y ~ normal(0, 2.5);
  y ~ normal(theta[movie], sigma_y);
}

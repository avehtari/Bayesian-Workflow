data {
  int<lower=0> N;                   // number of observations
  int<lower=0> K;                   // number of groups
  array[N] int<lower=1, upper=K> x; // discrete group indicators
  vector[N] y;                      // real valued observations
}
parameters {
  real mu0;                         // prior mean
  real<lower=0> sigma0;             // prior sd
  vector[K] z;                      // latent variable
  real<lower=0> sigma;              // common sd
}
transformed parameters {
  vector[K] mu = mu0 + sigma0 * z;  // group means
}
model {
  mu0 ~ normal(10, 10);             // weakly informative prior
  sigma0 ~ normal(0, 10);           // weakly informative prior
  z ~ normal(0, 1);                 // unit normal
  sigma ~ lognormal(0, .5);         // weakly informative prior
  y ~ normal(mu[x], sigma);         // observation model
}

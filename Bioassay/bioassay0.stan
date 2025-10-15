data {
  int<lower=0> J;  
  vector[J] x;
  array[J] int<lower=0> N;
  array[J] int<lower=0, upper=N> y;
}
parameters {
  real a;
  real b;
}
model {
  y ~ binomial_logit(N, a + b * x);
}

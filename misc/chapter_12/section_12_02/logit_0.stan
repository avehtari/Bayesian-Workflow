data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> y;
}
parameters {
  real a;
}
model {
  y ~ bernoulli_logit(a);
}


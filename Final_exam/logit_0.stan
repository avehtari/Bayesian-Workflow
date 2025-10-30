data {
  int J;
  array[J] int<lower=0, upper=1> y;
  vector[J] x;
}
parameters {
  real a, b;
}
model {
  y ~ bernoulli_logit(a + b*x);
}

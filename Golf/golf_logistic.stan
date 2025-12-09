data {
  int J;
  vector[J] x;
  array[J] int n;
  array[J] int y;
}
parameters {
  real a;
  real b;
}
model {
  y ~ binomial_logit(n, a + b*x);
}

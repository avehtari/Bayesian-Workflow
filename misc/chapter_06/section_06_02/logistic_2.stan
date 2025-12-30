data {
  int J;  
  vector[J] x;
  array[J] int n, y;
  real mu_a, mu_b, sigma_a, sigma_b;
}
parameters {
  real a, b;
}
model {
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ binomial_logit(n, a + b*x);
}

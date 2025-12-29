data {
  int N;
  vector[N] x;
  vector[N] y;
  int N_tilde;
  vector[N_tilde] x_tilde;
}
parameters {
  real a, b;
  real<lower=0> sigma;
}
model {
  y ~ normal(a + b*x, sigma); 
}
generated quantities {
  array[N_tilde] real y_tilde = normal_rng(a + b*x_tilde, sigma);
}

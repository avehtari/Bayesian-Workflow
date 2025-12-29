functions {
  vector Ey(vector b, vector x, vector z) {
    return inv_logit(b[1] + b[2]*x + b[3]*z + b[4]*x.*z);
  }
}
data {
  int<lower=0> N, K;
  vector[N] x, z;
  array[N] int y;
  int N_pop;
  vector[N_pop] x_pop, n_pop;
}
parameters {
  vector[K] b;
}
model {
  y ~ bernoulli(Ey(b, x, z));
}
generated quantities{
  real SATE = mean(Ey(b, x, rep_vector(1, N)) - Ey(b, x, rep_vector(0, N)));
  real PATE = sum(n_pop .* (Ey(b, x_pop, rep_vector(1, N_pop)) - Ey(b, x_pop, rep_vector(0, N_pop)))) / sum(n_pop);
}

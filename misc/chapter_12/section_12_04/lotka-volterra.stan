functions {
  vector dz_dt(real t,             // time
               vector z,           // system state {prey, predator}
               real alpha,
               real beta,
               real gamma,
               real delta) {
    real u = z[1];
    real v = z[2];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;
    return to_vector({ du_dt, dv_dt });
  }
}
data {
  int<lower = 0> N;           // number of measurement times
  array[N] real ts;                 // measurement times > 0
  array[2] real y_init;             // initial measured populations
  array[N,2] real<lower = 0> y;    // measured populations
}
parameters {
  array[4] real<lower = 0> theta;   // { alpha, beta, gamma, delta }
  array[2] real<lower = 0> z_init;  // initial population
  array[2] real<lower = 0> sigma;   // measurement errors
}
transformed parameters {
  array[N] vector[2] z
    = ode_rk45_tol(dz_dt, to_vector(z_init), 0, ts,
                   1e-5, 1e-3, 500,
                   theta[1], theta[2], theta[3], theta[4]);
}
model {
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(10), 1);
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(z_init[k]), sigma[k]);
    for (n in 1:N)
      y[n, k] ~ lognormal(log(z[n][k]), sigma[k]);
  }
}
generated quantities {
  array[2] real y_init_rep;
  array[N, 2] real y_rep;
  for (k in 1:2) {
    y_init_rep[k] = lognormal_rng(log(z_init[k]), sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(z[n][k]), sigma[k]);
  }
}

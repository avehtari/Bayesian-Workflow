data {
  int<lower=0> N;
  array[N] int y;
}
parameters {
  ordered[2] mu;
  real<lower=0, upper=1> theta;
}
model {
  for(n in 1:N) {
    target += log_mix(theta,
                      poisson_log_lpmf(y[n] | mu[1]),
                      poisson_log_lpmf(y[n] | mu[2]));
  }
  target += normal_lpdf(mu | 3, 1);
}


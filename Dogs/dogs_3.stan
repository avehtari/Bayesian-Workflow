data {
  int<lower=0> J;
  int<lower=0> T;
  array[J, T] int<lower=0, upper=1> y;
}
parameters {
  vector[J] logit_a;
  real mu_logit_a;
  real<lower=0> sigma_logit_a;
}
transformed parameters {
  vector[J] a = inv_logit(logit_a);
  matrix<lower=0, upper=1>[J, T] p;
  for (j in 1:J) {
    for (t in 1:T) {
      p[j, t] = a[j]^(t-1);
    }
  }
}
model {
  for (j in 1:J) {
    for (t in 1:T) {
      y[j, t] ~ bernoulli(p[j, t]);
    }
  }
  logit_a ~ normal(mu_logit_a, sigma_logit_a);
}
generated quantities {
  array[J, T] int<lower=0, upper=1> y_rep;
  for (j in 1:J) {
    for (t in 1:T) {
      y_rep[j, t] = bernoulli_rng(p[j, t]);
    }
  }
}


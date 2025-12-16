data {
  int<lower=0> J;
  int<lower=0> T;
  array[J, T] int<lower=0, upper=1> y;
}
parameters {
  real<lower=0, upper=1> a;
}
transformed parameters {
  matrix<lower=0, upper=1>[J, T] p;
  for (j in 1:J) {
    for (t in 1:T) {
      p[j, t] = a^(t-1);
    }
  }
}
model {
  for (j in 1:J) {
    for (t in 1:T) {
      y[j, t] ~ bernoulli(p[j, t]);
    }
  }
}
generated quantities {
  array[J, T] int<lower=0, upper=1> y_rep;
  for (j in 1:J) {
    for (t in 1:T) {
      y_rep[j, t] = bernoulli_rng(p[j, t]);
    }
  }
}


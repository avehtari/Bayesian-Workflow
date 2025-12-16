data {
  int<lower=0> J;
  int<lower=0> T;
  array[J, T] int<lower=0, upper=1> y;
}
transformed data {
  matrix<lower=0>[J, T] prev_shock;
  matrix<lower=0>[J, T] prev_avoid;
  for (j in 1:J) {
    prev_shock[j, 1] = 0;
    prev_avoid[j, 1] = 0;
    for (t in 2:T) {
      prev_shock[j, t] = prev_shock[j, t-1] + y[j, t-1];
      prev_avoid[j, t] = prev_avoid[j, t-1] + 1 - y[j, t-1];
    }
  }
}
parameters {
  matrix[J, 2] logit_ab;
  vector[2] mu_logit_ab;
  cov_matrix[2] Sigma_logit_ab;
}
transformed parameters {
  vector[J] a = inv_logit(logit_ab[, 1]);
  vector[J] b = inv_logit(logit_ab[, 2]);
}
model {
  for (j in 1:J) {
    for (t in 1:T) {
      real p = a[j]^prev_shock[j, t] * b[j]^prev_avoid[j, t];
      y[j, t] ~ bernoulli(p);
    }
    logit_ab[j, ] ~ multi_normal(mu_logit_ab, Sigma_logit_ab);
  }
}
generated quantities {
  array[J, T] int<lower=0, upper=1> y_rep;
  {
    real prev_shock_rep;
    real prev_avoid_rep;
    real p_rep;
    for (j in 1:J) {
      prev_shock_rep = 0;
      prev_avoid_rep = 0;
      y_rep[j, 1] = 1;
      for (t in 2:T) {
        prev_shock_rep = prev_shock_rep + y_rep[j, t-1];
        prev_avoid_rep = prev_avoid_rep + 1 - y_rep[j, t-1];
        p_rep = a[j]^prev_shock_rep * b[j]^prev_avoid_rep;
        y_rep[j, t] = bernoulli_rng(p_rep);
      }
    }
  }
}


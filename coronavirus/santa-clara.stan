data {
  int<lower = 0> y_sample, n_sample, y_spec, n_spec, y_sens, n_sens;
}
parameters {
  real<lower=0, upper = 1> p, spec, sens;
}
transformed parameters {
  real p_sample = p * sens + (1 - p) * (1 - spec);
}
model {
  y_sample ~ binomial(n_sample, p_sample);
  y_spec ~ binomial(n_spec, spec);
  y_sens ~ binomial(n_sens, sens);
}
generated quantities {
  real log_lik = binomial_lpmf(y_sample | n_sample, p_sample);
  real lprior = binomial_lpmf(y_spec | n_spec, spec) +
       	        binomial_lpmf(y_sens | n_sens, sens);   
}

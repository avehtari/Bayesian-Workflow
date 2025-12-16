functions {
  real integrand(real z, real notused, array[] real theta,
               array[] real X_j, array[] int data_j) {
    real sigma_epsilon = theta[1];
    real p_j = theta[2];
    int n_j = data_j[1];
    int y_j = data_j[2];
    real p = (z<1) ? exp(exponential_lpdf(z | 1/sigma_epsilon) + binomial_lpmf(y_j | n_j, p_j*(1-z))) : 0;
    return (is_inf(p) || is_nan(p)) ? 0 : p;
  }
}
data {
  int J;
  array[J] int n;
  array[J] int y;
}
parameters {
  real<lower=0> sigma_epsilon;
  vector[J] p_angle;
  vector[J] p_distance;
}
generated quantities {
  vector[J] log_lik;
  for (j in 1:J) {
    real p = p_angle[j] .* p_distance[j];
    real ps = p * (1-sigma_epsilon);
    log_lik[j] = log(integrate_1d(integrand,
                              0,
                              positive_infinity(), // this works better than 1!
                              append_array({sigma_epsilon}, {p}),
			      {0}, // not used, but an empty array not allowed
			      append_array({n[j]}, {y[j]}),
			      1e-6));
  }
}

functions {
  real integrand(real z, real notused, array[] real theta,
               array[] real X_j, array[] int data_j) {
    real sigma_epsilon = theta[1];
    real p_j = theta[2];
    int n_j = data_j[1];
    int y_j = data_j[2];
    real p = exp(exponential_lpdf(z | 1/sigma_epsilon) + binomial_lpmf(y_j | n_j, p_j*(1-z)));
    return (is_inf(p) || is_nan(p)) ? 0 : p;
  }
}
data {
  int J;
  array[J] int n;
  vector[J] x;
  array[J] int y;
  real r;
  real R;
}
transformed data {
  vector[J] threshold_angle = asin((R-r) ./ x);
  vector[J] raw_proportion = to_vector(y) ./ to_vector(n);
}
parameters {
  real<lower=0> sigma_angle;
  real<lower=0> sigma_distance;
  real<lower=0> sigma_epsilon;
  vector<lower=0, upper=1>[J] epsilon;
  real<lower=0> distance_tolerance;
  real<lower=0> overshot;
}
transformed parameters {
  vector[J] p_angle = 2*Phi(threshold_angle / sigma_angle) - 1;
  vector[J] p_distance = Phi((distance_tolerance - overshot) ./ ((x + overshot)*sigma_distance)) -
               Phi((- overshot) ./ ((x + overshot)*sigma_distance));
  vector[J] p = p_angle .* p_distance .* (1 - epsilon);
}
model {
  [sigma_angle, sigma_distance] ~ normal(0, 1);
  distance_tolerance ~ lognormal(log(3), .3);
  overshot ~ lognormal(log(1), .3);
  epsilon ~ exponential(1/sigma_epsilon);
  y ~ binomial(n, p);
}
generated quantities {
  real sigma_degrees = sigma_angle * 180 / pi();
  vector[J] residual = raw_proportion - p_angle .* p_distance .* (1 - sigma_epsilon);
  vector[J] log_lik;
  for (j in 1:J) {
    log_lik[j] = log(integrate_1d(integrand,
                              0.0,
                              positive_infinity(),
                              append_array({sigma_epsilon},{p_angle[j] .* p_distance[j]}),
			      {0}, // not used, but an empty array not allowed
			      append_array({n[j]},{y[j]}),
			      1e-3));
  }
}

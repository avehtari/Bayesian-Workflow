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
               Phi(-overshot ./ ((x + overshot)*sigma_distance));
  vector[J] p = p_angle .* p_distance .* (1 - epsilon);
}
model {
  [sigma_angle, sigma_distance] ~ normal(0, 1);
  distance_tolerance ~ lognormal(log(3), .2);
  overshot ~ lognormal(log(1), .2);
  epsilon ~ exponential(1/sigma_epsilon);
  y ~ binomial(n, p);
}
generated quantities {
  real sigma_degrees = sigma_angle * 180 / pi();
  vector[J] residual = raw_proportion - p_angle .* p_distance;
}

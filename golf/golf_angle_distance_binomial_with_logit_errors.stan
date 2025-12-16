data {
  int J;
  array[J] int n;
  vector[J] x;
  array[J] int y;
  real r;
  real R;
  real overshot;
  real distance_tolerance;
}
transformed data {
  vector[J] threshold_angle = asin((R-r) ./ x);
}
parameters {
  real<lower=0> sigma_angle;
  real<lower=0> sigma_distance;
  real<lower=0> sigma_eta;
  vector[J] eta;
}
model {
  eta ~ normal(0, 1);
  [sigma_angle, sigma_distance, sigma_eta] ~ normal(0, 1);
  vector[J] p_angle = 2*Phi(threshold_angle / sigma_angle) - 1;
  vector[J] p_distance = Phi((distance_tolerance - overshot) ./ ((x + overshot)*sigma_distance)) -
               Phi(-overshot ./ ((x + overshot)*sigma_distance));
  vector[J] p = inv_logit(logit(p_angle .* p_distance) + sigma_eta*eta);
  y ~ binomial(n, p);
}
generated quantities {
  real sigma_degrees = sigma_angle * 180 / pi();
}

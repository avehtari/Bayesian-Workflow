data {
  int K;
  int N;
  array[N] real y;
  array[K] real mu;
}
parameters {
  simplex[K] theta;
  real<lower=0> sigma;
}
model {
  array[K] real ps;
  sigma ~ exponential(0.1);
  for (n in 1:N) {
    for (k in 1:K) {
      ps[k] = log(theta[k]) + normal_lpdf(y[n] | mu[k], sigma);
    }
    target += log_sum_exp(ps);
  }
}

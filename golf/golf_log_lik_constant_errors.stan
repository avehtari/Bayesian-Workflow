data {
  int J;
  array[J] int n;
  array[J] int y;
}
parameters {
  vector[J] p;
}
generated quantities {
  vector[J] log_lik;
  for (j in 1:J) {
    log_lik[j] = binomial_lpmf(y[j] | n[j], p[j]);
  }
}

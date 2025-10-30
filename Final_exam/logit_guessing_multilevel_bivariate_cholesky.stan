data {
  int N;   // number of observations
  int J;   // number of students
  int K;   // number of items on exam
  array[N] int<lower=0, upper=J> student;
  array[N] int<lower=0, upper=K> item;
  array[N] int<lower=0, upper=1> y;
  vector[J] x;
  vector[2] mu_mu_ab;
  vector<lower=0>[2] sigma_mu_ab, mu_sigma_ab;
}
transformed data {
  vector[J] x_adj = (x - mean(x))/sd(x);
}
parameters {
  vector[2] mu_ab;
  vector<lower=0>[2] sigma_ab;
  array[K] vector[2] e_ab;
  cholesky_factor_corr[2] L_ab;
}
transformed parameters {
  vector[K] a;
  vector[K] b;
  for (k in 1:K) {
    a[k] = mu_ab[1] + sigma_ab[1] * e_ab[k][1]; 
    b[k] = mu_ab[2] + sigma_ab[2] * e_ab[k][2];
  }
}
model {
  e_ab ~ multi_normal_cholesky([0,0], L_ab);
  mu_ab ~ normal(mu_mu_ab, sigma_mu_ab);
  sigma_ab ~ exponential(1/mu_sigma_ab);
  L_ab ~ lkj_corr_cholesky(1);
  y ~ bernoulli(0.25 + 0.75*inv_logit(a[item] + b[item] .* x_adj[student]));
}
generated quantities {
  corr_matrix[2] Omega_ab = multiply_lower_tri_self_transpose(L_ab);
}

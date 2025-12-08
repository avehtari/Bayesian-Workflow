// Poisson regression with hierarchical intercept (varying effect)
functions {
  real integrand(real z, real notused, array[] real theta,
               array[] real X_i, array[] int y_i) {
    real sigmaz = theta[1];
    real mu_i = theta[2];
    real p = exp(normal_lpdf(z | 0, sigmaz) + poisson_log_lpmf(y_i | z + mu_i));
    return (is_inf(p) || is_nan(p)) ? 0 : p;
  }
}
data {
  int<lower=0> N;           // number of data points
  int<lower=0> P;           // number of covariates
  matrix[N,P] X;            // covariates
  array[N] int<lower=0> y;  // target
  vector[N] offsett;        // offset (offset variable name is reserved)
  real integrate_1d_reltol;
}
parameters {
  real alpha;               // intercept
  vector[P] beta;           // slope
  vector[N] z;              // individual intercept (varying effect)
  real<lower=0> sigmaz;     // prior scale for z
}
model {
  // priors
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  z ~ normal(0, sigmaz);
  sigmaz ~ normal(0, 1);
  // observation model
  y ~ poisson_log_glm(X, z+offsett+alpha, beta); 
}
generated quantities {
  // log_lik for PSIS-LOO
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] y_loorep;
  for (i in 1:N) {
    // z as posterior draws, this would be challenging for PSIS-LOO (and WAIC)
    // log_lik[i] = poisson_log_glm_lpmf({y[i]} | X[i,], z[i]+offsett[i]+alpha, beta);
    // posterior predictive replicates conditional on p(z[i] | y[i])
    // in theory these could be used with above (non-integrated) log_lik and PSIS-LOO
    // to get draws from LOO predictive distribution, but the distribution of the
    // ratios from above lok_lik is bad
    real mu_i = offsett[i] + alpha + X[i,]*beta;
    y_rep[i] = poisson_log_rng(z[i] + mu_i);
    // we can integrate each z[i] out with 1D adaptive quadrature to get more
    // stable log_lik and corresponding importance ratios
    log_lik[i] = log(integrate_1d(integrand,
                              negative_infinity(),
                              positive_infinity(),
                              append_array({sigmaz}, {mu_i}),
			      {0}, // not used, but an empty array not allowed
			      {y[i]},
			      integrate_1d_reltol));
    // conditional LOO predictive replicates conditional on p(z[i] | sigmaz, mu_i)
    // these combined with integrated log_lik and PSIS-LOO provide
    // more stable LOO predictive distributions
    y_loorep[i] = poisson_log_rng(normal_rng(0, sigmaz) + mu_i);
  }
}

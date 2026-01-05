// Gaussian linear model with adjustable priors
data {
  int<lower=0> N; // number of data points
  vector[N] x;    // covariate / predictor
  vector[N] y;    // target
  real pmualpha_c;// prior mean for alpha_c
  real psalpha;   // prior std for alpha
  real pmubeta;   // prior mean for beta
  real psbeta;    // prior std for beta
  real pssigma;   // prior std for half-normal prior for sigma
}
transformed data {
  // centering the predictor makes the posterior easier to sample
  real xmean = mean(x);
  vector[N] x_c = x - xmean;
}
parameters {
  real alpha_c;        // intercept for centered x
  real beta;           // slope
  real<lower=0> sigma; // standard deviation is constrained to be positive
}
model {
  alpha_c ~ normal(pmualpha_c, psalpha);   // prior
  beta ~ normal(pmubeta, psbeta);          // prior
  sigma ~ normal(0, pssigma);              // as sigma is constrained to be positive,
                                           // this is same as half-normal prior
  y ~ normal(alpha_c + beta*x_c, sigma); // observation model / likelihood
}

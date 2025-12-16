/* The error in the transformed data block has been fixed */
functions {
  real integrand(real z, real notused, array[] real theta,
               array[] real X_i, array[] int y_i) {
    real sigma_y = theta[1];
    real a_diff = theta[2];
    real sqrt_z = 2 * (step(z) - 0.5) * sqrt(abs(z));
    real p = exp(normal_lpdf(sqrt_z | a_diff, sigma_y) - log(abs(z))/2 -log(2));
    return (is_inf(p) || is_nan(p)) ? 0 : p;
  }
}

data {
  int N_teams;
  int N_games;
  vector[N_teams] prior_score;
  array[N_games] int team_1;
  array[N_games] int team_2;
  vector[N_games] score_1;
  vector[N_games] score_2;
  real df;
}
transformed data {
  vector[N_games] dif = score_1 - score_2;
  vector[N_games] sqrt_dif;
  for (i in 1:N_games){
    sqrt_dif[i] = 2 * (step(dif[i]) - 0.5) * sqrt(abs(dif[i]));
  }
}
parameters {
  vector[N_teams] alpha;
  real b;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[N_teams] a = b * prior_score + sigma_a * alpha;
}
model {
  alpha ~ normal(0, 1);
  b ~ normal(0, 1);
  sigma_a ~ normal(0, 1);
  sigma_y ~ normal(0, 1);
  for (n in 1:N_games) {
    target += (dif[n] == 0)
    ? log(integrate_1d(integrand,
                              -0.5,
                              0,
                              append_array(append_array({sigma_y}, {a[team_1[n]] - a[team_2[n]]}), {df}),
			      {0}, // not used, but an empty array not allowed
			      {0}, // not used, but an empty array not allowed
			      1e-5) +
	  integrate_1d(integrand,
                              0,
                              0.5,
                              append_array(append_array({sigma_y}, {a[team_1[n]] - a[team_2[n]]}), {df}),
			      {0}, // not used, but an empty array not allowed
			      {0}, // not used, but an empty array not allowed
			      1e-5))
    : log(integrate_1d(integrand,
                              dif[n]-0.5,
                              dif[n]+0.5,
                              append_array(append_array({sigma_y}, {a[team_1[n]] - a[team_2[n]]}), {df}),
			      {0}, // not used, but an empty array not allowed
			      {0}, // not used, but an empty array not allowed
			      1e-5));
//			      integrate_1d_reltol));
  }
  
}
generated quantities {
  vector[N_games] y_rep;
  vector[N_games] y_rep_original_scale;
  vector[N_games] log_lik;
  real lprior;
  for (n in 1:N_games) {
    y_rep[n] = normal_rng(a[team_1[n]] - a[team_2[n]], sigma_y);
    // log_lik[n] = student_t_lpdf(sqrt_dif[n] | df, a[team_1[n]] - a[team_2[n]], sigma_y) + (dif[n]==0 ? 0 : -log(abs(dif[n]))/2 -log(2));
    log_lik[n] = (dif[n] == 0)
    ? log(integrate_1d(integrand,
                              -0.5,
                              0,
                              append_array({sigma_y}, {a[team_1[n]] - a[team_2[n]]}),
			      {0}, // not used, but an empty array not allowed
			      {0}, // not used, but an empty array not allowed
			      1e-6) +
	  integrate_1d(integrand,
                              0,
                              0.5,
                              append_array({sigma_y}, {a[team_1[n]] - a[team_2[n]]}),
			      {0}, // not used, but an empty array not allowed
			      {0}, // not used, but an empty array not allowed
			      1e-6))
    : log(integrate_1d(integrand,
                              dif[n]-0.5,
                              dif[n]+0.5,
                              append_array({sigma_y}, {a[team_1[n]] - a[team_2[n]]}),
			      {0}, // not used, but an empty array not allowed
			      {0}, // not used, but an empty array not allowed
			      1e-6));
//			      integrate_1d_reltol));
  }
  y_rep_original_scale = y_rep .* abs(y_rep);
  lprior = normal_lpdf(b | 0, 1) +
           normal_lpdf(sigma_a | 0, 1) +
           normal_lpdf(sigma_y | 0, 1);
}

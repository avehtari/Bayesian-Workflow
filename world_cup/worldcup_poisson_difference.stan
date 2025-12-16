data {
  int N_teams;
  int N_games;
  vector[N_teams] prior_score;
  array[N_games] int team_1;
  array[N_games] int team_2;
  vector[N_games] score_1;
  vector[N_games] score_2;
}
transformed data {
  vector[N_games] dif = score_1 - score_2;
}
parameters {
  vector[N_teams] alpha;
  real b;
  real<lower=0> sigma_a;
}
transformed parameters {
  vector[N_teams] a = b * prior_score + sigma_a * alpha;
}
model {
  alpha ~ normal(0, 1);
  b ~ normal(0, 1);
  sigma_a ~ normal(0, 1);
  for (n in 1:N_games) {
    // log probability using the Poisson difference lpmf
    // e.g. Karlis and Ntzoufras (2003). Analysis of sports data by using bivariate Poisson models
    real l1 = exp(a[team_1[n]] - a[team_2[n]]);
    real l2 = exp(a[team_2[n]] - a[team_1[n]]);
    target += -(l1 + l2) + (log(l1) - log(l2))*(dif[n]/2.0) +
      log(modified_bessel_first_kind(to_int(dif[n]), 2.0*sqrt(l1*l2)));
  }
}
generated quantities {
  vector[N_games] y_rep;
  vector[N_games] log_lik;
  for (n in 1:N_games) {
    int score_1_rep = poisson_log_rng(a[team_1[n]] - a[team_2[n]]);
    int score_2_rep = poisson_log_rng(a[team_2[n]] - a[team_1[n]]);
    y_rep[n] = score_1_rep - score_2_rep;
    // log probability using the Poisson difference lpmf
    // e.g. Karlis and Ntzoufras (2003). Analysis of sports data by using bivariate Poisson models
    real l1 = exp(a[team_1[n]] - a[team_2[n]]);
    real l2 = exp(a[team_2[n]] - a[team_1[n]]);
    log_lik[n] = -(l1 + l2) + (log(l1) - log(l2))*(dif[n]/2.0) +
      log(modified_bessel_first_kind(to_int(dif[n]), 2.0*sqrt(l1*l2)));
  }
}

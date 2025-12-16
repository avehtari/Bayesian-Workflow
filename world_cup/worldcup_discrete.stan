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
  real<lower=0> sigma_z;
}
transformed parameters {
  vector[N_teams] a = b * prior_score + sigma_a * alpha;
}
model {
  alpha ~ normal(0, 1);
  b ~ normal(0, 1);
  sigma_a ~ normal(0, 1);
  sigma_z ~ normal(0, 1);
  for (n in 1:N_games) {
    target += log_diff_exp(normal_lcdf(dif[n]+0.5 | a[team_1[n]] - a[team_2[n]], sigma_z),
                           normal_lcdf(dif[n]-0.5 | a[team_1[n]] - a[team_2[n]], sigma_z));
  }
}
generated quantities {
  vector[N_games] y_rep;
  vector[N_games] z_rep;
  vector[N_games] log_lik;
  for (n in 1:N_games) {
    z_rep[n] = normal_rng(a[team_1[n]] - a[team_2[n]], sigma_z);
    y_rep[n] = round(z_rep[n]);
    log_lik[n] = log_diff_exp(normal_lcdf(dif[n]+0.5 | a[team_1[n]] - a[team_2[n]], sigma_z), normal_lcdf(dif[n]-0.5 | a[team_1[n]] - a[team_2[n]], sigma_z));
  }
}

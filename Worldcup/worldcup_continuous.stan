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
  vector[N_teams] alpha;     // latent team abilities
  real b;                    // prior weight
  real<lower=0> sigma_a;     // scale of ability prior distribution
  real<lower=0> sigma_y;     // scale of score difference distribution
}
transformed parameters {
  // team abilities as weighted sum of prior score and latent team ability
  vector[N_teams] a = b * prior_score + sigma_a * alpha;
}
model {
  alpha ~ normal(0, 1);
  b ~ normal(0, 1);
  sigma_a ~ normal(0, 1);
  sigma_y ~ normal(0, 1);
  dif ~ normal(a[team_1] - a[team_2], sigma_y);
}
generated quantities {
  vector[N_games] y_rep;
  vector[N_games] log_lik;
  for (n in 1:N_games) {
    y_rep[n] = normal_rng(a[team_1[n]] - a[team_2[n]], sigma_y);
    log_lik[n] = log_diff_exp(normal_lcdf(dif[n]+0.5 | a[team_1[n]] - a[team_2[n]], sigma_y), normal_lcdf(dif[n]-0.5 | a[team_1[n]] - a[team_2[n]], sigma_y));
  }
}

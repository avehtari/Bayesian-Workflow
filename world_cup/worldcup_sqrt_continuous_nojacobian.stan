/* The error in the transformed data block has been fixed */
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
  sqrt_dif ~ normal(a[team_1] - a[team_2], sigma_y);
}
generated quantities {
  vector[N_games] y_rep;
  vector[N_games] y_rep_original_scale;
  vector[N_games] log_lik;
  real lprior;
  for (n in 1:N_games) {
    y_rep[n] = normal_rng(a[team_1[n]] - a[team_2[n]], sigma_y);
    log_lik[n] = normal_lpdf(sqrt_dif[n] | a[team_1[n]] - a[team_2[n]], sigma_y);
  }
  y_rep_original_scale = y_rep .* abs(y_rep);
  lprior = normal_lpdf(b | 0, 1) +
           normal_lpdf(sigma_a | 0, 1) +
           normal_lpdf(sigma_y | 0, 1);
  real round15 = round(1.5);
  real round25 = round(2.5);
}

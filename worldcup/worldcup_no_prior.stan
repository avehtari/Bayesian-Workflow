data {
  int N_teams;
  int N_games;
  vector[N_teams] prior_score;
  array[N_games] int team_1;
  array[N_games] int team_2;
  vector[N_games] score_1;
  vector[N_games] score_2;
  real df;
  real b;  // To remove the prior score in the model, just set b=0 when running this program.
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
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[N_teams] a = b * prior_score + sigma_a * alpha;
}
model {
  alpha ~ normal(0, 1);
  sigma_a ~ normal(0, 1);
  sigma_y ~ normal(0, 1);
  sqrt_dif ~ student_t(df, a[team_1] - a[team_2], sigma_y);
}
generated quantities {
  vector[N_games] log_lik;
  for (n in 1:N_games) {
    log_lik[n] = student_t_lpdf(sqrt_dif[n] | df, a[team_1[n]] - a[team_2[n]], sigma_y);
  }
}

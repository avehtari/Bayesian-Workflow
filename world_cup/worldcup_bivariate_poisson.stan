functions {
 real bivariate_poisson_log_lpmf_(int score1, int score2,
                                  real log_lambda_1, real log_lambda_2, real log_lambda_3) {
    // bivariate Poisson with lambda_3 presenting the correlation
    // e.g. Karlis and Ntzoufras (2003). Analysis of sports data by using bivariate Poisson models
    int m = min(score1, score2);
    vector[m + 1] log_terms;
    for (z in 0 : m) {
      int r1 = score1 - z;
      int r2 = score2 - z;
      log_terms[z + 1] = r1 * log_lambda_1 - lgamma(r1 + 1)
                         + r2 * log_lambda_2 - lgamma(r2 + 1)
                         + z * log_lambda_3 - lgamma(z + 1);
    }
    return -exp(log_lambda_1) - exp(log_lambda_2) - exp(log_lambda_3) + log_sum_exp(log_terms);
}
  real bivariate_poisson_log_lpmf(array[,] int score,
                                  vector log_lambda_1, vector log_lambda_2, real log_lambda_3) {
    // bivariate Poisson lpmf for array[2,N] scores
    array[2] int score_dims = dims(score);
    real total = 0;
    for (i in 1:score_dims[2]) {
      total += bivariate_poisson_log_lpmf_(score[1,i], score[2,i],
                                           log_lambda_1[i], log_lambda_2[i], log_lambda_3);
    }
    return total;
  }
  real bivariate_poisson_log_lpmf(array[] int score,
                                  real log_lambda_1, real log_lambda_2, real log_lambda_3) {
    // bivariate Poisson lpmf for array[2] scores
    return bivariate_poisson_log_lpmf_(score[1], score[2],
                                       log_lambda_1, log_lambda_2, log_lambda_3);
  }
}
data {
  int N_teams;
  int N_games;
  vector[N_teams] prior_score;
  array[N_games] int team_1;
  array[N_games] int team_2;
  array[N_games] int score_1;
  array[N_games] int score_2;
}
parameters {
  vector[N_teams] z_o;
  vector[N_teams] z_d;
  real a;
  real b_o;
  real b_d;
  real c;
  real<lower=0> sigma_o;
  real<lower=0> sigma_d;
}
transformed parameters {
  vector[N_teams] o = b_o * prior_score + sigma_o * z_o;
  vector[N_teams] d = b_d * prior_score + sigma_d * z_d;
}
model {
  z_o ~ normal(0, 1);
  z_d ~ normal(0, 1);
  b_o ~ normal(0, 1);
  b_d ~ normal(0, 1);
  a ~ normal(0, 1);
  // prior centered on 1991–1992 Italian serie A data estimate as reported by
  // Karlis and Ntzoufras (2003) and very wide prior around.
  c ~ normal(-1.5, 20);
  sigma_o ~ normal(0, 1);
  sigma_d ~ normal(0, 1);
  {score_1, score_2} ~ bivariate_poisson_log(a + o[team_1] + d[team_2],
                                             a + o[team_2] + d[team_1],
					     c);
}
generated quantities {
  array[N_games] int y_rep;
  array[N_games] int score_1_rep;
  array[N_games] int score_2_rep;
  vector[N_games] log_lik;
  vector[N_games] log_lik2;
  for (n in 1:N_games) {
    int score_c = poisson_log_rng(c);
    score_1_rep[n] = poisson_log_rng(a + o[team_1[n]] + d[team_2[n]]) + score_c;
    score_2_rep[n] = poisson_log_rng(a + o[team_2[n]] + d[team_1[n]]) + score_c;
    y_rep[n] = score_1_rep[n] - score_2_rep[n];
    // log probability using the bivariate Poisson lpmf
    log_lik2[n] = bivariate_poisson_log_lpmf({score_1[n], score_2[n]} |
                                        a + o[team_1[n]] + d[team_2[n]],
                                        a + o[team_2[n]] + d[team_1[n]],
					c);
    // log probability using the Poisson difference lpmf
    // note that lambda3 cancels out in the difference
    int dif = score_1[n] - score_2[n];
    real l1 = exp(o[team_1[n]] + d[team_2[n]]);
    real l2 = exp(o[team_2[n]] + d[team_1[n]]);
    log_lik[n] = -(l1 + l2) + (log(l1) - log(l2))*(dif/2.0) +
                 log(modified_bessel_first_kind(dif, 2.0*sqrt(l1*l2)));
  }
}

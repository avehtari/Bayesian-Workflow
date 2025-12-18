#' ---
#' title: "Bayesian workflow book - World Cup continuous vs discrete"
#' author: "Andrew Gelman and Aki Vehtari"
#' date: 2021-01-12
#' date-modified: today
#' date-format: iso
#' format:
#'   html:
#'     toc: true
#'     toc-location: left
#'     toc-depth: 2
#'     number-sections: true
#'     smooth-scroll: true
#'     theme: readable
#'     code-copy: true
#'     code-download: true
#'     code-tools: true
#'     embed-resources: true
#'     anchor-sections: true
#'     html-math-method: katex
#' bibliography: ../casestudies.bib
#' ---
#' 
#' # World Cup 2014 team performance analysis
#'
#' We use an item-response model such as described in the previous
#' chapter to fit a model to estimate the abilities of the teams in
#' the 2014 football World Cup.  We use score differentials as data
#' (ignoring the shoot-outs).
#'
#' The forthcoming Bayesian workflow book has more information about
#' the modeling task, but this version of the code was made available
#' early, as it includes illustration of predictive performance
#' comparison of continuous and discrete models.
#' 
#' Game scores could be modelled with a discrete bivariate Poisson
#' model, and score differences with Poisson difference models (Karlis
#' and Ntzoufras, 2003), which are commonly used in analysis of sports
#' data. But we'll start with continuous models and discretized
#' continuous models to illustrate how those can be compared, and
#' eventually build also models using naturally discrete data models.
#'
#' -------------
#' 

#+ setup, include=FALSE
knitr::opts_chunk$set(message=FALSE, error=FALSE, warning=FALSE, comment=NA, cache=FALSE)
# switch this to TRUE to save figures in separate files
savefigs <- FALSE

#' **Load packages**
#| cache: FALSE
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
options(mc.cores = 4)
library(posterior)
options(
  tibble.print_max = 35,
  pillar.neg = FALSE,
  pillar.subtle = FALSE,
  pillar.sigfig = 2
)
library(bayesplot)
library(ggplot2)
theme_set(bayesplot::theme_default(base_family = "sans", base_size=14))
library(arm)
library(loo)
library(dplyr)
library(readr)

#' # Data
#'
#' Data include 64 game results from World Cup 2014 and Soccer Power
#' Index that was available on the internet a month before the
#' tournament (Silver, 2014).  We took the rankings, with Brazil at
#' the top (getting a score of 32) and Australia at the bottom (with a
#' score of 1), and then for simplicity in interpretation of the
#' parameters we rescaled these to have mean 0 and standard deviation
#' 1/2, to get ``prior scores'' that ranged from $-0.83$ to $+0.83$.
powerindex <- read_csv(root("world_cup", "data","soccerpowerindex.csv")) %>%
  mutate(prior_score = as.vector(scale(rev(index))/2))
teamnames <- powerindex$team
worldcup2014 <- read_csv(root("world_cup", "data", "worldcup2014.csv")) %>%
  mutate(team_1 = match(team1, teamnames),
         team_2 = match(team2, teamnames))
N_games <- nrow(worldcup2014)
gamenames <- with(worldcup2014, rev(paste(teamnames[team_1], "vs.", teamnames[team_2])))

#' ## Data for Stan
stan_data <- with(
  worldcup2014,
  list(
    N_teams = nrow(powerindex),
    N_games = N_games,
    team_1 = team_1,
    score_1 = score1,
    team_2 = team_2,
    score_2 = score2,
    prior_score = powerindex$prior_score,
    df = 7
  )
)

#' # First model
#'
#' The first model uses Student's $t$ for square root of score
#' difference. It's not particularly well justified model, but that is
#' what Andrew first thought and it turned out to be useful for
#' illustration.
#'
#' The model structure is as follows: if game $i$ has teams $j_1$ and
#' team $j_2$ playing, and they score $z_1$ and $z_2$ goals,
#' respectively, then the data point for this game is $y_i =
#' \mbox{sign}(z_1-z_2)*\sqrt{|z_1-z_2|}$, and the data model is: $y_i
#' \sim \mbox{normal}(a_{j_1[i]}-a_{j_2[i]}, \sigma_y)$, where
#' $a_{j_1}$ and $a_{j_2}$ are the ability parameters (to use
#' psychometrics jargon) for the two teams and $\sigma_y$ is a scale
#' parameter estimated from the data.
#' 
#' $$
#' y_i \sim \mbox{t}_{\nu}(a_{j_1[i]}-a_{j_2[i]}, \sigma_y),
#' $$
#' setting the degrees of freedom to $\nu=7$.
#'
#' 
#' 
#' Fit the model and show the results
#| label: fit_1
#| results: hide
#| cache: false
model_1 <- cmdstan_model(stan_file = root("world_cup", "worldcup_first_try.stan"))
fit_1 <- model_1$sample(data = stan_data, refresh = 0)

#| label: fit_1_summary
fit_1$summary(c("a", "b", "sigma_a", "sigma_y"))

#| label: fig-worldcup-mcmc-intervals-fit_1
#| fig-height: 7
#| fig-width: 8
fit_1$draws("a") |>
  mcmc_intervals(prob = 0) +
  scale_y_discrete(labels = rev(teamnames), limits = rev) +
  labs(x = "Team quality estimate with 90% intervals")
#if (savefigs) ggsave(root("world_cup/figs","worldcup_1.pdf", height=7, width=8))

#' ## Check fit of the first model
#| label: fit_1_rep
#| results: hide
#| cache: false
model_1_rep <- cmdstan_model(stan_file = root("world_cup", "worldcup_with_replication.stan"))
fit_1_rep <- model_1_rep$sample(data = stan_data, refresh = 0)

#| label: fig-worldcup-ppc-intervals-fit_1
#| fig-height: 10
#| fig-width: 8
ppc_intervals(
  y = stan_data$score_1 - stan_data$score_2,
  yrep = fit_1_rep$draws("y_rep_original_scale", format = "draws_matrix"),
  fatten = 0,
  prob = 1e-12
) +
  scale_x_continuous(labels = gamenames, breaks=1:64)+
  labs(y="Game score differentials\ncompared to 90% predictive interval from model", x="") +
  coord_flip() 
#if (savefigs) ggsave(root("world_cup/figs","worldcup_3.pdf", height=10, width=8))

#' # Second model without sqrt transformation
#| label: fit_2
#| results: hide
#| cache: false
model_2 <- cmdstan_model(stan_file = root("world_cup", "worldcup_no_sqrt.stan"))
fit_2 <- model_2$sample(data = stan_data, refresh = 0)

#| label: fit_2_summary
fit_2$summary(c("a", "b", "sigma_a", "sigma_y"))

#| label: fig-worldcup-mcmc-intervals-fit_2
#| fig-height: 7
#| fig-width: 8
fit_2$draws("a") |>
  mcmc_intervals(prob = 0) +
  scale_y_discrete(labels = rev(teamnames), limits = rev) +
  labs(x = "Team quality estimate with 90% intervals\n(model with no square root)")
#if (savefigs) ggsave(root("world_cup/figs","worldcup_4.pdf", height=7, width=8))

#' ## Check fit of the second model to data
#| label: fig-worldcup-ppc-intervals-fit_2
#| fig-height: 10
#| fig-width: 8
ppc_intervals(
  y = stan_data$score_1 - stan_data$score_2,
  yrep = fit_2$draws("y_rep", format = "draws_matrix"),
  fatten = 0,
  prob = 1e-12
) +
  scale_x_continuous(labels = gamenames, breaks = 1:64)+
  labs(y = "Game score differentials ncompared to 90% predictive interval\n(model with no square root)", x = "") +
  coord_flip() 
#if (savefigs) ggsave(root("world_cup/figs","worldcup_5.pdf", height=10, width=8))

#' # Fix the first model
#| label: fit_3
#| results: hide
#| cache: false
model_3 <- cmdstan_model(stan_file = root("world_cup", "worldcup_fixed.stan"))
fit_3 <- model_3$sample(data = stan_data, refresh = 0)

#| label: fit_3_summary
fit_3$summary(c("a", "b", "sigma_a", "sigma_y"))

#| label: fig-worldcup-mcmc-intervals-fit_3
#| fig-height: 10
#| fig-width: 8
fit_3$draws("a") |>
  mcmc_intervals(prob = 0) +
  scale_y_discrete(labels = rev(teamnames), limits = rev) +
  labs(x = "Team quality estimate with 90% intervals\n(corrected model)")
#if (savefigs) ggsave(root("world_cup/figs","worldcup_6.pdf", height=7, width=8))


#' ## Check the fit of the fixed first model
#| label: fig-worldcup-ppc-intervals-fit_3
#| fig-height: 10
#| fig-width: 8
ppc_intervals(
  y = stan_data$score_1 - stan_data$score_2,
  yrep = fit_3$draws("y_rep_original_scale", format = "draws_matrix"),
  fatten = 0,
  prob = 1e-12
) +
  scale_x_continuous(labels = gamenames, breaks = 1:64) +
  labs(y = "Game score differentials ncompared to 90% predictive interval\n(corrected model with square root)", x = "") +
  coord_flip() 
#if (savefigs) ggsave(root("world_cup/figs","worldcup_7.pdf", height=10, width=8))

#' ## Fit the same model without the powerindex prior
#' set b=0 in the data
#| label: fit_3_no_prior
#| results: hide
#| cache: false
model_3_no_prior <- cmdstan_model(stan_file = root("world_cup", "worldcup_no_prior.stan"))
fit_3_no_prior <- model_3_no_prior$sample(data = c(stan_data, b = 0), refresh = 0)

#| label: fig-worldcup-mcmc-intervals-fit_3_no_prior
#| fig-height: 7
#| fig-width: 8
fit_3_no_prior$draws("a") |>
  mcmc_intervals(prob = 0) +
  scale_y_discrete(labels = rev(teamnames), limits = rev) +
  labs(x = "Team quality estimate with 90% intervals\nModel without prior rankings")
#if (savefigs) ggsave(root("world_cup/figs","worldcup_2.pdf", height=7, width=8))


#' # Discrete models and LOO-CV comparison

#' ## Discrete model with explicit latent z sampled
#| label: fit_discr_z
#| results: hide
#| cache: false
model_discr_z <- cmdstan_model(stan_file = root("world_cup", "worldcup_discrete_z.stan"))
fit_discr_z <- model_discr_z$sample(data = stan_data, refresh = 0)
loo_discr_z <- fit_discr_z$loo()

#| label: fit_discr_z_summary
fit_discr_z$summary(c("b", "sigma_a", "sigma_z"))

#' ## Discrete model with latent z integrated out
#| label: fit_discr
#| results: hide
#| cache: false
model_discr <- cmdstan_model(stan_file = root("world_cup", "worldcup_discrete.stan"))
fit_discr <- model_discr$sample(data = stan_data, refresh = 0)
loo_discr <- fit_discr$loo()

#| label: fit_discr_summary
fit_discr$summary(c("a[1]", "a[32]", "b", "sigma_a", "sigma_z"))

#' ## Discrete model with no power score
#| label: fit_discr_nopwer
#| results: hide
#| cache: false
model_discr_nopower <- cmdstan_model(stan_file = root("world_cup", "worldcup_discrete_nopower.stan"))
fit_discr_nopower <- model_discr_nopower$sample(data = stan_data, refresh = 0)
loo_discr_nopower <- fit_discr_nopower$loo()

#| label: fit_discr_nopower_summary
fit_discr_nopower$summary(c("a[1]", "a[32]", "sigma_a", "sigma_z"))


#' ## Discrete model with power score only
#| label: fit_discr_poweronly
#| results: hide
#| cache: false
model_discr_poweronly <- cmdstan_model(stan_file = root("world_cup", "worldcup_discrete_poweronly.stan"))
fit_discr_poweronly <- model_discr_poweronly$sample(data = stan_data, refresh = 0)
loo_discr_poweronly <- fit_discr_poweronly$loo()

#| label: fit_discr_poweronly_summary
fit_discr_poweronly$summary(c("b0", "b", "sigma_z"))

#' ## Discrete model with no power score and pooled effect
#| label: fit_discr_pooled
#| results: hide
#| cache: false
model_discr_pool <- cmdstan_model(stan_file = root("world_cup", "worldcup_discrete_pooled.stan"))
fit_discr_pool <- model_discr_pool$sample(data = stan_data, refresh = 0)
loo_discr_pool <- fit_discr_pool$loo()

#| label: fit_discr_pool_summary
fit_discr_pool$summary(c("mu_z","sigma_z"))

#' ## Model comparison
#' 
#' Now that we have implemented discrete model, we compare various
#' models with or without different components. Hierarchical model
#' with power score is the best, but not much better than power score
#' only model or hierarchical model without power score. Clearly the
#' power score and the score differences are providing similar
#' information. Looking at the posterior we do see that using them
#' both, does decrease posterior uncertainty, but as the match
#' outcomes still have significant randomness the difference in
#' predictive performance is small.
loo_compare(
  list(
    "Hier. w power score" = loo_discr,
    "Power score only" = loo_discr_poweronly,
    "Hier. w/o power score" = loo_discr_nopower,
    "Pooled w/o power score" = loo_discr_pool
  )
)

#' # Discretizing continuous models
#' 
#' ## Continuous model with log predictive probability using midpoint rule
#' 
#| label: fit_cont_midp
#| results: hide
#| cache: false
model_cont_midp <- cmdstan_model(stan_file = root("world_cup", "worldcup_continuous_midpoint_ll.stan"))
fit_cont_midp <- model_cont_midp$sample(data = stan_data, refresh = 0)
loo_cont_midp <- fit_cont_midp$loo()

#| label: fit_cont_midp_summary
fit_cont_midp$summary(c("a[1]", "a[32]", "b", "sigma_a", "sigma_y"))

#' ## Continuous model with log predictive probability using exact integration
#' 
#| label: fit_cont
#| results: hide
#| cache: false
model_cont <- cmdstan_model(stan_file = root("world_cup", "worldcup_continuous.stan"))
fit_cont <- model_cont$sample(data = stan_data, refresh = 0)
loo_cont <- fit_cont$loo()

#| label: fit_cont_summary
fit_cont$summary(c("a[1]", "a[32]", "b", "sigma_a", "sigma_y"))

#' ## Model comparison
#' 
#' The elpd_diff's between models are small enough to be caused
#' by Monte Carlo variation
loo_compare(
  list(
    "Discrete model" = loo_discr,
    "Continuous + midpoint log_lik" = loo_cont_midp,
    "Continuous model" = loo_cont
  )
)

#' ## More discretized continuous models
#'
#' Continuous sqrt model with log predictive probability using midpoint rule and no Jacobian
#| label: fit_sqrt_cont_noj
#| results: hide
#| cache: false
model_sqrt_cont_noj <- cmdstan_model(stan_file = root("world_cup", "worldcup_sqrt_continuous_nojacobian.stan"))
fit_sqrt_cont_noj <- model_sqrt_cont_noj$sample(data = stan_data, refresh = 0)
loo_sqrt_cont_noj <- fit_sqrt_cont_noj$loo()

#| label: fit_sqrt_cont_noj_summary
fit_sqrt_cont_noj$summary(c("a[1]", "a[32]", "b", "sigma_a", "sigma_y"))

#' Continuous sqrt model with log predictive probability using
#' midpoint rule and Jacobian is not possible as Jacobian is infinite
#' for midpoint 0

#' Continuous sqrt model with log predictive probability using
#' Jacobian and quadrature integration
#| label: fit_sqrt_cont
#| results: hide
#| cache: false
model_sqrt_cont <- cmdstan_model(stan_file = root("world_cup", "worldcup_sqrt_continuous.stan"))
fit_sqrt_cont <- model_sqrt_cont$sample(data = stan_data, refresh = 0)
loo_sqrt_cont <- fit_sqrt_cont$loo()

#| label: fit_sqrt_cont_summary
fit_sqrt_cont$summary(c("a[1]", "a[32]", "b", "sigma_a", "sigma_y"))

#' Discrete sqrt model with log predictive probability using Jacobian
#' and quadrature integration. The sampling is slow as the quadrature
#' integration is done at each HMC/NUTS leapfrog step.
#| label: fit_sqrt_discr
#| results: hide
#| cache: false
model_sqrt_discr <- cmdstan_model(stan_file = root("world_cup", "worldcup_sqrt_discrete.stan"))
fit_sqrt_discr <- model_sqrt_discr$sample(data = stan_data, refresh = 0)
loo_sqrt_discr <- fit_sqrt_discr$loo()

#| label: fit_sqrt_discr_summary
fit_sqrt_discr$summary(c("a[1]", "a[32]", "b", "sigma_a", "sigma_y"))

#' ## Model comparison
#' 
#' Without taking Jacobian into account, the continuous model for
#' square root of score difference looks like the best
loo_compare(
  list(
    "Discrete" = loo_discr,
    "Discrete sqrt" = loo_sqrt_discr,
    "Cont-sqrt + midp, -Jacobian" = loo_sqrt_cont_noj
  )
)

#' Taking Jacobian into account, the square root models are worse, and
#' the difference between continuous and discrete square root model is
#' small enough to be explained by Monte Carlo variation.
loo_compare(
  list(
    "Discrete" = loo_discr,
    "Discrete sqrt" = loo_sqrt_discr,
    "Continuous sqrt +Jacobian" = loo_sqrt_cont
  )
)

#' # Posterior predictive checking
#' 
#' Posterior predictive checking for the discrete model looks fine
ppc_pit_ecdf(
  y = stan_data$score_1 - stan_data$score_2,
  yrep = fit_discr$draws(format = "matrix", variables = "y_rep")
)

#' Posterior predictive checking for the continuous model indicates
#' slight miscalibration with too many low PIT values (left tail of
#' the predictive disttribution is shorter than expected)
ppc_pit_ecdf(
  y = stan_data$score_1 - stan_data$score_2,
  yrep = fit_sqrt_cont$draws(format = "matrix", variables = "y_rep")
)

#' # Bivariate Poisson and Poisson difference models
#' 
#' Bivariate Poisson model is commonly used for football scores
#' [@Karlis-Ntzoufras:2003]. If we care only about the score
#' difference we can also use Poisson difference model
#' [@Karlis-Ntzoufras:2003].
#'
#| label: fit_bipois
#| results: hide
#| cache: false
model_bipois <- cmdstan_model(stan_file = root("world_cup", "worldcup_bivariate_poisson.stan"))
fit_bipois <- model_bipois$sample(data = stan_data, refresh = 0, adapt_delta = 0.95)

#| label: fit_bipois_summary
fit_bipois$summary(c("a","o[1]","o[32]","d[1]","d[32]","b_o","b_d","sigma_o","sigma_d"))
(loo_bipois <- fit_bipois$loo())

#| label: fit_poisdif
#| results: hide
#| cache: false
model_poisdif <- cmdstan_model(stan_file = root("world_cup", "worldcup_poisson_difference.stan"))
fit_poisdif <- model_poisdif$sample(data = stan_data, refresh = 0)

#| label: fit_poisdif_summary
fit_poisdif$summary(c("a[1]", "a[32]", "b", "sigma_a"))
(loo_poisdif <- fit_poisdif$loo())

#' The set in this case study is small, and we don't see practical
#' difference in the predictive performance.
loo_compare(
  list(
    "Discrete" = loo_discr,
    "Bivariate Poisson" = loo_bipois,
    "Poisson difference" = loo_poisdif
  )
)

#' In this case study, we used many other models for illustration, but
#' for real football score modeling, it is a good idea to start with
#' the bivariate Poisson model.

#' 
#' <br />
#' 
#' # References {.unnumbered}
#' 
#' <div id="refs"></div>
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2021-2025, Andrew Gelman, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2021-2025, Andrew Gelman, Aki Vehtari, licensed under CC-BY-NC 4.0.
#' 
#' # Original Computing Environment {.unnumbered}
#' 
sessionInfo()

#' 
#' <br />
#' 

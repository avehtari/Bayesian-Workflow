#' ---
#' title: "Illustration of simple problematic posteriors"
#' author: "Aki Vehtari"
#' date: 2021-06-10
#' date-modified: today
#' date-format: iso
#' format:
#'   html:
#'     number-sections: true
#'     code-copy: true
#'     code-download: true
#'     code-tools: true
#' bibliography: ../casestudies.bib
#' ---
#'
#' This notebook includes the code for Bayesian Workflow book Section
#' 12.3 *Failure modes and steps forward*.
#'
#' # Introduction
#'
#' This case study demonstrates using simple examples the most common
#' failure modes in Markov chain Monte Carlo based Bayesian inference,
#' how to recognize these using the diagnostics, and how to fix the
#' problems.
#'
#+ setup, include=FALSE
knitr::opts_chunk$set(
  cache = FALSE,
  message = FALSE,
  error = FALSE,
  warning = FALSE,
  comment = NA,
  out.width = "90%"
)

#' **Load packages**
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr) 
library(posterior)
options(pillar.neg = FALSE,
        pillar.subtle = FALSE,
        pillar.sigfig = 2)
options(width = 90)
library(tidyr) 
library(dplyr) 
library(ggplot2)
library(bayesplot)
library(RColorBrewer)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))
library(latex2exp)
library(patchwork)
set1 <- RColorBrewer::brewer.pal(7, "Set1")
SEED <- 48927 # set random seed for reproducibility

print_stan_file <- function(file) {
  code <- readLines(file)
  if (isTRUE(getOption("knitr.in.progress")) &
        identical(knitr::opts_current$get("results"), "asis")) {
    # In render: emit as-is so Pandoc/Quarto does syntax highlighting
    block <- paste0("```stan", "\n", paste(code, collapse = "\n"), "\n", "```")
    knitr::asis_output(block)
  } else {
    writeLines(code)
  }
}

#'
#' # Improper posterior
#'
#' An unbounded likelihood without a proper prior can lead to an
#' improper posterior. We recommend to always use proper priors
#' (integral over a proper distribution is finite) to guarantee proper
#' posteriors.
#'
#' A commonly used model that can have unbounded likelihood is
#' logistic regression with complete separation in data.
#'
#' ## Data
#'
#' Univariate continous predictor $x$, binary target $y$, and the two
#' classes are completely separable, which leads to unbounded
#' likelihood.
set.seed(SEED + 4)
M <- 1
N <- 10
x <- matrix(sort(rnorm(N)), ncol = M)
y <- rep(c(0, 1), each = N / 2)
data_logit <- list(M = M, N = N, x = x, y = y)
#| label: fig-separable_data
#| fig-height: 4
#| fig-width: 6
data.frame(data_logit) |>
  ggplot(aes(x, y)) +
  geom_point(size = 3, shape = 1, alpha = 0.6) +
  scale_y_continuous(breaks = c(0, 1))

#'
#' ## Model
#'
#' We use the following Stan logistic regression model, where we have
#' ``forgot'' to include prior for the coefficient `beta`.
code_logit <- root("problems", "logit_glm.stan")
#| output: asis
print_stan_file(code_logit)

#' Sample
#| label: fit_logit
#| results: hide
mod_logit <- cmdstan_model(stan_file = code_logit)
fit_logit <- mod_logit$sample(data = data_logit, seed = SEED, refresh = 0)

#'
#' ## Convergence diagnostics
#'
#' When running Stan, we get warnings. We can also
#' explicitly check the inference diagnostics:
fit_logit$diagnostic_summary()

#' We can also check $\widehat{R}$ end effective sample size (ESS) diagnostics
#' [@Vehtari-Gelman-Simpson-etal:2021]
draws <- as_draws_rvars(fit_logit$draws())
summarize_draws(draws)

#' We see that $\widehat{R}$ for both \texttt{alpha} and \texttt{beta}
#' are about 3 and Bulk-ESS is about 4, which indicate that the chains
#' are not mixing at all.
#'
#' The above diagnostics refer to a documentation
#' ([https://mc-stan.org/misc/warnings](https://mc-stan.org/misc/warnings))
#' that mentions possibility to adjust the sampling algorithm options
#' (e.g., increasing `adapt_delta` and `max_treedepth`), but it is
#' better first to investigate the posterior.
#'
#' The following Figure shows the posterior draws as marginal
#' histograms and joint scatterplots. The range of the values is huge,
#' which is typical for improper posterior, but the values of `alpha`
#' and `beta` in any practical application are likely to have much
#' smaller magnitude. In this case, increasing `adapt_delta` and
#' `max_treedepth` would not have solved the problem, and would have
#' just caused waste of modeler and computation time.
#'
#| label: fig-separable_pairs
#| fig-height: 4
#| fig-width: 6
mcmc_pairs(as_draws_array(draws), pars = c("alpha", "beta"),
           off_diag_args = list(alpha = 0.2))

#'
#' ## Stan compiler pedantic check
#'
#' The above diagnostics are applicable with any probabilistic
#' programming framework.  Stan compiler can also recognize some
#' common problems. By default the pedantic mode is not enabled, but
#' we can use option `pedantic = TRUE` at compilation time, or after
#' compilation with the `check_syntax` method.
#| results: hide
mod_logit$check_syntax(pedantic = TRUE)

#' The pedantic check correctly warns that `alpha` and `beta` don't
#' have priors.
#'
#' ## A fixed model with proper priors
#'
#' We add proper weak priors and rerun inference.
code_logit2 <- root("problems", "logit_glm2.stan")
#| output: asis
print_stan_file(code_logit2)
#' Sample
#| label: fit_logit2
#| results: hide
mod_logit2 <- cmdstan_model(stan_file = code_logit2)
fit_logit2 <- mod_logit2$sample(data = data_logit, seed = SEED, refresh = 0)

#'
#' ## Convergence diagnostics
#'
#' There were no convergence warnings. We can also
#' explicitly check the inference diagnostics:
fit_logit2$diagnostic_summary()

#' We check $\widehat{R}$ end ESS values, which in this case all look good.
draws <- as_draws_rvars(fit_logit2$draws())
summarize_draws(draws)

#' The following figure shows the more reasonable marginal histograms
#' and joint scatterplots of the posterior sample.
#| label: fig-separable_prior_pairs
#| fig-height: 4
#| fig-width: 6
mcmc_pairs(as_draws_array(draws), pars = c("alpha", "beta"),
           off_diag_args = list(alpha = 0.2))

#'
#' # A model with unused parameter
#'
#' When writing and editing models, a common mistake is to declare a
#' parameter, but not use it in the model. If the parameter is not
#' used at all, it doesn't have proper prior and the likelihood
#' doesn't provide information about that parameter, and thus the
#' posterior along that parameter is improper. We use the previous
#' logistic regression model with proper priors on `alpha` and
#' `beta`, but include extra parameter declaration `real
#' gamma`.
#'
#' ## Model
code_logit3 <- root("problems", "logit_glm3.stan")
#| output: asis
print_stan_file(code_logit3)
#' Sample
#| label: fit_logit3
#| results: hide
mod_logit3 <- cmdstan_model(stan_file = code_logit3)
fit_logit3 <- mod_logit3$sample(data = data_logit, seed = SEED, refresh = 0)

#'
#' ## Convergence diagnostics
#'
#' There is sampler warning. We can also explicitly call inference
#' diagnostics:
fit_logit3$diagnostic_summary()

#' Instead of increasing `max_treedepth`, we check the other convergence diagnostics. 
draws <- as_draws_rvars(fit_logit3$draws())
summarize_draws(draws)

#' $\widehat{R}$, Bulk-ESS, and Tail-ESS look good for `alpha` and
#' `beta, but really bad for `gamma`, clearly pointing where to look
#' for problems in the model code. The histogram of `gamma` posterior
#' draws show huge magnitude of values (values larger than $10^{20}$)
#' indicating improper posterior.
#| label: fig-unusedparam_pairs
#| fig-height: 4
#| fig-width: 6
mcmc_pairs(as_draws_array(draws), pars = c("alpha", "beta", "gamma"),
           off_diag_args = list(alpha = 0.2))

#' Non-mixing is well diagnosed by $\widehat{R}$ and ESS, but the
#' following Figure shows one of the rare cases where trace plots are
#' useful to illustrate the type of non-mixing in case of improper
#' uniform posterior for one the parameters.
#| label: fig-unusedparam_trace
#| fig-height: 4
#| fig-width: 7
mcmc_trace(as_draws_array(draws), pars = c("gamma"))

#' ## Stan compiler pedantic check
#'
#' Stan compiler pedantic check also recognizes that parameter `gamma` was
#' declared but was not used in the density calculation.
#| results: hide
mod_logit3$check_syntax(pedantic = TRUE)

#'
#' # A posterior with two parameters competing
#'
#' Sometimes the models have two or more parameters that have similar
#' or exactly the same role. We illustrate this by adding an extra
#' column to the previous data matrix. Sometimes the data matrix is
#' augmented with a column of 1’s to present the intercept effect. In
#' this case that is redundant as our model has the explicit intercept
#' term `alpha`, and this redundancy will lead to problems.
#'
#' ## Data
M <- 2
N <- 1000
x <- matrix(c(rep(1, N), sort(rnorm(N))), ncol = M)
y <- ((x[, 1] + rnorm(N) / 2) > 0) + 0
data_logit4 <- list(M = M, N = N, x = x, y = y)

#'
#' ## Model
#'
#' We use the previous logistic regression model with proper priors
#' (and no extra `gamma`).
#'
code_logit2 <- root("problems", "logit_glm2.stan")
#| output: asis
print_stan_file(code_logit2)
#' Sample
#| label: fit_logit4
#| results: hide
mod_logit4 <- cmdstan_model(stan_file = code_logit2)
fit_logit4 <- mod_logit4$sample(data = data_logit4, seed = SEED, refresh = 0)

#' The Stan sampling time per chain with the original data matrix was
#' less than 0.1s per chain. Now the Stan sampling time per chain is
#' several seconds, which is suspicious. There are no automatic
#' convergence diagnostic warnings and checking other diagnostics
#' don't show anything really bad.
#'
#' ## Convergence diagnostics
fit_logit4$diagnostic_summary()

draws <- as_draws_rvars(fit_logit4$draws())
summarize_draws(draws)
#' ESS estimates are above the recommended diagnostic thresholds
#' [@Vehtari-Gelman-Simpson-etal:2021], but
#' lower than what we would expect in general from Stan for such a
#' lower dimensional problem.
#'
#' The following figure shows marginal histograms and joint
#' scatterplots, and we can see that `alpha` and `beta[1]` are highly
#' correlated. 
#| label: fig-competing_params_pairs
#| fig-height: 4
#| fig-width: 6
mcmc_pairs(as_draws_array(draws), pars = c("alpha", "beta[1]", "beta[2]"),
           off_diag_args = list(alpha = 0.2))


#' We can compute the correlation.
cor(as_draws_matrix(draws)[, c("alpha", "beta[1]")])[1, 2]

#' The numerical value for the correlation is $-0.999$. The
#' correlation close to 1 can happen also from other reasons (see the
#' next example), but one possibility is that parameters have similar
#' role in the model. Here the reason is the constant column in $x$,
#' which we put there for the demonstration purposes. We may a have
#' constant column, for example, if the predictor matrix is augmented
#' with the intercept predictor, or if the observed data or subdata
#' used in the specific analysis happens to have only one unique
#' value.
#'
#' ## Stan compiler pedantic check
#'
#' Stan compiler pedantic check examining the code can’t
#' recognize this issue, as the problem depends also on the data.
#| results: hide
mod_logit4$check_syntax(pedantic = TRUE)

#'
#' # A posterior with very high correlation
#'
#' In the previous example the two parameters had the same role in the
#' model, leading to high posterior correlation. High posterior
#' correlations are common also in linear models when the predictor
#' values are far from 0. We illustrate this with a linear regression
#' model for the summer temperature in Kilpisjärvi, Finland,
#' 1952--2013. We use the year as the covariate $x$ without centering
#' it.
#'
#' ## Data
#'
#' The data are Kilpisjärvi summer month temperatures 1952-2013
#' measured by Finnish Meteorological Institute.
data_kilpis <- read.delim(root("problems/data", "kilpisjarvi-summer-temp.csv"), sep = ";")
data_lin <- list(M = 1,
                 N = nrow(data_kilpis),
                 x = matrix(data_kilpis$year, ncol = 1),
                 y = data_kilpis[, 5])

#| label: fig-kilpisjarvi_data
#| fig-height: 4
#| fig-width: 6
data.frame(data_lin) |>
  ggplot(aes(x, y)) +
  geom_point(size = 1) +
  labs(y = 'Summer temp. @Kilpisjärvi', x = "Year") +
  guides(linetype = "none")

#' ## Model
#'
#' We use the following Stan linear regression model
code_lin <- root("problems", "linear_glm_kilpis.stan")
#| output: asis
print_stan_file(code_lin)

#| label: fit_lin_kilpis
#| results: hide
mod_lin <- cmdstan_model(stan_file = code_lin)
fit_lin <- mod_lin$sample(data = data_lin, seed = SEED, refresh = 0)

#' ## Convergence diagnostics
#'
#' Stan gives a warning: There were X transitions after warmup that
#' exceeded the maximum treedepth. As in the previous example, there
#' are no other warnings.
fit_lin$diagnostic_summary()

draws <- as_draws_rvars(fit_lin$draws())
summarize_draws(draws)

#' ESS estimates are above the diagnostic threshold, but lower than we
#' would expect for such a low dimensional model, unless there are
#' strong posterior correlations. The following Figure shows the
#' marginal histograms and joint scatterplot for `alpha` and
#' `beta[1]`, which shows they are very highly correlated.
#| label: fig-correlating_params_pairs
#| fig-height: 4
#| fig-width: 6
mcmc_pairs(as_draws_array(draws), pars = c("alpha", "beta"),
           off_diag_args = list(alpha = 0.2))

#' Here the reason is that the $x$ values are in the range 1952--2013,
#' and the intercept `alpha` denotes the temperature at year 0, which
#' is very far away from the range of observed $x$. If the intercept
#' `alpha` changes, the slope `beta` needs to change too. The high
#' correlation makes the inference slower, and we can make it faster
#' by centering $x$. Here we simply subtract 1982.5 from the predictor
#' `year`, so that the mean of $x$ is 0. We could also include the
#' centering and back transformation to Stan code.
#'
#' ## Centered data
#'
data_lin <- list(
  M = 1,
  N = nrow(data_kilpis),
  x = matrix(data_kilpis$year - 1982.5, ncol = 1),
  y = data_kilpis[, 5]
)

#+ message=FALSE, error=FALSE, warning=FALSE
#| label: fit_lin_kilpis_centered_data
#| results: hide
fit_lin <- mod_lin$sample(data = data_lin, seed = SEED, refresh = 0)

#' ## Convergence diagnostics
#'
#' We check the diagnostics
fit_lin$diagnostic_summary()
draws <- as_draws_rvars(fit_lin$draws())
summarize_draws(draws)

#' The following figure shows the scatter plot.
#| label: fig-uncorrelating_params_pairs
#| fig-height: 4
#| fig-width: 6
mcmc_pairs(as_draws_array(draws), pars = c("alpha", "beta"),
           off_diag_args = list(alpha = 0.2))

#'
#' With this change, there is no posterior correlation, Bulk-ESS
#' estimates are 3 times bigger, and the mean time per chain goes from
#' 1.3s to less than 0.05s; that is, we get 2 orders of magnitude
#' faster inference. In a bigger problems this could correspond to
#' reduction of computation time from 24 hours to less than 20
#' minutes.
#'

#'
#' # A bimodal posterior
#'
#' Bimodal distributions can arise from many reasons as in mixture
#' models or models with non-log-concave likelihoods or priors (that
#' is, with distributions with thick tails). We illustrate the
#' diagnostics revealing the multimodal posterior. We use a simple toy
#' problem with $t$ model and data that is not from a $t$
#' distribution, but from a mixture of two normal distributions
#'
#' ## Data
#'
#' Bimodally distributed data
N <- 20
y <- c(rnorm(N / 2, mean = -5, sd = 1), rnorm(N / 2, mean = 5, sd = 1))
data_tt <- list(N = N, y = y)

#' ## Model
#'
#' Unimodal Student's $t$ model:
code_tt <- root("problems", "student.stan")
#| output: asis
print_stan_file(code_tt)
#' Sample
#| label: fit_tt_hard
#| results: hide
mod_tt <- cmdstan_model(stan_file = code_tt)
fit_tt <- mod_tt$sample(data = data_tt, seed = SEED, refresh = 0)

#' ## Convergence diagnostics
#'
#' We check the diagnostics
fit_tt$diagnostic_summary()
draws <- as_draws_rvars(fit_tt$draws())
summarize_draws(draws)

#' The $\widehat{R}$ values for `mu` are large and ESS values for `mu`
#' are small indicating convergence problems. The following figure
#' shows the histogram and trace plots of the posterior draws, clearly
#' showing the bimodality and that chains are not mixing between the
#' modes.
#'
#| label: fig-bimodal1_hist
#| fig-height: 4
#| fig-width: 6
mcmc_hist(as_draws_array(draws), pars = c("mu"))

#' In this toy example, with random initialization each chains has
#' 50\% probability of ending in either mode. We used Stan's default
#' of 4 chains, and when random initialization is used, there is 6\%
#' chance that when running Stan once, we would miss the
#' multimodality. If the attraction areas within the random
#' initialization range are not equal, the probability of missing one
#' mode is even higher. There is a tradeoff between the default
#' computation cost and cost of having higher probability of finding
#' multiple modes. If there is a reason to suspect multimodality, it
#' is useful to run more chains. Running more chains helps to diagnose
#' the multimodality, but the probability of chains ending in
#' different modes can be different from the relative probability mass
#' of each mode, and running more chains doesn't fix this. Other means
#' are needed to improve mixing between the modes (e.g. Yao et al.,
#' 2020) or to approximately weight the chains (e.g. Yao et al.,
#' 2022).
#'
#' # Easy bimodal posterior
#'
#' If the modes in the bimodal distribution are not strongly
#' separated, MCMC can jump from one mode to another and there are no
#' convergence issues.
N <- 20
y <- c(rnorm(N / 2, mean = -3, sd = 1), rnorm(N / 2, mean = 3, sd = 1))
data_tt <- list(N = N, y = y)

#| label: fit_tt_easy
#| results: hide
fit_tt <- mod_tt$sample(data = data_tt, seed = SEED, refresh = 0)

#' ## Convergence diagnostics
#'
#' We check the diagnostics
fit_tt$diagnostic_summary()
draws <- as_draws_rvars(fit_tt$draws())
summarize_draws(draws)

#' Two modes are visible.
#| label: fig-bimodal2_hist
#| fig-height: 4
#| fig-width: 6
mcmc_hist(as_draws_array(draws), pars = c("mu"))

#' Trace plot is not very useful. It shows the chains are jumping
#' between modes, but it's difficult to see whether the jumps happen
#' often enough and chains are mixing well.
#| label: fig-bimodal2_trace
#| fig-height: 4
#| fig-width: 7
mcmc_trace(as_draws_array(draws), pars = c("mu"))

#' Rank ECDF plot [@Sailynoja+etal:2022:PIT-ECDF] indicates good
#' mixing as all chains have their lines inside the envelope (the
#' envelope assumes no autocorrelation, which is the reason to thin
#' the draws here)
#| label: fig-bimodal2_rank_ecdf_diff
#| fig-height: 4
#| fig-width: 6
draws |> thin_draws(ndraws(draws) / ess_basic(draws$mu)) |>
  mcmc_rank_ecdf(pars = c("mu"), plot_diff = TRUE)

#' # Initial value issues
#'
#' MCMC requires some initial values. By default, Stan generates them
#' randomly from [-2,2] in unconstrained space (constraints on
#' parameters are achieved by transformations). Sometimes these
#' initial values can be bad and cause numerical issues. Computers,
#' (in general) use finite number of bits to present numbers and with
#' very small or large numbers, there can be problems of presenting
#' them or there can be significant loss of accuracy.
#'
#' The data is generated from a Poisson regression model. The Poisson
#' intensity parameter has to be positive and usually the latent
#' linear predictor is exponentiated to be positive (the
#' exponentiation can also be justified by multiplicative effects on
#' Poisson intensity).
#'
#' We intentionally generate the data so that there are initialization
#' problems, but the same problem is common with real data when the
#' scale of the predictors is large or small compared to the unit
#' scale. The following figure shows the simulated data.
#'
#' ## Data
set.seed(SEED)
M <- 1
N <- 20
x <- 1e3 * matrix(c(sort(rnorm(N))), ncol = M)
y <- rpois(N, exp(1e-3 * x[, 1]))
data_pois <- list(M = M, N = N, x = x, y = y)

#| label: fig-poisson_data
#| fig-height: 4
#| fig-width: 6
data.frame(data_pois) |>
  ggplot(aes(x, y)) +
  geom_point(size = 3)

#'
#' ## Model
#'
#' We use a Poisson regression model with proper priors. The line
#' `poisson_log_glm(x, alpha, beta)` corresponds to a distribution in
#' which the log intensity of the Poisson distribution is modeled with
#' `alpha + beta * x` but is implemented with better computational
#' efficiency.
code_pois <- root("problems", "pois_glm.stan")
#| output: asis
print_stan_file(code_pois)
#' Sample
#| label: fit_pois
#| results: hide
mod_pois <- cmdstan_model(stan_file = code_pois)
fit_pois <- mod_pois$sample(data = data_pois, seed = SEED, refresh = 0)

#' We get a lot of warnings:
#'
#'```
#' Chain 4 Rejecting initial value:
#' Chain 4   Log probability evaluates to log(0), i.e. negative infinity.
#' Chain 4   Stan can't start sampling from this initial value.
#'```
#'
#' ## Convergence diagnostics
#'
#' We check the diagnostics:
fit_pois$diagnostic_summary()
draws <- as_draws_rvars(fit_pois$draws())
summarize_draws(draws)

#' $\widehat{R}$ values are large and ESS values are small, indicating
#' bad mixing. Marginal histograms and joint scatterplots of the
#' posterior draws in the figure below clearly show that two
#' chains have been stuck away from two others.
#| label: fig-poisson_pairs
#| fig-height: 4
#| fig-width: 6
mcmc_pairs(as_draws_array(draws), pars = c("alpha", "beta"),
           off_diag_args = list(alpha = 0.2))

#' The reason for the issue is that the initial values for
#' `beta` are sampled from $(-2, 2)$ and `x` has some
#' large values. If the initial value for `beta` is higher than
#' about 0.3 or lower than $-0.4$, some of the values of
#' `exp(alpha + beta * x)` will overflow to floating point
#' infinity (`Inf`). 
#'
#' ## Scaled data
#'
#' Sometimes an easy option is to change the initialization range. For
#' example, in this the sampling succeeds if the initial values are
#' drawn from the range $(-0.001, 0.001)$. Alternatively we can scale
#' `x` to have scale close to unit scale. After this scaling, the
#' computation is fast and all convergence diagnostics look good.
data_pois <- list(M = M, N = N, x = x / 1e3, y = y)
#| label: fig-poisson_data2
#| fig-height: 4
#| fig-width: 6
data.frame(data_pois) |>
  ggplot(aes(x, y)) +
  geom_point(size = 3)

#| label: fit_pois_scaled_data
#| results: hide
mod_pois <- cmdstan_model(stan_file = code_pois)
fit_pois <- mod_pois$sample(data = data_pois, seed = SEED, refresh = 0)

#'
#' ## Convergence diagnostics
#'
#' We check the diagnostics:
fit_pois$diagnostic_summary()
draws <- as_draws_rvars(fit_pois$draws())
summarize_draws(draws)

#'
#' If the initial value warning comes only once, it is possible that
#' MCMC was able to escape the bad region and rest of the inference is
#' ok.
#'


#' # Thick tailed posterior
#'
#' We return to the logistic regression example with separable
#' data. Now we use proper, but thick tailed Cauchy prior.
#'
#' ## Model
code_logit_glm4 <- root("problems", "logit_glm4.stan")
#| output: asis
print_stan_file(code_logit_glm4)
#' Sample
#| label: fit_logit_glm4
#| results: hide
mod_logit_glm4 <- cmdstan_model(stan_file = code_logit_glm4)
fit_logit_glm4 <- mod_logit_glm4$sample(data = data_logit, seed = SEED, refresh = 0)

#'
#' ## Convergence diagnostics
#'
#' We check diagnostics
fit_logit_glm4$diagnostic_summary()
draws <- as_draws_rvars(fit_logit_glm4$draws())
summarize_draws(draws)

#' The rounded $\widehat{R}$ values look good, ESS values are
#' low. Looking at the marginal histograms and joint scatterplots of
#' the posterior draws in the following figure show a thick tail.
#| label: fig-thick_tail_pairs
#| fig-height: 4
#| fig-width: 6
mcmc_pairs(as_draws_array(draws), pars = c("alpha", "beta"),
           off_diag_args = list(alpha = 0.2))

#'
#' The dynamic HMC algorithm used by Stan, along with many other MCMC
#' methods, have problems with such thick tails and mixing is
#' slow.
#'
#' Rank ECDF plot indicates good mixing as all chains have their lines
#' inside the envelope (the envelope assumes no autocorrelation, which
#' is the reason to thin the draws here)
#| label: fig-thick_tail_rank_ecdf_diff
#| fig-height: 4
#| fig-width: 6
draws |> thin_draws(ndraws(draws) / ess_bulk(draws$alpha)) |>
  mcmc_rank_ecdf(pars = c("alpha"), plot_diff = TRUE)

#' More iterations confirm a reasonable mixing.
#| results: hide
fit_logit_glm4 <- mod_logit_glm4$sample(data = data_logit, seed = SEED, refresh = 0,
                                        iter_sampling = 4000)

draws <- as_draws_rvars(fit_logit_glm4$draws())
summarize_draws(draws)

#| label: fig-thick_tail_rank_ecdf_diff_more
#| fig-height: 4
#| fig-width: 6
draws |> thin_draws(ndraws(draws) / ess_bulk(draws$alpha)) |>
  mcmc_rank_ecdf(pars = c("alpha"), plot_diff = TRUE)

#' # Funnel
#'
#' A special case of varying curvature is known as the funnel, based
#' on the shape of the typical set of the distribution. Consider a
#' hierarchical model, $y_i\sim\normal(\mu_{k[i]},\sigma)$ for
#' $i=1,\dots,N$ and group membership variable $k[i]$ which takes on
#' values from 1 through $K$.  We shall assume the prior distribution,
#' $\mu_k\sim\normal(\mu_0,\sigma_0),j=1,\dots,J$.  When plotted on
#' the scale $\mu_1,\dots,\mu_K, \log\sigma_0$, this prior can be
#' visualized as a having shape of a funnel
#'
#' If the funnel-shaped prior is combined with a weak likelihood, the
#' posterior is also funnel shaped. As a toy example, we use the
#' Kilpisjärvi temperature data, with each group being one year, with
#' three summer month temperatures per year. With only three
#' observations per group, the likelihood is weak for each $\mu_k$ and
#' the prior is likely to dominate the posterior shape. The number of
#' groups is 71, and this high dimensionality makes the funnel
#' challenging.
data_kilpis <- read.delim(root("problems/data", "kilpisjarvi-summer-temp.csv"), sep = ";")
data_grpy <-list(N = length(data_kilpis$year)*ncol(data_kilpis[,2:4]),
             K = length(data_kilpis$year),
             x = rep(1:length(data_kilpis$year), ncol(data_kilpis[,2:4])),
             y = c(t(t(data_kilpis[,2:4]))))

#' Here is a direct implementation of the hierarchical model in
#' Stan. The parameterization used is also known as centered
#' parameterization.
code_hier_cp <- root("problems", "hier_cp.stan")
#| output: asis
print_stan_file(code_hier_cp)

#' We first try running Stan with its default settings.
#| label: fit_hier_cp
#| results: hide
SEED <- 48929
mod_hier_cp <- cmdstan_model(stan_file = code_hier_cp)
fit_hier_cp <- mod_hier_cp$sample(data = data_grpy, seed = SEED, refresh = 0)

#' We get a warning that some transitions ended with a divergence. The
#' convergence diagnostics $\widehat{R}$, bulk-ESS, and tail-ESS
#' reveal that the chains are not mixing well:
fit_hier_cp

#' Plot scatter plot of $\mu_1$ vs $\log\sigma_0$ with divergences shown in red (bayesplot)
np <- fit_hier_cp$sampler_diagnostics(format="df")|>
  mutate(Chain=.chain,Iteration=.iteration,Parameter="divergent__",Value=divergent__)|>
  select(Chain,Iteration,Parameter,Value)
fit_hier_cp$draws(format="df")|>
  mutate(log_sigma0=log(sigma0))|>
  mcmc_scatter(pars=c("mu[1]","sigma0"),transform=list(sigma0="log"), alpha=0.1, shape=20,
               np=np, np_style = scatter_style_np(div_shape = 18, div_size = 3)) +
  labs(x=TeX(r"($\mu_1$)"), y=TeX(r"($\log\,\sigma_0$)"), title="a) Centered param.")

#' Plot scatter plot of $\mu_1$ vs $\log\sigma_0$ with divergences shown in red (ggplot)
drws <- bind_draws(fit_hier_cp$draws(format="df"), fit_hier_cp$sampler_diagnostics(format="df")) |>
  mutate(log_sigma0=log(sigma0))
p1 <- ggplot(data = NULL, aes(x=`mu[1]`,`log_sigma0`)) +
  geom_point(data = drws |> filter(divergent__==0), shape = 20, color = bayesplot:::get_color("m"), fill = bayesplot:::get_color("lh"), alpha = 0.1, size = 2) +
  geom_point(data = drws |> filter(divergent__==1), shape = 23, fill = "red", color = "white", size = 2) +
  labs(x=TeX(r"($\mu_1$)"), y=TeX(r"($\log\,\sigma_0$)"), title="c) Centered param.")
p1

#' We change the adapt_delta tuning parameter of the NUTS algorithm to
#' force a smaller step size.
#| label: fit_hier_cp_999
#| results: hide
fit_hier_cp_999 <- mod_hier_cp$sample(data = data_grpy, seed = SEED, refresh=0, adapt_delta=0.999)

#' However, the convergence diagnostics still indicate serious mixing problems:
fit_hier_cp_999

#' Plot scatter plot of $\mu_1$ vs $\log\sigma_0$ with divergences shown in red (bayesplot)
np <- fit_hier_cp_999$sampler_diagnostics(format="df")|>
  mutate(Chain=.chain,Iteration=.iteration,Parameter="divergent__",Value=divergent__)|>
  select(Chain,Iteration,Parameter,Value)
fit_hier_cp_999$draws(format="df")|>
  mutate(log_sigma0=log(sigma0))|>
  mcmc_scatter(pars=c("mu[1]","sigma0"),transform=list(sigma0="log"), alpha=0.1, shape=20,
               np=np, np_style = scatter_style_np(div_shape = 18, div_size = 3)) +
  labs(x=TeX(r"($\mu_1$)"), y=TeX(r"($\log\,\sigma_0$)"), title="a) Centered param. + higher adapt_delta")

#' Plot scatter plot of $\mu_1$ vs $\log\sigma_0$ with divergences shown in red (ggplot)
drws <- bind_draws(fit_hier_cp_999$draws(format="df"), fit_hier_cp_999$sampler_diagnostics(format="df")) |>
  mutate(log_sigma0=log(sigma0))
p2 <- ggplot(data = NULL, aes(x=`mu[1]`,`log_sigma0`)) +
  geom_point(data = drws |> filter(divergent__==0), shape = 20, color = bayesplot:::get_color("m"), fill = bayesplot:::get_color("lh"), alpha = 0.1, size = 2) +
  geom_point(data = drws |> filter(divergent__==1), shape = 23, fill = "red", color = "white", size = 2) +
  labs(x=TeX(r"($\mu_1$)"), y=TeX(r"($\log\,\sigma_0$)"), title="b) Centered param. + higher adapt_delta")
p2

#' The usual approach to resolve the funnel problem is to change how
#' the model is parameterized. The so-called non-centered
#' parameterization provides the same model, but the sampling happens
#' in a transformed space that does not have the difficult funnel
#' geometry. 
code_hier_ncp <- root("problems", "hier_ncp.stan")
#| output: asis
print_stan_file(code_hier_ncp)

#' We run Stan with its default settings.
#| label: fit_hier_ncp
#| results: hide
mod_hier_ncp <- cmdstan_model(stan_file = code_hier_ncp)
fit_hier_ncp <- mod_hier_ncp$sample(data = data_grpy, seed = SEED, refresh=0)

#' The convergence diagnostics $\widehat{R}$, bulk-ESS, and tail-ESS
#' look good now.
fit_hier_ncp

#' Plot scatter plot of $\mu_1$ vs $\log\sigma_0$ with divergences shown in red (bayesplot)
np <- fit_hier_ncp$sampler_diagnostics(format="df")|>
  mutate(Chain=.chain,Iteration=.iteration,Parameter="divergent__",Value=divergent__)|>
  select(Chain,Iteration,Parameter,Value)
fit_hier_ncp$draws(format="df")|>
  mutate(log_sigma0=log(sigma0))|>
  mcmc_scatter(pars=c("mu[1]","sigma0"),transform=list(sigma0="log"), alpha=0.1, shape=20,
               np=np, np_style = scatter_style_np(div_shape = 18, div_size = 3)) +
  labs(x=TeX(r"($\mu_1$)"), y=TeX(r"($\log\,\sigma_0$)"), title="a) Non-centered param.")

#' Plot scatter plot of $\mu_1$ vs $\log\sigma_0$ with divergences shown in red (ggplot)
drws <- bind_draws(fit_hier_ncp$draws(format="df"), fit_hier_ncp$sampler_diagnostics(format="df")) |>
  mutate(log_sigma0=log(sigma0))
p3 <- ggplot(data = NULL, aes(x=`mu[1]`,`log_sigma0`)) +
  geom_point(data = drws |> filter(divergent__==0), shape = 20, color = bayesplot:::get_color("m"), fill = bayesplot:::get_color("lh"), alpha = 0.1, size = 2) +
  geom_point(data = drws |> filter(divergent__==1), shape = 23, fill = "red", color = "white", size = 2) +
  labs(x=TeX(r"($\mu_1$)"), y=TeX(r"($\log\,\sigma_0$)"), title="c) Non-centered param.")
p3

#' If we compare the scatter plots side by side, we clearly see that
#' increasing `adapt_delta` and getting rid of divergences did not
#' solve the funnel problem and the posterior estimates with centered
#' parameterization would be biased.
#| label: fig-kilpis-funnel
#| fig-width: 14
#| fig-height: 4.5
#| out-width: 100%
(p1 + p2 + p3) *
  scale_y_continuous(lim=c(-8,0.3)) *
  scale_x_continuous(lim=c(7.2,11.8)) *
  theme(plot.title = element_text(size=16)) +
  plot_layout(axis_titles="collect_y") 

#'
#' # Variance parameter that is not constrained to be positive
#'
#' Demonstration what happens if we forget to constrain a parameter
#' that has to be positive. In Stan the constraint can be added when
#' declaring the parameter as `real<lower=0> sigma;`
#'
#' ## Data
#'
#' We simulated x and y independently from independently from
#' normal(0,1) and normal(0,0.1) respectively. As $N=8$ is small,
#' there will be a lot of uncertainty about the parameters including
#' the scale sigma.
M <- 1
N <- 8
set.seed(SEED)
x <- matrix(rnorm(N), ncol = M)
y <- rnorm(N) / 10
data_lin <- list(M = M, N = N, x = x, y = y)

#' ## Model
#'
#' We use linear regression model with proper priors.
code_lin <- root("problems", "linear_glm.stan")
#| output: asis
print_stan_file(code_lin)
#' Sample
#| label: fit_lin
#| results: hide
mod_lin <- cmdstan_model(stan_file = code_lin)
fit_lin <- mod_lin$sample(data = data_lin, seed = SEED, refresh = 0)

#' We get many times the following warnings
#'```
#' Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#' Chain 4 Exception: normal_id_glm_lpdf: Scale vector is -0.747476, but must be positive finite! (in '/tmp/RtmprEP4gg/model-7caa12ce8e405.stan', line 16, column 2 to column 43)
#' Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#' Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#'```
#'
#' Sometimes these warnings appear in early phase of the sampling,
#' even if the model has been correctly defined. Now we have too many
#' of them, which indicates the samples is trying to jump to
#' infeasible values, which here means the negative scale parameter
#' values. Many rejections may lead to biased estimates.
#'
#' There are some divergences reported, which is also indication that
#' there might be some problem (as divergence diagnostic has an ad hoc
#' diagnostic threshold, there can also be false positive
#' warnings). Other convergence diagnostics are good, but due to many
#' rejection warnings, it is good to check the model code and
#' numerical accuracy of the computations.
#'
#' ## Convergence diagnostics
#'
#' We check diagnostics
fit_lin$diagnostic_summary()
draws <- as_draws_rvars(fit_lin$draws())
summarize_draws(draws)

#' ## Stan compiler pedantic check
#'
#' Stan compiler pedantic check can recognize that `A normal_id_glm
#' distribution is given parameter sigma as a scale parameter
#' (argument 4), but sigma was not constrained to be strictly
#' positive. The pedantic check is also warning about the very wide
#' priors.
#| results: hide
mod_lin$check_syntax(pedantic = TRUE)

#' After fixing the model with proper parameter constraint, MCMC runs
#' without warnings and the sampling efficiency is better. In this
#' specific case, the bias is negligible when running MCMC with the
#' model code without the constraint, but it is difficult to diagnose
#' without running the fixed model.
#'
#' Fixed model includes <lower=0> constraint for sigma.
code_lin2 <- root("problems", "linear_glm2.stan")
#| output: asis
print_stan_file(code_lin2)
#' Sample
#| label: fit_lin2
#| results: hide
mod_lin2 <- cmdstan_model(stan_file = code_lin2)
fit_lin2 <- mod_lin2$sample(data = data_lin, seed = SEED, refresh = 0)

#' We check diagnostics
draws2 <- as_draws_rvars(fit_lin2$draws())
summarize_draws(draws2)

#' In this specific case, the bias is negligible when running MCMC
#' with the model code without the constraint, but it is difficult to
#' diagnose without running the fixed model.
#'
#' <br />
#'
#' # References {.unnumbered}
#'
#' ::: {#refs}
#' :::
#'
#' # Licenses {.unnumbered}
#'
#' * Code &copy; 2021--2025, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2021--2025, Aki Vehtari, licensed under CC-BY-NC 4.0.

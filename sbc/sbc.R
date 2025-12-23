#' ---
#' title: "Simulation based calibration checking workflow"
#' author: "Martin Modrák"
#' date: 2022-12-16
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

#' # Introduction
#'
#' This notebook includes the code for the Bayesian Workflow book
#' Chapter 31 *Simulation-based calibration checking in model
#' development workflow*.
#' 
#' Here we describe a complete process to iteratively build and
#' validate the _implementation_ of a non-trivial, but still
#' relatively small model using simulation based calibration checking
#' (SBC, @sbc2023).  This is not a full Bayesian Workflow, instead the process
#' described here can be thought of as a subroutine in the full
#' workflow: here we take a relatively precise description of a model
#' as input and try to produce a Stan program that implements this
#' model. Once we have a Stan program we trust, it is still necessary
#' to validate its fit to actual data and other properties, which may
#' trigger a need to change the model. At this point you may want to
#' go back to simulations and make sure the modified model is
#' implemented correctly.
#' 
#' The workflow described here focuses on small models.  "Small" means
#' that the model is relatively fast to fit and we don't have to worry
#' about computation too much.  Still many of the approaches here
#' also apply to complex models (especially starting small and
#' building smaller submodels separately), and with proper separation
#' of the model into submodels, one can validate big chunks of Stan
#' code while working with small models only.
#' 
#' We expect the reader to be familiar with basics of the package. If not,
#' check out the [*Getting Started with SBC* vignette](https://hyunjimoon.github.io/SBC/articles/SBC.html).
#'
#' ## Example model
#' 
#' The example we'll investigate is building a two-component Poisson mixture,
#' where the mixing ratio is allowed to vary with some predictors while the means
#' of the components are the same for all observations.
#' A somewhat contrived real world situation where this could be a useful model:
#' there are two sub-species of an animal that are hard to observe directly, but leave
#' droppings (poop) behind, that we can find. Further, we know the subspecies differ in the 
#' average number of droppings they leave at one place. So we can take the number of droppings as a noisy
#' information about which subspecies was present at given location.
#' We observe the number of droppings at multiple locations and record some environmental covariates about the locations
#' (e.g. temperature, altitude) and want to learn something about the association between
#' those covariates and the prevalence of either subspecies.
#' 
#' A mathematical description would be:
#' 
#' $$
#' \begin{align*}
#'     y_i &\sim \mathrm{Poisson}(\mu_{z_i}) \\
#'     z_i &\sim \mathrm{Bernoulli}(\mathrm{logit}^{-1}(\theta_i)) \\
#'     \theta_i &= \alpha + \sum_k \beta_k \mathbf{X}_{i,k} \\
#'     \log \mu_{\{1,2\}} &\sim \mathrm{N(3, 1)} \\
#'     \alpha &\sim N(0, 2) \\
#'     \beta_k &\sim N(0, 1) 
#' \end{align*}
#' $$
#' 
#' This model naturally decomposes into two submodels: 
#' 
#' 1) the mixture submodel where the mixing ratio is the same
#' for all observations
#' 
#' 2) a beta regression where we take covariates and make a prediction of a probability,
#' assuming we (noisily) observe the probability.
#' 
#' It is good practice to start small and implement and validate each
#' of those submodels separately and then put them together and
#' validate the bigger model.  This makes is substantially easier to
#' locate bugs.  You'll notice that the process ends up involving a
#' lot of steps, but the fact is that we still ignore all the
#' completely invalid models we created while writing this vignette
#' (typos, compile errors, dimension mismatches, ...). Developing
#' models you can trust is hard work. More experienced users can
#' definitely make bigger steps at once, but we strongly discourage
#' anyone from writing a big model in one go.
#' 
#' ## Setup
#+ setup, include = FALSE
knitr::opts_chunk$set(
  cache = FALSE,
  message = FALSE,
  error = FALSE,
  warning = FALSE,
  comment = NA,
  out.width = "95%"
)
#'
#' Let's setup and get our hands dirty.
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
# remotes::install_github("hyunjimoon/SBC")
library(SBC)
options(SBC.min_chunk_size = 5)
library(cmdstanr)
library(posterior)
options(pillar.neg = FALSE,
        pillar.subtle = FALSE,
        pillar.sigfig = 2)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(patchwork) # Needed only for saving plots for the book
library(future)
plan(multisession) 

# Setup caching of results
cache_dir <- root("sbc", "_cache")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

#' # Mixture submodel
#' 
#' There is a good [guide to mixtures](https://mc-stan.org/docs/2_27/stan-users-guide/mixture-modeling-chapter.html) in the Stan user's guide.
#' Following the user's guide would save us from a lot of mistakes,
#' but for the sake of example, we will pretend we didn't really read
#' it - and we'll see the problems can be discovered via simulations.
#' 
#' So this is our first try at implementing the mixture submodel:
cat(readLines(root("sbc", "models/mixture_first.stan")), sep = "\n")
#| label: model_first
model_first <- cmdstan_model(root("sbc", "models/mixture_first.stan"))
backend_first <- SBC_backend_cmdstan_sample(model_first) 

#' And this is our code to simulate data for this model:
generator_func_first <- function(N) {
  mu1 <- rnorm(1, 3, 1)
  mu2 <- rnorm(1, 3, 1)
  theta <- runif(1)
  
  y <- numeric(N)
  for (n in 1:N) {
    if (runif(1) < theta) {
      y[n] <- rpois(1, exp(mu1))
    } else {
      y[n] <- rpois(1, exp(mu2))
    }
  }
  list(variables = list(mu1 = mu1,
                        mu2 = mu2,
                        theta = theta),
       generated = list(N = N,
                        y = y))
}
generator_first <- SBC_generator_function(generator_func_first, N = 50)

#' Let's start with just a single simulation:
set.seed(68455554)
datasets_first <- generate_datasets(generator_first, 1)
results_first <- compute_SBC(datasets_first, backend_first, 
                             cache_mode = "results", 
                             cache_location = file.path(cache_dir, "mixture_first"))

#' Let's examine the MCMC pairs plots
#| label: fig-sbcworkflow_mixture_first_pairs
# Fixing the condition for above/over diagonal chains, in a minority
# of runs the plot shows the problem less clearly, as discussed at
# https://github.com/stan-dev/bayesplot/issues/132
p_cond <- pairs_condition(chains = list(c(1, 3), c(2, 4)))
mixture_first_pairs <- mcmc_pairs(results_first$fits[[1]]$draws(),
                                  condition = p_cond,
                                  np = nuts_params(results_first$fits[[1]]))
mixture_first_pairs

#' There are divergences. One thing that stands out is that either
#' `mu1` is tightly determined and `mu2` is allowed the full prior
#' range or the other way around. We also don't learn anything about
#' theta.
#' 
#' This might be puzzling but relates to bad usage of `log_mix`. The
#' thing is that `poisson_log_lpmf(y | mu1)` returns a single number -
#' the total log likelihood of all elements of `y` given `mu1`. And
#' thus we are building a mixture where either all observations are
#' from the first component or all are from the second component. To
#' implement mixture where each observation is allowed to come from a
#' different component, we need to loop over observations and do a
#' separate `log_mix` call for each.
#' 
#' More details on the mathematical background are explained in the
#' ["Vectorizing mixtures"](https://mc-stan.org/docs/2_27/stan-users-guide/vectorizing-mixtures.html)
#' section of Stan User's guide.
#' 
#' ## Fixed mixture model
#' 
#' We make a new model fixing the `log_mix` problem.
#| label: model_fixed_log_mix
cat(readLines(root("sbc", "models/mixture_fixed_log_mix.stan")), sep = "\n")
model_fixed_log_mix <- cmdstan_model(root("sbc", "models/mixture_fixed_log_mix.stan"))
backend_fixed_log_mix <- SBC_backend_cmdstan_sample(model_fixed_log_mix)

#' So let's try once again with the same single simulation:
results_fixed_log_mix <- compute_SBC(datasets_first,
                                     backend_fixed_log_mix, 
                                     cache_mode = "results", 
                                     cache_location = file.path(cache_dir, "mixture_fixed_log_mix"))

#' No warnings this time. We look at the stats:
print(results_fixed_log_mix$stats, digits = 2)

#' We see nothing obviously wrong, the posterior means are relatively
#' close to simulated values (as summarised by the z-scores) - no
#' variable is clearly ridiculously misfit. So let's run a few more
#' iterations.
set.seed(8314566)
datasets_first_10 <- generate_datasets(generator_first, 10)

results_fixed_log_mix_2 <- compute_SBC(datasets_first_10,
                                       backend_fixed_log_mix, 
                                       cache_mode = "results", 
                                       cache_location = file.path(cache_dir, "mixture_fixed_log_mix_2"))

#' So there are some problems - we have quite a bunch of high R-hat
#' and low ESS values. This is the distribution of all rhats:
hist(results_fixed_log_mix_2$stats$rhat)

#' Let's examine a single pairs plot:
#| label: fig-sbcworkflow_mixture_fixed_log_mix_pairs
mixture_fixed_log_mix_pairs <- mcmc_pairs(results_fixed_log_mix_2$fits[[1]]$draws())
mixture_fixed_log_mix_pairs

#' We clearly see two modes in the posterior. And upon reflection, we
#' can see why: swapping `mu1` with `mu2` while also changing `theta`
#' for `1 - theta` gives _exactly_ the same likelihood - because the
#' ordering does not matter. A more detailed explanation of this type
#' of problem is at
#' https://betanalpha.github.io/assets/case_studies/identifying_mixture_models.html
#' 
#' ## Fixed parameter ordering
#' 
#' We can easily fix the ordering of the `mu`s by using the `ordered` built-in type.
#| label: model_fixed_ordered
cat(readLines(root("sbc", "models/mixture_fixed_ordered.stan")), sep = "\n")
model_fixed_ordered <- cmdstan_model(root("sbc", "models/mixture_fixed_ordered.stan"))
backend_fixed_ordered <- SBC_backend_cmdstan_sample(model_fixed_ordered) 

#' We also need to update the generator to match the new names and ordering constant:
generator_func_ordered <- function(N) {
  # If the priors for all components of an ordered vector are the same
  # then just sorting the result of a generator is enough to create
  # a valid draw from the ordered vector prior
  mu <- sort(rnorm(2, 3, 1)) 
  theta <- runif(1)
  y <- numeric(N)
  for (n in 1:N) {
    if (runif(1) < theta) {
      y[n] <- rpois(1, exp(mu[1]))
    } else {
      y[n] <- rpois(1, exp(mu[2]))
    }
  }
  list(variables = list(mu = mu, theta = theta),
       generated = list(N = N, y = y))
}
generator_ordered <- SBC_generator_function(generator_func_ordered, N = 50)

#' We are kind of confident (and the model fits quickly), so we'll already start with 10 simulations.
set.seed(3785432)
datasets_ordered_10 <- generate_datasets(generator_ordered, 10)

results_fixed_ordered <- compute_SBC(datasets_ordered_10,
                                     backend_fixed_ordered, 
                                     cache_mode = "results", 
                                     cache_location = file.path(cache_dir, "mixture_fixed_ordered"))

#' Now some fits still produce problematic Rhats or divergent
#' transitions, let's browse the `$backend_diagnostics` (which contain
#' Stan-specific diagnostic values) to see which simulations are
#' causing problems:
print(results_fixed_ordered$backend_diagnostics, digits = 2)

#' One of the fits has quite a lot of divergent transitions. Let's
#' look at the pairs plot for the model:
#| label: sbcworkflow_mixture_fixed_ordered_pairs
problematic_fit_id <- 2
problematic_fit <- results_fixed_ordered$fits[[problematic_fit_id]]
mixture_fixed_ordered_pairs <- mcmc_pairs(problematic_fit$draws(),
                                          np = nuts_params(problematic_fit))
mixture_fixed_ordered_pairs

#' There is a lot of ugly stuff going on. Notably, one can notice that
#' the posterior of theta is bimodal, preferring either almost 0 or
#' almost 1 - and when that happens, the mean of one of the components
#' is almost unconstrained. Why does that happen? The key to the
#' answer is in the simulated values for the component means:
subset_draws(datasets_ordered_10$variables, draw = problematic_fit_id)

#' We were unlucky enough to simulate data where both components have
#' almost the same mean and thus we are actually looking at data that
#' is not really a mixture. Mixture models can misbehave badly in such
#' cases (see once again the [case study by Michael
#' Betancourt](https://betanalpha.github.io/assets/case_studies/identifying_mixture_models.html#5_singular_components_and_computational_issues)
#' for a bit more detailed dive into this particular problem).
#' 
#' ## Fixing degenerate components?
#' 
#' What to do about this? Fixing the model to handle such cases
#' gracefully is hard. But the problem is basically our prior - we
#' want to express that (since we are fitting a two component model),
#' we don't expect the means to be too similar. So if we can change
#' our simulation to avoid this, we'll be able to proceed with SBC. If
#' such a pattern appeared in real data, we would still have a
#' problem, but we would notice thanks to the diagnostics.
#' 
#' This can definitely be done. But another way is to just ignore the
#' simulations that had divergences for SBC calculations. It turns out
#' that if we remove simulations in a way that only depends on the
#' observed data (and not on unobserved variables), the SBC identity
#' is preserved and we can use SBC without modifications. The
#' resulting check is however telling us something only for data that
#' were not rejected. In this case this is not a big issue: if a fit
#' had divergent transitions, we would not trust it anyway, so
#' removing fits with divergent transitions is not such a big deal.
#' 
#' For more details see the
#' [`rejection_sampling`](https://hyunjimoon.github.io/SBC/articles/rejection_sampling.html)
#' vignette.
#' 
#' So let us subset the results to avoid divergences:
sim_ids_to_keep <- 
  results_fixed_ordered$backend_diagnostics$sim_id[
    results_fixed_ordered$backend_diagnostics$n_divergent == 0]
# Equivalent tidy version if you prefer
# sim_ids_to_keep <- results_fixed_ordered$backend_diagnostics %>% 
#   dplyr::filter(n_divergent == 0) %>%
#   dplyr::pull(sim_id)
results_fixed_ordered_subset <- results_fixed_ordered[sim_ids_to_keep]
summary(results_fixed_ordered_subset)

#' This gives us no obvious problems.
#| label: fig-rank_hist-fixed_ordered_subset
plot_rank_hist(results_fixed_ordered_subset)
#| label: fig-ecdf_hist-fixed_ordered_subset
plot_ecdf_diff(results_fixed_ordered_subset)

#' Since we now have only `r length(results_fixed_ordered_subset)`
#' simulations, it is not surprising that we are still left with a
#' huge uncertainty about the actual coverage of our posterior
#' intervals - we can see that in a plot:
#| label: fig-plot_coverage-fixed_ordered_subset
plot_coverage(results_fixed_ordered_subset)

#' The coverage plot shows the observed coverage of central posterior
#' intervals of varying width and the associated uncertainty (black +
#' grey), the blue line represents perfect calibration.
#'
#' Or investigate numerically.
coverage <- empirical_coverage(results_fixed_ordered_subset$stats,
                               width = c(0.5, 0.9, 0.95))
coverage

theta_90_coverage_string <- paste0(round(100 * as.numeric(
  coverage[coverage$variable == "theta" & coverage$width == 0.9, c("ci_low", "ci_high")])),
  "%",
  collapse = " - ")


#' We can clearly see that while there are no terrible errors, a quite
#' big miscalibration is still consistent with the SBC results so far,
#' for example the 90% posterior interval for `theta` could (as far as
#' we know) contain `r theta_90_coverage_string` of the true
#' values. That's not very reassuring.
#' 
#' 
#' So we can run for more iterations - to reduce memory consumption,
#' we set `keep_fits = FALSE`.  You generally don't want to do this
#' unless you are really short on memory, as it makes you unable to
#' inspect any problems in your fits:
set.seed(54987622)
datasets_ordered_100 <- generate_datasets(generator_ordered, 100)
results_fixed_ordered_100 <- compute_SBC(datasets_ordered_100, backend_fixed_ordered, 
                    keep_fits = FALSE, cache_mode = "results", 
                    cache_location = file.path(cache_dir, "mixture_fixed_ordered_100"))

#' Once again we subset to keep only non-divergent fits - this also
#' removes all the problematic Rhats and ESS.
sim_ids_to_keep <- 
  results_fixed_ordered_100$backend_diagnostics$sim_id[
    results_fixed_ordered_100$backend_diagnostics$n_divergent == 0]
# Equivalent tidy version
# sim_ids_to_keep <- results_fixed_ordered_100$backend_diagnostics %>% 
#   dplyr::filter(n_divergent == 0) %>%
#   dplyr::pull(sim_id)

results_fixed_ordered_100_subset <- results_fixed_ordered_100[sim_ids_to_keep]
summary(results_fixed_ordered_100_subset)

#' And we can use `bind_results` to combine the new results with the
#' previous fits to not waste our computational effort.
results_fixed_ordered_combined <-
  bind_results(results_fixed_ordered_subset,
               results_fixed_ordered_100_subset)

#| label: fig-sbcworkflow_results_fixed_ordered_combined_results
ordered_combined_rank_hist <- plot_rank_hist(results_fixed_ordered_combined)
ordered_combined_ecdf_diff <- plot_ecdf_diff(results_fixed_ordered_combined)
ordered_combined_rank_hist / ordered_combined_ecdf_diff

#' Seems fairly well within the expected bounds. We could definitely
#' run more iterations if we wanted to have a more strict check, but
#' for now, we are happy and the remaining uncertainty about the
#' coverage of our posterior intervals is no longer huge, so it is
#' highly unlikely there is some big bug lurking down there. While we
#' see a potential problem where the coverage for `mu[1]` and `mu[2]`
#' is no longer consistent with perfect calibration, the `ecdf_diff`
#' plot takes precedence as the uncertainty in the coverage plot is
#' only approximate and we thus cannot take it too seriously (see
#' `help("empirical_coverage")` for some more details).
#| label: fig-plot_coverage-fixed_ordered_combined
plot_coverage(results_fixed_ordered_combined)

#' Note: it turns out that extending the model to more components
#' becomes somewhat tricky as the model can become sensitive to
#' initialization. Also the problems with data that can be explained
#' by fewer components than the model assumes become more prevalent.
#'
#' # Logistic regression submodel
#'
#' Let's move to the logistic regression submodel of our model.
#| label: model_logistic_first
cat(readLines(root("sbc", "models/logistic_first.stan")), sep = "\n")
model_logistic_first <- cmdstan_model(root("sbc", "models/logistic_first.stan"))
backend_logistic_first <- SBC_backend_cmdstan_sample(model_logistic_first) 

#' If you are good at reading code, you may notice there is a fatal
#' bug and fix it, if you do not, we can use simulations to discover
#' the bug.
#'
#' We also write a matching generator to check that the matrix algebra
#' in the model is correct.
generator_func_logistic_first <- function(N_obs, N_predictors) {
  alpha <- rnorm(1, 0, 2)
  beta <- rnorm(N_predictors, 0, 1)
  X <- matrix(rnorm(N_predictors * N_obs, 0, 1), nrow = N_obs, ncol = N_predictors)
  linpred <- array(alpha, N_obs)
  for (p in 1:N_predictors) {
    linpred <- linpred + X[, p] * beta[p]
  }
  y <- rbinom(N_obs, size = 1, prob = plogis(linpred))
  list(variables = list(alpha = alpha, beta = beta),
       generated = list(N_obs = N_obs, N_predictors = N_predictors, y = y, X = X)
  )
}
generator_logistic_first <- SBC_generator_function(generator_func_logistic_first, N_obs = 50, N_predictors = 2)

#' We'll start with 10 simulations once again.
set.seed(2345688)
datasets_logistic_first <- generate_datasets(generator_logistic_first, 10)

results_logistic_first_10 <- compute_SBC(datasets_logistic_first,
                                         backend_logistic_first, 
                                         cache_mode = "results", 
                                         cache_location = file.path(cache_dir, "logistic_first_10"))

#' Already with 10 datasets we are likely to see suspicious rank/ECDF plots:
#| label: fig-sbcworkflow_logistic_first_results
logistic_first_ranks <- plot_rank_hist(results_logistic_first_10)
logistic_first_ecdf <- plot_ecdf_diff(results_logistic_first_10)
logistic_first_ranks / logistic_first_ecdf

#' At this point, we could use more simulations to see if the
#' discrepancy is real. But this is also a good opportunity to
#' introduce additional ways to diagnose mismatches between the Stan
#' code and the simulator. A simple diagnostic is to plot the
#' simulated values against posterior estimates, which can be done via
#' the [plot_sim_estimated()] function.
#| label: fig-sbcworkflow_logistic_first_sim_estimated
logistic_first_sim_estimated <- plot_sim_estimated(results_logistic_first_10) + 
  labs(x = "Simulated value", y = "Mean, 95% CI")
logistic_first_sim_estimated

#' One thing that immediately stands out is that the posterior
#' inferences for `alpha` seem to be independent of the simulated
#' value.  Indeed, in our Stan code, `alpha` never enters the
#' likelihood part of the model: where we have ``X * beta` we should
#' have had `alpha + X * beta`.
#'
#' This problem manifests as suspicious rank/ECDF plots for the `beta`
#' parameters. In fact, using this model, SBC will never show a
#' failure for `alpha`, as just sampling from the prior (and ignoring
#' likelihood) will always satisfy the SBC equality.
#'
#' An additional useful diagnostic in this case is running SBC for a
#' derived quantity: the SBC equality has to be satisfied not only for
#' parameters of the model, but also for all quantities derived from
#' parameters and data. Quite often, adding the model likelihood as an
#' additional quantity increases sensitivity of SBC checks, as the
#' likelihood is a complex function of all parameters. So lets do just
#' that. We don't need to refit the model, we just call
#' `recompute_SBC_statistics` to compute the derived quantities.
logistic_loglik_dq <- derived_quantities(
  log_lik =  sum(dbinom(y,
                        size = 1,
                        prob = plogis(alpha + X %*% beta),
                        log = TRUE)))
results_logistic_first_10_dq <-
  recompute_SBC_statistics(results_logistic_first_10,
                           datasets_logistic_first,
                           backend_logistic_first,
                           dquants = logistic_loglik_dq)

#' The rank and ECDF plots are shown below.
#| label: fig-sbcworkflow_logistic_first_results_dq
logistic_first_ranks <- plot_rank_hist(results_logistic_first_10_dq)
logistic_first_ranks <- plot_ecdf_diff(results_logistic_first_10_dq)
logistic_first_ranks / logistic_first_ecdf

#' While the failures for the `beta` parameters are barely visible
#' with 10 simulations, `log_lik` signals a clear failure.
#| label: fig-sbcworkflow_results_dq_loglik_only
results_dq_loglik_only <- plot_rank_hist(results_logistic_first_10_dq, variables = "log_lik") /
  plot_ecdf_diff(results_logistic_first_10_dq, variables = "log_lik")
results_dq_loglik_only

#' ## Merging the intercept with predictors
#'
#' One way to resolve the problem -- and simplify our Stan code --- is
#' by treating the intercept as just another predictor which happens
#' to have all 1's in its column of \texttt{X} and has a different
#' prior. This is also how most common regression modelling packages
#' handle the situation. We thus modify our Stan code to:
cat(readLines(root("sbc", "models/logistic_merged_intercept.stan")), sep = "\n")

#' This looks cleaner, but you may notice one additional issue that we
#' created during the rewrite. We will see that it will quickly
#' manifest if we use SBC. We also need to modify our simulation code
#' in R to include intercept in the design matrix, but in R we will
#' keep the explicit loop to decrease chances of having the same
#' problem in both R and Stan.
#| label: model_logistic_merged_intercept
model_logistic_merged_intercept <- cmdstan_model(root("sbc", "models/logistic_merged_intercept.stan"))
backend_logistic_merged_intercept <- SBC_backend_cmdstan_sample(model_logistic_merged_intercept, chains = 2) 

#' We now update the generator code to match:
generator_func_logistic_merged_intercept <- function(N_obs, N_predictors) {
  beta <- c(rnorm(1, 0, 2), rnorm(N_predictors - 1, 0, 1))
  X <- matrix(rnorm(N_predictors * N_obs, 0, 1), nrow = N_obs, ncol = N_predictors)
  X[, 1] <- 1 # Intercept
  y <- array(NA_real_, N_obs)
  linpred <- array(0, N_obs)
  for (p in 1:N_predictors) {
    linpred <- linpred + X[, p] * beta[p]
  }
  y <- rbinom(N_obs, size = 1, prob = plogis(linpred))
  list(
    variables = list(beta = beta),
    generated = list(N_obs = N_obs,
                     N_predictors = N_predictors,
                     y = y,
                     X = X))
}
generator_logistic_merged_intercept <- SBC_generator_function(generator_func_logistic_merged_intercept, N_obs = 50, N_predictors = 3)

set.seed(125488)
datasets_logistic_merged_intercept_10 <- generate_datasets(generator_logistic_merged_intercept, 10)

#' And we also need to update the definition of our derived quantity to reflect the renaming.
logistic_merged_intercept_loglik_dq <- derived_quantities(log_lik = sum(dbinom(y, size = 1, prob = plogis(X %*% beta), log = TRUE)))

results_logistic_merged_intercept_10 <- compute_SBC(
  datasets_logistic_merged_intercept_10,
  backend_logistic_merged_intercept,
  dquants = logistic_merged_intercept_loglik_dq,
  cache_mode = "results",
  cache_location = file.path(cache_dir, "logistic_merged_intercept_10"))

#' The results for 10 simulations are:
#| label: sbcworkflow_logistic_merged_intercept_results
logistic_merged_intercept_ranks <-  plot_rank_hist(results_logistic_merged_intercept_10)
logistic_merged_intercept_ecdf <- plot_ecdf_diff(results_logistic_merged_intercept_10)
logistic_merged_intercept_ranks / logistic_merged_intercept_ecdf

#' This signals a potential problem with `beta[1]` (the
#' intercept). The reason for the failure is that we have included two
#' separate statements for prior for \texttt{beta[1]}.
#'
#' This example also shows that the `log_lik` term is not magic, as it
#' does not signal this failure earlier than `beta[1]`.
#'
#' ## Fixing prior definition
#'
#' To avoid declaring two priors for `beta[1]` we need to modify the
#' last line of the \texttt{model} block to
#' ```
#' target += normal_lpdf(beta[2:N_predictors] | 0, 1);    
#' ```
#' so the full model now is:
#| label: model_logistic_fixed_prior
cat(readLines(root("sbc", "models/logistic_fixed_prior.stan")), sep = "\n")
model_logistic_fixed_prior <- cmdstan_model(root("sbc", "models/logistic_fixed_prior.stan"))
backend_logistic_fixed_prior <- SBC_backend_cmdstan_sample(model_logistic_fixed_prior, chains = 2) 

#' The results for ten simulations are:
results_logistic_fixed_prior_10 <- compute_SBC(
  datasets_logistic_merged_intercept_10,
  backend_logistic_fixed_prior,
  dquants = logistic_merged_intercept_loglik_dq,
  cache_mode = "results",
  cache_location = file.path(cache_dir, "logistic_fixed_prior_10"))

#| label: fig-sbcworkflow_logistic_fixed_prior_results_10
logistic_fixed_prior_ranks <- plot_rank_hist(results_logistic_fixed_prior_10)
logistic_fixed_prior_ecdf <- plot_ecdf_diff(results_logistic_fixed_prior_10)
logistic_fixed_prior_ranks / logistic_fixed_prior_ecdf

#' No obvious problem, so let's add 200 additional simulations:
set.seed(32464655)
datasets_logistic_merged_intercept_200 <- generate_datasets(generator_logistic_merged_intercept, 200)
results_logistic_fixed_prior_200 <- compute_SBC(
  datasets_logistic_merged_intercept_200,
  backend_logistic_fixed_prior,
  keep_fits = FALSE,
  dquants = logistic_merged_intercept_loglik_dq,
  cache_mode = "results",
  cache_location = file.path(cache_dir, "logistic_fixed_prior_200"))

#| label: fig-sbcworkflow_logistic_fixed_prior_results_200
logistic_fixed_prior_ranks <- plot_rank_hist(results_logistic_fixed_prior_200)
logistic_fixed_prior_ecdf <- plot_ecdf_diff(results_logistic_fixed_prior_200)
logistic_fixed_prior_ranks / logistic_fixed_prior_ecdf

#' Looking good! We can also check that we are indeed able to learn
#' the model parameters with reasonable precision from the data.
#| label: fig-sbcworkflow_logistic_fixed_prior_sim_estimated
plot_sim_estimated(results_logistic_fixed_prior_200)

#' # Full model
#'
#' We are finally ready to make a first attempt at the full model:
#| label: model_combined_first
cat(readLines(root("sbc", "models/combined_first.stan")), sep = "\n")
model_combined <- cmdstan_model(root("sbc", "models/combined_first.stan"))
backend_combined <- SBC_backend_cmdstan_sample(model_combined)

#' And this is our generator for the full model:
generator_func_combined <- function(N_obs, N_predictors) {
  # If the priors for all components of an ordered vector are the same
  # then just sorting the result of a generator is enough to create
  # a valid draw from the ordered vector prior
  mu <- sort(rnorm(2, 3, 1)) 
  beta <- c(rnorm(1, 0, 2), rnorm(N_predictors - 1, 0, 1))
  X <- matrix(rnorm(N_predictors * N_obs, 0, 1), nrow = N_obs, ncol = N_predictors)
  X[, 1] <- 1 # Intercept
  y <- array(NA_real_, N_obs)
  for (n in 1:N_obs) {
    linpred <- 0
    for (p in 1:N_predictors) {
      linpred <- linpred + X[n, p] * beta[p]
    }
    theta <- plogis(linpred)
    if (runif(1) < theta) {
      y[n] <- rpois(1, exp(mu[1]))
    } else {
      y[n] <- rpois(1, exp(mu[2]))
    }
  }
  list(variables = list(beta = beta,
                        mu = mu),
    generated = list(N_obs = N_obs,
                     N_predictors = N_predictors,
                     y = y,
                     X = X))
}
generator_combined <- SBC_generator_function(generator_func_combined, N_obs = 50, N_predictors = 3)

#' We are confident (and the fits are fast anyway), so we start with 200 simulations:
set.seed(5749955)
dataset_combined <- generate_datasets(generator_combined, 200)

results_combined <- compute_SBC(dataset_combined,
                                backend_combined, 
                                keep_fits = FALSE,
                                cache_mode = "results", 
                                cache_location = file.path(cache_dir, "combined"))

#' We get some amount of divergent transitions, but the ranks look pretty good:
#| label: fig-sbcworkflow_combined_results
combined_ranks <- plot_rank_hist(results_combined)
combined_ecdf <- plot_ecdf_diff(results_combined)
combined_ranks / combined_ecdf

#' Indeed it seems the model works pretty well.
#'
#' ## Adding rejection sampling
#'
#' As done previously, we could just exclude the fits that had
#' divergences, but just to complete our tour of possibilities, we'll
#' show one more option to dealing with this type of problem.
#'
#' The general idea is that although we might not want to/be able to
#' express our prior belief about the model (here that the two mixture
#' components are distinct) by priors on model parameters, we still
#' may be able to express our prior belief about the data itself.
#'
#' And it turns out that if we remove simulations that don't meet a
#' certain condition imposed on the observed data, the implied prior
#' on parameters becomes an additive constant and we can use exactly
#' the same model to fit only the non-rejected simulations. Note that
#' this does not hold if we rejected simulations based on some
#' unobserved variables - for more details see the
#' [`rejection_sampling`](https://hyunjimoon.github.io/SBC/articles/rejection_sampling.html)
#' vignette.
#'
#' The main advantage is that if we can do this, we can avoid wasting
#' computation on fitting data that would likely produce divergences
#' anyway. The downside is that it means we no longer have a guarantee
#' the model works for non-rejected data, so we need to check if the
#' data we want to analyze would not be rejected by our criterion.
#'
#' How to build such a criterion here? We'll note that for
#' Poisson-distributed variables the ratio of mean to variance (a.k.a
#' the Fano factor) is always 1. So if the components are too similar,
#' the data should resemble a Poisson distribution and have Fano
#' factor of 1, while if the components are distinct the Fano factor
#' will be larger.
#'
#' All the divergence are for low fano factors - this is the histogram
#' of Fano factor for diverging fits:
#| label: fig-sbcworkflow_fanos
fanos <- vapply(dataset_combined$generated, 
                function(dataset) { var(dataset$y) / mean(dataset$y) }, 
                FUN.VALUE = 0)
fanos_df <- data.frame(fano = fanos, 
                       type = ifelse(results_combined$backend_diagnostics$n_divergent > 0,
                                      "Has divergent transitions", "No divergent transitions"))
fano_threshold <- 1.8
fanos_plot <-  ggplot(fanos_df, aes(x = fano)) + 
  geom_histogram(binwidth = 0.1) + 
  geom_vline(xintercept = fano_threshold, color = "blue", size = 2) +
  scale_x_log10("Fano factor") + facet_wrap(~type, ncol = 1)
fanos_plot

#' So what we'll do is that we'll reject any simulation where the
#' observed data have Fano factor < `r fano_threshold`. In practice a
#' simple way to implement this is to wrap our generator code in a
#' loop and break from the loop only when the generated data meet our
#' criteria (i.e. is not rejected). This is our code:
generator_func_combined_reject <- function(N_obs, N_predictors) {
  if (N_obs < 5) {
    stop("Too low N_obs for this simulator")
  }
  repeat {
    # If the priors for all components of an ordered vector are the same
    # then just sorting the result of a generator is enough to create
    # a valid draw from the ordered vector prior
    mu <- sort(rnorm(2, 3, 1)) 
    beta <- c(rnorm(1, 0, 2), rnorm(N_predictors - 1, 0, 1))
    X <- matrix(rnorm(N_predictors * N_obs, 0, 1), nrow = N_obs, ncol = N_predictors)
    X[, 1] <- 1 # Intercept
    y <- array(NA_real_, N_obs)
    for (n in 1:N_obs) {
      linpred <- 0
      for (p in 1:N_predictors) {
        linpred <- linpred + X[n, p] * beta[p]
      }
      theta <- plogis(linpred)
      if (runif(1) < theta) {
        y[n] <- rpois(1, exp(mu[1]))
      } else {
        y[n] <- rpois(1, exp(mu[2]))
      }
      
    }
    if (var(y) / mean(y) > fano_threshold) {
      break;
    }
  }
  list(variables = list(beta = beta,
                        mu = mu),
       generated = list(N_obs = N_obs,
                        N_predictors = N_predictors,
                        y = y,
                        X = X))
}
generator_combined_reject <- 
  SBC_generator_function(generator_func_combined_reject, N_obs = 50, N_predictors = 3)

#' We'll once again fit our model to 200 simulations:
set.seed(44685226)
dataset_combined_reject <- generate_datasets(generator_combined_reject, 200)

results_combined_reject <- compute_SBC(dataset_combined_reject,
                                       backend_combined, 
                                       keep_fits = FALSE,
                                       cache_mode = "results", 
                                       cache_location = file.path(cache_dir, "combined_reject"))

#' No more divergences! And the ranks look nice.
#| label: fig-sbcworkflow_combined_reject_results
combined_ranks <- plot_rank_hist(results_combined_reject)
combined_ecdf <- plot_ecdf_diff(results_combined_reject)
combined_ranks / combined_ecdf

#' And our coverage is pretty tight:
#| label: fig-sbcworkflow_combined_reject_coverage
plot_coverage(results_combined_reject)

#' Below we show the uncertainty for two variables and some widths of
#' central posterior intervals numerically:
stats_subset <- results_combined_reject$stats[
  results_combined_reject$stats$variable %in% c("beta[1]", "mu[1]"), ]
empirical_coverage(stats_subset, c(0.25, 0.5, 0.9, 0.95))

#' Maybe we think the remaining uncertainty is too big, so we'll run
#' 300 more simulations, just to be sure:
set.seed(1395367854)
dataset_combined_reject_more <- generate_datasets(generator_combined_reject, 300) 
results_combined_reject_more <- bind_results(
  results_combined_reject,
  compute_SBC(dataset_combined_reject_more,
              backend_combined, 
              keep_fits = FALSE,
              cache_mode = "results", 
              cache_location = file.path(cache_dir, "combined_reject_more")))

#' We get some very small number of problematic fits, which we can
#' ignore in this volume (but probably more aggressive rejection
#' sampling would remove those as well).
#'
#' Our plots and coverage are now pretty decent:
#| label: fig-sbcworkflow_combined_reject_more_results
combined_ranks <- plot_rank_hist(results_combined_reject_more)
combined_ecdf <- plot_ecdf_diff(results_combined_reject_more)
combined_ranks / combined_ecdf

#| label: fig-sbcworkflow_combined_reject_coverage_more
combined_reject_coverage <- plot_coverage_diff(results_combined_reject_more)
combined_reject_coverage

stats_subset <- results_combined_reject_more$stats[
  results_combined_reject_more$stats$variable %in% c("beta[1]", "mu[2]"), ]
empirical_coverage(stats_subset, c(0.25, 0.5, 0.9, 0.95))

#' This actually shows a limitation of the coverage results - for
#' `mu[1]` the approximate CI for coverage excludes exact calibration
#' for a bunch of intervals, but above we see that the more
#' trustworthy `plot_ecdf_diff` is not showing a problem (although
#' there is some tendency towards slight underdispersion).
#'
#' Still, this might warrant further investigation if small
#' discrepancies in `mu` are considered important, if we are
#' interested only in the `beta` coefficients, we can stay assured
#' that their calibration is pretty good. We give you our word that we
#' ran additional simulations and the discrepancy disappears.
#'
#' Finally, we can also use this simulation exercise to understand
#' what would we be likely to learn from an experiment matching the
#' simulations (50 observations, 3 predictors) and plot the true
#' values (simulated by the generator) against estimated mean + 90%
#' posterior credible interval:
#| label: fig-sbcworkflow_sim_estimated_final
sim_estimated_final <- plot_sim_estimated(results_combined_reject_more, alpha = 0.2)
sim_estimated_final

#' We see that we get very precise information about `mu` and a decent
#' picture about all `beta` elements, but the remaining uncertainty is
#' large. We could for example compute the probability that the
#' posterior 90% interval for `beta[2]` excludes zero, i.e. that we
#' learn something about the sign of the association with a continuous
#' predictor:
stats_beta2 <- 
  results_combined_reject_more$stats[
    results_combined_reject_more$stats$variable == "beta[2]",]
mean(sign(stats_beta2$q5) == sign(stats_beta2$q95))

#' Turns out the probability is only around 50%. Depending on your
#' aims, this might be a reason to plan for a larger sample size!
#'
#' # Take home message
#'
#' There are couple lessons I hope this exercise showed: First,
#' building models you can trust is hard work and it is very easy to
#' make mistakes. Despite the models presented here being relatively
#' simple, diagnosing the problems in them was not straightforward and
#' required non-trivial background knowledge. For this reason, moving
#' in small steps during model development is crucial and can save you
#' time as diagnosing the same problems in a 300-line Stan model with
#' 50 parameters can be basically impossible.
#'
#' We also hope we convinced you that the SBC package lets you get
#' high-quality information from your simulation efforts and not only
#' diagnose problems but also get some sort of assurance in the end
#' that your model is at least pretty close to your simulator.
#'
#' 
#' ## References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' ## Original Computing Environment {.unnumbered}
#' 
sessionInfo()

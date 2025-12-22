#' ---
#' title: "Roaches cross-validation case study"
#' author: "Aki Vehtari"
#' date: 2017-01-10
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
#' # Setup  {.unnumbered}
#' 
#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=FALSE, comment=NA, out.width='95%')

#' 
#' **Load packages**
#| cache: FALSE
library(loo)
library(rstantools)
library(brms)
library(cmdstanr)
options(mc.cores = 4)
library(ggplot2)
library(ggdist)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size=16))
library(posterior)
options(posterior.num_args=list(digits=2))
library(priorsense)
library(dplyr)
library(tibble)
library(reliabilitydiag)

#' 
#' # Introduction
#' 
#' This case study demonstrates cross-validation model comparison, and
#' posterior and cross-validation predictive checking of
#' models. Furthermore the notebook demonstrates how to use integrated
#' PSIS-LOO with varying intercept ("random effect") models.
#' 
#' The roaches data example comes from Chapter 8.3 of [Gelman and Hill
#' (2007)](http://www.stat.columbia.edu/~gelman/arm/).
#' 
#' > the treatment and control were applied to 160 and 104 apartments,
#' respectively, and the outcome measurement $y_i$ in each apartment
#' $i$ was the number of roaches caught in a set of traps. Different
#' apartments had traps for different numbers of days
#' 
#' In addition to an intercept, the regression predictors for the
#' model are the pre-treatment number of roaches `roach1`, the
#' treatment indicator `treatment`, and a variable indicating whether
#' the apartment is in a building restricted to elderly residents
#' `senior`. The distribution of `roach1` is very skewed and we take a
#' square root of it. Because the number of days for which the roach
#' traps were used is not the same for all apartments in the sample,
#' we include it as an `exposure2` by adding $\ln(u_i)$) to the linear
#' predictor $\eta_i$ and it can be specified using the `offset`
#' argument in `brms`.
#' 
#' Load data
data(roaches, package = "rstanarm")
roaches$sqrt_roach1 <- sqrt(roaches$roach1)
head(roaches) |> print(digits=2)

#' # Poisson model
#' 
#' Fit a Poisson regression model with `brms`
#| results: hide
#| cache: true
fit_p <- brm(y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)),
                data = roaches,
                family = poisson,
                prior = prior(normal(0,1), class = b),
                refresh = 0)

#' Plot posterior
#| label: fig-posterior-poisson
#| fig-height: 3
#| fig-width: 6
mcmc_areas(fit_p, prob_outer = .999,
           regex_pars = c("sqrt_roach1", "treatment", "senior")) +
  coord_cartesian(xlim = c(-0.65, 0.25))

#' 
#' All marginal posteriors are clearly away from zero. We need to do
#' some model checking before trusting these.
#' 
#' ## Posterior predictive checking
#' 
#' Posterior predictive checking can often detect problems and also
#' provide more information about the reason. As the range of counts
#' is large, we can use kernel density estimate plot
#' [@Sailynoja-Johnson-Martin-etal:2025].
#| label: fig-ppc_dens-poisson
#| fig-height: 4
#| fig-width: 7
pp_check(fit_p, type = "dens_overlay", ndraws = 20) +
  scale_x_sqrt(breaks = c(0, 1, 3, 10, 30, 100, 300),
               lim = c(0, 400))

#' We see that the marginal distribution of model replicated data is
#' clearly different from the observed data which are more
#' dispersed. Posterior predictive checking can sometimes have
#' significant bias due to double use of the data, but in this case
#' the discrepancy is so big that further checks are not needed.
#' 
#' Although in this case the model misspecification is obvious with
#' kernel density plot, too, for count the often a better choice is a
#' rootogram variant proposed by @Sailynoja-Johnson-Martin-etal:2025.
#| label: fig-ppc_rootogram-poisson
#| fig-height: 4
#| fig-width: 7
#| out-width: 100%
pp_check(fit_p, type = "rootogram", style = "discrete") +
  scale_x_sqrt() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))

#'
#' ## Cross-validation model comparison
#'
#' For demonstration, we show what would happen in cross validation
#' based model comparison if we would ignore posterior predictive
#' checking result.
#' 
#' We use Pareto-smoothed importance sampling leave-one-out (PSIS-LOO)
#' cross-validation [@Vehtari+etal:PSIS-LOO:2017] as it is very fast to compute.
fit_p <- add_criterion(fit_p, criterion="loo")
# store the result for later use
loo_p1 <- loo(fit_p)

#' We get warning about high Pareto-$\hat{k}$'s in PSIS-LOO
#' computation and recompute using moment matching
#' [@Paananen+etal:2021:implicit] which translates and scales the
#' posterior draws to better match the moments of the leave-one-out
#' posteriors (done only for the folds with high Pareto-$\hat{k}$)
fit_p <- add_criterion(
  fit_p,
  criterion = "loo",
  moment_match = TRUE,
  overwrite = TRUE
)

#' Now all Pareto-$\hat{k}$'s are ok, and we examine LOO results.
loo(fit_p)

#' `p_loo` is about 278, which is much higher than the number of
#' parameters (4), which indicates bad misspecification, which we did
#' already see also with posterior predictive checking. Many high
#' Pareto-\(\hat{k}\)'s in PSIS-LOO without moment matching were
#' likely also caused by model misspecification.
#' 
#' For demonstration, we show what would happen if we would try to use cross validation
#' based model comparison to assess predictive relevance of the covariates. We later
#' compare these to results when using better models.
#'
#' We form 3 models by dropping each of the covariates out at a time. 
#| results: hide
#| cache: true
fit_p_m1 <- update(fit_p, formula = y ~ treatment + senior) |>
  add_criterion(criterion = "loo", moment_match = TRUE)
fit_p_m2 <- update(fit_p, formula = y ~ sqrt_roach1 + senior)  |>
  add_criterion(criterion = "loo", moment_match = TRUE)
fit_p_m3 <- update(fit_p, formula = y ~ sqrt_roach1 + treatment) |>
  add_criterion(criterion = "loo", moment_match = TRUE)

#' 
#' Moment matching is able to assist PSIS-LOO computation and all
#' Pareto $\hat{k}$ values are good.
loo_compare(fit_p,
            fit_p_m1,
            fit_p_m2,
            fit_p_m3,
            model_names=c("Poisson full model",
                          "Poisson w/o sqrt(roach1)",
                          "Poisson w/o treatment",
                          "Poisson w/o senior"))

#' Based on `elpd_diff` and `se_diff` the roaches covariate would be
#' relevant, but although dropping treatment or senior covariate will
#' make a large change to elpd, the uncertainty is also very large and
#' cross-validation indicates that these covariates are not
#' necessarily relevant. The posterior marginals are conditional on
#' the model, but cross-validation is more cautious by not using any
#' model for the future data distribution. The column `p_worse`
#' provides a normal approximation base probability that the model is
#' wrose than the model with the best performance. The column
#' `diag_diff` indicates that the distributions of the pointwise
#' performance differences are likely to have so thick tails that the
#' normal approximation based on `elpd_diff` and `diff_se` is not well
#' calibrated, and we can't trust the `p_worse` values.  In this case,
#' the thick tails in the pointwise performance differences are likely
#' due to Poisson model being underdispersed compared to the data.
#' 
#' 
#' # Negative binomial model
#' 
#' We change the Poisson model to a more robust negative binomial
#' model. Often it would be sensible to start with negative binomial
#' model for counts and skip the Poisson model.  The negative-binomial
#' shape parameter has the `brms` default prior, which is
#' inverse-gamma$(.4, .3)$ [@Vehtari:2024].
#' 
#| results: hide
#| cache: true
fit_nb <- update(fit_p, family = negbinomial)

#' 
#' Plot posterior
#| label: fig-posterior-nb
#| fig-height: 3
#| fig-width: 6
mcmc_areas(fit_nb, prob_outer = .999,
           regex_pars = c("sqrt_roach1", "treatment", "senior"))

#' 
#' Treatment effect is much closer to zero, and senior effect has lot
#' of probability mass on both sides of 0. So it matters, which model
#' we use, and we should trust posteriors only if the model passes
#' predictive checking.
#' 
#' 
#' ## Posterior and LOO predictive checking
#' 
#' We use posterior predictive checking to compare marginal data and
#' predictive distributions. 
#| label: fig-ppc_dens-nb
#| fig-height: 4 
#| fig-width: 7
pp_check(fit_nb, type = "dens_overlay", ndraws = 20) +
  scale_x_sqrt(breaks = c(0, 1, 3, 10, 30, 100, 300),
               lim = c(0, 400))

#' We see that the negative-binomial model is much better although it
#' seems that the model predictive distribution has more mass for
#' small counts than the real data. This discrepancy can also be an
#' artefact from the kernel density estimate, and it is better to
#' examine discrete rootogram
#| label: fig-ppc_rootogram-nb
#| fig-height: 4
#| fig-width: 7
#| out-width: 100%
pp_check(fit_nb, type = "rootogram", style = "discrete") +
  scale_x_sqrt() +
  coord_cartesian(xlim = c(0, 400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))

#' We see that the negative-binomial model is quite good, and there is
#' no clear discrepancy based on the predictive intervals and data.
#'
#' Instead of looking at the marginal distribution, we can look at the
#' pointwise predictive distribution and corresponding observations.
#| label: fig-ppc_intervals-nb
#| fig-height: 4
#| fig-width: 6
pp_check(fit_nb, "intervals") +
  scale_y_sqrt() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.9))

#'
#' If the model predictive distributions are well calibrated, then the
#' observations should look like randomly drawn from the predictive
#' distributions. We can use a `pit_ecdf` to make probability integral
#' transformation (PIT) plot to make the comparison of pointwise
#' posterior predictive distributions and data. One PIT value is the
#' cumulative density of one observation given the corresponding
#' conditional predictive distribution. If the pointwise predictive
#' distributions are well calibrated, PIT values are (almost)
#' uniformly distributed. PIT-ECDF plot compares observed PIT values
#' to uniform distribution.
#| label: fig-ppc_pit-nb
#| fig-height: 4
#| fig-width: 5
pp_check(fit_nb, type = "pit_ecdf")

#' Now that posterior predictive check looks quite good, it is useful
#' to be more careful and use LOO predictive checking, too. We first
#' run PSIS-LOO computation to check it works.
fit_nb <- add_criterion(fit_nb, criterion = "loo")
# store the result for later use
loo_nb1 <- loo(fit_nb)

#' We get warning about high Pareto-$\hat{k}$ and to improve
#' computation accuracy, we re-run with moment matching. We later need
#' some Pareto smoothed importance sampling intermediate computation
#' results and save them, too.
fit_nb <- add_criterion(
  fit_nb,
  criterion = "loo",
  moment_match = TRUE,
  save_psis = TRUE,
  overwrite = TRUE
)

#' Let's look at the LOO results.
(loo_nb <- loo(fit_nb))

#' All Pareto-$\hat{k}$ are good indicating PSIS-LOO computation with
#' moment matching works well.  `p_loo` is closer to the actual number
#' of parameters, but still slightly larger than the total number of
#' parameters 5 which is slightly suspicious. `p_loo` is small
#' compared to the number of observations and we may expect pointwise
#' LOO predictive intervals to be similar to pointwise posterior
#' predictive intervals, and LOO-PIT-ECDF to look similar to posterior
#' PIT-ECDF. Although we had used moment matching above, some
#' intermediate results are not stored as they may take a lot of
#' memory, and we use moment matching argument again when computing
#' LOO intervals and LOO-PIT's.  
#| label: fig-ppc_loo_intervals-nb 
#| fig-height: 4
#| fig-width: 5
pp_check(fit_nb, "loo_intervals", moment_match = TRUE) +
  scale_y_sqrt() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.9))

#| label: fig-ppc_loopit-nb 
#| fig-height: 4
#| fig-width: 5
pp_check(fit_nb, type = "loo_pit_ecdf", moment_match = TRUE)
#'
#' As we guessed, there is not much difference between LOO and
#' posterior inervals, and LOO-PIT-ECDF and PIT-ECDF.
#'
loo_compare(fit_p,
            fit_nb,
            model_names=c("Poisson",
                          "Neg-bin"))

#' In LOO model comparison the negative binomial model is much better
#' than the Poisson model. The diagnostic `diag_diff` indicates the
#' pointwise performance differences are likely to have thick tails,
#' which is again likely due to the Poisson model being
#' underdispersed. As we have made the model checking, we know that
#' the Poisson model is misspecified and could be dropped just based
#' on that fact. As the difference is very big and many times bigger
#' than `diff_se`, we can also be certain that negative binomiasl
#' model has better predictive performance.
#'
#' The next comparison shows the result without moment matching
#' improved computation.
loo_compare(list(`Poisson`=loo_p1,
                 `Neg-bin`= loo_nb1))

#' The difference is now a bit smaller, but still the difference is
#' huge. The column `diag_elpd` is reminding us that there were high
#' Pareto k values in PSIS-LOO computation. This illustrates that we
#' don't always need to improve PSIS-LOO computation to get rid of all
#' high Pareto-$\hat{k}$'s, if the difference between the models is
#' clearly bigger than possible bias.
#'
#' 
#' As Poisson is a special case of negative-Binomial, instead of using
#' cross validation model comparison, we could have also seen that
#' Poisson is not a good model by looking at the posterior of the
#' over-dispersion parameter (which gets very small values), and there
#' would not have been need to fit Poisson model at all.
#| label: fig-posterior-nb_dispersion
#| fig-height: 2
#| fig-width: 6
mcmc_areas(fit_nb, prob_outer = .999, pars = "shape")

#' Posterior predictive rootogram and LOO-PIT-ECDF looked good, but
#' for discrete models where some discrete values have relatively high
#' probability it may miss things and it is good to examine more
#' carefully how well calibrated predictive probabilities are. In
#' roaches, 64% of observations are 0, and it makes sense to look how
#' well we predict zeros or non-zeros. To check calibration of binary
#' target predictive probabilities, we use PAV-adjusted reliability
#' diagram [@Dimitriadis-Gneiting-Jordan:2021,@Sailynoja-Johnson-Martin-etal:2025].
#' Although in this case, the difference is small, for demonstration
#' we use LOO predictive probabilities instead of posterior predictive
#' probabilities for non-zeros. In case of good calibration, the
#' calibration line would stay mostly inside of the envelope.
#| label: fig-pava-nb
#| fig-height: 4
#| fig-width: 5
rd <- reliabilitydiag(
  EMOS = E_loo((posterior_predict(fit_nb) > 0) + 0, loo_nb$psis_object)$value, 
  y = as.numeric(roaches$y > 0)
)
autoplot(rd) +
  labs(x = "Predicted probability of non-zero", 
       y = "Conditional event probabilities") +
  bayesplot::theme_default(base_family = "sans", base_size = 16)

#' There is a slight miscalibration indicated by red curve being quite
#' much outside of the blue envelope. We later build a zero-inflated
#' model to test whether that would improve calibration and predictive
#' performance.
#'
#' # Poisson model with varying intercepts
#' 
#' Sometimes overdispersion is modelled by adding varying intercepts
#' ("random effects") for each individual. Negative-binomial model can
#' be considered as mixture of Poissons with mixing distribution being
#' a gamma distribution. It is common to use normal prior for varying
#' intercepts in log intensity scale. This kind of varying intercept
#' model may fit the data better than negative binomial model for two
#' reasons. First, it is possible that the normal variation in log
#' intensity scale is closer to actual variation. Second, explicitly
#' presenting the latent values and using a prior makes the model
#' slightly more flexible as the posterior of the varying intercepts
#' can be different from the prior. However, in this case with just
#' one observation per varying intercept, it is likely that the
#' posterior is close to normal as likelihood contribution from one
#' observation for each varying intercept is weak.
#' The following example demonstrates computational challenges with
#' varying intercept approach.
#'
#' We add varying intercept for each observation with normal prior
#' using the formula term `(1 | id)`.  We use run the sampling four times
#' longer to get high enough effective sample sizes.
#' 
#| results: hide
#| cache: true
roaches$id <- 1:dim(roaches)[1]
fit_pvi <- brm(y ~ sqrt_roach1 + treatment + senior + (1 | id) + offset(log(exposure2)),
                  data = roaches, family = poisson, 
                  prior = prior(normal(0, 1), class = b),
                  warmup = 1000, iter = 5000, thin = 4,
                  refresh = 0)

#' 
#' ## Analyse posterior
#' 
#' Plot posterior
#| label: fig-posterior-poisson_varying
#| fig-height: 3
#| fig-width: 6
mcmc_areas(as.matrix(fit_pvi), prob_outer = .999,
    regex_pars = c("sqrt_roach1", "treatment", "senior"))

#' 
#' The marginals are similar as with negative-binomial model, but
#' slightly closer to 0.
#' 
#' ## Cross-validation checking
#' 
#' Compute LOO using PSIS-LOO.
fit_pvi <- add_criterion(fit_pvi, criterion="loo")
loo(fit_pvi)

#' `p_loo` is about 164, which is less than the number of parameters
#' 267, but it is relatively large compared to to the number of
#' observations (`p_loo >>N/5`), which indicates very flexible
#' model. In this case, this is due to having an intercept parameter for
#' each observation. Removing one observation changes the posterior
#' for that intercept so much that importance sampling fails (even
#' with Pareto smoothing).  Also the very large number of high $k$
#' values is probaly due to having very flexible model.
#' We can try to improve computation with moment matching.
#' By default `brm()` does not store the varying coefficients and
#' we need to re-run with `save_pars = save_pars(all = TRUE)` before
#' using moment matching.
#| results: hide
#| cache: true
fit_pvi <- update(fit_pvi, save_pars = save_pars(all = TRUE))
fit_pvi <- add_criterion(fit_pvi, criterion="loo", moment_match = TRUE, overwrite = TRUE)
loo(fit_pvi)

#' Moment matching is able to reduce the number of high
#' Pareto-$\hat{k}$'s from 204 to 46, which is still a lot. Varying
#' coefficient models is the one special case where moment matching is
#' likely to not be able to help enough as the posterior of the group
#' specific parameter is changing too much and moment matching for
#' high dimensional non-normal posterior is not able to help
#' enough. `brms` supports running MCMC for the LOO folds with high
#' Pareto-$\hat{k}$'s with `reloo = TRUE`, but that would in this case
#' require 46 refits.  We can use $K$-fold-CV instead to re-fit the
#' model 10 times, each time leaving out 10% of the observations. This
#' shows that cross-validation itself is not infeasible for varying
#' parameter models.
#| results: hide
#| cache: true
(kcvpvi <- kfold(fit_pvi, K = 10))

#' 
#' loo package allows comparing PSIS-LOO and $K$-fold-CV results
#| cache: false
loo_compare(list(`Neg-bin`=loo(fit_nb),
                 `Poisson var. int.`=kcvpvi))

#' 
#' There is not much difference, and the uncertainty is big. The
#' column `diag_diff` shows that the distribution of the differences
#' have thick tails, which indicates that `elpd_diff`, `diff_se`, and
#' `p_worse` can not be fully trusted, and the conclusion is that we can't
#' know which model has better predictive performance.
#'
#' To verify, that we can compare PSIS-LOO and $K$-fold-CV results, we
#' can run $K$-fold-CV also for negative-binomial model.
#| cache: true
(kcvnb <- kfold(fit_nb, K=10))
#| cache: false
loo_compare(list(`Neg-bin`=kcvnb,
                 `Poisson var. int.`=kcvpvi))

#' 
#' The difference in predictive performance is very similar, and the
#' small change can be partially explained by the higher variability
#' due to thick tailed distribution of the pointwise differences.
#' 
#' Now that we've seen that based on robust $K$-fold-CV there is not
#' much difference between negative-binomial and
#' varying-intercept-Poisson models, we can also check how bad the
#' comparison would have been with PSIS-LOO without using moment
#' matching and having many high Pareto-$\hat{k}$ warnings.
loo_compare(fit_nb, fit_pvi)

#' 
#' If we would have ignored Pareto-$k$ warnings, we would have
#' mistakenly assumed that varying intercept model is much better.
#' The column `diag_elpd` is reminding that there are many high
#' Pareto-$k$ values in the elpd computation.
#'
#' Note that WAIC is (as usual) even worse (see also
#' @Vehtari+etal:PSIS-LOO:2017)
loo_compare(waic(fit_nb), waic(fit_pvi))

#' 
#' ## Posterior predictive checking
#' 
#' We do posterior predictive checking for varying intercept Poisson  model.
#| label: fig-ppc_dens-poisson_varying
#| fig-height: 4
#| fig-width: 7
pp_check(fit_pvi, type = "dens_overlay", ndraws = 20) +
  scale_x_sqrt(breaks = c(0, 1, 3, 10, 30, 100, 300),
               lim = c(0, 400))
#' The match looks perfect, but that can be explained with having one
#' parameter for each observation and kernel density estimate hiding something.
#'
#' Looking at the discrete rootogram looks also fine.
#' 
#| label: fig-ppc_rootogram-pvi
#| fig-height: 4
#| fig-width: 7
#| out-width: 100%
pp_check(fit_pvi, type = "rootogram", style = "discrete") +
  scale_x_sqrt() +
  coord_cartesian(xlim = c(0, 400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))

#' PIT-ECDF plot seems to indicate problems
#| label: fig-ppc_pit-poisson_varying
#| fig-height: 4
#| fig-width: 5
pp_check(fit_pvi, type = "pit_ecdf")

#' There are too many PIT values near 0.5. If we look at the
#' predictive intervals and observations, we see that many the
#' observations are in the middle of the posterior predictive
#' interval, which can be explained by having very flexible model with
#' one parameter for each observation. Posterior predictive checking
#' is likely to fail with flexible models (having big `p_loo`).
#| label: fig-ppc_intervals-poisson_varying
#| fig-height: 4
#| fig-width: 6
pp_check(fit_pvi, "intervals") +
  scale_y_sqrt() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.9))

#' LOO-PIT's can take into account the model flexibility, if the
#' computation works.  In this case LOO-PIT's look slightly better,
#' but still showing problems, but this is because PSIS-LOO fails (as
#' discussed above)
#| label: fig-ppc_loopit-poisson_varying
#| fig-height: 4
#| fig-width: 5
pp_check(fit_pvi, type = "loo_pit_ecdf")

#' # Poisson model with varying intercept and integrated LOO
#' 
#' Removing one observation changes the posterior for that intercept
#' so much that importance sampling fails (even with Pareto
#' smoothing). We can get improved stability by integrating that
#' intercept out with something more accurate than the importance
#' sampling [@Vehtari+etal:2016:LOO_for_GLVM].  If there is only one
#' group or individual specific parameter then we can integrate that
#' out easily with 1D adaptive quadrature. 2 parameters can be handled
#' with nested 1D quadratures, but for more parameters nested
#' quadrature is likely to be too slow and other integration methods
#' are needed.
#'
#' As we can easily integrate out only 1 (or 2) parameters (per group
#' / individual), it's not easy to make a generic approach for `rstanarm`
#' and `brms`, and thus we illustrate the approach with direct Stan model
#' code and `cmdstanr`.
#' 
#' The following Stan model uses individual specific intercept term
#' `z[i]` (with common prior with scale `sigmaz`). In the generated
#' quantities, the usual `log_lik` computation is replaced with
#' integrated approach. We also generate `y_loorep` which is the LOO
#' predictive distribution given other parameters than `z`. This is
#' needed to get the correct LOO predictive distributions when
#' combined with integrated PSIS-LOO.
poisson_vi_int <- "poisson_vi_integrate.stan"
writeLines(readLines(poisson_vi_int))

#' We could also move the integrated likelihood to the model block and
#' not use MCMC to sample the varying intercepts at all. This would
#' make marginal posterior easier for MCMC, but using quadrature
#' integration $N$ times for each leapfrog step in Hamiltonian Monte
#' Carlo sampling will increase the sampling time more than what would
#' be the benefit from the simpler marginal posterior. When we do the
#' integration in the generetad quantities, the quadrature is computed
#' only for each saved iteration making the computation faster.
#' 
#' Next we compile the model, prepare the data, and sample. `integrate_1d_reltol`
#' sets the relative tolerance for the adaptive 1D quadrature function
#' `integrate_1d` (if with other model and data you see messages about
#' error estimate of integral exceeding the given relative tolerance
#' times norm of integral, you will get NaNs and need to increase the
#' relative tolerance).
#'
#' We increase the number of sampling iterations to improve effective
#' sample size needed for better LOO-PIT plot.
#| results: hide
#| cache: false
mod_p_vi <- cmdstan_model(stan_file = poisson_vi_int)
datap <- list(N = dim(roaches)[1],
              P = 3,
              offsett = log(roaches$exposure2),
              X = roaches[,c("sqrt_roach1", "treatment", "senior")],
              y = roaches$y,
              integrate_1d_reltol = 1e-6)
fit_p_vi <- mod_p_vi$sample(
  data = datap,
  refresh = 0,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 8000,
  thin = 8
)

#' 
#' The posterior is similar as above for varying intercept Poisson model,
#' as it should be, as the generated quantities is not affecting
#' the posterior.
#| label: fig-posterior-poisson_var_int
#| fig-height: 3
#| fig-width: 6
mcmc_areas(as_draws_matrix(fit_p_vi$draws(variables = c("beta", "sigmaz"))),
           prob_outer = .999)

#' Now the PSIS-LOO doesn't give warnings and the result is close to K-fold-CV.
(loo_p_vi <- fit_p_vi$loo(save_psis = TRUE))

loo_compare(list(`Poisson var. int. int-LOO` = loo_p_vi,
            `Neg-bin` = loo_nb))

#' Comparing to the negative binomial model there is not much
#' difference, and the difference in `elpd_diff` compared to
#' $K$-fold-CV may be partially explained by high variability due to
#' thick tails of the pointwise elpd differences (there may be yhash
#' warning, which in this case can be ignored).
#'
#' LOO-PIT plot looks good, although a bit more variation, which is
#' probably due to lower effective sample sizes due to more
#' challenging posterior than for negative binomial model.
#| fig-height: 4
#| fig-width: 5
#| label: fig-ppc_loopit-poisson_varying_int
ppc_loo_pit_ecdf(
  y = roaches$y,
  yrep = fit_p_vi$draws(variables = "y_loorep", format = "matrix"),
  psis_object = loo_p_vi$psis_object
)

#' We check the calibration of predictive probabilities for zeros vs
#' non-zeros. The calibration plot looks better than with negative
#' binomial model. The predicted probabilities have wider range and
#' the calibration curve stays better inside the envelope. Not shown
#' here, but the varying intercepts are quite close to normally
#' distributed (as we expected), and thus the difference to the
#' negative binomial would be mostly the different distribution for
#' the individual variation.
#| label: fig-pava-p_vi
#| fig-height: 4
#| fig-width: 5
rd <- reliabilitydiag(
  EMOS = E_loo((fit_p_vi$draws(variables="y_loorep", format = "matrix")>0)+0, 
               loo_p_vi$psis_object)$value,
  y = as.numeric(roaches$y > 0)
)
autoplot(rd) +
  labs(x = "Predicted probability of non-zero", 
       y = "Conditional event probabilities") +
  bayesplot::theme_default(base_family = "sans", base_size = 16)

#' # Zero-inflated negative-binomial model
#' 
#' As the proportion of zeros is quite high in the data and the
#' calibration plot for negative binomial model indicated slight
#' miscalibration in prediction of zeros and non-zeros, it is
#' worthwhile to test also a zero-inflated negative-binomial model,
#' which is a mixture of two models
#' - logistic regression to model the proportion of extra zero counts
#' - negative-binomial model
#' 
#' Fir zero-inflated negative-binomial model.
#| results: hide
#| cache: true
fit_zinb <-
  brm(bf(y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)),
         zi ~ sqrt_roach1 + treatment + senior + offset(log(exposure2))),
      family = zero_inflated_negbinomial(), data = roaches, 
      prior = c(prior(normal(0, 1), class = "b"), 
                prior(normal(0, 1), class = "b", dpar = "zi"), 
                prior(normal(0, 1), class = "Intercept", dpar = "zi")), 
      seed = 1704009, refresh = 1000)

#' Based on PSIS-LOO, zero-inflated negative-binomial is clearly better.
fit_zinb <- add_criterion(fit_zinb, criterion = "loo", save_psis = TRUE)

#' We get warning about one high Pareto-$\hat{k}$'s in PSIS-LOO, which
#' we can fix with moment matching.
fit_zinb <- add_criterion(fit_zinb, criterion = "loo", save_psis = TRUE,
                          moment_match = TRUE, overwrite = TRUE)

#' `p_loo` is close to the total number of parameters in the model
#' (9), which is a good sign.
(loozinb <- loo(fit_zinb))

#' Zero-inflated negative binomial (ZINB) model has clearly better predictive
#' performance than the negative binomial (and all diagnostics are good).
loo_compare(fit_nb, fit_zinb,
            model_names=c("Neg-bin",
                          "ZINB"))

#' Posterior predictive checking looks good, but there is no clear difference
#' when looking at marginal predictive distributions or LOO-PIT.
#| label: fig-ppc_dens-zinb
#| fig-height: 4
#| fig-width: 7
pp_check(fit_zinb, type='dens_overlay', ndraws = 20) +
  scale_x_sqrt(breaks = c(0, 1, 3, 10, 30, 100, 300),
               lim = c(0, 400))

#| label: fig-ppc_rootogram-zinb
#| fig-height: 4
#| fig-width: 7
#| out-width: 100%
pp_check(fit_zinb, type = "rootogram", style = "discrete") +
  scale_x_sqrt() +
  coord_cartesian(xlim=c(0, 400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))

#| label: fig-ppc_pit-zinb
#| fig-height: 4
#| fig-width: 5
pp_check(fit_zinb, type= "pit_ecdf")

#' LOO-PIT-ECDF
#| label: fig-ppc_loopit-zinb
#| fig-height: 4
#| fig-width: 5
pp_check(fit_zinb, type= "loo_pit_ecdf", moment_match = TRUE)

#' Reliability diagram assessing the calibration of predicted
#' probabilities of zero vs non-zero looks clearly better than the one
#' for negative binomial model.
#| label: fig-pava-zinb
#| fig-height: 4
#| fig-width: 5
rd <- reliabilitydiag(
  EMOS = E_loo((posterior_predict(fit_zinb)>0)+0, loozinb$psis_object)$value,
  y = as.numeric(roaches$y > 0)
)
autoplot(rd) +
  labs(x = "Predicted probability of non-zero", 
       y = "Conditional event probabilities") +
  bayesplot::theme_default(base_family = "sans", base_size=16)

#' Although the models are different, with finite data and wide LOO
#' predictive distributions, there is a limit in which differences can
#' be see in LOO-PIT values. Both negative-binomial and zero-inflated
#' negative binomial are close enough the LOO-PIT can't see
#' discrepancy from the data, but elpd_loo and calibration plot were
#' able to to show that zero-inflation component improves the
#' predictive accuracy and calibration.
#' 
#' ## Analyse posterior
#' 
#' Plot posterior
#| label: fig-posterior-zinb
#| fig-height: 4
#| fig-width: 8
mcmc_areas(as.matrix(fit_zinb)[, 3:8], prob_outer = .999)

#' 
#' The posterior marginals for negative-binomial part are similar to marginals
#' in the plain negative-binomial model. The marginal effects for the
#' logistic part have opposite sign as the logistic part is modelling
#' the extra zeros.
#'
#' The treatment effect is now divided between negative-binomial and
#' logistic part. We can use the model to make predictions for the
#' expected number of roaches given treatment and no-treatment to get
#' one marginal posterior to examine.
#'

#' Expectations of posterior predictive distributions given
#' treatment=0 and treatment=1
pred_zinb <- posterior_epred(fit_zinb,
                           newdata = rbind(mutate(roaches, treatment = 0),
                                           mutate(roaches, treatment = 1)))
#' Ratio of expected number of roaches with vs without treatment
#| label: fig-effect-zinb
#| fig-height: 3
#| fig-width: 6
ratio_zinb <- array(rowMeans(pred_zinb[, 263:524] / pred_zinb[, 1:262]), 
                    c(1000, 4, 1)) |>
  as_draws_df() |>
  set_variables(variables= "ratio")
ratio_zinb |>
  ggplot(aes(x = ratio)) +
  stat_dots(quantiles = 100) +
  stat_slab(fill=NA, color="gray") +
  labs(x= "Ratio of roaches with vs without treatment", y = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line.y = element_blank(), strip.text.y = element_blank()) +
  xlim(c(0, 1)) +
  geom_vline(xintercept=1, linetype = "dotted")

#' The treatment clearly reduces the expected number of roaches.
#'
#' Assuming our model is causally sensible, we can trust the posterior
#' more if the model has well calibrated predictive distributions
#' (passes posterior and LOO predictive, and calibration checking).
#' For illustration purposes, we compare the posteriors using Poisson,
#' negative binomial and zero-inflated negative binomial.
#'
pred_p <- posterior_epred(fit_p,
                           newdata = rbind(mutate(roaches, treatment = 0),
                                           mutate(roaches, treatment = 1)))
ratio_p <- array(rowMeans(pred_p[, 263:524] / pred_p[, 1:262]), 
                 c(1000, 4, 1)) |>
  as_draws_df() |>
  set_variables(variables= "ratio")

pred_nb <- posterior_epred(fit_nb,
                           newdata = rbind(mutate(roaches, treatment = 0),
                                           mutate(roaches, treatment = 1)))
ratio_nb <- array(rowMeans(pred_nb[, 263:524] / pred_nb[, 1:262]), 
                  c(1000, 4, 1)) |>
  as_draws_df() |>
  set_variables(variables= "ratio")

#| label: fig-effect-p-nb-zinb
#| fig-height: 3
#| fig-width: 6
clr <- khroma::colour("bright", names = FALSE)(7)
ratio_zinb |>
  ggplot(aes(x = ratio)) +
  stat_slab(data = ratio_p, fill = NA, color = clr[1], alpha = 0.6) +
  stat_slab(data = ratio_nb, fill = NA, color = clr[2], alpha = 0.6) +
  stat_slab(fill = NA, color = clr[3], alpha = 0.6) +
  labs(x = "Ratio of roaches with vs without treatment", y = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line.y = element_blank(), strip.text.y = element_blank()) +
  xlim(c(0, 1)) +
  geom_vline(xintercept = 1, linetype = "dotted") +
  annotate(geom = "text", label = "Poisson", x = 0.58, y = 0.9,
           hjust = 0, color = clr[1], size = 5) +
  annotate(geom = "text", label = "NB", x = 0.33, y = 0.9, 
           hjust = 1, color = clr[2], size = 5) +
  annotate(geom = "text", label = "ZINB", x = 0.44, y = 0.9, 
           hjust = 0, color = clr[3], size = 5)

#' All models show benefit of treatment, but Poisson model is
#' overconfident with too narrow posterior. Posteriors using negative
#' binomial (NB) and zero-inflated negative binomial (ZINB) are quite
#' similar and in this case it would not probably matter for decision
#' making even if the slightly worse negative binomial model would
#' have been used. 
#'
#' 
#' ## Prior sensitivity analysis
#' 
#' Finally we make prior sensitivity analysis by power-scaling both
#' prior and likelihood. As there is posterior dependency between the
#' negative binomial and zero-inflated coefficients, it is not that
#' useful to look at the prior sensitivity for the parameters.  We
#' focus on the actual quantity of interest, that is, the ratio of
#' expected number of roaches with vs without treatment.
#| label: fig-priorsense-zinb
powerscale_sensitivity(fit_zinb, prediction = \(x, ...) ratio_zinb) |>
                         filter(variable== "ratio") |>
                         mutate(across(where(is.double),  ~num(.x, digits = 2)))
#'
#' There is no prior sensitivity and the likelihood is informative.
#'

#' ## Predictive relevance of covariates
#' 
#' Let's finally check cross-validation model comparison to see
#' whether improved model has effect on the predictive performance
#' comparison.
#| results: hide
#| cache: true
fit_zinb_m2 <-
  update(fit_zinb,
         formula = bf(y ~ sqrt_roach1 + senior + offset(log(exposure2)),
                      zi ~ sqrt_roach1 + senior + offset(log(exposure2)))) |>
  add_criterion(criterion = "loo", moment_match = TRUE)
#+
loo_compare(fit_zinb,
            fit_zinb_m2,
            model_names = c("ZINB full model",
                            "ZINB w/o treatment"))

#' Treatment effect improves the predictive performance with 96%
#' probability. However, as the distribution of the pointwise
#' differences seem to have a thick tail, the estimated probability
#' using the normal approximation can be sensitive to a few biggest
#' pointwise differences, which are due to the a few biggest numbers
#' of roaches observed. Interestingly in this case, both models pass
#' model checking, seem to be well specified, the distributions of
#' pointwise elpds do not have thick tails, and we are able to
#' reliably estimate `elpd_loo` for each model. This indicates that
#' `elpd_diff` is likely to be well estimated, but the distribution of
#' pointwise differences having thicker tail make `se_diff` and the
#' normal approximation based `p_worse` unreliable. The conclusion is
#' that the full model has slightly better predictive performance.
#'
#' As the variation in the number of roaches is high in both treatment
#' groups, it is difficult to predict the number of roaches per
#' apartment making the difference in predictive performance also
#' small, even we can clearly see clear difference in the expected
#' number of roaches in the groups.
#' 

#' 
#' <br />
#' 
#' # References {.unnumbered}
#' 
#' <div id="refs"></div>
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2017-2025, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2017-2025, Aki Vehtari, licensed under CC-BY-NC 4.0.
#' 
#' # Original Computing Environment {.unnumbered}
#' 
sessionInfo()

#' 
#' <br />
#' 

#' ---
#' title: "Kilpisjärvi PIT demo for Workflow book"
#' author: "Aki Vehtari"
#' date: 2024-06-27
#' date-modified: today
#' date-format: iso
#' format:
#'   html:
#'     number-sections: true
#'     code-copy: true
#'     code-download: true
#'     code-tools: true
#' bibliography: ../../../casestudies.bib
#' ---
#' 
#' This notebook includes part of the code for the Bayesian Workflow
#' book Section 8.2 about posterior predictive checking.

#+ setup, include=FALSE
knitr::opts_chunk$set(
  cache = FALSE,
  message = FALSE,
  error = FALSE,
  warning = FALSE,
  comment = NA,
  out.width = "95%"
)

#' 
#' **Load packages**
#| cache: FALSE
library(loo)
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 1)
devtools::load_all("~/proj/bayesplot")
## library(bayesplot)
ggplot2::theme_set(bayesplot::theme_default(base_family = "sans", base_size=14))
SEED <- 48927 # set random seed for reproducibility

#' # Linear Gaussian model
#' 
#' Use the Kilpisjärvi summer month temperatures 1952--2022 data from `aaltobda` package. We can read data directly from a URL.
load(url("https://github.com/avehtari/BDA_course_Aalto/raw/master/rpackage/data/kilpisjarvi2022.rda"))
data_lin <- data.frame(year = kilpisjarvi2022$year,
                       temp = kilpisjarvi2022$temp.summer)

#' To analyse whether there has been change in the average summer month
#' temperature we use a linear model with Gaussian model for the
#' unexplained variation. By default brms uses uniform prior for the
#' coefficients.
#' 
#' Formula `temp ~ year` corresponds to model $\mathrm{temp} ~
#' \mathrm{normal}(\alpha + \beta \times \mathrm{temp}, \sigma)$.  The
#' model could also be defined as `temp ~ 1 + year` which explicitly
#' shows the intercept ($\alpha$) part. Using the variable names
#' `brms` uses the model can be written also as
#' $$
#' \mathrm{temp} \sim \mathrm{normal}(\mathrm{b\_Intercept}*1 + \mathrm{b\_year}*\mathrm{year}, \mathrm{sigma})
#' $$

#' In this case, we are happy with the default prior for the
#' intercept. In this specific case, the flat prior on coefficient is
#' also fine, but we add an weakly informative prior just for the
#' illustration. Let's assume we expect the temperature to change less
#' than 1°C in 10 years. With `student_t(3, 0, 0.03)` about 95% prior
#' mass has less than 0.1°C change in year, and with low degrees of
#' freedom (3) we have thick tails making the likelihood dominate in
#' case of prior-data conflict. In real life, we do have much more
#' information about the temperature change, and naturally a
#' hierarchical spatio-temporal model with all temperature measurement
#' locations would be even better.
#| results: hide
fit_lin <- brm(temp ~ year, data = data_lin, family = gaussian(),
               prior = prior(student_t(3, 0, 0.03), class="b"),
               seed = SEED, refresh = 0)

#' Check the summary of the posterior and inference diagnostics.
fit_lin

#' Posterior predictive check with density overlays examines the whole
#' temperature distribution. We generate replicate data using 20 different
#' posterior draws (with argument `ndraws`).
#| label: fig-kilpis_ppc_dens_overlay
#| fig-height: 4
#| fig-width: 4
pp_check(fit_lin, type="dens_overlay", ndraws=20)

#' In a posterior predictive check which compares pointwise posterior
#' predictive intervals and the observations, the proportion of the
#' observations within the shown 50% and 90% intervals should be
#' approximately 50% and 90%.
#| label: fig-kilpis_ppc_intervals
#| fig-height: 4
#| fig-width: 4
pp_check(fit_lin, type="intervals")

#' Instead of counting the proportion of observations only in some
#' interval lengths, we can emaine all cumulative densities. If the
#' predictive distributions match well the data (are well calibrated)
#' then the pointwise cumulative densities have distribution which is
#' close to uniform. The cumulative densities are also known as
#' probability integral transformations (PIT). We can compare
#' empirical CDF of observed PITs to uniform distribution and to
#' associated simulatanous confidence interval shown as an
#' envelope. If the PIT ECDF stays inside the envelope, the predictive
#' distributions are well calibrated.
#| label: fig-kilpis_ppc_pit_ecdf
#| fig-height: 4
#| fig-width: 4
pp_check(fit_lin, type="pit_ecdf", method = "correlated")

#' There can be a lot of white space in an ECDF plot, and sometimes we
#' prefer to plot a difference from the expected ECDF to improve the
#' dynamic range of the plot.
#| label: fig-kilpis_ppc_pit_ecdf_diff
#| fig-height: 4
#| fig-width: 4
pp_check(fit_lin, type="pit_ecdf", plot_diff=TRUE, method = "correlated")

#' In this example, the number of observations is much higher than the
#' number of parameters in the model and there is not much difference
#' in whether we use posterior predictive or cross-validation
#' predictive distributions in the model checking
#' [@Tesso-Vehtari:2026]. In case of more flexible models it is better
#' to use LOO predictive intervals and LOO-PIT [@Tesso-Vehtari:2026].
#| label: fig-kilpis_ppc_loo_intervals
#| fig-height: 4
#| fig-width: 4
pp_check(fit_lin, type="loo_intervals")

#| label: fig-kilpis_ppc_loo_pit_ecdf_diff
#| fig-height: 4
#| fig-width: 4
pp_check(fit_lin, type="loo_pit_ecdf", plot_diff=TRUE, method = "correlated")

#' ---
#' title: "How many iterations to run and how many digits to report"
#' author: "Aki Vehtari"
#' date: 2020-12-05
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
#' This notebook includes the code for Bayesian Workflow book Sections
#' 11.4 *How many parallel chains and how many iterations* and 11.6
#' *How many digits to report based on posterior uncertainty*.
#' 
#' # Introduction
#'
#' When running any iterative algorithm we need to decide when to stop
#' or how many iterations and how many sequences to run. In the
#' context of posterior inference, the answer can be divided in two
#' parts: (1) how many chains and iterations do we need to assess that
#' the algorithm is sampling from something close to the target
#' distribution, and (2) how many draws do we need so that the
#' reported posterior summaries would not change in important ways if
#' the inference would be repeated? This notebook discusses the second
#' part.
#' 
#' Before we can answer how many chains and iterations we need to run,
#' we need to know how many significant digits we want to report. Too
#' often, we see tables filled with numbers like $1.7705$. It's very
#' unlikely that all these digits are accurately estimated, and also
#' very unlikely that the accuracy in magnitude of one in ten thousand
#' would be needed for any practical purpose. Reporting too many
#' digits makes it more difficult to read summary tables. Thus, before
#' considering how many iterations we need it's useful to consider how
#' many digits it is sensible to report.
#' 
#' Without any additional information, we may assume that relative
#' errors smaller than 1\% are typically negligible and thus most of
#' the time reporting at most 2 significant digits would be
#' sufficient. Accuracy to just one significant digit can be
#' sufficient in the early stages of analysis and can be sensible and
#' convenient even in final reporting, as discussed in the example
#' below.
#' 
#' MCMC and other Monte Carlo methods are stochastic and if the
#' inference would be repeated with a different random seed (random
#' number generators in computers are usually producing deterministic
#' pseudo-random sequences), the estimates would vary. The amount of
#' variation reduces when more iterations are used. For example, if
#' the posterior would be close to normal with standard deviation 1,
#' $\mathrm{normal}(\mu, 1)$, then 2000 independent draws from the
#' posterior would provide enough accuracy that the second significant
#' digit of posterior mean would only sometimes vary to one smaller or
#' larger value. On the other hand, for 1 significant digit accuracy,
#' 100 independent draws would be often sufficient, but reliable
#' convergence diagnostics may need more iterations than 100 (see,
#' e.g., @Vehtari-Gelman-Simpson-etal:2021). The above thumb-rules are useful
#' for the necessary number of independent draws for a posterior mean
#' of a distribution which is relatively close to normal and in other
#' cases and, for exmaple, for posterior quantiles more draws may be
#' needed.
#' 
#' MCMC in general doesn't produce independent draws and the effect of
#' dependency affects how many draws are needed to estimate different
#' expectations. As in general, we don't know beforhand how MCMC will
#' perform for a new posterior, and we don't know what is the scale of
#' that posterior beforehand, we need to start with some initial guess
#' of number of iterations to run.
#' 
#' # Summary
#' 
#' Summary of workflow for how many digits to report
#'
#' 1. Run inference with some default number of iterations
#' 2. Check convergence diagnostics for all parameters
#' 3. Check that ESS is big enough for reliable convergence
#'    diagnostics for all quantities of interest
#' 4. Look at the posterior for quantities of interest and decide how
#'    many significant digits is reasonable taking into account the
#'    posterior uncertainty (using SD, MAD, or tail quantiles)
#' 5. Check that MCSE is small enough for the desired accuracy of
#'    reporting the posterior summaries for the quantities of
#'    interest.
#' 
#'  - If the accuracy is not sufficient, report less digits or run
#'    more iterations.
#'  - Halving MCSE requires quadrupling the number of iterations
#'    (if CLT holds).
#'  - Different quantities of interest have different MCSE and may
#'    require different number of iterations for the desired accuracy.
#'  - Some quantities of interest may have posterior distribution with
#'    infinite variance, and then the ESS and MCSE are not defined for
#'    the expectation. In such cases use, for example, median instead
#'    of mean and mean absolute deviation (MAD) instead of standard
#'    deviation. ESS and MCSE for (non-extreme) quantiles can be
#'    derived from the (non-extreme) cumulative probabilities that
#'    always have finite mean and variance.
#' 
#' # Example
#'
#' As an example, we analyse the trend in summer months average
#' temperature 1952--2013 at Kilpisjärvi in northwestern Finnish
#' Lapland. Summer months are June, July, and August, and we analyse
#' the average temperature over these months in each year.
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
library(tidyr)
library(dplyr)
library(cmdstanr)
library(posterior)
options(pillar.neg = FALSE,
        pillar.subtle = FALSE,
        pillar.sigfig = 2)
library(tinytable)
options(tinytable_format_num_fmt = "significant_cell",
        tinytable_format_digits = 2,
        tinytable_tt_digits = 2)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
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

#' ## Kilpisjärvi data and model
#' 
#' Load Kilpisjärvi summer month average temperatures 1952-2013
#' recorded by Finnish Meteorological Institute (CC BY
#' 4.0). Kilpisjärvi is in northwestern Finnish Lapland (about
#' 69°03'N, 20°50'E).
data_kilpis <- read.delim(root("digits/data", "kilpisjarvi-summer-temp.csv"), sep = ";")
data_lin <- list(N = nrow(data_kilpis),
                 x = data_kilpis$year,
                 xpred = 2016,
                 y = data_kilpis[, 5])
#' 
#' Plot the data
ggplot() +
  geom_point(aes(x, y), data = data.frame(data_lin), size = 1) +
  labs(y = "Summer temperature\n at Kilpisjärvi", x = "Year")

#' We use a simple linear model with normal observation model, weakly
#' informative normal prior, and the predictor (time) centered to have
#' 0 mean. We assume a priori that increase or decrease in temperature
#' is equally likely, and that it is unlikely that temperature would
#' change more 10 degrees in 100 years. We also assume that yearly
#' variation in summer temperatures is likely to be less than 3
#' degrees.
#' 
#' The following Stan code centers the covariate to reduce posterior
#' dependency of slope and coefficient parameters. It also makes it
#' easier to define the prior on average temperature in the center of
#' the time range (instead defining prior for temperature at year 0).
code_lin <- root("digits", "linear.stan")
#| output: asis
print_stan_file(code_lin)

#' Prior parameter values for weakly informative priors
data_lin_priors <- c(list(
    pmualpha_c = 10,     # prior mean for average temperature
    psalpha = 10,        # weakly informative
    pmubeta = 0,         # a priori incr. and decr. as likely
    psbeta = 0.1/3,      # avg temp prob does does not incr. more than a degree per 10 years:  setting this to +/-3 sd's
    pssigma = 1),        # setting sd of total variation in summer average temperatures to 1 degree implies that +/- 3 sd's is +/-3 degrees: 
  data_lin)

#' ## Run inference for some number of iterations
#'
#' We run MCMC with Stan's current default settings. All convergence
#' diagnostics pass and all effective sample size are big enough that
#' we can assume that the convergence diagnostics are reliable.
#| label: fit_lin
#| results: hide
mod_lin <- cmdstan_model(stan_file = code_lin)
fit_lin <- mod_lin$sample(data = data_lin_priors, seed = SEED, refresh = 0)

#' ## Run convergence diagnostics
#' 
#' There were no convergence issues reported by sampling. We can also
#' explicitly call CmdStan inference diagnostics:
fit_lin$cmdstan_diagnose()

#' We check $\widehat{R}$ end ESS values, which in this case all look good.
draws <- as_draws_rvars(fit_lin$draws())
summarize_draws(draws) |> tt()

#' Compute posterior draws for the linear fit
draws$mu <- draws$alpha_c + draws$beta * (data_lin$x - mean(data_lin$x))
#' Plot the linear fit with 90% posterior interval
#| label: fig-kilpisjarvi_fit
#| fig-height: 3
#| fig-width: 6
data.frame(x = data_lin$x,
           y = data_lin$y,
           Emu = mean(draws$mu),
           q05 = as.vector(quantile(draws$mu, 0.05)),
           q95 = as.vector(quantile(draws$mu, 0.95))) |>
  ggplot() +
  geom_ribbon(aes(x = x, ymin = q05, ymax = q95), fill = "grey90") +
  geom_line(aes(x = x, y = Emu, )) +
  geom_point(aes(x, y), size = 1) +
  labs(y = "Summer temperature\n at Kilpisjärvi (°C)", x = "Year") +
  guides(linetype = "none")

#' At this point it is sufficient that diagnostics are OK
#' (cmdstan_diagnose says "no problems detected") and effective sample
#' sizes are large enough (>400) that we can assume the diagnostics to
#' be reliable.
#'
#' We could also do posterior predictive checking and residual plots,
#' but now we focus on checking how many MCMC iterations are needed
#' and how many digits to report in posterior summary results.
#' 

#' ## How many digits to report based on posterior uncertainty
#'
#' We want to report posterior summaries for the slope, that is the
#' increase in average summer temperature, and for the probability
#' that the slope is positive. We first consider what is the
#' reasonable reporting accuracy if the posterior inference would be
#' exact (ie ignoring first the Monet Carlo variability).
#' 
#' We start looking at the mean and 90% interval for the slope parameter beta
draws |>
  subset_draws("beta") |>
  summarize_draws(mean, ~quantile(.x, probs = c(0.05, 0.95))) |>
  tt()

#' These values correspond to the temperature increase per year, but
#' to improve readability we switch looking at the temperature
#' increase per 100 years. 
draws <- draws |>
  mutate_variables(beta100 = 100*beta)

#' Let's look at the mean and 90% interval for the expected temperature
#' increase per 100 years.
draws |>
  subset_draws("beta100") |>
  summarize_draws(mean, ~quantile(.x, probs = c(0.05, 0.95))) |>
  tt()

#' The number of digits shown in R varies depending on the print
#' method and options for the specific object. Above we have used
#' `tinytable::tt()` to make pretty tables. In the beginning of the code
#' we did set options
#' ```
#' options(tinytable_format_num_fmt = "significant_cell",
#'         tinytable_format_digits = 2,
#'         tinytable_tt_digits = 2)
#' ```
#' which make each cell to show 2 significant digits, which is most
#' often enough. Without `tt()` the output of `summarize_draws()`
#' is a tibble, which looks like this
draws |>
  subset_draws("beta100") |>
  summarize_draws(mean, ~quantile(.x, probs = c(0.05, 0.95))) |>
  tt()
#' The printing of tibbles is handled by `pillar`, which by default
#' prints 3 significant digits, but we have set in the beginning an
#' option to show 2 significant digits. 
#' ```
#' options(pillar.sigfig = 2)
#' ```
#' 
#' Depending on the print method, we may also see many more digits. For
#' example, the default number of digits for R numeric output is 7.
mean(draws$beta)
quantile(draws$beta, probs = c(0.05, 0.95))
#' We could use `options(digits = 2)` to show only 2 decimal digits,
#' but that is not as useful as showing 2 significant digits.
#' 
#' The more digits are shown, the more of them are likely to be
#' unnecessary clutter distracting the reader from the important
#' message.  Considering the width of the 90% interval, practically
#' meaningful accuracy would be here to report posterior mean as 2.0
#' and the posterior interval as [0.7 , 3.2]. Depending on the
#' context, it might be even better to round more, and report that the
#' increase is estimated to be 1 to 3 degrees per century (81%
#' probability, or rounded to 80% probability), or 0 to 4 degrees per
#' century (99% probability). There is no need to stick to reporting
#' 90% or 95% interval.
#'
#' The number of significant digits needed for reporting can be often
#' determined also with rough posterior approximations that indicate
#' the order of magnitude of the posterior mean and scale, and then
#' more iterations may be needed to get more accurate estimate for
#' those digits.
#'
#' ## How many digits to report based on Monte Carlo standard error
#'
#' Now that we have an estimate for the posterior uncertainty of the
#' slope parameter, we can check whether we have enough many
#' iterations for the desired reporting accuracy. Markov chain Monte
#' Carlo method used by Stan to make the posterior estimates is
#' stochastic. If we repeat the computation with different random
#' number generator seed, we get slightly different estimates.
#| results: hide
estimates <- t(sapply(1:10, function(i) {
  mod_lin$sample(data = data_lin_priors, seed = SEED + i, refresh = 0,
                 show_messages = FALSE)$draws() |>
    mutate_variables(beta100 = 100*beta) |>
    subset_draws("beta100") |>
    summarize_draws(mean, ~quantile(.x, probs = c(0.05, 0.95))) |>
    select(-variable)})) |> unlist() |> matrix(nrow = 10)
colnames(estimates) <- c("mean", "5%", "95%")

as_tibble(estimates) |> tt()

#' We see that for mean the third digit is varying and the rounded
#' value is between 1.9 and 2.0. For the 5% quantile even the first
#' significant digit is sometimes varying and the rounded value would
#' vary between 0.6 and 0.7.  For the 95% quantile the second digit is
#' varying and the rounded values would vary between 3.2 and
#' 3.3. Based on this, it would be OK to report the mean as 1.9 or 2
#' and 90% interval as [0.7, 3.2] as based on the first Monte Carlo
#' estimate. Considering the scale, the minor variation in the last
#' digit is not affecting the interpretation of the
#' results. Alternatively the courser ranges [1, 3] or [0, 4] could be
#' reported as discussed above.
#' 
#' Instead of repeating the estimation many times we can estimate the
#' accuracy of the original sampling by computing Monte Carlo standard
#' error (MCSE) that takes into account the quantity of interest, and the
#' effective sample size of MCMC draws (see, e.g., @Vehtari-Gelman-Simpson-etal:2021).
#'
draws |>
  subset_draws("beta100") |>
  summarize_draws(mcse_mean, ~mcse_quantile(.x, probs = c(0.05, 0.95))) |>
  tt()

#' We show also the ESS values for mean and quantiles, to illustrate
#' that the ESS values can be different for different quantities, but
#' also that the same ESS doesn't lead to the same MCSE, and similar
#' ESS doesn't lead to similar MCSE as MCSE depends also on the
#' quantity. Here, for example, although ESSs are similar for all
#' quantities, 5\% and 95\% quantiles have clearly higher MCSE than
#' the mean.
#' 
draws |>
  subset_draws("beta100") |>
  summarize_draws(ess_mean, ~ess_quantile(.x, probs = c(0.05, 0.95))) |>
  tt()

#' The MCSE for the mean estimate is about 0.01 and for 5% and 95%
#' quantiles about 0.03. If we multiply these by 2, the likely range
#' of variation due to Monte Carlo is $\pm 0.02$ for mean and $\pm
#' 0.07$ for 5% and 95% quantiles. From this we can interpret that
#' it's unlikely there would be variation in the reported estimate for
#' the mean, if it is reported as 2.0. For 5% and 95% quantiles
#' there can be variation in the first decimal digit, but that
#' difference would not be meaningful in most cases. We could run more
#' iterations to get the MCSE for quantiles down to something like
#' 0.01, which would require about 10 times more iterations. Assuming
#' the posterior has finite variance and MCMC is mixing well, we
#' expect that running four times more iterations would halve the
#' MCSEs.
#'
#' As discussed before, it might be anyway better to report the
#' posterior uncertainty interval with less digits, and if we reported
#' either 80% interval $[1, 3]$ or 99% interval $[0, 4]$, we would
#' already have sufficient accuracy for the number of shown
#' digits. These MCSE estimates illustrate also the fact that usually
#' tail quantiles have lower accuracy than the posterior mean.
#'
#' We can also report the probability that the temperature change is
#' positive.
draws |>
  mutate_variables(beta0p = beta100 > 0) |>
  subset_draws("beta0p") |>
  summarize_draws("mean", mcse = mcse_mean) |>
  tt()

#' The probability is simply estimated as a posterior mean of an
#' indicator function and the usual MCSE for mean estimate is used.
#' The MCSE indicates that we have enough MCMC iterations for
#' practically meaningful reporting that the probability that the
#' temperature is increasing is larger than 99%. There is not much
#' practical difference to reporting that the probability is 99.3%
#' and to estimate that third digit accurately would require 64 times
#' more iterations. For this simple problem, sampling that many
#' iterations would not be time consuming, but we might also instead
#' consider to obtain more data to verify that the summer temperature
#' in northern Finland has been increasing since 1952.
#'

#' We additionaly compute probabilities that the temperature increase
#' is more than 1, 2, 3 or 4 degrees, and corresponding MCSEs and ESSs.
draws |>
  subset_draws("beta100") |>
  mutate_variables(beta1p = beta100 > 1,
                   beta2p = beta100 > 2,
                   beta3p = beta100 > 3,
                   beta4p = beta100 > 4) |>
  subset_draws("beta[1-4]p", regex = TRUE) |>
  summarize_draws("mean", mcse = mcse_mean, ESS = ess_mean) |>
  tt()

#' Taking into account MCSEs given the current posterior sample, we
#' can summarise these as p(beta100>1) = 88%--91%, p(beta100>2) =
#' 46%--51%, p(beta100>3) = 7%--10%, p(beta100>4) = 0.2%--1%. To get
#' these probabilities estimated with 2 digit accuracy would again
#' require more iterations (16-300 times more iterations depending on
#' the quantity), but the added iterations would not change the
#' conclusion radically.
#' 

#' # MCSE computation details
#' 
#' The details of how MCSE is estimated for posterior expectations and
#' quantiles are provided by @Vehtari-Gelman-Simpson-etal:2021, and
#' implementations are available, for example in the `posterior` R
#' package used in this notebook and in `ArviZ` Python package.
#'
#' # Rough estimates for how many iterations to run initially
#' 
#' Before running any MCMC, we may obtain some rough estimates of how
#' many iterations we would need for certain relative accuracy for
#' the posterior mean (if MCMC runs well).
#'
#' If we assume that the posterior of parameter $\theta$ has finite
#' mean and variance, and we have many independent Monte Carlo draws
#' $S$, we can use central limit theorem to justify that the uncertainty
#' related to the estimated expectation can be approximated with
#' normal distribution, $\mathrm{normal}(\hat{\theta}, \hat{\sigma}_\theta /
#' \sqrt{S})$, where $\hat{\theta}$ is the estimate and
#' $\hat{\sigma}_\theta$ is the estimated posterior standard deviation
#' of $\theta$.
#'
#' We can see that given any $\hat{\theta}$ and $\hat{\sigma}_\theta$,
#' if $S=100$, the uncertainty is 10\% of the posterior scale which
#' would often be sufficient for initial experiments, and if $S=10000$,
#' the uncertainty is 1\% of the posterior scale which would often be
#' sufficient for two significant digit accuracy.
#'
#' 10\% relative accuracy means that despite the magnitude of
#' $\hat{\sigma}_\theta$ we get accuracy approximately a worth at
#' least one digit. Let's first assume $\hat{\theta} = 1.234$. We can
#' simplify the analysis by noting that if examine the behavior given
#' different $\hat{\sigma}_\theta$ values as, e.g., $10^{-1}$,
#' $10^{0}$, $10^{1}$, and correspondingly present the $\hat{\theta}$
#' as $12.34 \times 10^{-1}$, $1.234 \times 10^{0}$, $0.1234 \times
#' 10^{1}$. We get equivalent, but simpler analysis if we set
#' $\hat{\sigma}_\theta=1$ and analyse the cases with different values
#' $\hat{\theta}$ (ie we change the coordinates to have unit scale
#' with respect to the posterior standard deviation).
#' 
#'   - If $S=100$ (independent draws) and $\hat{\sigma}_\theta=1$,
#'     MCSE is $0.1$ and with 99\% probability that variation in
#'     E[$\theta$] is $\pm 0.3$. Thus, if the estimation would be
#'     repeated the first significant digit would stay the same or
#'     have minor variability to one smaller or one lager digit. Thus
#'     we have approximately one significant digit accuracy.
#'   - If $S=2000$ (independent draws) and $\hat{\sigma}_\theta=1$,
#'     MCSE is $0.02$ and with 99\% probability that variation in
#'     E[$\theta$] is $\pm 0.07$. Thus, if the estimation would be
#'     repeated the first significant digit would stay the same, and
#'     the second significant digit would have have minor variability
#'     to one smaller or one lager digit. Thus we have approximately
#'     two significant digit accuracy. With a larger $S$, there is
#'     less variability in the second significant digit.
#'
#' Dynamic Hamiltonian Monte Carlo in Stan is often so efficient that
#' ESS>S/2. Thus running with the default options 4 chains with 1000
#' iterations after warmup is likely to give near two significant
#' digit accuracy for the posterior mean. The accuracy for 5\% and
#' 95\% quantiles would be between one and two significant digits.
#'
#' The above analysis shows the benefit of interpreting ESS as a scale
#' free diagnostic whether we are likely to have enough iterations. A
#' scale free here means we don't need to compare ESS values to
#' posterior standard deviations or to domain knowledge of the
#' quantity of interest, making it faster to check that we have high
#' enough ESS for many quantities of interest. However, high ESS is
#' not sufficient to guarantee certain accuracy, as MCSE depends also
#' on the quantity of interest and thus in the end it is useful to
#' check MCSEs for the values to be reported. For example, above, the
#' estimate for whether the temperature increase is larger than 4
#' degrees per century has high ESS, but the indicator variable
#' contains less information (than continuous values) and thus much
#' higher ESS would be needed for two significant digit accuracy.
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
#' * Code &copy; 2020--2025, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2020--2025, Aki Vehtari, licensed under CC-BY-NC 4.0.

#' ---
#' title: "No vehicles in the park"
#' author: "Andrew Gelman and Aki Vehtari"
#' date: 2025-09-16
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
#'     css: ../_styles.css
#'     code-copy: true
#'     code-download: true
#'     code-tools: true
#'     embed-resources: true
#'     anchor-sections: true
#'     html-math-method: katex
#' bibliography: ../casestudies.bib
#' ---
#' 
#' This notebook includes the code for Bayesian Workflow book Chapter
#' 29 *Sampling problems with latent variables: No vehicles in the park*.
#'
#' # Introduction
#'
#' It can be hard to sample from a posterior of even a simple
#' multilevel model. We demonstrate with an example that began with a
#' post from [@luu:2024], which pointed to an online quiz from
#' [@Turner:2024] adapted from [@Hart:1958]
#' (see also @Schlag:1999)
#'
#+ setup, include=FALSE
knitr::opts_chunk$set(
  cache = FALSE,
  message = FALSE,
  error = FALSE,
  warning = FALSE,
  comment = NA,
  out.width = '95%'
)

#' 
#' **Load packages**
#| cache: FALSE
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(arm)
library(posterior)
options(posterior.num_args = list(digits = 2), digits = 2, width = 90)
library(lme4)
library(cmdstanr)
options(mc.cores = 4)
# CmdStanR output directory makes Quarto cache to work
dir.create(root("park_rule", "stan_output"))
options(cmdstanr_output_dir = root("park_rule", "stan_output"))
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))
library(ggdist)
library(patchwork)
library(tictoc)
mytoc <- \() {
  toc(func.toc = \(tic, toc, msg) {
    sprintf("%s took %s sec", msg, as.character(signif(toc - tic, 2)))
  })}

print_stan_code <- function(code) {
  if (isTRUE(getOption("knitr.in.progress")) &
        identical(knitr::opts_current$get("results"), "asis")) {
    # In render: emit as-is so Pandoc/Quarto does syntax highlighting
    block <- paste0("```stan", "\n", paste(code, collapse = "\n"), "\n", "```")
    knitr::asis_output(block)
  } else {
    writeLines(code)
  }
}

#' # Data
#' 
#' Data has 51953 responses from 2409 people giving yes/no answers to
#' 27 questions about no vehicles in the park rule. Each row of the
#' data file includes indexes for the respondent and the question,
#' indicators if the person described in the survey question has a
#' male-sounding or white-sounding name, and the outcome: 1 for the
#' response, ``Yes, this violates the rule'' and 0 for the response,
#' ``No, this does not violate the rule.''
#'

park <- read.csv(root("park_rule", "data", "park.csv"))
N <- nrow(park)
respondents <- sort(unique(park$submission_id))
J <- length(respondents)
respondent <- rep(NA, N)
for (j in 1:J) {
  respondent[park$submission_id == respondents[j]] <- j
}
items <- sort(unique(park$question_id))
K <- length(items)
item <- park$question_id
y <- park$answer
n_responses <- rep(NA, J)
for (j in 1:J) {
  n_responses[j] <- sum(respondent == j)
}
male_name <- park$male
white_name <- park$white
#' add number of answered questions as a predictor
n_responses_full <- n_responses[respondent]

#' # lme4 model
#'
#' `lme4` performs fast fitting of simple multilevel models using
#' approximate marginal maximum likelihood for the variance
#' parameters. We start with lme4 for quick analysis to get
#' some insights to the data.

data_park <- data.frame(y,
                        respondent,
                        item,
                        male_name,
                        white_name,
                        n_responses_full)
tic("lme4 fit 1")
fit_lme4 <- glmer(
  y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_responses_full,
  family = binomial(link = "logit"),
  data = data_park
)
mytoc()
display(fit_lme4)

#' ## Refit the model
#'
#' Code the last predictor in terms of the number of items skipped so
#' that zero is a reasonable baseline.
n_skipped <- K - n_responses
n_skipped_full <- n_skipped[respondent]
data_park <- data.frame(y,
                        respondent,
                        item,
                        male_name,
                        white_name,
                        n_skipped_full)

tic("lme4 fit 2")
fit_lme4 <- glmer(
  y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_skipped_full,
  family = binomial(link = "logit"),
  data = data_park
)
mytoc()
display(fit_lme4)

mean(y)
mean(y[male_name == 0 & white_name == 0 & n_skipped_full == 0])
sum(male_name == 0 & white_name == 0 & n_skipped_full == 0)
invlogit(-2.42 + 0.07*mean(male_name) + 0.08*mean(white_name) + 0.03*mean(n_skipped_full))
mean(predict(fit_lme4, type = "response"))

#' ## Simulate data and refit
#'
#' Check with simulated data that inference works.
set.seed(123)
a_respondent_sim <- rnorm(J, 0, sqrt(VarCorr(fit_lme4)$respondent))
a_item_sim <- rnorm(K, 0, sqrt(VarCorr(fit_lme4)$item))
b_sim <- fixef(fit_lme4)
X <- cbind(1, male_name, white_name, n_skipped_full)
p_sim <- invlogit(a_respondent_sim[respondent] + a_item_sim[item] + X %*% b_sim)
y_sim <- rbinom(N, 1, p_sim)
data_sim <- data.frame(data_park, y_sim)

tic("lme4 fit sim")
fit_lme4_sim <- glmer(
  y_sim ~ (1 | item) + (1 | respondent) + male_name + white_name + n_responses_full,
  family = binomial(link = "logit"),
  data = data_sim
)
mytoc()
display(fit_lme4_sim)

a_item_hat <- ranef(fit_lme4)$item
print(a_item_hat, digits = 1)

#' ## Plots

wordings <- read.csv(root("park_rule", "data", "park.txt"), header = FALSE)$V2
wordings <- substr(wordings, 2, nchar(wordings) - 1)
a_item_hat <- unlist(a_item_hat)
names(a_item_hat) <- wordings
print(sort(a_item_hat), digits=1)

#| label: fig-park-lme4-1
item_avg <- rep(NA, K)
for (k in 1:K) {
  item_avg[k] <- mean(y[item==k])
}
plot(item_avg, a_item_hat, pch=20)

#| label: fig-park-lme4-2
plot(logit(item_avg), a_item_hat, type = "n")
text(logit(item_avg), a_item_hat, names(a_item_hat), cex = .5)

#| label: fig-park-lme4-3
a_respondent_hat <- unlist(ranef(fit_lme4)$respondent)
respondent_avg <- rep(NA, J)
for (j in 1:J) {
  respondent_avg[j] <- mean(y[respondent == j])
}
plot(respondent_avg, a_respondent_hat, pch = 20, cex = .4)

#' # Stan models
#'
#' In general we prefer full Bayesian inference, even if it may take
#' more time.
#' 
#' ## Stan data
#'
#' CmdStanR needs the data in list format.
X <- cbind(male_name, white_name, n_skipped_full)
stan_data <- list(
  N = N,
  J = J,
  K = K,
  L = ncol(X),
  y = y,
  respondent = respondent,
  item = item,
  X = X
)

#' ## Hierarchical logistic regression with non-centered parameterization
#'
#' Hierarchical models of often benefit from non-centered
#' parameterization, so we start with that. We use `<multiplier=...>`
#' declaration to implement the non-centered parameterization. In the
#' model block, we use `bernoulli_logit_glm()`, which is more
#' efficient than `bernoulli_logit()` and can be used when the latent
#' model can be presented as a linear model. We use weak priors for
#' the coefficients and varying effect population scales.
park_1 <- cmdstan_model(root("park_rule", "park_1.stan"))
#| output: asis
print_stan_code(park_1$code())

#' When using the default sampling options and working interactively, we
#' quickly see very slow sampling. To investigate, we switch to use
#' one fifth of the iterations, and sampling takes about 6
#' minutes. Furthermore, we use option `init = 0.1` to initialize the
#' unconstrained parameters with random uniform values from range
#' [-0.1,0.1], which is often better initialization for varying
#' effects than the default range [-2,2].
#| label: fit_1
#| results: hide
#| cache: true
tic("Stan sampling model 1")
fit_1 <- park_1$sample(data = stan_data, init = 0.1,
                       iter_warmup = 200, iter_sampling = 200)
#'
#| cache: true
mytoc()

print(fit_1)

#' We see some suspiciously high $\widehat{R}$ values and low ESSs. We
#' didn't get divergence warnings, so the reason is unlikely to be a
#' funnel shaped posterior. We also didn't get maximum treedepth
#' exceedences warnings, so the there is no obvious high
#' correlations. We can further examine the sampler diagnostics:
fit_1$sampler_diagnostics() |> as_draws_rvars()

#' We are specifically interested in how efficient each Hamiltonian
#' Monte Carlo iteration is. This can be measured by the the number of
#' leapfrog steps `n_leapfrog__`, which is close to the number of log
#' density and gradient evaluations. Instead of examining
#' `n_leapfrog__` directly, it is common to examine `treedepth__` as
#' it scales logarithmically with respect to `n_leapfrog__`. More
#' specifically, $\text{treedepth\_\_}=\log_2(\text{treedepth\_\_}+1)$.
#' Average `treedepth__` is about 7 which is not high for hierarchical
#' model posteriors, and the variation measured by standard deviation
#' is low which indicates that the posterior curvature is not highly
#' varying and thus there is no strong funnel shape.

#' As ESSs are low for the parameter `a`, we check the trace plot.
draws_1 <- fit_1$draws(format = "df")
#| label: fig-fit_1-trace-a
#| out-width: 80%
draws_1 |>
  mcmc_trace(pars = "a") +
  labs(x = "Iteration")

#' There is clearly high auto-correlation. As we have `r N`
#' observations, we would expect the posterior for `a` to be narrow,
#' but now the posterior standard deviation is `r round(sd(draws_1$a),2)`.
#'
#' Examining the model code, we see
#' ```stan
#'   y ~ bernoulli_logit_glm(X, a + a_respondent[respondent] + a_item[item], b);
#' ```
#' and remember the discussion in Section 12.3 about parameters with
#' similar roles and identifiability.  Here all `a`, `a_respondent`
#' and `a_item` influence the total intercept. If value of `a`
#' increases, the total intercept stays the same if values of
#' `a_respondent` and `a_item` get lower at the same time. These
#' parameters are not well identified alone.
#'
#' Let's check the diagnostics for `a_respondent` and `a_item`, too.
draws_1 |>
  subset_draws(variable = "a_respondent") |>
  summarize_draws()
draws_1 |>
  subset_draws(variable = "a_item") |>
  summarize_draws()

#' The ESSs for `a_respondent` are high and ESSs for `a_item` are
#' low. As each `a_respondent` depends on a smaller number of
#' observations, the related uncertainty swamps the autocorrelation
#' due to the dependency.
#'
#' We can examine the scatter plot of `a` and sum of `a_item` to see
#' the strong correlation.
#| label: fig-fit_1-scatter-a-a_item
#| out-width: 80%
draws_rvars(a = as_draws_rvars(draws_1)$a,
            sum_a_item = as_draws_rvars(draws_1)$a_item |> rvar_sum()) |>
  mcmc_scatter(alpha = 0.5) +
  labs(y = "sum(a_item)")

#' ## Hierarchical logistic regression with zero sum-to-zero parameterization
#'
#' We can remove the identifiability problem by constraining
#' `a_respondent` and `a_item` to have sum zero. This can be easily
#' achieved in Stan with `sum_to_zero_vector` data type. In addition
#' of potentially improving sampling speed, sum-to-zero constraint
#' reducing posterior dependencies improves also interpretability of
#' the marginal posteriors of the parameters.  As `sum_to_zero_vector`
#' data type does not allow `multiplier`, we need to change how we
#' implement the non-centered parameterization.
park_2 <- cmdstan_model(root("park_rule", "park_2.stan"))
#| output: asis
print_stan_code(park_2$code())
#| label: fit_2
#| results: hide
#| cache: true
tic("Stan sampling model 2")
fit_2 <- park_2$sample(data = stan_data, init = 0.1,
                       iter_warmup = 200, iter_sampling = 200)
#'
#| cache: true
mytoc()

draws_2 <- fit_2$draws(format = "df")
print(fit_2)

#' The sampling time is reduced about 10%, and we get a big
#' improvement in $\widehat{R}$ and ESSs for `a`. The posterior
#' standard deviation of `a` is `r round(sd(draws_2$a),2)`, which is
#' much more sensible considering the data size.
#'
#' Now $\widehat{R}$ and ESSs indicate problems with
#' `sigma_item`. Examining the sampler diagnostics, shows that
#' `treedepth` is lower, indicating smaller posterior correlations,
#' and still with low standard deviation, indicating that curvature is
#' not highly varying.
fit_2$sampler_diagnostics() |> as_draws_rvars()

#' We examine the convergence diagnostics for `a_item`.
draws_2 |>
  subset_draws(variable = "a_item") |>
  summarize_draws()

#' The usual way to investigate funnels in hierarchical models, is to
#' look at the scatter plot of one of the variables and corresponding
#' prior scale. We examine `a_item` and `sigma_item`. Instead of
#' showing plot for all variables in vector `a_item`, we show here the
#' scatter plot for `a_item[1]` and `sigma_item´.
#| label: fig-fit_2-scatter-a_item-sigma_item
#| out-width: 80%
draws_2 |>
  subset_draws(variable = c("a_item[1]","sigma_item")) |>
  mcmc_scatter(alpha = 0.5)

#' There is no indication of funnel. This is hiding the issue, as the
#' sampling was actually done using parameters `z_respondent` and
#' `z_item`!
#' ```stan
#'   sum_to_zero_vector[J] z_respondent;
#'   sum_to_zero_vector[K] z_item;
#' ```
#' Thus, we should investigate the sampling performance for `z_item`:
draws_2 |>
  subset_draws(variable = "z_item") |>
  summarize_draws()

#' Now $\widehat{R}$ and ESSs indicate problems with `z_item`.
#' We examine the scatter plot for `z_item[1]` and `sigma_item`:
#| label: fig-fit_2-scatter-z_item-sigma_item
#| out-width: 80%
draws_2 |>
  subset_draws(variable = c("z_item[1]","sigma_item")) |>
  mcmc_scatter(alpha = 0.5)

#' There is a strong correlation. There is also slight banana shape,
#' but as `sigma_item` is constrained to be positive, the sampling is
#' done in log space, and we should use that also for the scatter
#' plot.
#| label: fig-fit_2-scatter-z_item-log_sigma_item
#| out-width: 80%
draws_2 |>
  subset_draws(variable = c("z_item[1]","sigma_item")) |>
  mcmc_scatter(transformations = list(sigma_item=log), alpha = 0.5) +
  labs(y = "log(sigma_item)")

#' The banana shape is weaker, and we have mostly linear dependency
#' which matches what inferred from the diagnostic values.
#' 
#' When non-centered parameterization is used, but the likelihood
#' contribution is strong, we get strong dependency between the latent
#' value `z` and `sigma`. In this case, we have a large number of
#' observations per each item, and the centered parameterization is
#' likely to be a better choice.
#'
#' In our first model, the non-centered parameterization was
#' implemented using vector data type with `<multiplier=...>`, which
#' hides the latent parameter, and it is easier to miss how to detect
#' the problem. With the explicit latent parameterization `z_item`,
#' the issue was easier to detect.
#' 
#' ## Hierarchical logistic regression with zero sum-to-zero parameterization
#'
#' We switch to using the centered parameterization for both
#' `a_respondent` and `a_item`.  We could still use non-centered for
#' `a_respondent`, but further experiments not shown here, indicate
#' that there is not much difference between the parameterizations for
#' `a_respondent` and thus we use the simpler form.
park_3 <- cmdstan_model(root("park_rule", "park_3.stan"))
#| output: asis
print_stan_code(park_3$code())
#| label: fit_3
#| results: hide
#| cache: true
tic("Stan sampling model 3")
fit_3 <- park_3$sample(data = stan_data, init = 0.1,
                       iter_warmup = 200, iter_sampling = 200)
#'
#| cache: true
mytoc()

draws_3 <- fit_3$draws(format = "df")
print(fit_3)

#' The sampling time is reduced about 10%, and we get a big
#' improvement in $\widehat{R}$ and ESSs for `sigma_item`.
#'
#' Examining the sampler diagnostics, shows that
#' `treedepth` is further reduced, indicating the posterior is easier than
#' with the first and second model parameterizations.
fit_3$sampler_diagnostics() |> as_draws_rvars()

#' ## Hierarchical logistic regression with sum-to-zero parameterization and centered predictors
#'
#' We can further reduce posterior dependencies by centering the
#' predictor values. We can do the centering in Stan code block
#' `transformed data`.
park_4 <- cmdstan_model(root("park_rule", "park_4.stan"))
#| output: asis
print_stan_code(park_4$code())
#| label: fit_4
#| results: hide
#| cache: true
tic("Stan sampling model 4")
fit_4 <- park_4$sample(data = stan_data, init = 0.1,
                       iter_warmup = 200, iter_sampling = 200)
#'
#| cache: true
mytoc()

draws_4 <- fit_4$draws(format = "df")
print(fit_4)

#' The sampling time is reduced by about 10%. Compared to the first
#' model the sampling time has reduced about 30%, all convergence
#' diagnostics look better and effective sample size per iteration is
#' much higher. The The posterior standard deviation of `a` is further
#' reduced to `r round(sd(draws_3$a),2)` as the correlation between
#' `a` and predictor coefficients has been removed by centering the
#' predictors.
#'
#' Looking at the sampler diagnostics, the `treedepth` is further
#' reduced compared to the previous model and posterior.
fit_4$sampler_diagnostics() |> as_draws_rvars()

#' We refit the final model using the default number of iterations.
#| label: fit_4-refit
#| results: hide
#| cache: true
tic("Stan sampling model 4")
fit_4 <- park_4$sample(data = stan_data, init = 0.1)
#'
#| cache: true
mytoc()

print(fit_4)

#' $\widehat{R}$ and ESS diagnostics look good. Although we are
#' running five times more iterations, the sampling time is only 50%
#' longer, which is due to 1) the most of the time spend in the warmup
#' before the adaptation is good, and 2) with longer adaptation the
#' actual sampling is faster.

#' Let's check the sampler diagnostics.
fit_4$sampler_diagnostics() |> as_draws_rvars()

#' The average `treedepth__` has further reduced, which means the
#' number of log density and gradient evaluations per iteration is
#' reduced approximately by 23% compared to running the algorithm with
#' fewer iteration. This further reduction comes from better
#' adaptation of the mass matrix and step size during the warmup. The
#' total sampling time is about 6 minutes with my laptop.
#'
#' ## Comparison of lme4 and Bayesian posterior estimates
#'
#' As there are quite many observations and lme4 is also using
#' marginalization over the latent values when estimating the standard
#' deviations of the varying intercepts, is there much difference in
#' the results?
#'
#' We compare the point estimates (conditional mode for lme4,
#' posterior mean for Bayes) and 90% intervals (normal approximation
#' for lme4, central posterior interval for Bayes) for `a_item` and
#' `a_respondent`.
#| label: fig-lme4-vs-stan-a_item
#| out-width: 80%
dr4 <- fit_4$draws() |> as_draws_rvars()
a_item_hat <- ranef(fit_lme4)$item
a_item_sd <- sqrt(as.numeric(attr(a_item_hat, "postVar")))
a_item_h <- as.numeric(a_item_hat$`(Intercept)`)
ggplot(data=NULL) +
  geom_abline(color="gray") + 
  stat_pointinterval(aes(x = a_item_h, ydist = dr4$a_item),
                     .width = c(0.90),
                     interval_size_range = c(0.4, 0.8),
                     alpha = 0.5) +
  geom_pointinterval(aes(y = mean(dr4$a_item), x = a_item_h,
                         xmin = a_item_h - 1.64*a_item_sd,
                         xmax = a_item_h + 1.64*a_item_sd),
                     interval_size_range = c(0.4, 0.8),
                     alpha = 0.5) +
  labs(x = "lme4", y = "Stan", title = "a_item")

#| label: fig-lme4-vs-stan-a_respondent
#| out-width: 80%
a_respondent_hat <- ranef(fit_lme4)$respondent
a_respondent_sd <- sqrt(as.numeric(attr(a_respondent_hat, "postVar")))
a_respondent_h <- as.numeric(a_respondent_hat$`(Intercept)`)
ggplot(data=NULL) +
  geom_abline(color="gray") + 
  stat_pointinterval(aes(x = a_respondent_h, ydist = dr4$a_respondent),
                     .width = c(0.90),
                     interval_size_range = c(0.4, 0.8),
                     alpha = 0.1) +
  geom_pointinterval(aes(y = mean(dr4$a_respondent), x = a_respondent_h,
                         xmin = a_respondent_h - 1.64*a_respondent_sd,
                         xmax = a_respondent_h + 1.64*a_respondent_sd),
                     interval_size_range = c(0.4, 0.8),
                     alpha = 0.1) +
  labs(x = "lme4", y = "Stan", title = "a_respondent")

#' There is not much difference, although Stan estimates have slightly
#' wider range.  For most purposes lme4 result would be just fine and
#' computationally faster than full Bayes for this kind of big data.
#'
#' The wider range of Stan estimates is likely due to effect of
#' integrating over the uncertainty in `sigma_item` and
#' `sigma_respondent`:
#| label: fig-lme4-vs-stan-sigmas
#| fig-height: 4
#| fig-width: 8
#| out-width: 80%
ggplot(data=NULL) +
  stat_slab(aes(xdist = dr4$sigma_item), density = "unbounded", trim = TRUE, fill = NA, color = "black") +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line.y = element_blank(), strip.text.y = element_blank()) +
  coord_cartesian(expand = FALSE) +
  ## xlim(c(0.005, 0.025)) +
  labs(x = "sigma_item", y = "") +
  geom_vline(xintercept = as.data.frame(VarCorr(fit_lme4))$sdcor[2], color = "black", linetype="dashed")+
  annotate(geom = "text", x = as.data.frame(VarCorr(fit_lme4))$sdcor[2]*1.02, y = .97, hjust=0, label = "lme4 estimate") +
ggplot(data=NULL) +
  stat_slab(aes(xdist = dr4$sigma_respondent), density = "unbounded", trim = TRUE, fill = NA, color = "black") +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line.y = element_blank(), strip.text.y = element_blank()) +
  coord_cartesian(expand = FALSE) +
  ## xlim(c(0.005, 0.025)) +
  labs(x = "sigma_respondent", y = "") +
  geom_vline(xintercept = as.data.frame(VarCorr(fit_lme4))$sdcor[1], color = "black", linetype="dashed") +
  annotate(geom = "text", x = as.data.frame(VarCorr(fit_lme4))$sdcor[1]*1.002, y = .97, hjust=0, label = "lme4 estimate")

#' 
#' <br />
#' 
#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2025--2026, Andrew Gelman and Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2025--2026, Andrew Gelman and Aki Vehtari, licensed under CC-BY-NC 4.0.

#' ---
#' title: "Sleepstudy: Prior Specification and model checking"
#' author: "Paul Bürkner and Aki Vehtari"
#' date: 2025-10-22
#' date-modified: today
#' date-format: iso
#' format:
#'   html:
#'     encoding: UTF-8
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
#' ---

#' # Setup  {.unnumbered}
#' 
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
library("rprojroot")
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(ggplot2)
library(patchwork)
library(dplyr)
library(loo)
library(brms)
options(mc.cores = 4)
dir.create(root("sleep_study","models/"))
BRMS_MODEL_DIR <- root("sleep_study", "models/")
options(future.globals.maxSize = 1e9)
library(priorsense)
options(priorsense.plot_help_text = FALSE)
theme_set(bayesplot::theme_default(base_family = "sans"))

#' # Main Story and Messages
#'
#' Prior distributions are at the heart of Bayesian statistics and are
#' mentioned as one of its defining features in almost all
#' introductions. Yet, in practice, specifying priors remains a highly
#' challenging and complex topic that tends to cause a lot of
#' confusion for people having to deal with it. In this case study, I
#' will try to clarify some of this confusion by explaining the
#' different purposes of priors and things that should be considered
#' when specifying them. Towards the end, I will highlight future
#' research directions to make prior specification easier for
#' everyone.
#' 
#' Main Messages:
#' 
#' - Priors specification can follow different purposes.
#' - Prior specification is hard for everyone, as is statistical modeling in general.
#' - We should aim to identify when thinking about the prior is important and when it is not.
#' - We need a lot of further research in the areas of prior specification and elicitation.
#' 
#' 
#' ## The Sleepstudy Data
#' 
#' Reasons for choosing the `sleepstudy` data set:
#' 
#' - Few variables all of which are easy to understand
#' - easy yet important multilevel structure 
#' - sensible to express with both linear and generalized linear models
#' - non-trivial error distributions
#' - independent priors are sensible(ish) due to the small number of parameters
#' - well known to a lot of R users
#' 
#' Introduce the research questions, gather variables, and available data
#'
#' Days 0-1 were adaptation and training (T1/T2), day 2 was baseline (B);
#' sleep deprivation started after day 2. We Drop days 0-1, and make the
#' baseline to be new 0.
data("sleepstudy", package = "lme4")
conditions <- make_conditions(sleepstudy, "Subject", incl_vars = FALSE)
sleepstudy <- sleepstudy |>
  filter(Days >= 2) |>
  mutate(Days = Days - 2)

#' Plot the data
#| label: fig-sleepstudy-data
#| fig-height: 4
#| fig-width: 8
sleepstudy |>
  ggplot(aes(Days, Reaction)) + 
  geom_point() +
  facet_wrap("Subject", ncol = 6) +
  scale_x_continuous(breaks = 0:7) +
  labs(y = "Reaction time (ms)")
#ggsave("plots/sleepstudy_data.pdf", height = 4, width = 8)

#' ## Simple Linear Model
#'
#' Prior base
prior_lin_base <- prior(normal(200, 100), class = b, coef = "Intercept") +
  prior(normal(0, 20), class = b, coef = "Days") +
  prior(exponential(0.02), class = sigma)
#'
#' Model base and sample from the posterior
#| results: hide
#| cache: true
fit_lin_base <- brm(
  Reaction ~ 0 + Intercept + Days, 
  data = sleepstudy,
  family = gaussian(),
  prior = prior_lin_base, 
  file = paste0(BRMS_MODEL_DIR, "fit_lin_base"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit_lin_base, priors = TRUE)

#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit_lin_base
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit_lin_base), points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")
  

#' ## Simple Linear Model (centered predictors)
#' 
#' Points to discuss:
#' 
#' - priors on original or centered intercept?
#' - dependency of the prior on marginal moments of the data?
#' - different qualitative options for priors on b and sigma
#'
#' Prior 1
prior1 <- prior(normal(250, 100), class = Intercept) +
  prior(normal(0, 20), class = b) +
  prior(exponential(0.02), class = sigma)

#' Sample from the prior 1
fit1_prior <- brm(
  Reaction ~ 1 + Days, 
  data = sleepstudy,
  family = gaussian(),
  prior = prior1, 
  sample_prior = "only",
  file = paste0(BRMS_MODEL_DIR, "fit1_prior"),
  file_refit = "on_change"
) 

#' Prior predictive checking
#| label: fig-sleepstudy-pp_check-fit1_prior
#| fig-height: 2.5
#| fig-width: 5
set.seed(652312)
pp_check(fit1_prior, ndraws = 100) + 
  ylim(c(0, 0.02)) +
  theme_sub_axis_y(line = element_blank())
#ggsave("plots/prior_pred_lin_sleep.pdf", height = 2.5, width = 5)

#' Model 1: sample from the posterior
fit1 <- brm(
  Reaction ~ 1 + Days, 
  data = sleepstudy,
  family = gaussian(),
  prior = prior1, 
  file = paste0(BRMS_MODEL_DIR, "fit1"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit1, priors = TRUE)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit1
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit1), points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")
#' Posterior predictive checking
#| label: fig-sleepstudy-pp_check-fit1
#| fig-height: 2.5
#| fig-width: 5
pp_check(fit1, ndraws = 50) +
  theme_sub_axis_y(line = element_blank())

#' Prior sensitivity analysis
#| label: fig-sleepstudy-priorsense-fit1
#| fig-height: 3
#| fig-width: 7
powerscale_plot_dens(fit1, variable = c("b_Intercept", "b_Days", "sigma"),
                     component = "prior")
## ggsave("plots/sleep_priorsense_weak.pdf", width = 7, height = 3)

#' ## Simple Linear Model (informative priors)
#' 
#' Points to discuss:
#' 
#' - Priors will be influencing the posterior if chosen to be informative enough
#' - For models that are simple relative to the amount of data, prior distributions
#'   are unlikely to affect the posterior strongly, unless prior are very informative
#' 
#' Prior 2
prior2 <- prior(normal(250, 100), class = Intercept) +
  prior(normal(0, 1), class = b) +
  prior(exponential(0.02), class = sigma)

#' Model 2: sample from the posterior
fit2 <- brm(
  Reaction ~ 1 + Days, 
  data = sleepstudy,
  family = gaussian(),
  prior = prior2, 
  file = paste0(BRMS_MODEL_DIR, "fit2"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit2)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit2
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit2), points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' Prior sensitivity analysis
#| label: fig-sleepstudy-priorsense-fit2
#| fig-height: 3
#| fig-width: 7
powerscale_plot_dens(fit2, variable = c("b_Intercept", "b_Days", "sigma"),
                     component = "prior")
## ggsave("plots/sleep_priorsense_strong.pdf", width = 7, height = 3)

#' ## Simple Linear Model (informative priors with fat tails)
#' 
#' Points to discuss:
#' - tails of the priors (normal vs. student-t)
#' 
#' Prior 2b
prior2b <- prior(normal(250, 100), class = Intercept) +
  prior(student_t(7, 0, 1), class = b) +
  prior(exponential(0.02), class = sigma)

#' Model 2b: sample from the posterior
fit2b <- brm(
  Reaction ~ 1 + Days, 
  data = sleepstudy,
  family = gaussian(),
  prior = prior2b, 
  file = paste0(BRMS_MODEL_DIR, "fit2b"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit2b)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit2b
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit2b), points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' Prior sensitivity analysis
#| label: fig-sleepstudy-priorsense-fit2b
#| fig-height: 3
#| fig-width: 7
powerscale_plot_dens(fit2b, variable = c("b_Intercept", "b_Days", "sigma"),
                     component = "prior")
## ggsave("plots/sleep_priorsense_student.pdf", width = 7, height = 3)

#' Illustrate difference between normal and Student-t prior:
#| label: fig-compare-normal-student-density
#| fig-height: 2
#| fig-width: 5
x <- seq(-4, 4, 0.01)
d1 <- dnorm(x, 0, 1)
d2 <- dstudent_t(x, 7, 0, 1)
data.frame(
  d = c(d1, d2), 
  x = rep(x, 2), 
  Prior = rep(c("normal(0, 1)", "Student-t(7, 0, 1)"), each = length(x))
) |>
  ggplot(aes(x, d, color = Prior)) +
  geom_line(size = 0.8) +
  xlab(expression(b[1])) +
  ylab("Density")
##ggsave("plots/normal_student_density.pdf", height = 2, width = 5)

#' Compute CI-bound for an exponential prior:
qexp(c(0.025, 0.975), 0.02)

#' ## Linear Varying Intercept Model
#' 
#' Points to discuss:
#' 
#' - How to represent unidimensional multilevel structures via priors
#' - priors on hyperparameters (SDs)
#' - shall the prior on sigma change now that we add more terms?
#'
#' Prior 3
prior3 <- prior(normal(250, 100), class = Intercept) +
  prior(normal(0, 20), class = b) +
  prior(exponential(0.02), class = sigma) +
  prior(exponential(0.02), class = sd)

#' Model 3: sample from the posterior
fit3 <- brm(
  Reaction ~ 1 + Days + (1 | Subject), 
  data = sleepstudy,
  family = gaussian(),
  prior = prior3, 
  file = paste0(BRMS_MODEL_DIR, "fit3"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit3)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-1-fit3
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit3), plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-2-fit3
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit3, conditions = conditions, re_formula = NULL),
     ncol = 6, points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' ## Linear Varying Intercept and Slope Model
#' 
#' Points to discuss:
#' 
#' - How to represent multidimensional multilevel structures via priors
#' - priors on hyperparameters (SDs and correlations)
#' - The implications of the LKJ prior for correlation matrices
#' 
#' Prior 4
prior4 <- prior(normal(250, 100), class = Intercept) +
  prior(normal(0, 20), class = b) +
  prior(exponential(0.04), class = sigma) +
  prior(exponential(0.04), class = sd, group = Subject, coef = Intercept) +
  prior("exponential(1.0/15)", class = sd, group = Subject, coef = Days) +
  prior(lkj(1), class = cor)

#' Model 4: sample from the posterior
fit4 <- brm(
  Reaction ~ 1 + Days + (1 + Days | Subject), 
  data = sleepstudy,
  family = gaussian(),
  prior = prior4, 
  file = paste0(BRMS_MODEL_DIR, "fit4"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit4)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-1-fit4
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit4), plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' Posterior predictive checking
pp_check(fit4) +
  theme_sub_axis_y(line = element_blank())

#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-2-fit4
#| fig-height: 4
#| fig-width: 7
plot(conditional_effects(fit4, conditions = conditions, re_formula = NULL),
     ncol = 6, points = TRUE, plot = FALSE)[[1]] +
  scale_x_continuous(breaks = 0:9) +
  labs(y = "Reaction time (ms)")
## ggsave("plots/sleep_multilevel_ceffects.pdf", height = 4, width = 7)

#' Illustrate the marginal LKJ(1) prior for different dimensions.
#| label: fig-sleepstudy-LKJ1
#| fig-height: 2.5
#| fig-width: 5
dmLKJ <- function(x, eta, d) {
  dbeta((x + 1) / 2, eta + (d - 2)/2, eta + (d - 2)/2)
}
data.frame(x = rep(seq(-0.999, 0.999, 0.001), 3)) |>
  mutate(
    d = rep(c(2, 5, 10), each = n()/3),
    dens = dmLKJ(x, 1, d),
    d = factor(d)
  ) |> 
  ggplot(aes(x, dens, color = d)) +
  geom_line(size = 0.8) +
  scale_color_viridis_d() +
  ylab("Density") +
  xlab(expression(rho))
## ggsave("plots/LKJ_1_density.pdf", height = 2.5, width = 5)

#' ## Log-Linear Prior-Only Model
#'
#' Reuse priors from the normal linear model: Prior ln1
prior_ln1 <- prior(normal(250, 100), class = Intercept) +
  prior(normal(0, 20), class = b) +
  prior(exponential(0.02), class = sigma)

#' Sample from the prior ln1
fit_ln1_prior <- brm(
  Reaction ~ 1 + Days, 
  data = sleepstudy,
  family = lognormal(),
  prior = prior_ln1, 
  sample_prior = "only",
  file = paste0(BRMS_MODEL_DIR, "fit_ln1_prior"),
  file_refit = "on_change"
) 

#' Prior predictive checking
#| label: fig-sleepstudy-pp_check-fit_ln1_prior
#| fig-height: 2.5
#| fig-width: 5
pp_check(fit_ln1_prior) +
  theme_sub_axis_y(line = element_blank())

#' Prior predictive checking
set.seed(652312)
prp_ln <- apply(posterior_predict(fit_ln1_prior), 2, function(x) mean(x[is.finite(x)]))
prp_ln_dat <- data.frame(y = sleepstudy$Reaction, yrep = prp_ln)
gg_ln1_prior <- ggplot(prp_ln_dat, aes(y, yrep)) +
  geom_point() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  ylab("Mean yrep")
# ggsave("plots/prior_pred_ln1_sleep.pdf", height = 2.5, width = 5)

#' Use more sensible priors: Prior ln2
prior_ln2 <- prior(normal(5, 0.55), class = Intercept) +
  prior(normal(0, 0.2), class = b) +
  prior(exponential(3), class = sigma)

#' Sample from the prior ln2
fit_ln2_prior <- brm(
  Reaction ~ 1 + Days, 
  data = sleepstudy,
  family = lognormal(),
  prior = prior_ln2, 
  sample_prior = "only",
  file = paste0(BRMS_MODEL_DIR, "fit_ln2_prior"),
  file_refit = "on_change"
) 

#' Prior predictive checking
#| label: fig-sleepstudy-pp_check-fit_ln2_prior
#| fig-height: 2.5
#| fig-width: 5
pp_check(fit_ln2_prior) +
  theme_sub_axis_y(line = element_blank())

#' Prior predictive checking
set.seed(652312)
prp_ln <- apply(posterior_predict(fit_ln2_prior), 2, function(x) mean(x[is.finite(x)]))
prp_ln_dat <- data.frame(y = sleepstudy$Reaction, yrep = prp_ln)
gg_ln2_prior <- ggplot(prp_ln_dat, aes(y, yrep)) +
  geom_point() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  ylab("Mean yrep")
# ggsave("plots/prior_pred_ln2_sleep.pdf", height = 2.5, width = 5)

#' Prior predictive checking comparing priors ln1 and ln2
#| label: fig-sleepstudy-pp_check-fit_ln1_ln2_prior
#| fig-height: 2.5
#| fig-width: 6
gg_ln1_prior + gg_ln2_prior
## ggsave("plots/prior_pred_ln_sleep.pdf", height = 2.5, width = 6)

#' Sample from the posterior
fit_ln2 <- brm(
  Reaction ~ 1 + Days, 
  data = sleepstudy,
  family = lognormal(),
  prior = prior_ln2, 
  file = paste0(BRMS_MODEL_DIR, "fit_ln2"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit_ln2)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-1-ln2
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit_ln2), plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' Prior sensitivity analysis
#| label: fig-sleepstudy-priorsense-fit_ln2
#| fig-height: 3
#| fig-width: 7
powerscale_plot_dens(fit_ln2, variable = c("b_Intercept", "b_Days", "sigma"),
                     component = "prior")
## ggsave("plots/sleep_priorsense_ln_reg.pdf", width = 7, height = 3)

#' ## Log-Linear Varying Intercept and Slope Model
#' 
#' Points to discuss:
#' 
#' - positive only family may be preferred theoretically but may not always
#'   be required
#' - How a non-identity link (log in this case) messes with our intuition
#'   about parameters and hence with prior specification
#' - how Jensen's inequality makes direct translations of prior difficult
#' - how Jacobian adjustment would be needed to ensure equivalence of two
#'   prior on different scales.

#' Prior 5
prior5 <- prior(normal(5, 0.55), class = Intercept) +
  prior(normal(0, 0.2), class = b) +
  prior(exponential(6), class = sigma) +
  prior(exponential(6), class = sd, group = Subject, coef = Intercept) +
  prior(exponential(10), class = sd, group = Subject, coef = Days) +
  prior(lkj(1), class = cor)

#' Model5: sample from the posterior
fit5 <- brm(
  Reaction ~ 1 + Days + (1 + Days | Subject), 
  data = sleepstudy,
  family = lognormal(),
  prior = prior5, 
  file = paste0(BRMS_MODEL_DIR, "fit5"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit5)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit5
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit5), plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' Posterior predictive checking
#| label: fig-sleepstudy-pp_check-fit5
#| fig-height: 2.5
#| fig-width: 5
pp_check(fit5) +
  theme_sub_axis_y(line = element_blank())

#' Posterior predictive checking comparing fit4 and fit5
#| label: fig-sleepstudy-pp_check-fit4-fit5
#| fig-height: 4
#| fig-width: 7
pp_check(fit4, type = "intervals") + labs(y = "Reaction time (ms)") +
  pp_check(fit5, type = "intervals") + 
  plot_layout(guides = "collect") 
## ggsave("plots/sleep_multilevel_ppc.pdf", height = 4, width = 7)

#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-plus-fit5
#| fig-height: 4
#| fig-width: 7
plot(conditional_effects(fit5, conditions = conditions, re_formula = NULL),
     ncol = 6, points = TRUE, plot = FALSE)[[1]] +
  scale_x_continuous(breaks = 0:9) +
  labs(y = "Reaction time (ms)")
## ggsave("plots/sleep_multilevel_ceffects_ln.pdf", height = 4, width = 7)

#' Prior sensitivity analysis
#| label: fig-sleepstudy-priorsense-fit5
#| fig-height: 4
#| fig-width: 13
vars5 <- c("b_Intercept", "b_Days", "sd_Subject__Intercept", 
           "sd_Subject__Days", "cor_Subject__Intercept__Days", "sigma")
powerscale_plot_dens(fit5, variable = vars5, component = "prior")
## ggsave("plots/sleep_priorsense_ln_mlm.pdf", width = 13, height = 4)

#' Compare with linear multilevel model which indicated the lognormal model
#' having a little better predictive performance.
loo(fit4, fit5)

#' ## Log-Linear Distributional Multilevel Model

#' ## Points to discuss:
#'
#' - Parameters for standard deviations on the log or log-log scale
#'   or hard to understand and hence set priors on.
#' - This is specifically true for standard deviations of standard deviatons
#'   on the log or log-log scale.
#' - In theory, we would not need the
#'   exponential link on sigma but then we had to care for the positivity of the
#'   varying intercepts on sigma and hence would have to specify, for example, an
#'   hierarchical Gamma rather than an hierarchical normal prior.
#' - Look at prior predictions for the correlations to demonstrate
#'   the effect of the LKJ prior for larger than 2x2 matrices
#' 
#' Prior 6
prior6 <- prior(normal(5.5, 0.55), class = Intercept) +
  prior(normal(0, 0.2), class = b) +
  prior(normal(0, 0.3), class = Intercept, dpar = sigma) +
  prior(exponential(3), class = sd, group = Subject, coef = Intercept) +
  prior(exponential(5), class = sd, group = Subject, coef = Days) +
  prior(exponential(3), class = sd, dpar = sigma, group = Subject) +
  prior(lkj(1), class = cor)

#' Model 6: sample from the posterior
fit6 <- brm(
  bf(Reaction ~ 1 + Days + (1 + Days |S| Subject),
     sigma ~ 1 + (1 |S| Subject)),
  data = sleepstudy,
  family = lognormal(),
  prior = prior6, 
  file = paste0(BRMS_MODEL_DIR, "fit6"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit6)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit6
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit6, conditions = conditions, re_formula = NULL),
     ncol = 6, points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

## ## Exgaussian Distributional Multilevel Model

#' ## Points to discuss:
#'
#' - All the models before are relatively robust to the choice of priors:
#'   Even completely flat priors work ok and don't produce much different results
#'   because the models are parsimonious relative to the amount of data and quality
#' - This robustness is no longer the case for the upcoming models
#' - Exgaussian is linear yet can handle strong right skewness at the expense
#'   of not strictly respecting the lower boundary of zero
#' - Avoids the log-log awkwardness of the lognormal and related models
#' - Interpretation of the beta (skewness) parameter at least somewhat
#'   intuitive
#'
#' Prior 7
prior7 <- prior(normal(250, 100), class = Intercept) +
  prior(normal(0, 20), class = b) +
  prior(normal(0, 5), class = Intercept, dpar = sigma) +
  prior(normal(0, 5), class = beta) +
  prior(exponential(0.04), class = sd, group = Subject, coef = Intercept) +
  prior(exponential(0.05), class = sd, group = Subject, coef = Days) +
  prior(exponential(0.2), class = sd, dpar = sigma, group = Subject) +
  prior(lkj(1), class = cor)

#' Model 7: sample from the posterior
fit7 <- brm(
  bf(Reaction ~ 1 + Days + (1 + Days |S| Subject),
     sigma ~ 1 + (1 |S| Subject)), 
  data = sleepstudy,
  family = exgaussian(),
  prior = prior7, 
  inits = 0,
  cores = 4,
  file = paste0(BRMS_MODEL_DIR, "fit7"),
  file_refit = "on_change"
) 

#' Posterior summary for model 7
summary(fit7)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit7
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit7, conditions = conditions, re_formula = NULL),
     ncol = 6, points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

## ## Exgaussian Distributional Multilevel Model (default priors)

#' Points to discuss:
#' 
#' - some parameters are better informed by the data than others
#' - using our custom priors and default priors barely matters for most
#'   except for the skewness parameter where the effect is clearly visible.
#' - default priors are already chosen in an effort to be sensible(ish)
#'
#' Model 8: sample from the posterior
fit8 <- brm(
  bf(Reaction ~ 1 + Days + (1 + Days |S| Subject),
     sigma ~ 1 + (1 |S| Subject)), 
  data = sleepstudy,
  family = exgaussian(),
  # no prior argument: default priors are used
  inits = 0,
  cores = 4,
  file = paste0(BRMS_MODEL_DIR, "fit8"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit8)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit8
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit8, conditions = conditions, re_formula = NULL),
     ncol = 6, points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' ## Exgaussian Distributional Multilevel Model (flat priors)
#' 
#' Points to discuss:
#'
#' - without reasonable(ish) priors, sampling may fall apart with some
#'   really long-running chains, divergent transitions and increased Rhats.
#' - priors matter, but if you don't specify quite informative priors
#'   or have small data relative to the model complexity, you will
#'   likely not see a difference compared to sensible default priors
#' 
#' Extract all default priors:
prior9 <- get_prior(
  bf(Reaction ~ 1 + Days + (1 + Days |S| Subject),
     sigma ~ 1 + (1 |S| Subject)), 
  data = sleepstudy,
  family = exgaussian()
)
#' Make all priors flat (except for the varying coefficients' priors):
prior9$prior <- ""

#' Model 9: sample from the posterior
fit9 <- brm(
  bf(Reaction ~ 1 + Days + (1 + Days |S| Subject),
     sigma ~ 1 + (1 |S| Subject)), 
  data = sleepstudy,
  family = exgaussian(),
  prior = prior9, 
  inits = 0,
  cores = 4,
  file = paste0(BRMS_MODEL_DIR, "fit9"),
  file_refit = "on_change"
) 

#' Posterior summary
summary(fit9)
#' Posterior conditional effects
#| label: fig-sleepstudy-conditional_effects-fit9
#| fig-height: 4
#| fig-width: 8
plot(conditional_effects(fit9, conditions = conditions, re_formula = NULL),
     ncol = 6, points = TRUE, plot = FALSE)[[1]] +
  labs(y = "Reaction time (ms)")

#' # Model comparison and checking
#'
#' Compute LOO-CV for models 3, 4, and 5
fit3 <- add_criterion(fit3, criterion = "loo", save_psis = TRUE, reloo = TRUE)
fit4 <- add_criterion(fit4, criterion = "loo", save_psis = TRUE, reloo = TRUE)
fit5 <- add_criterion(fit5, criterion = "loo", save_psis = TRUE, reloo = TRUE)

#' Compare varying intercept and slope models with normal and log-normal data models
loo_compare(fit4, fit5)

#' Log-normal is not different from normal model, which makes sense as
#' all reaction times are fra from 0.
#' 
#' Posterior predictive checking of fit4 using LOO predictive intervals
#| label: fig-sleepstudy-ppc-loo_intervals-fit4
#| fig-height: 4
#| fig-width: 11
pp_check(fit4, type = "loo_intervals") +
  labs(y = "Reaction time (ms)")
## ggsave(file = "sleepstudy_fit4_loo_intervals.pdf", width = 11, height = 4)

#' There are clearly some outliers.
#'
#' Create varying intercept (fit3t) and varying intercept and slope
#' (fit4t) models with Student's t data model
#| results: hide
#| cache: true
fit3t <- update(fit3, family = student())
fit4t <- update(fit4, family = student())

#' Compute LOO-CV for models 3, 4, and 5
fit3t <- add_criterion(fit3t, criterion = "loo", save_psis = TRUE)
fit4t <- add_criterion(fit4t, criterion = "loo", save_psis = TRUE)

#' Posterior predictive checking of fit4t using LOO predictive intervals
#| label: fig-sleepstudy-ppc-loo_intervals-fit4t
#| fig-height: 4
#| fig-width: 11
pp_check(fit4t, type = "loo_intervals") +
  labs(y = "Reaction time (ms)")

#' The LOO predictive intervals look better now.

#| label: fig-sleepstudy-ppc-loo_intervals-fit4-fit5-fit4t
#| fig-height: 3.5
#| fig-width: 11
pp_check(fit4, type = "loo_pit_ecdf", ndraws = 4000) +
  pp_check(fit5, type = "loo_pit_ecdf", ndraws = 4000) +
  pp_check(fit4t, type = "loo_pit_ecdf", ndraws = 4000) 
## ggsave(file = "sleepstudy_loo_pit_ecdf.pdf", width = 11, height = 3.5)

#' We see that normal and log-normal models have too wide predictive
#' distribution for most observations, which is due to a few outliers
#' inflating the residual scale. LOO-PIT plot for Student's t model
#' looks better.

#' Compare normal and Student's t models
loo_compare(fit4, fit4t)

#' Student's t model has much better predictive performance.
#' 
#' Examine how much adding varying slope improved predictive performance in case of normal data model:
loo_compare(fit3, fit4)

#' Examine how much adding varying slope improved predictive
#' performance in case of Student's t data model:
loo_compare(fit3t, fit4t)

#' When using Student's t model, the predictive performance difference
#' between not using or using varying slope is bigger.
#' 

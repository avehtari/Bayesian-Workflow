#' ---
#' title: "Uncertainty in Bayesian LOO-CV Model Comparison"
#' author: "Aki Vehtari"
#' date: 2025-06-23
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
#' This notebook includes the code for the Bayesian Workflow book
#' Section 9.4 *Model selection using predictive checking and predictive
#' performance* (adapted from Section 5 of paper by @Sivula-Magnusson-Vehtari:2025).
#'
#' # Introduction
#' 
#' We demonstrate the use of uncertainty quantification of the
#' predictive performance difference with three real-data examples. We
#' assume that the true data generating processes are more complex
#' than the models used. We cover all three scenarios
#'  1. very similar predictions,
#'  2. model misspecification,
#'  3. and small data,
#' that can affect how well calibrated the normal approximation.
#'
#' Inference was made using MCMC with 4 chains with 1000 warmup and
#' 1000 sampling iterations. Convergence diagnostics
#' [@Vehtari-Gelman-Simpson-etal:2021], using the
#' `posterior` package, [@Burkner-Gabry-Kay-etal:2024]
#' indicated reliable posterior inference. For LOO-CV we used the
#' `loo` package [@vehtari2022loopkg], which uses fast
#' PSIS-LOO [@Vehtari+Gelman+Gabry:2017_practical] for
#' computation.
#'
#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=TRUE, comment=NA, out.width='95%')
#' **Load packages**
#| cache: FALSE
library(tidyr)
library(dplyr)
library(tibble)
library(loo)
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 1)
library(rstanarm)
library(posterior)
options(posterior.num_args=list(digits=2))
library(pillar)
options(pillar.negative = FALSE)
library(tinytable)
options(tinytable_format_num_fmt = "significant_cell", tinytable_format_digits = 2, tinytable_tt_digits=2)
library(ggplot2)
library(ggdist)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))

#' Define LaTeX commands
#' \newcommand*{\Mk}{{{\mathrm{M}_k}}}
#' \newcommand*{\Md}{{{\mathrm{M}_a,\mathrm{M}_b}}}
#' \newcommand*{\yobs}{{y}}
#' \newcommand*{\elpdHat}[2]{{\widehat{\mathrm{elpd}}_\mathrm{\scriptsize LOO}\bigr(#1 \mid #2\bigl)}}
#' \newcommand*{\elpdHati}[3]{{\widehat{\mathrm{elpd}}_{\mathrm{\scriptsize LOO},\, #3}\bigr(#1 \mid #2\bigl)}}
#' \newcommand*{\seHat}[2]{{\widehat{\mathrm{SE}}_\mathrm{\scriptsize LOO}\bigr(#1 \mid #2\bigl)}}
#' 
#'
#' # Primate milk
#'
#' @mcelreath2020statistical describes the primate milk data: *``A
#' popular hypothesis has it that primates with larger brains produce
#' more energetic milk, so that brains can grow quickly... The
#' question here is to what extent energy content of milk, measured
#' here by kilocalories, is related to the percent of the brain mass
#' that is neocortex... We'll end up needing female body mass as well,
#' to see the masking that hides the relationships among the
#' variables.''* The data include 17 different primate species. The
#' target variable is the energy content of milk (kcal.per.g) and the
#' covariates are the percent of the brain mass that is neocortex
#' (neocortex) and the logarithm of female body mass (log(mass)). The
#' covariates and target are centered and scaled to have unit
#' variance.
data(milk)
milk <- milk |>
  drop_na() |>
  mutate(neocortex = neocortex.perc / 100)
head(milk)
#'
#' We use the following four models with weakly
#' informative normal(0, 1) priors for the coefficients and an
#' exponential(1) prior for the residual scale:
#' \begin{align*}
#'   \mathrm{M}_1: \quad & \mathrm{kcal.per.g} \sim \operatorname{normal}(\alpha, \sigma) \\
#'   \mathrm{M}_2: \quad& \mathrm{kcal.per.g} \sim \operatorname{normal}(\alpha + \beta_1\times\mathrm{neocortex}, \sigma) \\
#'   \mathrm{M}_3: \quad & \mathrm{kcal.per.g} \sim \operatorname{normal}(\alpha + \beta_2\times\mathrm{log(mass)}, \sigma) \\
#'   \mathrm{M}_4: \quad & \mathrm{kcal.per.g} \sim \operatorname{normal}(\alpha + \beta_1\times\mathrm{neocortex} + \beta_2\times\mathrm{log(mass)}, \sigma).
#' \end{align*}
#'
#' We fit the models with the `rstanarm` package
#' [@Goodrich-Gabry-Ali-etal:2024]
#| label: milk_models
#| results: hide
#| cache: true
M_1 <- stan_glm(kcal.per.g ~ 1, data = milk, seed = 2030,
                 prior = normal(0, 1, autoscale = TRUE), refresh = 0)
M_2 <- update(M_1, formula = kcal.per.g ~ neocortex)
M_3 <- update(M_1, formula = kcal.per.g ~ log(mass))
M_4 <- update(M_1, formula = kcal.per.g ~ neocortex + log(mass))

#' We compute LOO-CV estimates using the fast PSIS-LOO method
#' [@Vehtari+Gelman+Gabry:2017_practical, @Vehtari+etal:2024:PSIS]
loo1 <- loo(M_1)
loo2 <- loo(M_2)
loo3 <- loo(M_3)
loo4 <- loo(M_4)

#' We compare the models $\mathrm{M}_1, \mathrm{M}_2,\mathrm{M}_3,\mathrm{M}_4$.
#'
#' The current version of `loo_compare()` shows `elpd_diff` and `diff_se`
(loo_comp <- loo_compare(loo1, loo2, loo3, loo4))

#' We add the probability that a model has worse predictive
#' performance than the model with the best predictive performance
#' using the normal approximation.
loo_comp |>
  as.data.frame() |>
  rownames_to_column("model") |>
  dplyr::select(model, elpd_diff, se_diff) |>
  mutate(p = ifelse(elpd_diff==0,NA,pnorm(0, elpd_diff, se_diff))) |>
  tt() |>
  format_tt(replace = "-")

#' In the paper, we the comparison was reported using the $M_1$ as the
#' baseline, and the probability that a model has better predictive
#' performance than the baseline.
loo_comp2 <- sapply(list(loo1, loo2, loo3, loo4), \(x) pointwise(x, "elpd_loo")) 
colnames(loo_comp2) <- c("M_1","M_2","M_3","M_4")
loo_comp2 <- loo_comp2 |>
  as_tibble() |>
  mutate(across(M_1:M_4, ~ .x - M_1))
tibble(model=colnames(loo_comp2),
       elpd_diff = apply(loo_comp2, 2, sum),
       se_diff = apply(loo_comp2, 2, \(x) sd(x) * sqrt(length(x)))) |>
  mutate(p = ifelse(elpd_diff==0,NA,pnorm(0, -elpd_diff, se_diff))) |>
  tt() |>
  format_tt(replace = "-")

#' 
#' Based on model checking and the distribution of pointwise
#' $\elpdHati{\Mk}{\yobs}{i}$, the models seem to be reasonably
#' specified and we are fine with respect to Scenario 2 (model
#' misspecification) . Models $\mathrm{M}_2$ and $\mathrm{M}_3$ have
#' very small $\elpdHat{\Md}{\yobs}$ compared to model
#' $\mathrm{M}_1$. The direct use of the normal approximation gives
#' probabilities 0.16 and 0.6 that these models have better predictive
#' performance than model $\mathrm{M}_1$. As $\elpdHat{\Md}{\yobs}$ is
#' small (Scenario 1) and the number of observations is small
#' (Scenario 3), we may assume $\seHat{\Md}{\yobs}$ to be
#' underestimated and the error distribution to be more skewed than
#' normal. However, since $\elpdHat{\Md}{\yobs}$ is small, we can
#' state that there is no practical or statistical difference in the
#' predictive performance.
#'
#' The direct use of $\elpdHat{{\mathrm{M}_4,\mathrm{M}_1}}{\yobs}$
#' and $\seHat{{\mathrm{M}_4,\mathrm{M}_1}}{\yobs}$ would give
#' probability $0.96$ that model $\mathrm{M}_4$ has better predictive
#' than model $\mathrm{M}_1$. This difference (4.2) is big enough that
#' we are fine with respect to Scenario 1, but the number of
#' observations is small (Scenario 3), and on expectation we may
#' assume $\seHat{{\mathrm{M}_4,\mathrm{M}_1}}{\yobs}$ to be
#' underestimated. If we multiply
#' $\seHat{{\mathrm{M}_4,\mathrm{M}_1}}{\yobs}$ by 2 (heuristic
#' based on the limit of equations by @Bengio-Grandvalet:2004) to
#' make a more conservative estimate, the probability that model
#' $\mathrm{M}_4$ has better predictive performance is bigger than
#' 0.81. Considering we have only 17 observations, this is quite
#' good. Collecting more data is, however, recommended.
#'
#' As the predictive distribution includes the aleatoric uncertainty
#' (modelled by the data model), there is often more uncertainty in
#' the predictive performance model comparison than in the posterior
#' distribution (see, e.g., @Wang-Gelman:2015). In simple
#' models, we can also look at the posterior for the quantities of
#' interest. With model $\mathrm{M}_4$, $95\%$ central posterior
#' intervals for $\beta_1$ and $\beta_2$ are $(1.1,3.7)$ and
#' $(-0.12,-0.04)$ respectively, which indicates data have information
#' about the parameters. The covariates neocortex and log(mass) are
#' collinear, which causes correlation in the posterior of the
#' coefficients, which could make the marginal posteriors overlap 0,
#' even if the joint posterior does not, in which case, looking at the
#' predictive performance is useful. In this case, although neocortex
#' and log(mass) are collinear, they don't have useful information
#' alone, and the useful predictive information is along the second
#' principal component of their joint distribution, which explains why
#' the models with only one of the covariates are not better than the
#' intercept-only model.
#'
#' # Sleep study
#'
#' @Belenky-Wesensten-Thorne-etal:2003] collected data on the effect
#' of chronic sleep restriction. We use a subset of data in the R
#' package `lme4` [@Bates-Maechler-Bolker-etal:2015]. The data
#' contains average reaction times (in milliseconds) for 18 subjects
#' with sleep restricted to 3 hours per night for 7 consecutive
#' nights (days 0 and 1 were adaptation and training and removed
#' from this analysis).
data(sleepstudy, package="lme4")
sleepstudy2 <- sleepstudy |>
  filter(Days >= 2)

#' The compared models are a linear model, a linear model with varying
#' intercept for each subject, and a linear model with varying
#' intercept and slope for each subject. All models use a normal data
#' model. The models were fitted using `brms` [@Burkner:2017],
#' and the default `brms` priors; prior for the coefficient for
#' Days is uniform, the prior for the varying intercept is normal with
#' unknown scale having a half-normal prior, and the prior for the
#' varying intercept and slope is bivariate normal with unknown scales
#' having half-normal priors and correlation having LKJ prior
#' [@Lewandowski-Kurowicka-Joe:2009].
#'
#' Using the `brms` formula notation, the compared models are
#' \begin{align*}
#'   \mathrm{M}_1: \quad & \mathrm{Reaction} \sim \mathrm{Days} \\
#'   \mathrm{M}_2: \quad & \mathrm{Reaction} \sim \mathrm{Days} + (1\,\mid\,\mathrm{Subject}) \\
#'   \mathrm{M}_3: \quad & \mathrm{Reaction} \sim \mathrm{Days} + (\mathrm{Days}\,\mid\, \mathrm{Subject}).
#' \end{align*}
#'
#' Based on the study design, $\mathrm{M}_3$ is the appropriate model
#' for the analysis, but comparing models is useful for assessing how
#' much information the data has about the varying intercepts and
#' slopes. For a few LOO-folds with high Pareto-$\hat{k}$ diagnostic
#' value ($>0.7$, @Vehtari+etal:2024:PSIS) we re-ran MCMC (with
#' `reloo=TRUE` in `brms`). We use `add_criterion()` to store the
#' loo object inside the brmsfit objects.
#'
#| label: sleep_models
#| results: hide
#| cache: true
M_1 <- brm(Reaction ~ Days,
           data = sleepstudy2,
           family = gaussian(),
           refresh=0) |>
  add_criterion(criterion="loo", save_psis=TRUE, reloo=TRUE)
M_2 <- brm(Reaction ~ Days + (1 | Subject),
           data = sleepstudy2,
           family = gaussian(),
           refresh=0) |>
  add_criterion(criterion="loo", save_psis=TRUE, reloo=TRUE)
M_3 <- brm(Reaction ~ Days + (Days | Subject),
           data = sleepstudy2,
           family = gaussian(),
           refresh=0) |>
  add_criterion(criterion="loo", save_psis=TRUE, reloo=TRUE)

#' We compare the models $\mathrm{M}_1, \mathrm{M}_2,\mathrm{M}_3,\mathrm{M}_4$
loo_compare(M_1, M_2, M_3) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  dplyr::select(model, elpd_diff, se_diff) |>
  mutate(p = ifelse(elpd_diff==0,NA,pnorm(0, elpd_diff, se_diff))) |>
  tt() |>
  format_tt(j=4, replace = "-") |>
  format_tt(i=3, j=4, digits=5)

#' Model $\mathrm{M}_3$ is estimated to have better predictive
#' performance, but only with 0.9 probability of having better
#' performance than model $\mathrm{M}_2$. Model-checking reveals that
#' two observations are clear outliers with respect to these models,
#' making the normal approximation likely to be poorly calibrated
#' (Scenario 3).
#'
#' The difference is slightly smaller as can be expected as usually elpd estimates
#' with high Pareto-k's are optimistic. Now the probability that model $\mathrm{M}_3$
#' has better predictive performance than $\mathrm{M}_2$ is
#' `r lc=loo_compare(M_2, M_3);round(pnorm(0,lc[2,"elpd_diff"],lc[2,"se_diff"]),2)`.
#'
#' Although we can get more accurate `elpd` and `elpd_diff` estimates
#' with more computation, the accuracy of the normal approximation for
#' the uncertainty in the difference can still be compromised by the
#' outliers.  Model-checking reveals that two observations are clear
#' outliers with respect to these models, making the normal
#' approximation likely to be poorly calibrated (Scenario 3).  We see
#' the outliers, for example, by comparing LOO predictive intervals to
#' observations.
pp_check(M_2, type="loo_intervals") +
  labs(y="Reaction time (ms)") +
  bayesplot::theme_default(base_family = "sans", base_size = 18)

pp_check(M_3, type="loo_intervals") +
  labs(y="Reaction time (ms)") +
  bayesplot::theme_default(base_family = "sans", base_size = 18)

#' If we examine the pointwise differences, we see that model 3 is
#' almost always better, but outliers cause couple values with big
#' magnitude, making it more likely that the normal approximation for
#' quantifying uncertainty in elpd_diff is not accurate.
ggplot(data=NULL, aes(x=1:nrow(sleepstudy2),
           pointwise(M_3$criteria$loo, "elpd_loo")-pointwise(M_2$criteria$loo, "elpd_loo"))) +
  geom_point() +
  hline_0(alpha=0.5) +
  labs(x="Data index", y="pointwise elpd_diff")

#' We also fitted models using a Student's $t$ model to create models
#' $\mathrm{M}_{1t}$, $\mathrm{M}_{2t}$, and $\mathrm{M}_{3t}$. Based
#' on model checking, there is no obvious model misspecification.
#'
#| results: hide
#| cache: true
M_1t <- brm(Reaction ~ Days,
            data = sleepstudy2,
            family = student(),
            refresh=0) |>
  add_criterion(criterion="loo", save_psis=TRUE, reloo=TRUE)
M_2t <- brm(Reaction ~ Days + (1 | Subject),
            data = sleepstudy2,
            family = student(),
            refresh=0) |>
  add_criterion(criterion="loo", save_psis=TRUE, reloo=TRUE)
M_3t <- brm(Reaction ~ Days + (Days | Subject),
            data = sleepstudy2,
            family = student(),
            refresh=0) |>
  add_criterion(criterion="loo", save_psis=TRUE, reloo=TRUE)

#' LOO predictive intervals look now better.
pp_check(M_2t, type="loo_intervals") +
  labs(y="Reaction time (ms)") +
  bayesplot::theme_default(base_family = "sans", base_size = 18)

pp_check(M_3t, type="loo_intervals") +
  labs(y="Reaction time (ms)") +
  bayesplot::theme_default(base_family = "sans", base_size = 18)

#' We first compare $\mathrm{M}_{3}$ and $\mathrm{M}_{3t}$ to see
#' whether a Student's $t$ model is more appropriate.
loo_compare(M_3, M_3t) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  dplyr::select(model, elpd_diff, se_diff) |>
  mutate(p = ifelse(elpd_diff==0,NA,pnorm(0, elpd_diff, se_diff))) |>
  tt() |>
  format_tt(j=4, replace = "-", digits=3)

#' Although in this comparison $\mathrm{M}_3$ is misspecified, the
#' better specified model $\mathrm{M}_{3t}$ shows much better
#' predictive performance, and as we can expect $\seHatPlain$ to be
#' inflated, the actual probability that $\mathrm{M}_{3t}$ is better
#' than $\mathrm{M}_{3}$ is likely to be bigger than $0.999$. We then
#' compare the three Student's $t$ models:
#' 
loo_compare(M_1t, M_2t, M_3t) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  dplyr::select(model, elpd_diff, se_diff) |>
  mutate(p = ifelse(elpd_diff==0,NA,pnorm(0, elpd_diff, se_diff))) |>
  tt() |>
  format_tt(j=4, replace = "-")

#' The probability that model $\mathrm{M}_{3t}$ is better than models
#' $\mathrm{M}_{1t}$ and $\mathrm{M}_{2t}$ is close to 1. The models
#' appear sufficiently well specified, the number of observations is
#' bigger than 100, and the differences are not small, so we can
#' assume that the normal approximation is well calibrated.
#' 
#' If we examine the pointwise differences, we see that model 3 is
#' almost always better, and as there are no outliers, the
#' distribution of pointwise elpd differences is such that the normal
#' approximation for quantifying uncertainty in elpd_diff is likely to
#' be accurate.
ggplot(data=NULL, aes(x=1:nrow(sleepstudy2),
           pointwise(M_3t$criteria$loo, "elpd_loo")-pointwise(M_2t$criteria$loo, "elpd_loo"))) +
  geom_point() +
  hline_0(alpha=0.5) +
  labs(x="Data index", y="pointwise elpd_diff")

#' In this case, the effect of days with sleep constrained to 3 hours
#' is so big that the main conclusion stays the same with all the
#' models. Still, for example, model $\mathrm{M}_{3t}$ does indicate
#' higher variation between subjects than model $\mathrm{M}_{3}$. As
#' $\mathrm{M}_{3t}$ passes the model checking and has higher
#' predictive performance, we should continue looking at the posterior
#' of model $\mathrm{M}_{3t}$.
#'
#' # Roaches
#'
#' @Gelman-Hill:2007 describe in Chapter 8.3 the roaches data as
#' follows: *``the treatment and control were applied to 160 and 104
#' apartments, respectively, and the outcome measurement $y_i$ in
#' each apartment $i$ was the number of roaches caught in a set of
#' traps. Different apartments had traps for different numbers of
#' days''.* The goal is to estimate the efficacy of a pest
#' management system at reducing the number of roaches.
data(roaches)
# Roach1 is very skewed and we take a square root
roaches$sqrt_roach1 <- sqrt(roaches$roach1)

#' The target is the number of roaches (y), and the covariates include
#' the square root of the pre-treatment number of roaches
#' (sqrt\_roach1), a treatment indicator variable (treatment), and a
#' variable indicating whether the apartment is in a building
#' restricted to elderly residents (senior). As the number of days for
#' which the roach traps were used is not the same for all apartments,
#' the offset argument includes the logarithm of the number of days
#' the traps were used (log(exposure2)). The latent regression model
#' presented with `brms` formula notation is:
#' 
#' $$
#' \mathrm{y} \sim
#'  \mathrm{sqrt\_roach1} + \mathrm{treatment} + \mathrm{senior} + \mathrm{offset(log(exposure2))}.
#' $$
#' 
#' We fit the following models using the `brms` package. 
#' \begin{align*}
#'   \mathrm{M}_1: \quad & \text{Poisson} \\
#'   \mathrm{M}_2: \quad & \text{Negative.binomial}
#'   \mathrm{M}_3: \quad & \text{Zero-inflated negative-binomial} \\
#' \end{align*}
#'
#' The zero-inflation is modelled using the same latent formula (with
#' its own parameters). All coefficients have $\operatorname{normal}(0,1)$ priors
#' and the negative-binomial shape parameter has the `brms`
#' default prior, which is inverse-gamma$(.4, .3)$ [@Vehtari:2024].
#'
#' For the Poisson model we re-ran MCMC for all LOO-folds with high
#' Pareto-$\hat{k}$ diagnostic value (>0.7) (with `reloo=TRUE` in
#' `brms`), and for negative-binomial and zero-inflated
#' negative-binomial we used moment matching
#' [@Paananen+etal:2021:implicit} for a few LOO-folds with high
#' Pareto-$\hat{k}$ diagnostic value (>0.7) (with `moment_match=TRUE`
#' in `brms`).
#| label: roaches_models
#| results: hide
#| cache: true
M_1 <- brm(y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)),
           family=poisson(), data=roaches,
           prior=c(prior(normal(0,1), class='b')),
           seed=1704009, refresh=0) |>
  add_criterion(criterion='loo', save_psis=TRUE, reloo=TRUE)
M_2 <- brm(y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)),
           family=negbinomial(), data=roaches,
           prior=c(prior(normal(0,1), class='b'),
                   prior(inv_gamma(0.4, 0.3), class='shape')),
           seed=1704009, refresh=0) |>
  add_criterion(criterion='loo', save_psis=TRUE, moment_match=TRUE)
M_3 <- brm(bf(y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)),
              zi ~ sqrt_roach1 + treatment + senior + offset(log(exposure2))),
           family=zero_inflated_negbinomial(), data=roaches,
           prior=c(prior(normal(0,1), class='b'),
                   prior(normal(0,1), class='b', dpar='zi'),
                   prior(normal(0,1), class='Intercept', dpar='zi'),
                   prior(inv_gamma(0.4, 0.3), class='shape')),
           seed=1704009, refresh=0) |>
  add_criterion(criterion='loo', save_psis=TRUE, moment_match=TRUE)

#' 
loo_compare(M_1, M_2, M_3) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  dplyr::select(model, elpd_diff, se_diff) |>
  mutate(p = ifelse(elpd_diff==0,NA,pnorm(0, elpd_diff, se_diff))) |>
  tt() |>
  format_tt(j=4, replace = "-", digits=4)

#' The zero-inflated negative-binomial model ($\mathrm{M}_3$) is
#' clearly the best. Based on model checking, the Poisson model
#' ($\mathrm{M}_1$) is underdispersed which indicates Scenario 2
#' (model misspecification), but the difference is so big that we can
#' be certain that the zero-inflated negative-binomial model is
#' better. As the number of observations is larger than 100, and the
#' difference to model $\mathrm{M}_2$ is not small, we may assume the
#' normal approximation is well calibrated.
#'
#' As we had used an ad-hoc square root transformation of
#' pre-treatment number of roaches, we fitted a model $\mathrm{M}_4$
#' replacing the latent linear term for the square root of
#' pre-treatment number of roaches with a spline.
#' 
#| results: hide
#| cache: true
M_4 <- brm(bf(y ~ s(sqrt_roach1) + treatment + senior + offset(log(exposure2)),
              zi ~ s(sqrt_roach1) + treatment + senior + offset(log(exposure2))),
           family=zero_inflated_negbinomial(), data=roaches,
           prior=c(prior(normal(0,1), class='b'),
                   prior(normal(0,1), class='b', dpar='zi'),
                   prior(normal(0,1), class='Intercept', dpar='zi'),
                   prior(inv_gamma(0.4, 0.3), class='shape')),
           save_pars = save_pars(all = TRUE),
           seed=1704009) |> 
  add_criterion(criterion='loo', save_psis=TRUE, moment_match=TRUE)

#' 
loo_compare(M_3, M_4) |>
  as.data.frame() |>
  rownames_to_column("model") |>
  dplyr::select(model, elpd_diff, se_diff) |>
  mutate(p = ifelse(elpd_diff==0,NA,pnorm(0, elpd_diff, se_diff))) |>
  tt() |>
  format_tt(j=4, replace = "-")


#' Model $\mathrm{M}_4$ (with spline) seems to be slightly better, but
#' now the difference is so small (Scenario 1) that the normal
#' approximation is likely to be not perfectly calibrated. As the
#' difference is small, we can proceed with either model.
#' 
#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2025--206, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2025--206, Aki Vehtari, licensed under CC-BY-NC 4.0.

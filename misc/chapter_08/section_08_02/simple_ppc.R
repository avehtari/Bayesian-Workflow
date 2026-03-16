#' ---
#' title: "Kilpisjärvi PIT demo for Bayesian Workflow book"
#' author: "Martin Modrák"
#' date: 2020-10-07
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
#' This notebook includes code for the Bayesian Workflow book Figure
#' 8.5 in Section 8.2 about posterior predictive checking.

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
library(dplyr)
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 4)
library(bayesplot)
library(ggplot2)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(patchwork)

#' # Simple posterior predictive checking example
set.seed(85524)
data_lnorm <- data.frame(y = rlnorm(15, 3, 1))
#| results: hide
fit_lnorm <- brm(y ~ 1,
                 family = gaussian(),
                 data = data_lnorm,
                 refresh = 0)
#' 
pp_lnorm <- pp_check(fit_lnorm, type = "dens_overlay", nsamples = 20)
pp_lnorm


set.seed(197854312)
data_betabinom <- data.frame(y = rbinom(15, 7, prob = rbeta(15, 3, 1)))
#| results: hide
fit_betabinom <- brm(y | trials(7) ~ 1,
                     family = binomial(),
                     data = data_betabinom,
                     refresh = 0)
#' 
pp_betabinom <- pp_check(fit_betabinom, type = "stat", stat = "sd", binwidth = 0.1)
pp_betabinom

set.seed(32148234)
data_groups <- data.frame(group = c("Low","High")[rbinom(70, 1, 0.5) + 1]) |>
  mutate(prob = if_else(group == "Low", 0.4, 0.6), y = rbinom(n(), size = 6, prob = prob))
#| results: hide
fit_groups <- brm(y | trials(6) ~ 1,
                  family = binomial(),
                  data = data_groups,
                  refresh = 0)
#' 
preds <- posterior_predict(fit_groups, nsamples = 1000)
preds[preds > max(data_groups$y)] <- max(data_groups$y) + 1
pp_groups_whole <- ppc_bars(data_groups$y, preds)
pp_groups_whole
pp_groups <- ppc_bars_grouped(data_groups$y, preds, group = data_groups$group)
pp_groups

#| label: fig-posterior_predictive_simple
#| fig-width: 7
#| fig-height: 4
pp_all <- ((pp_lnorm + labs(tag = "A")) + (pp_groups_whole  + labs(tag = "C")) + (pp_betabinom  + labs(tag = "B")) + (patchwork::free(pp_groups +  labs(tag = "D")))) + plot_layout(nrow = 2) 
pp_all
#ggsave("posterior_predictive_simple.pdf", pp_all, width = 7, height = 4)

#' ---
#' title: "Mixture model for time series competition"
#' author: "Andrew Gelman"
#' date: 2022-08-18
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
#' This notebook includes the code for the Bayesian Workflow book
#' Chapter 20 *Using a fitted model for decision analysis: Mixture
#' model for time series competition*.
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
library(scales)
library(arm)
library(posterior)
library(cmdstanr)
options(mc.cores = 4)
set.seed(123)

#' # Data
series <- matrix(
  scan(root("timeseries", "data", "Series1000.txt")),
  nrow = 1000,
  ncol = 135,
  byrow = TRUE
)
T <- 135
N <- 1000

#| label: fig-series_1
#| fig-height: 3.4
#| fig-width: 5.5
par(mar = c(3, 3, 2, 0), tck = -.01, mgp = c(1.5, .5, 0))
plot(c(1, T), range(series), 
     xaxs = "i", bty = "l", type = "n", 
     xlab = "Time", ylab = "y")
for (n in 1:N){
  lines(1:T, series[n,], lwd = .5, col = alpha("black", 0.2))
}

#' # Maximum likelihood fits
slope <- rep(NA, N)
se <- rep(NA, N)
for (n in 1:N){
  data <- series[n,]
  time <- 1:T
  fit <- lm(data ~ time)
  slope[n] <- 100*coef(fit)["time"]
  se[n] <- 100*se.coef(fit)["time"]
}

#| label: fig-series_2
#| fig-height: 3.6
#| fig-width: 5.5
par(mar = c(3, 3, 2, 0), tck = -.01, mgp = c(1.5, .5, 0))
plot(slope, se, 
     ylim = c(0, 1.05 * max(se)), yaxs = "i", bty = "l", 
     xlab = "Estimated slope", ylab = "SE", pch = 20, cex = .5)

#| label: fig-series_3
#| fig-height: 3.6
#| fig-width: 5.5
par(mar = c(3, 3, 2, 0), tck = -.01, mgp = c(1.5, .5, 0))
hist(
  slope,
  xlab = "Estimated slope",
  breaks = seq(floor(10 * min(slope)), ceiling(10 * max(slope))) / 10,
  main = ""
)

#' # Stan mixture model 1
y <- slope
K <- 3
mu <- c(0, -1, 1)
data <- list(y = y, K = K, N = N, mu = mu)
mod <- cmdstan_model(root("timeseries", "mixture.stan"))
mod
#| label: fit_mix
#| results: hide
fit_mix <- mod$sample(data = data, refresh = 0)
#'
print(fit_mix)

#' # Stan mixture model 2
mod2 <- cmdstan_model(root("timeseries", "mixture_2.stan"))
mod2
#| label: fit_mix2
#| results: hide
fit_mix_2 <- mod2$sample(data = data, refresh = 0)
#'

prob_sims <- as_draws_rvars(fit_mix_2$draws())
prob <- mean(prob_sims$p)
max_prob <- apply(prob, 1, max)
choice <- apply(prob, 1, which.max)
expected_correct <- sum(max_prob)
sd_correct <- sqrt(sum(max_prob * (1 - max_prob)))
print(table(choice))

pfround(c(expected_correct, sd_correct), 1)

pnorm(expected_correct, 900, sd_correct)
1/pnorm(expected_correct, 900, sd_correct)
1/pnorm(expected_correct, 899.5, sd_correct)

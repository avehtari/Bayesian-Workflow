#' ---
#' title: "Incremental development and testing: Black cat adoptions"
#' author: "Richard McElreath and Aki Vehtari"
#' date: 2025-12-27
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
#' Chapter 22 *Incremental development and testing: Black cat adoptions*.
#'
#' # Introduction
#'
#' Intro text
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
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(scales)
library(cmdstanr)
# CmdStanR output directory makes Quarto cache to work
dir.create(root("cat_adoptions", "stan_output"), showWarnings = FALSE)
options(cmdstanr_output_dir = root("cat_adoptions", "stan_output"))
library(posterior)
library(survival)
library(rethinking)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))
library(ggdist)
library(ggsurvfit)

#' Utility functions for sampling using cmdstanr and plotting
cstan <- function(f, data = list(), seed = 123, chains = 4) {
  model <- cmdstan_model(f)
  fit <- model$sample(
    data = data,
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    refresh = 0
  )
  return(fit)
}
dens <- function(x, adj = 0.5, norm.comp = FALSE, main = "",
                 show.HPDI = FALSE, add = FALSE, ...) {
  thed <- density(x, adjust = adj)
  if (add == FALSE) {
    plot(thed, main = main, ...)
  } else {
    lines(thed$x, thed$y, ...)
  }
}
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

#' # Data
#'
#' Cat adoptions data is available in rethinking package. We read the
#' data from github, so there is no need to install the rethinking
#' package.
urlfile <- "https://raw.githubusercontent.com/rmcelreath/rethinking/master/data/AustinCats.csv"
d <- read_delim(urlfile, delim = ";")
glimpse(d)
d <- d |> mutate(days = days_to_event,
                 adopted = ifelse(out_event == "Adoption", 1, 0),
                 color = ifelse(color == "Black", 1, 2))

#' Prepare data for Stan models
dat <- list(N = nrow(d),
            days = d$days,
            adopted = d$adopted,
            color = d$color)

#' Plot individual cats as lines
## command below makes the plot window with right aspect ratio
# rethinking::blank2(w=1.2)
#| label: fig-cats-data
n <- 100
idx <- sample(1:dat$N, size = n)
ymax <- max(dat$days[idx])
plot(NULL, xlim = c(0, ymax), ylim = c(1, n),
     xlab = "Days observed", ylab = "Cat")
for (i in 1:n) {
  j <- idx[i]
  cat_color <- ifelse(dat$color[j] == 1, "black", "orange")
  lines(c(0, dat$days[j]), c(i, i), lwd = 4, col = cat_color)
  if (dat$adopted[j] == 1) {
    points(dat$days[j], i, pch = 16, cex = 1.5, col = cat_color)
  }
}

#' ggplot version
#| label: fig-gg-cats-data
d[idx, ] |>
  ggplot(aes(y = seq_along(days), x = days, color = factor(color))) +
  geom_segment(aes(yend = seq_along(days), xend = 0), size = 1) +
  geom_point(aes(shape = factor(adopted)), size = 3) +
  scale_shape_manual(values = c(1, 16), labels = c("Other", "Adopted")) +
  scale_color_manual(values = c("1" = "black", "2" = "orange"),
                     labels = c("Black", "Other")) +
  labs(y = "Cat", x = "Days observed", color = "Color", shape = "Event") +
  coord_cartesian(expand = c(left = FALSE)) +
  theme(legend.position = "inside",
        legend.justification.inside = c(1, 0.5))

#' # Generative models
#'
#' How should we model these data? Think about how they were
#' generated.  We start with process of adoption and then add
#' observation (censoring) process.
cat_adopt <- function(day, prob) {
  if (day > 1000 ) return(day)
  if (runif(1) > prob) {
    # keep waiting...
    day <- cat_adopt(day + 1, prob)
  }
  day # adopted
}
sim_cats1 <- function(n = 10, p = c(0.1, 0.2)) {
  color <- rep(NA, n)
  days <- rep(NA, n)
  for (i in 1:n) {
    color[i] <- sample(c(1, 2), size = 1, replace = TRUE)
    days[i] <- cat_adopt(1, p[color[i]])
  }
  return(list(N = n, days = days, color = color, adopted = rep(1, n)))
}

# version using rgeom - just for comparison
#sim_cats1 <- function(n=1e3,p=c(0.1,0.2)) {
#    color <- sample(c(1,2),size=n,replace=TRUE)
#    days <- rgeom( n , p[color] ) + 1
#    return(list(N=n,days=days,color=color,adopted=rep(1,n)))
#}

#' Simulate using the generative model
synth_cats <- sim_cats1(1e3)

#' Plot empirical K-M curves
#| label: fig-cats-km-curves
sfit <- survfit(Surv(days, adopted) ~ color, data = synth_cats)
plot(sfit, lty = 1, lwd = 3, col = c("black", "orange"),
     xlab = "Days", ylab = "Proportion un-adopted")

#' Plot empirical K-M curves using ggplot
#| label: fig-gg-cats-km-curves
#| fig-height: 3.5
#| fig-width: 6
synth_cats |>
  survfit2(formula = Surv(days, adopted) ~ color, data = _) |> # ggsurvfit version
  tidy_survfit() |>
  ggplot(aes(x = time, y = estimate, color = strata)) +
  geom_step(linewidth = 1) +
  xlim(c(0, 50)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(values = c("1" = "black", "2" = "orange"),
                     labels = c("Black", "Other")) +
  labs(x = "Days", y = "Proportion un-adopted", color = "Color") +
  theme(legend.position = "none")


#' ## First Stan model
cat_code1 <- root("cat_adoptions", "adoptions_observed.stan")
#| output: asis
print_stan_file(cat_code1)

#' Prior predictive simulation
n <- 12
sim_prior <- replicate(n, rbeta(2, 1, 10))
#| label: fig-prior-predictive-1
# rethinking::blank2(w=1.2)
plot(NULL, xlab = "Days", ylab = "Proportion un-adopted",
     xlim = c(0, 50), ylim = c(0, 1))
mtext("Prior predictive distribution")
for (i in 1:n) {
  days_rep <- sim_cats1(n = 1e3, p = sim_prior[, i])
  xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
  lines(xfit, lwd = 2, col = c("black", "orange"))
}

#' Prior predictive simulation with ggplot
#| label: fig-gg-prior-predictive-1
#| fig-height: 3.5
#| fig-width: 6
lapply(1:n, \(i) sim_cats1(n = 1e3, p = sim_prior[, i]) |>
         as.data.frame() |> mutate(sim = i)) |>
  bind_rows() |>
  survfit2(formula = Surv(days, adopted) ~ color + sim, data = _) |> # ggsurvfit version
  tidy_survfit() |>
  mutate(color = str_split_i(strata, ", ", 1),
         sim = str_split_i(strata, ", ", 2)) |>
  ggplot(aes(x = time, y = estimate, color = color, group = interaction(color, sim))) +
  geom_step() +
  xlim(c(0, 50)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(values = c("1" = "black", "2" = "orange"),
                     labels = c("Black", "Other")) +
  labs(x = "Days", y = "Proportion un-adopted", color = "Color")

#' Test the first model code using simulated data
p <- c(0.1, 0.15)
sim_dat <- sim_cats1(n = 1000, p = p)
#| results: hide
#| cache: true
fit1s <- cstan(cat_code1, data = sim_dat)

#' Posterior summary
print(fit1s)

#' Posterior with simulated data
#| label: fig-post1-sim
post1s <- fit1s$draws(format = "df")
plot(density(post1s$`p[1]`), lwd = 3, xlab = "Probability of adoption",
     xlim = c(0.07, 0.2), main = "")
k <- density(post1s$`p[2]`)
lines(k$x, k$y, lwd = 3, col = "orange")
abline(v = p[1], lwd = 2)
abline(v = p[2], lwd = 2, col = "orange")

#' Posterior with simulated data with ggplot
#| label: fig-gg-post1-sim
#| fig-height: 3.5
#| fig-width: 6
post1s |>
  ggplot() +
  stat_slab(aes(x = `p[1]`), density = "unbounded", trim = FALSE, fill = NA, color = "black") +
  stat_slab(aes(x = `p[2]`), density = "unbounded", trim = FALSE, fill = NA, color = "orange") +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line.y = element_blank(), strip.text.y = element_blank()) +
  coord_cartesian(expand = FALSE) +
  xlim(c(0.07, 0.2)) +
  labs(x = "Probability of adoption", y = "") +
  geom_vline(xintercept = p, color = c("black", "orange"))

#' Sample from the posterior using the real data
#| results: hide
#| cache: true
fit1 <- cstan(cat_code1, data = dat)

#' Posterior summary
print(fit1)

#' Kaplan-Meier posterior simulations
#| label: fig-post1-km
post1 <- fit1$draws(format = "df")
# rethinking::blank2(w=1.2)
plot(NULL, xlab = "Days", ylab = "Proportion un-adopted", xlim = c(0, 50), ylim = c(0, 1))
mtext("Posterior predictive distribution (1000 cats)")
n <- 12
for (i in 1:n) {
  days_rep <- sim_cats1(n = 1e3, p = post1[i, c("p[1]", "p[2]")])
  xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
  lines(xfit, lwd = 2, col = alpha(c("black", "orange"), 0.5))
}

#' Posterior Kaplan-Meier with ggplot
#| label: fig-gg-post1-km
#| fig-height: 3.5
#| fig-width: 6
lapply(1:n, \(i) sim_cats1(n = 1e3, p = post1[i, c("p[1]", "p[2]")]) |>
         as.data.frame() |> mutate(sim = i)) |>
  bind_rows() |>
  survfit2(formula = Surv(days, adopted) ~ color + sim, data = _) |> # ggsurvfit version
  tidy_survfit() |>
  mutate(color = str_split_i(strata, ", ", 1),
         sim = str_split_i(strata, ", ", 2)) |>
  ggplot(aes(x = time, y = estimate, color = color, group = interaction(color, sim))) +
  geom_step(alpha = 0.5) +
  xlim(c(0, 50)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(values = c("1" = "black", "2" = "orange"),
                     labels = c("Black", "Other")) +
  labs(x = "Days", y = "Proportion un-adopted", color = "Color") +
  theme(legend.position = "none")

#' ## Add observation (censoring) model
#'
#' Simulate from the generative model
sim_cats2 <- function(n = 10, p = c(0.1, 0.2), cens = 50) {
  color <- rep(NA, n)
  days <- rep(NA, n)
  for (i in 1:n) {
    color[i] <- sample(c(1, 2), size = 1, replace = TRUE)
    days[i] <- cat_adopt(1, p[color[i]])
  }
  adopted <- ifelse(days < cens, 1, 0)
  days <- ifelse(adopted == 1, days, cens)
  return(list(N = n, days = days, color = color, adopted = adopted))
}

cat_code2 <- root("cat_adoptions", "adoptions_censored.stan")
#| output: asis
print_stan_file(cat_code2)

#' Test censoring model using simulated data
sim_dat <- sim_cats2(n = 1e3, p = c(0.01, 0.02))
#| results: hide
#| cache: true
fit2s <- cstan(cat_code2, data = sim_dat)

#' Posterior summary
print(fit2s)

#| label: fig-post2-sim
post2s <- fit2s$draws(format = "df")
dens(post2s$`p[1]`, lwd = 3, xlab = "Probability of adoption", xlim = c(0, 0.03), ylim = c(0, 600))
dens(post2s$`p[2]`, add = TRUE, lwd = 3, col = "orange")
abline(v = 0.01, lwd = 2); abline(v = 0.02, lwd = 2, col = "orange")

#' Posterior with simulated data with ggplot
#| label: fig-gg-post2-sim
#| fig-height: 3.5
#| fig-width: 6
post2s |>
  ggplot() +
  stat_slab(aes(x = `p[1]`), density = "unbounded", trim = FALSE, fill = NA, color = "black") +
  stat_slab(aes(x = `p[2]`), density = "unbounded", trim = FALSE, fill = NA, color = "orange") +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line.y = element_blank(), strip.text.y = element_blank()) +
  coord_cartesian(expand = FALSE) +
  xlim(c(0.005, 0.025)) +
  labs(x = "Probability of adoption", y = "") +
  geom_vline(xintercept = c(0.01, 0.02), color = c("black", "orange"))

#' Test previous model with new censored data
#| results: hide
#| cache: true
fit1s <- cstan(cat_code1, data = sim_dat)

print(fit1s)

#| label: fig-post1-sim2
post1s <- fit1s$draws(format = "df")
dens(post1s$`p[1]`, lwd = 3, xlab = "Probability of adoption",
     xlim = c(0, 0.06), ylim = c(0, 170))
dens(post1s$`p[2]`, add = TRUE, lwd = 3, col = "orange")
abline(v = 0.01, lwd = 2); abline(v = 0.02, lwd = 2, col = "orange")

#' Posterior with simulated data with ggplot
#| label: fig-gg-post1-sim2
#| fig-height: 3.5
#| fig-width: 6
post1s |>
  ggplot() +
  stat_slab(aes(x = `p[1]`), density = "unbounded", trim = FALSE, fill = NA, color = "black") +
  stat_slab(aes(x = `p[2]`), density = "unbounded", trim = FALSE, fill = NA, color = "orange") +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line.y = element_blank(), strip.text.y = element_blank()) +
  coord_cartesian(expand = FALSE) +
  xlim(c(0.008, 0.062)) +
  labs(x = "Probability of adoption", y = "") +
  geom_vline(xintercept = c(0.01, 0.02), color = c("black", "orange"))

#' Sample using real data
#| results: hide
#| cache: true
fit1 <- cstan(cat_code1, data = dat)
fit2 <- cstan(cat_code2, data = dat)

#' Kaplan-meier posterior simulations
#| label: fig-post1-post2-km
# rethinking::blank2(w=1.2)
post1 <- fit1$draws(format = "df")
post2 <- fit2$draws(format = "df")
plot(NULL, xlab = "Days", ylab = "Proportion un-adopted", xlim = c(0, 50), ylim = c(0, 1))
mtext("Posterior predictive distribution (1000 cats)")
n <- 12
# New estimates
for (i in 1:n) {
  days_rep <- sim_cats1(n = 1e3, p = post2[i, c("p[1]", "p[2]")])
  xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
  lines(xfit, lwd = 2, col = alpha(c("black", "orange"), 0.5))
}
# Add a few simulations from first model, to show impact of censoring
n <- 1
for (i in 1:n) {
  days_rep <- sim_cats1(n = 1e4, p = post1[i, c("p[1]", "p[2]")])
  xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
  lines(xfit, lwd = 4, col = alpha(c("black", "orange"), 0.5))
}

#| label: fig-gg-post2-km
#| fig-height: 3.5
#| fig-width: 6
n <- 12
lapply(1:n, \(i) sim_cats1(n = 1e3, p = post2[i, c("p[1]", "p[2]")]) |>
         as.data.frame() |> mutate(sim = i)) |>
  bind_rows() |>
  survfit2(formula = Surv(days, adopted) ~ color + sim, data = _) |> # ggsurvfit version
  tidy_survfit() |>
  mutate(color = str_split_i(strata, ", ", 1),
         sim = str_split_i(strata, ", ", 2)) |>
  ggplot(aes(x = time, y = estimate, color = color, group = interaction(color, sim))) +
  geom_step(alpha = 0.5) +
  xlim(c(0, 50)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_color_manual(values = c("1" = "black", "2" = "orange"),
                     labels = c("Black", "Other")) +
  labs(x = "Days", y = "Proportion un-adopted", color = "Color") +
  theme(legend.position = "none")
  

#' ## Model that uses parameters for censored observations
cat_code3 <- root("cat_adoptions", "adoptions_imputation.stan")
#| output: asis
print_stan_file(cat_code3)

#'
#| results: hide
#| cache: true
fit3s <- cstan(cat_code3, data = sim_dat)

#' Posterior summary
#| echo: false
oldo <- options(width = 90)
#| echo: true
print(fit3s)
#| echo: false
options(oldo)

#' ## Poisson model
#'
#' Model that uses Poisson count outcomes instead of duration outcomes
#' to handle censoring. This depends upon constant hazard function
#' though?
#'

cat_code4 <- root("cat_adoptions", "adoptions_poisson.stan")
#| output: asis
print_stan_file(cat_code4)

#| results: hide
#| cache: true
fit4s <- cstan(cat_code4, data = sim_dat)

#' Posterior summary
print(fit4s)

#' ## Varying effects model

#' Simulate background traits that differentiate cats of same color.
sim_cats3 <- function(n = 10, p = c(0.1, 0.2), cens = 50, xsd = c(0.1, 0.2)) {
  color <- rep(NA, n)
  days <- rep(NA, n)
  for (i in 1:n) {
    color[i] <- sample(c(1, 2), size = 1, replace = TRUE)
    z <- rnorm(1, 0, xsd[color[i]])
    pp <- inv_logit(logit(p[color[i]]) + z)
    days[i] <- cat_adopt(1, pp)
  }
  adopted <- ifelse(days < cens, 1, 0)
  days <- ifelse(adopted == 1, days, cens)
  return(list(N = n, days = days, color = color, adopted = adopted))
}
sim_dat <- sim_cats3(n = 1000, p = c(0.2, 0.1), xsd = c(0.1, 0.1))

#' Varying effects model
cat_code5 <- root("cat_adoptions", "adoptions_varying.stan")
#| output: asis
print_stan_file(cat_code5)

#| results: hide
#| cache: true
fit2s <- cstan(cat_code2, data = sim_dat)
fit5s <- cstan(cat_code5, data = sim_dat)


#' # Workflow

#' ## Prior predictive distributuon
#'
#' Repeatedly sample from prior, simulate observations

#' Prior draws
n <- 100
p1 <- rbeta(n, 1, 10)
p2 <- rbeta(n, 1, 10)
sim_cats2 <- function(n = 1e3, p = c(0.01, 0.02), cens = 50) {
  color <- sample(c(1, 2), size = n, replace = TRUE)
  days <- rgeom(n, p[color]) + 1
  adopted <- ifelse(days < cens, 1, 0)
  days <- ifelse(adopted == 1, days, cens)
  return(list(N = n, days = days, color = color, adopted = adopted))
}
prior_days <- sapply(1:n, function(i) sim_cats2(1, p = c(p1[i], p2[i]))$days)
plot(prior_days, xlab = "simulated cat", ylab = "days",
     pch = ifelse(prior_days == 50, 1, 16))

#' ## Posterior predictive distribution
#'
#' Sample from posterior, simulate observations.
#' Problem with this example: need to impute censored values
#' so we'll simulate Kaplan-Meier curves to compare to empirical curve.

#' Plot empirical K-M curves
#| label: fig-km-2
sfit <- survfit(Surv(days, adopted) ~ color, data = dat)
plot(sfit, lty = 1, lwd = 0.1, col = c("black", "orange"), xlim = c(0, 90),
     xlab = "Days", ylab = "Proportion un-adopted")

#' Simulate and draw
# n <- 12
# for (i in 1:n) {
#   days_rep <- sim_cats2(n = 1e3, p = post2$p[i, ], cens = 200)
#   xfit <- survfit(Surv(days, adopted) ~ color, data = days_rep)
#   lines(xfit, lwd = 1, col = alpha(c("black", "orange"), 0.5))
# }

#' Overlay empirical curves
#| label: fig-km-3
# lines(sfit, lwd = 5, col = c("white", "white"))
# lines(sfit, lwd = 3, col = c("black", "orange"))

#' # Licenses {.unnumbered}
#'
#' * Code &copy; 2025, Richard McElreath, licensed under BSD-3.
#' * Text &copy; 2025, Richard McElreath, licensed under CC-BY-NC 4.0.

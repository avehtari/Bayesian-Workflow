#' ---
#' title: "Bioassay case study"
#' author: "Andrew Gelman and Aki Vehtari"
#' date: 2025-10-15
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
#' This notebook includes the code for Bayesian Workflow book Section
#' 3.5 *A simple example of probabilistic programming*.
#'
#' # Introduction
#'
#' We demonstrate the role of statistical programming environments and
#' probabilistic programming by going through the steps of data
#' analysis and computation in the context of a simple example.
#'
#+ setup, include=FALSE
knitr::opts_chunk$set(message=FALSE, error=FALSE, warning=FALSE, comment=NA, cache=FALSE)

#' *Load packages*
#| cache: FALSE
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 4)
library(posterior)
library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 16))
library(marginaleffects)

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

#' # Bioassay data
#'
#' Read the data from csv
# df_bioassay <- read.csv("bioassay/data/bioassay.csv")
#'
#' Data as data frame
df_bioassay <- data.frame(
  dose = c(-0.86, -0.30, -0.05, 0.73),
  batch_size = c(5, 5, 5, 5),
  deaths = c(0, 1, 3, 5)
)

#' Data plotted with base R
#| label: fig-bioassay-data-plot
#| fig-height: 4
#| fig-width: 7
with(df_bioassay,
     plot(dose, deaths, xlab = "Dose log(g/ml)", ylab = "# of deaths",
          pch = 19, cex = 1.5, bty = "l"))

#' Data plotted with ggplot
#| label: fig-bioassay-data-ggplot
#| fig-height: 4
#| fig-width: 8
df_bioassay |>
  ggplot(aes(x = dose, y = deaths)) +
  geom_point(size = 3) +
    labs(x = "Dose log(g/ml)", y = "# of deaths")

#' Data as list for CmdStanR
bioassay_data <- with(df_bioassay,
                      list(J = nrow(df_bioassay),
                           x = dose,
                           N = batch_size,
                           y = deaths))

#' # Stan models and inference
#'
#' Stan model 0 (without priors)
bioassay_stan_file <- root("bioassay","bioassay0.stan")
#| output: asis
print_stan_file(bioassay_stan_file)

#' Compile the Stan model code using pedantic mode
#| label: mod0
mod0 <- cmdstan_model(bioassay_stan_file, pedantic = TRUE)

#' Stan model 1 (with priors)
bioassay_stan_file <- root("bioassay","bioassay1.stan")
#| output: asis
print_stan_file(bioassay_stan_file)

#' Compile the updated Stan model code using pedantic mode
#| label: mod1
mod1 <- cmdstan_model(bioassay_stan_file, pedantic = TRUE)

#' Sample with default settings (expect no progress output)
#| label: fit1
#| results: hide
fit1 <- mod1$sample(data = bioassay_data, refresh = 0)

#' # Posterior summary
#'
#' Posterior summary and MCMC diagnostics
fit1

#' Posterior draws as data frame
draws1 <- fit1$draws(format="df")

#' Plot with base plot: data, 20 logistic curves given 20 posterior
#' draws, and one logistic given the posterior mean
#| label: fig-bioassay-posterior-plot
#| fig-height: 3
#| fig-width: 5
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
with(df_bioassay,
     plot(dose, deaths/batch_size, xlab = "Dose log(g/ml)", ylab = "Pr (death)",
          pch=19, cex=1.5, bty="l"))
invlogit <- plogis
for (s in sample(nrow(draws1), 20)) {
  curve(invlogit(draws1$a[s] + draws1$b[s] * x),
        col = "red", lwd = 0.5, add = TRUE)
}
curve(invlogit(mean(draws1$a) + mean(draws1$b) * x),
      col = "blue", lwd = 2, add = TRUE)

#' Plot with ggplot: data, 20 logistic curves given 20 posterior
#' draws, and one logistic given the posterior mean
#| label: fig-bioassay-posterior-ggplot
#| fig-height: 4
#| fig-width: 8
draws1 |>
  resample_draws(ndraws = 20) |>
  expand_grid(x=seq(-1, 1, length = 100)) |>
  mutate(y = plogis(a + b * x)) |>
  ggplot() +
  geom_point(data = df_bioassay, aes(x = dose, y = deaths / batch_size), size = 3) +
  geom_line(aes(x = x, y = y, group = .draw), alpha = .5, color = "red") +
  geom_function(fun = \(x) plogis(mean(draws1$a) + mean(draws1$b)*x),
                color = "blue",
                linewidth = 1) +
  labs(x = "Dose log(g/ml)", y = "Pr(death)")

#' # Derived quantities
#'
#' Compute posterior draws for LD50 in log(g/ml) and in mg/ml
draws1 <- draws1 |>
  mutate_variables(LD50_log_g_ml = -a / b,
                   LD50_mg_ml = 1000 * exp(LD50_log_g_ml))

#' Alternatively using base R (as posterior package draws object is
#' more than just plain data frame, it is better to use
#' `mutate_variables()`, which will also work for all other posterior
#' draws objects (`draws_list`, `draws_array`, `draws_rvars`))
#'
#| eval: false
draws1$LD50_log_g_ml <- -draws1$a /  draws1$b
draws1$LD50_mg_ml <- 1000 * exp(draws1$LD50_log_g_ml)

#' Summarise LD50 posterior
draws1 |>
  subset_draws(variable = "LD50_log_g_ml") |>
  summarize_draws()

#' Common part for the next two LD50 plots (without data)
ggplot_LD50_mg_ml <-
  ggplot(mapping = aes(x = LD50_mg_ml)) +
  scale_x_log10(limits = c(375,2700)) +
  labs(x = "LD50 mg/ml") +
  geom_vline(xintercept = c(500, 2000), alpha=0.1) +
  annotate(geom = "text", x = 500*0.96, y = .97, hjust = 1, label = "Category 3") +
  annotate(geom = "text", x = 1000, y = .97, label = "Category 4") +
  annotate(geom = "text", x = 2000*1.04, y = .97, hjust = 0, label = "Category 5") +
  theme_sub_axis_y(text = element_blank(),
                                        line = element_blank(),
                                        ticks = element_blank(),
                                        title = element_blank())

#' Quantile dot plot of the LD50 (mg/ml) posterior
#| label: fig-bioassay-LD50-quantile-dots
#| fig-height: 4
#| fig-width: 8
ggplot_LD50_mg_ml +
  stat_dots(data = draws1, quantiles = 100) +
  coord_cartesian(expand = c(bottom = FALSE))

#' Kernel density plot of the LD50 (mg/ml) posterior
#| label: fig-bioassay-LD50-KDE
#| fig-height: 4
#| fig-width: 8
ggplot_LD50_mg_ml +
  stat_slab(data = draws1, color = "gray", fill = NA) +
  coord_cartesian(expand = c(bottom = FALSE))

#' # brms model and inference
#'
#' The same model with `brms`
#| label: bfit1
#| results: hide
bfit1 <- brm(deaths | trials(batch_size) ~ dose,
             family = binomial(),
             prior = c(prior(normal(0, 5), class = Intercept),
                       prior(normal(0, 5), lb = 0, class = b)),
             data = df_bioassay,
             refresh = 0)

#' Summary of the inference
bfit1

#' Posterior draws as data frame
bdraws1 <- as_draws_df(bfit1)

#' Plot data, 20 logistic curves given 20 posterior draws, and one
#' logistic given the posterior mean
#| label: fig-bioassay-posterior-ggplot-brms
#| fig-height: 4
#| fig-width: 8
bdraws1 |>
  resample_draws(ndraws = 20) |>
  expand_grid(x=seq(-1, 1, length = 100)) |>
  mutate(y = plogis(b_Intercept + b_dose * x)) |>
  ggplot() +
  geom_line(aes(x = x, y = y, group = .draw), alpha = .5, color = "red") +
  geom_point(data = df_bioassay, aes(x = dose, y = deaths / batch_size), size = 3) +
  geom_function(fun = \(x) plogis(mean(bdraws1$b_Intercept) + mean(bdraws1$b_dose)*x),
                color = "blue",
                linewidth = 1) +
  labs(x = "Dose log(g/ml)", y = "Pr(death)")

#' `brms` provides also shortcut for plotting the posterior mean and
#' uncertainty. We use `plot=FALSE` and `[[1]]` to return a ggplot
#' object, so that we can modify the axes labels.
#| label: fig-bioassay-posterior-conditional_effects-brms
#| fig-height: 4
#| fig-width: 8
p1 <- plot(conditional_effects(bfit1), plot=FALSE)[[1]] +
  labs(x = "Dose log(g/ml)", y = "Pr(death)")
p1

#' We can add data points with `ggplot::geom_point()`.
#| label: fig-bioassay-posterior-conditional_effects-data-brms
#| fig-height: 4
#| fig-width: 8
p1 +
  geom_point(data = df_bioassay,
             inherit.aes = FALSE,
             aes(x = dose, y = deaths / batch_size),
             size = 3)

#' We can draw also individual posterior curves with `spaghetti=TRUE,
#' ndraws=20`, but then the posterior mean would be computed only from
#' 20 curves. Increasing `ndraws` would make the plotted posterior
#' mean more accurate, but would make the spaghetti plot more messy.
#| label: fig-bioassay-posterior-conditional_effects-spaghetti-brms
#| fig-height: 4
#| fig-width: 8
p1 <- plot(conditional_effects(bfit1, spaghetti=TRUE, ndraws=20),
           plot=FALSE)[[1]] +
  labs(x = "Dose log(g/ml)", y = "Pr(death)")
p1

#' `marginaleffects` package also provides function for plotting the
#' posterior prediction and uncertainty. By default it predicts on the
#' scale of actual outcome, but as all batch sizes are equal we can
#' transform to the interval from 0 to 1. We need to add the observed
#' proportions of deaths with `geom_point()`.
#| label: fig-bioassay-posterior-marginaleffects-brms
#| fig-height: 4
#| fig-width: 8
p2 <- marginaleffects::plot_predictions(bfit1,
                                        condition = "dose",
                                        transform = \(x) x/5) +
  labs(x = "Dose log(g/ml)", y = "Pr(death)") +
  geom_point(data = df_bioassay,
             inherit.aes = FALSE,
             aes(x = dose, y = deaths / batch_size),
             size = 3)
p2

#' The `brms` and `marginaleffects` plotting functions demonstrate that
#' often shortcut functions can provide quick plot or summary, but to
#' have more control you may still need write more code to get what
#' you really want.
#'
#' # LD50 posterior
#'
#' Compute posterior draws for lethal dose 50 (LD50) in log(g/ml) and
#' in mg/ml
bdraws1 <- bdraws1 |>
  mutate_variables(LD50_log_g_ml = -b_Intercept/b_dose,
                   LD50_mg_ml = 1000 * exp(LD50_log_g_ml))

#' Posterior summary
bdraws1 |>
  subset_draws(variable = "LD50_log_g_ml") |>
  summarize_draws()

#' Quantile dot plot of the LD50 (mg/ml) posterior.
#| label: fig-bioassay-LD50-quantile-dots-brms
#| fig-height: 4
#| fig-width: 8
ggplot_LD50_mg_ml +
  stat_dots(data = bdraws1, quantiles = 100) +
  coord_cartesian(expand = c(bottom = FALSE))

#' Kernel density plot of the LD50 (mg/ml) posterior.
#| label: fig-bioassay-LD50-KDE-brms
#| fig-height: 4
#| fig-width: 8
ggplot_LD50_mg_ml +
  stat_slab(data = bdraws1, color = "gray", fill = NA) +
  coord_cartesian(expand = c(bottom = FALSE))

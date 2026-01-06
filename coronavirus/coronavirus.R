#' ---
#' title: "Sensitivity and specificity in coronavirus testing"
#' author: "Andrew Gelman and Aki Vehtari"
#' date: 2022-08-22
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
#' Chapter 19 *Building up to a hierarchical model: Coronavirus
#' testing*.
#'
#' # Introduction
#'
#' In early April, 2020, @Bendavid-Mulaney-Sood-etal:2020a recruited
#' 3330 residents of Santa Clara County, California and tested them
#' for SARS-CoV-2 antibodies. We re-analyze the data.
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
library(cmdstanr)
options(mc.cores = 4)
library(posterior)
library(priorsense)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))
library(dplyr)
library(ggh4x)

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

#' Define function to compute shortest posterior interval (SPIN; @Liu-Gelman-Zheng:2015)
spin <- function(x, lower=NULL, upper=NULL, conf=0.95) {
  x <- sort(as.vector(x))
  if (!is.null(lower)) {
    if (lower > min(x)) stop("lower bound is not lower than all the data")
    else x <- c(lower, x)
  }
  if (!is.null(upper)) {
    if (upper < max(x)) stop("upper bound is not higher than all the data")
    else x <- c(x, upper)
  }
  n <- length(x)
  gap <- round(conf*n)
  width <- x[(gap+1):n] - x[1:(n-gap)]
  index <- min(which(width==min(width)))
  x[c(index, index + gap)]
}

#' # Simple model fit using pooled specificity and sensitivity

code <- root("coronavirus", "santa-clara.stan")
#| output: asis
print_stan_file(code)
#'
sc_model <- cmdstan_model(code)

#' Compute posterior with data from @Bendavid-Mulaney-Sood-etal:2020a, 11 April 2020
#| label: fit_1
#| results: hide
fit_1 <- sc_model$sample(
  data = list(
    y_sample = 50,
    n_sample = 3330,
    y_spec = 369 + 30,
    n_spec = 371 + 30,
    y_sens = 25 + 78,
    n_sens = 37 + 85
  ),
  refresh = 0,
  iter_warmup = 1e4,
  iter_sampling = 1e4
)
print(fit_1, digits = 3)

#' Check prior-likelihood sensitivity using powerscaling
powerscale_sensitivity(fit_1)
#| label: fig-specificity_priorsense_1
#| fig-height: 6.5
#| fig-width: 12
powerscale_plot_dens(
  fit_1,
  variables = c("p", "spec", "sens"),
  lower_alpha = 0.5,
  help_text = FALSE
)

#' Inference for the population prevalence
draws_1 <- fit_1$draws()
subset <- sample(1e4, 1e3)
x <- as.vector(draws_1[subset, , "spec"])
y <- as.vector(draws_1[subset, , "p"])
#| label: fig-scatter
#| fig-height: 3.5
#| fig-width: 4.5
#| out-width: 70%
par(mar = c(3, 3, 0, 1), mgp = c(2, .7, 0), tck = -.02)
plot(x, y, 
     xlim = c(min(x), 1), ylim = c(0, max(y)), 
     xaxs = "i", yaxs = "i", 
     xlab = expression(paste("Specificity, ", gamma)), 
     ylab = expression(paste("Prevalence, ", pi)), bty = "l", 
     pch=20, cex=.3)

#' ggplot version
#| label: fig-gg-scatter
#| fig-height: 3.5
#| fig-width: 4.5
#| out-width: 70%
fit_1$draws() |>
  resample_draws(ndraws = 1e3) |>
  mcmc_scatter(pars = c("spec", "p"), size = 1, alpha = 0.5) +
  lims(x = c(NA, 1), y = c(0, NA)) +
  coord_cartesian(expand = c(bottom = FALSE)) +
  labs(x = expression(paste("Specificity, ", gamma)), 
       y = expression(paste("Prevalence, ", pi)))

#| label: fig-hist
#| fig-height: 3.5
#| fig-width: 5.5
#| out-width: 70%
par(mar = c(3, 3, 0, 1), mgp = c(2, .7, 0), tck = -.02)
hist(y, 
     yaxt = "n", yaxs = "i", 
     xlab = expression(paste("Prevalence, ", pi)), 
     ylab = "", main = "")

#| label: fig-gg-hist
#| fig-height: 3.5
#| fig-width: 5.5
#| out-width: 70%
fit_1$draws() |>
  mcmc_hist(pars = c("p"), breaks = seq(0, 0.03, length.out = 30)) +
  theme(axis.line.y = element_blank()) + 
  labs(x = expression(paste("Prevalence, ", pi)))

#' Use the shortest posterior interval, which makes more sense than a central 
#' interval because of the skewness of the posterior and the hard boundary at 0.
print(spin(draws_1[, , "p"], lower = 0, upper = 1, conf = 0.95), digits = 2)

#' Compute posterior with data from @Bendavid-Mulaney-Sood-etal:2020b, 27 April 2020
#'
#| label: fit_2
#| results: hide
fit_2 <- sc_model$sample(
  data = list(
    y_sample = 50,
    n_sample = 3330,
    y_spec = 3308,
    n_spec = 3324,
    y_sens = 130,
    n_sens = 157
  ),
  refresh = 0,
  iter_warmup = 1e4,
  iter_sampling = 1e4
)

print(fit_2, digits = 3)
draws_2 <- fit_2$draws()
print(spin(draws_2[, , "p"], lower = 0, upper = 1, conf = 0.95), digits = 2)


#' Check prior-likelihood sensitivity using powerscaling
powerscale_sensitivity(fit_2)
#| label: fig-specificity_priorsense_x
powerscale_plot_dens(
  fit_2,
  variables = c("p", "spec", "sens"),
  lower_alpha = 0.5,
  help_text = FALSE
)

#' # Hierarchical model for sensitivity and specificity
#'
#' Hierarchical model that allows sensitivity and specificity to vary across studies.
code_hierarchical <- root("coronavirus", "santa-clara-hierarchical.stan")
#| output: asis
print_stan_file(code_hierarchical)
#'
sc_model_hierarchical <- cmdstan_model(code_hierarchical)

#' Compute posterior using data from @Bendavid-Mulaney-Sood-etal:2020b 27 April 2020
#| label: fit_3a
#| results: hide
santaclara_data <- list(
  y_sample = 50,
  n_sample = 3330,
  J_spec = 14,
  y_spec = c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50),
  n_spec = c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52),
  J_sens = 4,
  y_sens = c(0, 78, 27, 25),
  n_sens = c(0, 85, 37, 35),
  logit_spec_prior_scale = 1,
  logit_sens_prior_scale = 1
)
fit_3a <- sc_model_hierarchical$sample(
  data = santaclara_data,
  refresh = 0,
  iter_warmup = 1e4,
  iter_sampling = 1e4, 
  adapt_delta = 0.99 # to reduce divergences 
)

print_variables <- c("p", "spec[1]", "sens[1]", 
                     "mu_logit_spec", "mu_logit_sens", 
                     "sigma_logit_spec", "sigma_logit_sens")
print(fit_3a, variables = print_variables, digits = 3)

draws_3a <- fit_3a$draws()
print(spin(draws_3a[, , "p"], lower = 0, upper = 1, conf = 0.95), digits = 2)

#' Check prior-likelihood sensitivity using powerscaling
powerscale_sensitivity(fit_3a, variable = c("p", "sens", "spec")) |>
  print(n = Inf)
#| label: fig-specificity_priorsense_2
#| fig-height: 6.5
#| fig-width: 12
powerscale_plot_dens(
  fit_3a,
  variables = c("p", "spec[1]", "sens[1]"),
  lower_alpha = 0.5,
  help_text = FALSE
) +
  ggh4x::facetted_pos_scales(x=rep(list(scale_x_continuous(limits=c(0, 0.35)),
                                        scale_x_continuous(limits=c(0, 1)),
                                        scale_x_continuous(limits=c(0.98, 1))), 2))

#' # Hierarchical model with stronger priors
#' 
#' Compute posterior using data from @Bendavid-Mulaney-Sood-etal:2020b
#' 27 April 2020 and stronger priors
#' 
santaclara_data$logit_spec_prior_scale <- 0.3
santaclara_data$logit_sens_prior_scale <- 0.3
#| label: fit_3b
#| results: hide
#' MCMC gets sometimes stuck on minor mode, so we initialize with Pathfinder
pth_3b <- sc_model_hierarchical$pathfinder(
  data = santaclara_data,
  refresh = 0,
  max_lbfgs_iters = 100,
  num_paths = 10
)
fit_3b <- sc_model_hierarchical$sample(
  data = santaclara_data,
  refresh = 0,
  iter_warmup = 1e4,
  iter_sampling = 1e4, 
  adapt_delta = 0.99, # not as necessary with the stronger priors, but include just in case
  init = pth_3b
)

print(fit_3b, variables = print_variables, digits = 3)

draws_3b <- fit_3b$draws()
print(spin(draws_3b[, , "p"], lower = 0, upper = 1, conf = 0.95), digits = 2)

print(spin(draws_3a[, , "p"], lower = 0, upper = 1, conf = 0.95), digits = 2)
print(spin(draws_3a[, , "spec[1]"], lower = 0, upper = 1, conf = 0.95), digits = 2)
print(spin(draws_3a[, , "sens[1]"], lower = 0, upper = 1, conf = 0.95), digits = 2)
print(spin(draws_3a[, , "mu_logit_spec"], conf = 0.95), digits = 2)
print(spin(draws_3a[, , "mu_logit_sens"], conf = 0.95), digits = 2)
print(spin(draws_3a[, , "sigma_logit_spec"], conf = 0.95), digits = 2)
print(spin(draws_3a[, , "sigma_logit_sens"], conf = 0.95), digits = 2)

print(spin(draws_3b[, , "p"], lower = 0, upper = 1, conf = 0.95), digits = 2)
print(spin(draws_3b[, , "spec[1]"], lower = 0, upper = 1, conf = 0.95), digits = 2)
print(spin(draws_3b[, , "sens[1]"], lower = 0, upper = 1, conf = 0.95), digits = 2)
print(spin(draws_3b[, , "mu_logit_spec"], conf = 0.95), digits = 2)
print(spin(draws_3b[, , "mu_logit_sens"], conf = 0.95), digits = 2)
print(spin(draws_3b[, , "sigma_logit_spec"], conf = 0.95), digits = 2)
print(spin(draws_3b[, , "sigma_logit_sens"], conf = 0.95), digits = 2)


#' # MRP model
#'
#' MRP model allowing prevalence to vary by sex, ethnicity, age
#' category, and zip code.  Model is set up to use the ethnicity, age,
#' and zip categories of @Bendavid-Mulaney-Sood-etal:2020b
code_model_hierarchical_mrp <- root("coronavirus", "santa-clara-hierarchical-mrp.stan")
#| output: asis
print_stan_file(code_model_hierarchical_mrp)
#'
sc_model_hierarchical_mrp <- cmdstan_model(code_model_hierarchical_mrp)

#' To fit the model, we need individual-level data.  
#' These data are not publicly available, so just to get the program running, 
#' we take the existing 50 positive tests and assign them at random to the 3330 people.
N <- 3330
y <- sample(rep(c(0, 1), c(3330 - 50, 50)))
n <- rep(1, 3330)

#' Here are the counts of each sex, ethnicity, and age from @Bendavid-Mulaney-Sood-etal:2020b.  
#' We don't have zip code distribution but we looked it up and there are 58 zip codes 
#' in Santa Clara County; for simplicity we assume all zip codes are equally likely.  
#' We then assign these traits to people at random.  
#' This is wrong--actually, these variable are correlated in various ways--but, again,
#' now we have fake data we can use to fit the model.
male <- sample(rep(c(0, 1), c(2101, 1229)))
eth <- sample(rep(1:4, c(2118, 623, 266, 306 + 17)))
age <- sample(rep(1:4, c(71, 550, 2542, 167)))
N_zip <- 58
zip <- sample(1:N_zip, 3330, replace = TRUE)

#' Setting up the zip code level predictor.  
#' In this case we will use a random number with mean 50 and standard deviation 20.  
#' These are arbitrary numbers that we chose just to be able to test the centering 
#' and scaling in the model.   
#' In real life we might use %Latino or average income in the zip code
x_zip <- rnorm(N_zip, 50, 20)

#' Setting up the poststratification table.  
#' For simplicity we assume there are 1000 people in each cell in the county.  
#' Actually we'd want data from the Census.
J <- 2 * 4 * 4 * N_zip
N_pop <- rep(NA, J)
count <- 1
for (i_zip in 1:N_zip){
  for (i_age in 1:4){
    for (i_eth in 1:4){
      for (i_male in 0:1){
        N_pop[count] <- 1000
        count = count + 1
      }
    }
  }
}

#' Put together the data and fit the model
#| label: fit_4
#| results: hide
santaclara_mrp_data <- list(
  N = N,
  y = y,
  male = male,
  eth = eth,
  age = age,
  zip = zip,
  N_zip = N_zip,
  x_zip = x_zip,
  J_spec = 14,
  y_spec = c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50),
  n_spec = c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52),
  J_sens = 4,
  y_sens = c(0, 78, 27, 25),
  n_sens = c(0, 85, 37, 35),
  logit_spec_prior_scale = 0.3,
  logit_sens_prior_scale = 0.3,
  coef_prior_scale = 0.5,
  J = J,
  N_pop = N_pop
)
fit_4 <- sc_model_hierarchical_mrp$sample(data = santaclara_mrp_data)

#' Show inferences for some model parameters. In addition to `p_avg`, the population
#' prevalence, we also look at the inferences for the first three
#' poststratification cells just to check that everything makes sense
print_variables <- c("p_avg", "b", "a_age", "a_eth", 
                     "sigma_eth", "sigma_age", "sigma_zip", 
                     "mu_logit_spec", "sigma_logit_spec",  
                     "mu_logit_sens", "sigma_logit_sens", 
                     "p_pop[1]", "p_pop[2]", "p_pop[3]")
print(fit_4, variables = print_variables, max_rows = 100)

draws_4 <- fit_4$draws()
print(spin(draws_4[, , "p_avg"], lower = 0, upper = 1, conf = 0.95), digits = 2)

#' # Additional prior sensitivity analysis
pos_tests <- c(78, 27, 25)
tests <- c(85, 37, 35)
sens_df <- data.frame(pos_tests, tests, sample_sens = pos_tests / tests)

neg_tests <- c(368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50)
tests <- c(371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52)
spec_df <- data.frame(neg_tests, tests, sample_spec = neg_tests / tests)

pos_tests <- 50
tests <- 3330
unk_df <- data.frame(pos_tests, tests, sample_prev = pos_tests / tests)

#' Stan model for prior sensitivity analysis
code_sens <- root("coronavirus", "prior-sensitivity.stan")
#| output: asis
print_stan_file(code_sens)
#'
model_sens <- cmdstan_model(code_sens)

#' Data
data <- list(
  K_pos = nrow(sens_df),
  N_pos = array(sens_df$tests),
  n_pos = array(sens_df$pos_tests),
  K_neg = nrow(spec_df),
  N_neg = array(spec_df$tests),
  n_neg = array(spec_df$neg_tests),
  K_unk = nrow(unk_df),
  N_unk = array(unk_df$tests),
  n_unk = array(unk_df$pos_tests)
)

#' Empty data frame to store posterior intervals with different priors
ribbon_df <- data.frame(
  sigma_sens = c(),
  sigma_spec = c(),
  prev05 = c(),
  prev50 = c(),
  prev95 = c()
)

#' Loop over different prior parameter values
sigma_senss <- c(0.01, 0.25, 0.5, 0.75, 1)
sigma_specss <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
for (sigma_sens in sigma_senss) {
  for (sigma_spec in sigma_specss) {
    data2 <- append(data, list(sigma_sigma_logit_sens = sigma_sens,
                               sigma_sigma_logit_spec = sigma_spec))
    fit <-  model_sens$sample(
      data = data2,
      iter_warmup = 5e4,
      iter_sampling = 5e4,
      seed = 1234,
      refresh = 0,
      adapt_delta = 0.95,
      show_messages = FALSE, 
      show_exceptions = FALSE
    )
    pis <- fit$draws()[, , "pi"]
    ribbon_df <- rbind(ribbon_df,
                       data.frame(sigma_sens = paste("sensitivity hyperprior", "=", sigma_sens),
                                  sigma_spec = sigma_spec,
                                  prev05 = quantile(pis, 0.05),
                                  prev50 = quantile(pis, 0.5),
                                  prev95 = quantile(pis, 0.95)))
  }
}
#| label: fig-prior-sensitivity-2
#| fig-height: 2.5
#| fig-width: 9
ggplot(ribbon_df, aes(x = sigma_spec)) +
  facet_wrap(~ sigma_sens, nrow = 1) +
  geom_ribbon(aes(ymin = prev05, ymax = prev95), fill = "gray95") +
  geom_line(aes(y = prev50), linewidth = 0.5) +
  geom_line(aes(y = prev05), color = "darkgray", linewidth = 0.25) +
  geom_line(aes(y = prev95), color = "darkgray", linewidth = 0.25) +
  scale_y_log10(limits = c(0.0009, 1.1), breaks = c(0.001, 0.01, 0.1, 1)) +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 1),
                     breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  ylab("prevalence") +
  xlab("specificity hyperprior") +
  theme_bw() +
  theme(panel.spacing = unit(0.25, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2022--2025, Andrew Gelman and Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2022--2025, Andrew Gelman and Aki Vehtari, licensed under CC-BY-NC 4.0.

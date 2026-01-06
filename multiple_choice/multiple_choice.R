#' ---
#' title: "Modeling performance on a multiple choice exam"
#' author: "Andrew Gelman"
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
#' Chapter 4 *Introduction to workflow: Modeling performance on a
#' multiple choice exam*.
#'
#' # Introduction
#'
#' We demonstrate with an example from Section 4.15 of
#' @Gelman-Vehtari:2024, assessing the grading of a multiple-choice test.
#' We analyze data from a 24-question final exam from a class of 32
#' students, where our applied goal is to check that the individual
#' test questions are doing a good job at discriminating between
#' poorly- and well-performing students.  No external data are
#' available on the students, so we assess their abilities using their
#' total score on the exam.
#'
#' Each item is a multiple choice question with 4 possible answers and
#' is scored either 1 (correct) or 0 (incorrect). The total scores
#' across the 24 question exam range from 12 to 21, with an average of
#' 16. Across the questions the hardest question has 4 of the 32
#' students answering it correctly, while all students manage to
#' answer the easiest question correctly.
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
library(arm)
set.seed(123)
# Bring in plotting functions from a separate file
source(root("multiple_choice", "plot_functions.R"))

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

#' # Stan models
#' 
#' Compile all Stan programs to use throughout the file
#| label: compile-stan-models
logit_0 <- cmdstan_model(root("multiple_choice", "logit_0.stan"))
logit_prior <- cmdstan_model(root("multiple_choice", "logit_prior.stan"))
logit_guessing <- cmdstan_model(root("multiple_choice", "logit_guessing.stan"))
logit_guessing_uncentered <- cmdstan_model(root("multiple_choice", "logit_guessing_uncentered.stan"))
logit_threshold <- cmdstan_model(root("multiple_choice", "logit_threshold.stan"))
logit_threshold_prior <- cmdstan_model(root("multiple_choice", "logit_threshold_prior.stan"))
logit_guessing_multilevel <- cmdstan_model(root("multiple_choice", "logit_guessing_multilevel.stan"))
logit_guessing_uncentered_multilevel <- cmdstan_model(root("multiple_choice", "logit_guessing_uncentered_multilevel.stan"))
logit_guessing_multilevel_bivariate <- cmdstan_model(root("multiple_choice", "logit_guessing_multilevel_bivariate.stan"))
logit_guessing_multilevel_bivariate_cholesky <- cmdstan_model(root("multiple_choice", "logit_guessing_multilevel_bivariate_cholesky.stan"))
irt_guessing <- cmdstan_model(root("multiple_choice", "irt_guessing.stan"))
irt_guessing_discrimination <- cmdstan_model(root("multiple_choice", "irt_guessing_discrimination.stan"))

#' # Data
#' 
#' Read in data and construct score for each student
responses <- read.csv(root("multiple_choice", "data", "final_exam_responses.csv"))
answers <- read.csv(root("multiple_choice", "data", "final_exam_answers.csv"))
J <- nrow(responses)  # number of students
K <- ncol(responses)  # number of items
correct <- array(NA, c(J,K))
for (k in 1:K){
  correct[,k] <- ifelse(responses[,k] == as.character(answers[k]), 1, 0)
}
score <- rowSums(correct)
item <- colSums(correct)
summary(score)
summary(item)

score_jitt <- score + jitter(rep(0, J), amount = 0.3)
score_adj <- (score - mean(score)) / sd(score)
score_adj_jitt <- (score_jitt - mean(score)) / sd(score)
data <- list(J = J, x = score, y = correct)

item_id_0 <- LETTERS[1:J]  # Only works here because J is no more than 26!

#' # Simple models
#' 
#' ## Base model logit_0
#| output: asis
print_stan_code(logit_0$code())

#| label: fig-final_exams_1
#| fig-height: 4
#| fig-width: 9
plot_logit(
  "Fit to item A on exam",
  "Score on exam",
  logit_0,
  list(J = data$J, x = data$x, y = data$y[, 1]),
  score_jitt,
  guessprob = 0
)

#' ## Add priors
#| output: asis
print_stan_code(logit_prior$code())

#| label: fig-final_exams_2
#| fig-height: 4
#| fig-width: 9
plot_logit(
  "Fit to item A:  rescaled predictor and weakly informative prior",
  "Standardized exam score",
  logit_prior,
  list(J = data$J, x = data$x, y = data$y[, 1], 
       mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5),
  score_adj_jitt,
  guessprob = 0
)

#| label: fig-final_exams_2b
#| fig-height: 6
#| fig-width: 11
plot_logit_grid(
  "Rescaled predictor and weakly informative prior",
  "Standardized exam score",
  logit_prior,
  c(data, list(mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5)),
  score_adj_jitt,
  item_id_0,
  guessprob = 0
)

#| label: fig-final_exams_2_challenge
#| fig-height: 4
#| fig-width: 9
plot_logit(
  "Fit to item G: rescaled predictor and weakly informative prior",
  "Standardized exam score",
  logit_prior,
  list(J = data$J, x = data$x, y = data$y[, 7], 
       mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5),
  score_adj_jitt,
  guessprob = 0
)

#' ## Fix the data problems
answers[c(5, 14, 17)] <- c("d", "d", "c")
correct <- array(NA, c(J,K))
for (k in 1:K){
  correct[,k] <- ifelse(responses[,k] == as.character(answers[k]), 1, 0)
}
score <- rowSums(correct)
score_jitt <- score + jitter(rep(0, J), amount = 0.3)
score_adj <- (score - mean(score)) / sd(score)
score_adj_jitt <- (score_jitt - mean(score)) / sd(score)
data <- list(J = J, x = score, y = correct)
item_id <- rank(colSums(correct), ties = "first")

#| label: fig-final_exams_3
#| fig-height: 6
#| fig-width: 11
plot_logit_grid(
  "After fixing the data problem",
  "Standardized exam score",
  logit_prior,
  c(data, list(mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5)),
  score_adj_jitt,
  item_id,
  guessprob = 0
)

#' ## Allow for guessing
#| output: asis
print_stan_code(logit_guessing$code())

#| label: fig-final_exams_4
#| fig-height: 6
#| fig-width: 11
plot_logit_grid(
  "Probabilities constrained to range from 0.25 to 1",
  "Standardized exam score",
  logit_guessing,
  c(data, list(mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5)),
  score_adj_jitt,
  item_id,
  guessprob = 0.25
)

#' # Multilevel models
#' 
#' In preparation for multilevel model, create long dataset
N <- J*K
y <- rep(NA, N)
student <- rep(NA, N)
item <- rep(NA, N)
count <- 0
for (j in 1:J){
  for (k in 1:K){
    count <- count + 1
    y[count] <- correct[j,k]
    student[count] <- j
    item[count] <- k
  }
}
longdata <- list(
  N = N, J = J, K = K, 
  student = student, 
  item = item, 
  y = y, 
  x = score
)

#' ## Multilevel model
#| output: asis
print_stan_code(logit_guessing_multilevel$code())

#| label: fig-final_exams_5
#| fig-height: 6
#| fig-width: 11
fit_5 <- plot_logit_grid_2(
  "Multilevel model, partially pooling across the 24 exam questions",
  "Standardized exam score",
  logit_guessing_multilevel,
  c(longdata, list(
    mu_mu_a = 0, sigma_mu_a = 5, 
    mu_mu_b = 0, sigma_mu_b = 5,
    mu_sigma_a = 5, mu_sigma_b = 5
  )),
  score_adj_jitt,
  item_id,
  guessprob = 0.25
)

print(fit_5, variables = c("mu_a", "sigma_a", "mu_b", "sigma_b"))

#' ## Multilevel model with correlation
#| output: asis
print_stan_code(logit_guessing_multilevel_bivariate$code())

#| label: fig-final_exams_6
#| fig-height: 6
#| fig-width: 11
fit_6 <- plot_logit_grid_2(
  "Multilevel model with correlation",
  "Standardized exam score",
  logit_guessing_multilevel_bivariate,
  c(longdata, list(
    mu_mu_ab = c(0, 0),
    sigma_mu_ab = c(5, 10),
    mu_sigma_ab = c(5, 10)
  )),
  score_adj_jitt,
  item_id,
  guessprob = 0.25
)

print(fit_6, variables = c("mu_ab", "sigma_ab", "Omega_ab"))

#' ## Multilevel model with correlation using Cholesky
#| output: asis
print_stan_code(logit_guessing_multilevel_bivariate_cholesky$code())

#| label: fig-final_exams_7
#| fig-height: 6
#| fig-width: 11
fit_7 <- plot_logit_grid_2(
  "Multilevel model with correlation:  Cholesky parameterization",
  "Standardized exam score",
  logit_guessing_multilevel_bivariate_cholesky,
  c(longdata, list(
    mu_mu_ab = c(0, 0),
    sigma_mu_ab = c(5, 10),
    mu_sigma_ab = c(5, 10)
  )),
  score_adj_jitt,
  item_id,
  guessprob = 0.25
)

print(fit_7, variables = c("mu_ab", "sigma_ab", "Omega_ab"))


#' # Item-response theory (IRT) models
#'
#' ## Item-response model
#| output: asis
print_stan_code(irt_guessing$code())

#| label: fig-final_exams_11
#| fig-height: 6
#| fig-width: 11
fit_11 <- plot_irt(
  "Item-response model",
  irt_guessing,
  c(longdata, list(
    mu_mu_beta = 0, sigma_mu_beta = 5,
    mu_sigma_alpha = 5, mu_sigma_beta = 5
  )),
  item_id,
  guessprob = 0.25
)

#' ## Item-response model with discrimination parameters
#| output: asis
print_stan_code(irt_guessing_discrimination$code())

#| label: fig-final_exams_12
#| fig-height: 6
#| fig-width: 11
fit_12 <- plot_irt(
  "Item-response model with discrimination parameters",
  irt_guessing_discrimination,
  c(longdata, list(
    mu_mu_beta = 0, sigma_mu_beta = 5,
    mu_sigma_alpha = 5, mu_sigma_beta = 5,
    guessprob = 0.25,
    mu_sigma_gamma = 0.5
  )),
  item_id,
  guessprob = 0.25
)

#' ## Item-response model with discrimination parameters with init
#| label: fig-final_exams_13
#| fig-height: 6
#| fig-width: 11
fit_13 <- plot_irt(
  "Item-response model with discrimination parameters",
  irt_guessing_discrimination,
  c(longdata, list(
    mu_mu_beta = 0, sigma_mu_beta = 5,
    mu_sigma_alpha = 5, mu_sigma_beta = 5,
    mu_sigma_gamma = 0.5
  )),
  item_id,
  guessprob = 0.25,
  init = 0.1
)

#' IRT plots
alpha_sims <- as.matrix(fit_13$draws("alpha", format = "df"))[, 1:J]
beta_sims <- as.matrix(fit_13$draws("beta", format = "df"))[, 1:K]
gamma_sims <- as.matrix(fit_13$draws("gamma", format = "df"))[, 1:K]
alpha_hat <- apply(alpha_sims, 2, median)
alpha_sd <- apply(alpha_sims, 2, mad)
beta_hat <- apply(beta_sims, 2, median)
beta_sd <- apply(beta_sims, 2, mad)
gamma_hat <- apply(gamma_sims, 2, median)
gamma_sd <- apply(gamma_sims, 2, mad)

#| label: fig-irt_displays_1
#| fig-height: 4
#| fig-width: 6
par(mar = c(3, 0, 0, 0), mgp = c(1.5, .2, 0), tck = -.01)
rng <- range(
  alpha_hat - 3*alpha_sd,
  beta_hat - 3*beta_sd,
  alpha_hat + 3*alpha_sd,
  beta_hat + 3*beta_sd
)
plot(
  x = rng, 
  y = c(-1, 1),  
  xlab = "Posterior distributions for student abilities (above) and item difficulties (below)",
  ylab = "", yaxt = "n",
  bty = "n", type = "n"
)
for (j in 1:J){
  curve(dnorm(x, alpha_hat[j], alpha_sd[j]), col = "red", add = TRUE)
}
for (k in 1:K){
  curve(-dnorm(x, beta_hat[k], beta_sd[k]), col = "red", add = TRUE)
}

#| label: fig-irt_displays_2
#| fig-height: 3.2
#| fig-width: 4
par(mar = c(2.5, 2.5, .5, .5), mgp = c(1.5, .2, 0), tck = -.01)
x_rng <- range(beta_hat - beta_sd, beta_hat + beta_sd)
y_rng <- range(gamma_hat - gamma_sd, gamma_hat + gamma_sd)
plot(
  x_rng,
  y_rng,
  xlab = expression(beta[k]),
  ylab = expression(gamma[k]),
  bty = "l", type = "n"
)
for (k in 1:K) {
  lines(
    beta_hat[k] + c(-1, 1) * 0,
    gamma_hat[k] + c(-1, 1) * gamma_sd[k],
    col = "red", lwd = .5
  )
  lines(
    beta_hat[k] + c(-1, 1) * beta_sd[k],
    gamma_hat[k] + c(-1, 1) * 0,
    col = "red", lwd = .5
  )
}
text(beta_hat, gamma_hat, item_id, col = "blue", cex = .9)

#' # Prior predictive simulations
prior_predictive <- function(x, x_jitt, mu_a, sigma_a, mu_b, sigma_b) {
  a <- rnorm(1, mu_a, sigma_a)
  b <- rnorm(1, mu_b, sigma_b)
  y <- rbinom(length(x), 1, invlogit(a + b*x))
  plot(
    range(x), c(0, 1),
    xlab = "x", ylab = "y",
    xaxt = "n", yaxt = "n",
    yaxs = "i", bty = "l", type = "n"
  )
  axis(1, seq(-2,2,1))
  axis(2, c(0, 1))
  points(x_jitt, 0.5 + 0.96 * (y - 0.5), cex = .7, pch = 20, col = "blue")
}

#| label: fig-multiplechoice_prior_predictive_1
#| fig-height: 2.5
#| fig-width: 7.5
par(oma = c(0, 0, 1.5, 0), mfrow = c(2, 5), mar = c(3, 3, 1, 1), 
    mgp = c(1.3, .2, 0), tck = -.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 0.5, 0, 0.5)
}
mtext("10 prior predictive simulations with a ~ normal(0, 0.5) and b ~ normal(0, 0.5)", 
      side = 3, line = .5, outer = TRUE, cex = .7)

#| label: fig-multiplechoice_prior_predictive_2
#| fig-height: 2.5
#| fig-width: 7.5
par(oma = c(0, 0, 1.5, 0), mfrow = c(2, 5), mar = c(3, 3, 1, 1), 
    mgp = c(1.3, .2, 0), tck = -.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 5, 0, 5)
}
mtext("10 prior predictive simulations with a ~ normal(0, 5) and b ~ normal(0, 5)", 
      side = 3, line = .5, outer = TRUE, cex = .7)

#| label: fig-multiplechoice_prior_predictive_3
#| fig-height: 2.5
#| fig-width: 7.5
par(oma = c(0, 0, 1.5, 0), mfrow = c(2, 5), mar = c(3, 3, 1, 1), 
    mgp = c(1.3, .2, 0), tck = -.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 50, 0, 50)
}
mtext("10 prior predictive simulations with a ~ normal(0, 50) and b ~ normal(0, 50)", 
      side = 3, line = .5, outer = TRUE, cex = .7)


#' # Breaking the model
#| output: asis
print_stan_code(logit_guessing_uncentered$code())

#' Simulate data
set.seed(123)
J <- 32
x <- runif(J, 10, 20)
a_ <- -6
b_ <- 0.4
y <- rbinom(J, 1, 0.25 + 0.75 * invlogit(a_ + b_ * x))
m_x <- mean(x)
s_x <- sd(x)
x_adj <- (x - m_x)/s_x

break_data <- list(
  J = J,
  x = x,
  y = y,
  mu_a = 0,
  sigma_a = 1000,
  mu_b = 0,
  sigma_b = 1000
)
#| results: hide
break_1_fit <- logit_guessing_uncentered$sample(data = break_data, refresh = 0)

print(break_1_fit)

a <- extract_variable(break_1_fit, "a")
b <- extract_variable(break_1_fit, "b")
n_sims <- length(a)

#| label: fig-break_1
#| fig-height: 3
#| fig-width: 4
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0), tck = -.01)
plot(
  x, y,
  xlab = "Exam score",
  ylab = "Pr (correct answer)",
  yaxs = "i", bty = "l", type = "n"
)
for (s in sample(n_sims, 20)) {
  curve(
    0.25 + 0.75 * invlogit(a[s] + b[s] * x),
    lwd = .5, col = "red", add = TRUE
  )
}
points(x, 0.5 + 0.985 * (y - 0.5), cex = .7, pch = 20)
curve(0.25 + 0.75 * invlogit(a_ + b_ * x), add = TRUE)

#| results: hide
break_2_fit <- logit_guessing$sample(data = break_data, refresh = 0)

print(break_2_fit)

a <- extract_variable(break_2_fit, "a")
b <- extract_variable(break_2_fit, "b")
n_sims <- length(a)

#| label: fig-break_2
#| fig-height: 3
#| fig-width: 4
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0), tck = -.01)
plot(
  x_adj, y,
  xlim = c(-2, 2),
  xlab = "Standardized exam score",
  ylab = "Pr (correct answer)",
  yaxs = "i", bty = "l", type = "n"
)
for (s in sample(n_sims, 20)) {
  curve(
    0.25 + 0.75 * invlogit(a[s] + b[s] * x),
    lwd = .5, col = "red", add = TRUE
  )
}
points(x_adj, 0.5 + 0.985 * (y - 0.5), cex = .7, pch = 20)
curve(0.25 + 0.75 * invlogit(a_ + b_ * (m_x + s_x * x)), add = TRUE)

#| label: fig-final_exams_break_1
#| fig-height: 6
#| fig-width: 11
plot_logit_grid_2(
  "Breaking the model",
  "Exam score",
  logit_guessing_multilevel,
  c(longdata, list(
    mu_mu_a = 0, sigma_mu_a = 5, 
    mu_mu_b = 0, sigma_mu_b = 5, 
    mu_sigma_a = 5, mu_sigma_b = 5
  )),
  score_jitt,
  item_id,
  guessprob = 0.25
)

#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2022--2025, Andrew Gelman, licensed under BSD-3.
#' * Text &copy; 2022--2025, Andrew Gelman, licensed under CC-BY-NC 4.0.

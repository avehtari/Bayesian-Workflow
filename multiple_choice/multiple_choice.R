library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
options(mc.cores = 4)
library(posterior)
library(arm)
set.seed(123)

# Bring in plotting functions from separate file
source(root("multiple_choice", "plot_functions.R"))

# Compile Stan programs to use throughout the file
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

# Read in data and construct score for each student
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

# The plots_logit function fits the model and makes plots
plots_logit(
  root("multiple_choice", "final_exams_1"),
  "Fit to item A on exam",
  "Score on exam",
  logit_0,
  list(J = data$J, x = data$x, y = data$y[, 1]),
  score_jitt,
  guessprob = 0
)

# Add priors
plots_logit(
  root("multiple_choice", "final_exams_2"),
  "Fit to item A:  rescaled predictor and weakly informative prior",
  "Standardized exam score",
  logit_prior,
  list(J = data$J, x = data$x, y = data$y[, 1], 
       mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5),
  score_adj_jitt,
  guessprob = 0
)
plots_logit_grid(
  root("multiple_choice", "final_exams_2"),
  "Rescaled predictor and weakly informative prior",
  "Standardized exam score",
  logit_prior,
  c(data, list(mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5)),
  score_adj_jitt,
  item_id_0,
  guessprob = 0
)

plots_logit(
  root("multiple_choice", "final_exams_2_challenge"),
  "Fit to item G:  rescaled predictor and weakly informative prior",
  "Standardized exam score",
  logit_prior,
  list(J = data$J, x = data$x, y = data$y[, 7], 
       mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5),
  score_adj_jitt,
  guessprob = 0
)

# Fix the data problems
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

plots_logit_grid(
  root("multiple_choice", "final_exams_3"),
  "After fixing the data problem",
  "Standardized exam score",
  logit_prior,
  c(data, list(mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5)),
  score_adj_jitt,
  item_id,
  guessprob = 0
)

# Allow for guessing
plots_logit_grid(
  root("multiple_choice", "final_exams_4"),
  "Probabilities constrained to range from 0.25 to 1",
  "Standardized exam score",
  logit_guessing,
  c(data, list(mu_a = 0, sigma_a = 5, mu_b = 0, sigma_b = 5)),
  score_adj_jitt,
  item_id,
  guessprob = 0.25
)

# In preparation for multilevel model, create long dataset
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


fit_5 <- plots_logit_grid_2(
  root("multiple_choice", "final_exams_5"),
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

fit_6 <- plots_logit_grid_2(
  root("multiple_choice", "final_exams_6"),
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

fit_7 <- plots_logit_grid_2(
  root("multiple_choice", "final_exams_7"),
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


# IRT models
fit_11 <- plots_irt(
  root("multiple_choice", "final_exams_11"),
  "Item-response model",
  irt_guessing,
  c(longdata, list(
    mu_mu_beta = 0, sigma_mu_beta = 5,
    mu_sigma_alpha = 5, mu_sigma_beta = 5
  )),
  item_id,
  guessprob = 0.25
)
fit_12 <- plots_irt(
  root("multiple_choice", "final_exams_12"),
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
fit_13 <- plots_irt(
  root("multiple_choice", "final_exams_13"),
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

alpha_sims <- as.matrix(fit_13$draws("alpha", format = "df"))[, 1:J]
beta_sims <- as.matrix(fit_13$draws("beta", format = "df"))[, 1:K]
gamma_sims <- as.matrix(fit_13$draws("gamma", format = "df"))[, 1:K]
alpha_hat <- apply(alpha_sims, 2, median)
alpha_sd <- apply(alpha_sims, 2, mad)
beta_hat <- apply(beta_sims, 2, median)
beta_sd <- apply(beta_sims, 2, mad)
gamma_hat <- apply(gamma_sims, 2, median)
gamma_sd <- apply(gamma_sims, 2, mad)

pdf(root("multiple_choice", "irt_displays_1.pdf"), height=4, width=6)
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
dev.off()

pdf(root("multiple_choice", "irt_displays_2.pdf"), height = 3.2, width = 4)
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
dev.off()

# Prior predictive simulations
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

pdf(root("multiple_choice", "multiplechoice_prior_predictive_1.pdf"), height = 2.5, width = 7.5)
par(oma = c(0, 0, 1.5, 0), mfrow = c(2, 5), mar = c(3, 3, 1, 1), 
    mgp = c(1.3, .2, 0), tck = -.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 0.5, 0, 0.5)
}
mtext("10 prior predictive simulations with a ~ normal(0, 0.5) and b ~ normal(0, 0.5)", 
      side = 3, line = .5, outer = TRUE, cex = .7)
dev.off()

pdf(root("multiple_choice", "multiplechoice_prior_predictive_2.pdf"), height = 2.5, width = 7.5)
par(oma = c(0, 0, 1.5, 0), mfrow = c(2, 5), mar = c(3, 3, 1, 1), 
    mgp = c(1.3, .2, 0), tck = -.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 5, 0, 5)
}
mtext("10 prior predictive simulations with a ~ normal(0, 5) and b ~ normal(0, 5)", 
      side = 3, line = .5, outer = TRUE, cex = .7)
dev.off()

pdf(root("multiple_choice", "multiplechoice_prior_predictive_3.pdf"), height = 2.5, width = 7.5)
par(oma = c(0, 0, 1.5, 0), mfrow = c(2, 5), mar = c(3, 3, 1, 1), 
    mgp = c(1.3, .2, 0), tck = -.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 50, 0, 50)
}
mtext("10 prior predictive simulations with a ~ normal(0, 50) and b ~ normal(0, 50)", 
      side = 3, line = .5, outer = TRUE, cex = .7)
dev.off()


# Breaking the model
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
break_1_fit <- logit_guessing_uncentered$sample(data = break_data, refresh = 0)
print(break_1_fit)
a <- extract_variable(break_1_fit, "a")
b <- extract_variable(break_1_fit, "b")
n_sims <- length(a)


pdf(root("multiple_choice", "break_1.pdf"), height = 3, width = 4)
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
dev.off()


break_2_fit <- logit_guessing$sample(data = break_data, refresh = 0)
print(break_2_fit)
a <- extract_variable(break_2_fit, "a")
b <- extract_variable(break_2_fit, "b")
n_sims <- length(a)



pdf(root("multiple_choice", "break_2.pdf"), height = 3, width = 4)
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
dev.off()


plots_logit_grid_2(
  root("multiple_choice", "final_exams_break_1"),
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


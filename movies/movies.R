#' ---
#' title: "Simulated data of movie ratings"
#' author: "Andrew Gelman"
#' date: 2022-08-15
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
#' Chapter 16 *Coding a series of models: Simulated data of movie
#' ratings*.
#' 
#' # Introduction
#'
#' Consider the following scenario.  You are considering which of two
#' movies to go see.  Both have average online ratings of 4 out of 5
#' stars, but one is based on 2 ratings and the other is based on 100.
#' Which movie should you choose?
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

#' **Load packages and set options**
#| cache: FALSE
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
options(mc.cores = 4)
library(posterior)
library(ggplot2)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(ggdist)
library(patchwork)
library(dplyr)
set.seed(1234)

#' # Model for two movies
y_1 <- c(3, 5)
y_2 <- rep(c(2, 3, 4, 5), c(10, 20, 30, 40))
y <- c(y_1, y_2)
N <- length(y)
movie <- rep(c(1, 2), c(length(y_1), length(y_2)))
movie_data <- list(y = y, N = N, movie = movie)
mod_1 <- cmdstan_model(root("movies", "ratings_1.stan"))
#' Stan model code
mod_1
#' Sample
#| label: fit_1
#| results: hide
fit_1 <- mod_1$sample(data = movie_data, refresh = 0)
#' Posterior summary
print(fit_1)

#' # Extending the model to J movies
J <- 40
N_ratings <- sample(0:100, J, replace = TRUE)
N <- sum(N_ratings)
movie <- rep(1:J, N_ratings)
theta <- rnorm(J, 3.0, 0.5)
y <- rnorm(N, theta[movie], 2.0)
movie_data <- list(y = y, N = N, J = J, movie = movie)
mod_2 <- cmdstan_model(root("movies", "ratings_2.stan"))
#' Stan model code
mod_2

#' Sample
#| label: fit_2
#| results: hide
fit_2 <- mod_2$sample(data = movie_data, refresh = 0)
#' Posterior summary
print(fit_2)

theta_post <- fit_2$draws("theta", format = "matrix")
theta_post_quants <- t(apply(theta_post, 2, function(x) 
  quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
))

#| label: fig-movies_1
#| fig-height: 4
#| fig-width: 5
par(mar = c(3, 3, 2, 1), mgp = c(1.7, 0.5, 0), tck = -0.02)
par(pty = "s")
rng <- range(theta_post_quants, theta)
plot(rng, rng, xlab = "Posterior median, 50%, and 95% interval",
     ylab = "True parameter value", bty = "l", type = "n")
abline(0, 1, col = "gray")
points(theta_post_quants[ , "50%"], theta, pch = 20)
for (j in 1:J) {
  lines(c(theta_post_quants[j, "25%"], theta_post_quants[j, "75%"]),
        rep(theta[j], 2), lwd = 2)
  lines(c(theta_post_quants[j, "2.5%"], theta_post_quants[j, "97.5%"]),
        rep(theta[j], 2), lwd = 0.5)
}
mtext(expression(paste("Comparing parameters ", theta[j], 
                       " to their posterior inferences")), side = 3)

#' ggplot version
#| label: fig-gg-movies_1
#| fig-height: 4.5
#| fig-width: 4.5
#| out-width: 70%
draws <- as_draws_rvars(fit_2$draws("theta"))
rng <- range(summarize_draws(
  draws, ~quantile(.x, probs = c(0.025, 0.975)))[,c("2.5%","97.5%")])
ggplot(data = NULL) +
  coord_fixed(xlim = rng, ylim = rng) +
  geom_abline(color = "gray") +
  stat_pointinterval(aes(y = theta, xdist = draws$theta),
                     .width = c(0.5, 0.95),
                     interval_size_range = c(0.4, 0.8),
                     alpha = 0.5) +
  labs(x = "Posterior median, 50%, and 95% interval",
       y = "True parameter value",
       subtitle = expression(paste("Comparing parameters ", theta[j], 
                                " to their posterior inferences")))

#| label: fig-movies_2
#| fig-height: 4
#| fig-width: 6
interval_width <- theta_post_quants[,"75%"] - theta_post_quants[,"25%"]
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(c(0, 1.02 * max(N_ratings)), c(0, 1.02 * max(interval_width)), 
     xlab = "Number of ratings", ylab = "Width of 50% posterior interval", 
     yaxs = "i", yaxs = "i", bty = "l", type = "n")
points(N_ratings, interval_width, pch = 20)
mtext("Where you have more data, you have less uncertainty", side = 3)

#' ggplot version
#| label: fig-gg-movies_2
#| fig-height: 4
#| fig-width: 6
bind_cols(tibble(N = N_ratings),
          summarize_draws(draws, ~quantile(.x, probs = c(0.025, 0.975)))) |>
  mutate(interval_width = `97.5%` - `2.5%`) |>
  ggplot(aes(x = N, y = interval_width)) +
  geom_point() +
  ylim(c(0, NA)) +
  labs(x = "Number of ratings",
       y = "Width of 50% posterior interval",
       title = "Where you have more data, you have less uncertainty")

#' # Item-response model with parameters for raters and for movies
#' 
#' ## Balanced data
J <- 40
K <- 100
N <- J * K
movie <- rep(1:J, rep(K, J))
rater <- rep(1:K, J)
mu <- 3
sigma_a <- 0.5
sigma_b <- 0.5
sigma_y <- 2
alpha <- rnorm(J, 0, 1)
beta <- rnorm(K, 0, 1)
y <- rnorm(N, mu + sigma_a * alpha[movie] - sigma_b * beta[rater], sigma_y)
data_3 <- list(N = N, J = J, K = K, movie = movie, rater = rater, y = y)
mod_3 <- cmdstan_model(root("movies", "ratings_3.stan"))
#' Stan model code
mod_3
#' Sample
#| label: fit_3
#| results: hide
fit_3 <- mod_3$sample(data = data_3, refresh = 0)
#' Posterior summary
print(fit_3, variables = c("mu", "sigma_a", "sigma_b", "sigma_y"))

alpha_post <- fit_3$draws("alpha", format = "matrix")
beta_post <- fit_3$draws("beta", format = "matrix")
quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)
alpha_post_quants <- t(apply(alpha_post, 2, function(x) quantile(x, probs = quants)))
beta_post_quants <- t(apply(beta_post, 2, function(x) quantile(x, probs = quants)))

#| label: fig-movies_3
#| fig-height: 4
#| fig-width: 9
par(mfrow = c(1, 2))
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -0.02)
par(pty = "s")
rng <- range(alpha_post_quants, alpha)
plot(rng, rng, 
     xlab = "Posterior median, 50%, and 95% interval", 
     ylab = "True parameter value", 
     bty = "l", type = "n")
abline(0, 1, col = "gray")
points(alpha_post_quants[, "50%"], alpha, pch = 20)
for (j in 1:J){
  lines(c(alpha_post_quants[j, "25%"], alpha_post_quants[j, "75%"]), 
        rep(alpha[j], 2), lwd = 2)
  lines(c(alpha_post_quants[j, "2.5%"], alpha_post_quants[j, "97.5%"]), 
        rep(alpha[j], 2), lwd = 0.5)
}
mtext(expression(paste("Checking the ", alpha[j], "'s")), side = 3)
rng <- range(beta_post_quants, beta)
plot(rng, rng, 
     xlab = "Posterior median, 50%, and 95% interval", 
     ylab = "True parameter value", 
     bty = "l", type = "n")
abline(0, 1, col = "gray")
points(beta_post_quants[, "50%"], beta, pch = 20)
for (k in 1:K){
  lines(c(beta_post_quants[k, "25%"], beta_post_quants[k, "75%"]), 
        rep(beta[k], 2), lwd = 2)
  lines(c(beta_post_quants[k, "2.5%"], beta_post_quants[k, "97.5%"]), 
        rep(beta[k], 2), lwd = 0.5)
}
mtext(expression(paste("Checking the ", beta[j], "'s")), side = 3)

#' ggplot version
#| label: fig-gg-movies_3
#| fig-height: 4
#| fig-width: 8
drawsa <- as_draws_rvars(fit_3$draws(c("alpha")))
rnga <- range(summarize_draws(
  drawsa, ~quantile(.x, probs = c(0.025, 0.975)))[,c("2.5%","97.5%")])
p3a <- ggplot(data = NULL) +
  coord_fixed(xlim = rnga, ylim = rnga) +
  geom_abline(color = "gray") +
  stat_pointinterval(aes(y = alpha, xdist = drawsa$alpha),
                     .width = c(0.5, 0.95),
                     interval_size_range = c(0.4, 0.8),
                     alpha = 0.5) +
  labs(x = "Posterior median, 50%, and 95% interval",
       y = "True parameter value",
       title = expression(paste("Checking the ", alpha[j], "'s")))
drawsb <- as_draws_rvars(fit_3$draws(c("beta")))
rngb <- range(summarize_draws(
  drawsb, ~quantile(.x, probs = c(0.025, 0.975)))[,c("2.5%","97.5%")])
p3b <- ggplot(data = NULL) +
  coord_fixed(xlim = rngb, ylim = rngb) +
  geom_abline(color = "gray") +
  stat_pointinterval(aes(y = beta, xdist = drawsb$beta),
                     .width = c(0.5, 0.95),
                     interval_size_range = c(0.4, 0.8),
                     alpha = 0.5) +
  labs(x = "Posterior median, 50%, and 95% interval",
       y = "True parameter value",
       title = expression(paste("Checking the ", beta[j], "'s")))
p3a + p3b

#' ## Unbalanced data
genre <- rep(c("romantic", "crime"), c(round(J / 2), J - round(J / 2)))
prob_of_rated <- ifelse(beta[rater] > 0,
                        ifelse(genre[movie] == "romantic", 0.2, 0.7),
                        ifelse(genre[movie] == "romantic", 0.7, 0.2))
rated <- rbinom(N, 1, prob_of_rated)  == 1 # TRUE if movie was rated, FALSE if not
data_3a <- list(N = sum(rated), J = J, K = K, 
                movie = movie[rated], rater = rater[rated], 
                y = y[rated])
#' Sample
#| label: fit_3a
#| results: hide
fit_3a <- mod_3$sample(data = data_3a, refresh = 0)
#' Posterior summary
print(fit_3a, variables = c("mu", "sigma_a", "sigma_b", "sigma_y"))

alpha_post <- fit_3a$draws("alpha", format = "matrix")
beta_post <- fit_3a$draws("beta", format = "matrix")
quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)
alpha_post_quants <- t(apply(alpha_post, 2, function(x) quantile(x, probs = quants)))
beta_post_quants <- t(apply(beta_post, 2, function(x) quantile(x, probs = quants)))

add_legend <- function(text, pch, range) {
  legend(0.6 * min(range) + 0.4 * max(range), 
         min(range) + 0.12 * (max(range) - min(range)), 
         text[1], pch = pch[1], cex = 0.8, bty = "n")
  legend(0.6 * min(range) + 0.4 * max(range), 
         min(range) + 0.06 * (max(range) - min(range)), 
         text[2], pch = pch[2], cex = 0.8, bty = "n")
}

#| label: fig-movies_4
#| fig-height: 4
#| fig-width: 9
par(mfrow = c(1, 2), oma = c(0, 0, 1, 0))
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
par(pty = "s")
rng <- range(alpha_post_quants, alpha)
plot(rng, rng, 
     xlab = "Posterior median, 50%, and 95% interval", 
     ylab = "True parameter value", bty = "l", type = "n")
abline(0, 1, col = "gray")
points(alpha_post_quants[genre=="romantic", "50%"], alpha[genre=="romantic"],
       pch = 1, cex = 0.9)
points(alpha_post_quants[genre == "crime", "50%"], alpha[genre == "crime"], 
       pch = 20)
for (j in 1:J){
  lines(c(alpha_post_quants[j, "25%"], alpha_post_quants[j, "75%"]), 
        rep(alpha[j], 2), lwd = 2)
  lines(c(alpha_post_quants[j, "2.5%"], alpha_post_quants[j, "97.5%"]), 
        rep(alpha[j], 2), lwd = 0.5)
}
add_legend(c("Romantic comedies", "Crime movies"), pch = c(1, 20), range = rng)
mtext(expression(paste("Checking the ", alpha[j], "'s")), side = 3)

rng <- range(beta_post_quants, beta)
plot(rng, rng, 
     xlab = "Posterior median, 50%, and 95% interval", 
     ylab = "True parameter value", 
     bty = "l", type = "n")
abline(0, 1, col = "gray")
points(beta_post_quants[beta < 0, "50%"], beta[beta < 0], pch = 1, cex = 0.9)
points(beta_post_quants[beta > 0, "50%"], beta[beta > 0], pch = 20)
for (k in 1:K){
  lines(c(beta_post_quants[k, "25%"], beta_post_quants[k, "75%"]), 
        rep(beta[k], 2), lwd = 2)
  lines(c(beta_post_quants[k, "2.5%"], beta_post_quants[k, "97.5%"]), 
        rep(beta[k], 2), lwd = 0.5)
}
add_legend(c("Nice raters", "Difficult raters"), pch = c(1, 20), range = rng)
mtext(expression(paste("Checking the ", beta[j], "'s")), side = 3)
mtext("Checking fits for model when difficult reviewers were more likely to rate certain genres",
      side = 3, outer = TRUE)


#' ggplot version
#| label: fig-gg-movies_4
#| fig-height: 4
#| fig-width: 8
drawsa <- as_draws_rvars(fit_3a$draws(c("alpha")))
rnga <- range(summarize_draws(
  drawsa, ~quantile(.x, probs = c(0.025, 0.975)))[,c("2.5%","97.5%")])
p4a <- ggplot(data = NULL) +
  coord_fixed(xlim = rnga, ylim = rnga) +
  geom_abline(color = "gray") +
  stat_pointinterval(aes(y = alpha, xdist = drawsa$alpha, shape = genre),
                     .width = c(0.5, 0.95),
                     interval_size_range = c(0.4, 0.8),
                     point_size = 2,
                     alpha = 0.5) +
  scale_shape_manual(values = c(1, 19), labels = c("Romantic comedies", "Crime movies")) +
  labs(x = "Posterior median, 50%, and 95% interval",
       y = "True parameter value",
       title = expression(paste("Checking the ", alpha[j], "'s"))) +
  guides(shape = guide_legend(position = "inside")) +
  theme(legend.justification.inside = c(0.5, 0),
        legend.title = element_blank())
drawsb <- as_draws_rvars(fit_3a$draws(c("beta")))
rngb <- range(summarize_draws(
  drawsb, ~quantile(.x, probs = c(0.025, 0.975)))[,c("2.5%","97.5%")])
p4b <- ggplot(data = NULL) +
  coord_fixed(xlim = rngb, ylim = rngb) +
  geom_abline(color = "gray") +
  stat_pointinterval(aes(y = beta, xdist = drawsb$beta, shape = factor(beta>0)),
                     .width = c(0.5, 0.95),
                     interval_size_range = c(0.4, 0.8),
                     point_size = 2,
                     alpha = 0.5) +
  scale_shape_manual(values = c(1, 19), labels = c("Nice raters", "Difficult raters")) +
  labs(x = "Posterior median, 50%, and 95% interval",
       y = "True parameter value",
       title = expression(paste("Checking the ", beta[j], "'s"))) +
  guides(shape = guide_legend(position = "inside")) +
  theme(legend.justification.inside = c(0.5, 0),
        legend.title = element_blank())
p4a + p4b

#' # Comparison to naive data averaging
ybar <- rep(NA, J)
for (j in 1:J) {
  ybar[j] <- mean(y[movie == j & rated])
}

a_true <- mu + sigma_a * alpha
a_post_median <- fit_3a$summary(variables =  "a")$median

#| label: fig-movies_5
#| fig-height: 4
#| fig-width: 9
par(mfrow = c(1, 2))
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -0.02)
par(pty = "s")
rng <- range(a_post_median, ybar, a_true)
plot(rng, rng, 
     xlab = "Raw average rating for movie j", 
     ylab = expression(paste("True ",  a[j])), bty = "l", type = "n")
abline(0, 1, col = "gray")
points(ybar[genre == "romantic"], a_true[genre == "romantic"], pch = 1, cex = 0.9)
points(ybar[genre == "crime"], a_true[genre == "crime"], pch = 20)
mtext("Problems with raw averaging", side = 3)
add_legend(c("Romantic comedies", "Crime movies"), pch = c(1, 20), range = rng)
plot(rng, rng, 
     xlab = "Posterior median estimate for movie j", 
     ylab = expression(paste("True ",  a[j])), bty = "l", type = "n")
abline(0, 1, col = "gray")
points(a_post_median[genre == "romantic"], a_true[genre == "romantic"], pch = 1, cex = 0.9)
points(a_post_median[genre == "crime"], a_true[genre == "crime"], pch = 20)
mtext("Model-based estimates do better", side = 3)
add_legend(c("Romantic comedies", "Crime movies"), pch = c(1, 20), range = rng)
mtext("Problems with raw averages when difficult reviewers were more likely to rate certain genres",
      side = 3, outer = TRUE)

#' ggplot version
#| label: fig-gg-movies_5
#| fig-height: 4
#| fig-width: 8
rng <- range(a_post_median, ybar, a_true)
p5a <- ggplot(data = NULL, aes(x = ybar, y = a_true, shape = genre)) +
  coord_fixed(xlim = rng, ylim = rng) +
  geom_abline(color = "gray") +
  geom_point(size = 2) +
  scale_shape_manual(values = c(1, 19), labels = c("Romantic comedies", "Crime movies")) +
  labs(x = "Raw average rating for movie j",
       y = expression(paste("True ",  a[j])),
       title = "Problems with raw averaging") +
  guides(shape = guide_legend(position = "inside")) +
  theme(legend.justification.inside = c(0.5, 0),
        legend.title = element_blank())
p5b <- ggplot(data = NULL, aes(x = a_post_median, y = a_true, shape = genre)) +
  coord_fixed(xlim = rng, ylim = rng) +
  geom_abline(color = "gray") +
  geom_point(size = 2) +
  scale_shape_manual(values = c(1, 19), labels = c("Romantic comedies", "Crime movies")) +
  labs(x = "Posterior median estimate for movie j",
       y = expression(paste("True ",  a[j])),
       title = "Model-based estimates do better") +
  guides(shape = guide_legend(position = "inside")) +
  theme(legend.justification.inside = c(0.5, 0),
        legend.title = element_blank())
p5a + p5b

#' ---
#' title: "Stochastic learning in dogs"
#' author: "Andrew Gelman"
#' date: 2022-07-16
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
#' This notebook includes the `CmdStanR` code for the Bayesian Workflow book
#' Chapter 18 *Posterior predictive checking: Stochastic learning in dogs*.
#'
#' # Introduction
#'
#' We analyse stochastic learning in dogs data by
#' @Bush+Mosteller:1955.
#' 
#+ setup, include=FALSE
knitr::opts_chunk$set(
  cache=FALSE,
  message=FALSE,
  error=FALSE,
  warning=TRUE,
  comment=NA,
  out.width='95%'
)

#' **Load packages**
#| cache: FALSE
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
options(mc.cores = 4)
library(posterior)
library(MASS)
library(arm)
set.seed(123)

#' # Data
dogs <- read.table(root("dogs", "data", "dogs.dat"), skip = 2)
shock <- ifelse(as.matrix(dogs[, 2:26]) == "S", 1, 0)
dogs_data <- list(y = shock, J = nrow(shock), T = ncol(shock))

#' # Models

dogs_0 <- cmdstan_model(root("dogs", "dogs_0.stan"))
#| results: hide
fit_0 <- dogs_0$sample(data = dogs_data, refresh = 0)
#'
print(fit_0)

dogs_1 <- cmdstan_model(root("dogs", "dogs_1.stan"))
#| results: hide
fit_1 <- dogs_1$sample(data = dogs_data, refresh = 0)
#'
print(fit_1)

dogs_2 <- cmdstan_model(root("dogs", "dogs_2.stan"))
#| results: hide
fit_2 <- dogs_2$sample(data = dogs_data, refresh = 0)
#'
print(fit_2)

dogs_3 <- cmdstan_model(root("dogs", "dogs_3.stan"))
#| results: hide
fit_3 <- dogs_3$sample(data = dogs_data, refresh = 0)
#'
print(fit_3, variables = c("mu_logit_a", "sigma_logit_a"))

dogs_4 <- cmdstan_model(root("dogs", "dogs_4.stan"))
#| results: hide
fit_4 <- dogs_4$sample(data = dogs_data, refresh = 0)
#'
print(fit_4, variables = c("mu_logit_ab", "Sigma_logit_ab"))

dogs_5 <- cmdstan_model(root("dogs", "dogs_5.stan"))
#| results: hide
fit_5 <- dogs_5$sample(data = dogs_data, refresh = 0)
#'
print(fit_5, variables = c("mu_logit_ab", "sigma_logit_ab", "Omega_logit_ab[1,2]",
                           "a[1]", "b[1]"))

#' # Plots
empty_plot <- function(a = "") {
  plot(0, 0, bty = "n", xaxt = "n", yaxt = "n", type = "n")
  text(0, 0, a, cex = .8)
}

plot_dogs <- function(y, ...) {
  J <- nrow(y)
  T <- ncol(y)
  max_y_times <- rep(NA, J)
  for (j in 1:J) {
    max_y_times[j] <- max((1:T)[y[j, ] == 1])
  }
  y_ordered <- y[rev(order(max_y_times)), ]
  image(t(y_ordered), bty = "n", xaxt = "n", yaxt = "n", ...)
}

plot_ppc <- function(fit, label){
  post <- as_draws_rvars(fit$draws())
  empty_plot(label)
  for (k in 1:3) {
    for (i in sample(1000, 2)) {
      rep <- sum(subset_draws(post, iter = i, chain = k)$y_rep)
      plot_dogs(rep)
    }
  }
}


#| label: fig-dogs_ppc
#| fig-height: 5.5
#| fig-width: 5.5
par(mfrow = c(7, 7), mar = c(.5, .5, .5, .5))
empty_plot("Real dogs")
plot_dogs(shock)
for (k in 1:5){
  empty_plot()
}
plot_ppc(fit_0, "PPsims from M0:\nlogit model")
plot_ppc(fit_1, "PPsims from M1:\n1-parameter\nlog model")
plot_ppc(fit_2, "PPsims from M2:\n2-parameter\nlog model")
plot_ppc(fit_3, "PPsims from M3:\nhier 1-par\nlog model")
plot_ppc(fit_4, "PPsims from M4:\nhier 2-par\nlog model")
plot_ppc(fit_5, "PPsims from M5:\nhier 2-par\nlog model\nwith prior")

#| label: fig-dogs_inference
#| fig-height: 3
#| fig-width: 7.5
post <- as_draws_rvars(fit_5$draws())
par(mfrow = c(2, 5), pty = "s", 
    mar = c(2.5, 2.5, 0.5, 0.5), mgp = c(1.5, 0.2, 0), 
    tck = -0.02, oma = c(0, 0, 1, 0))
for (k in 1:2){
  index <- sample(1000, 5)
  for (i in 1:5) {
    a_sim <- sum(subset_draws(post, iter = index[i], chain = k)$a)
    b_sim <- sum(subset_draws(post, iter = index[i], chain = k)$b)
    plot(c(0.55, 1), c(0.55, 1), 
         xlab= if (k == 2) "a" else "", ylab = if (i == 1) "b" else "",
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", type = "n")
    if (k==2) axis(1, c(0.6, 0.8, 1), c("0.6", "0.8", "1")) else axis(1, c(0.6, 0.8, 1), c("", "", "")) 
    if (i==1) axis(2, c(0.6, 0.8, 1), c("0.6", "0.8", "1")) else axis(1, c(0.6, 0.8, 1), c("", "", "")) 
    abline(0, 1, lwd = .5, col = "gray")
    points(a_sim, b_sim, pch = 20, cex = .6)
    mtext("10 posterior simulations of the parameters of the 30 dogs", 
          3, 0, cex = 0.8, outer = TRUE)
  }
}

#| label: fig-dogs_point_estimate
#| fig-height: 4
#| fig-width: 4
post <- as_draws_rvars(fit_5$draws())
par(pty = "s", mar = c(3, 3.5, 2, 1), mgp = c(2, .5, 0), tck = -.01)
plot(median(post$a), median(post$b), 
     xlim = c(0.55, 1), ylim = c(0.55, 1), xaxs = "i", yaxs = "i", 
     xlab = expression(hat(a)), ylab = expression(hat(b)), pch = 20, cex = 0.6, 
     main = "Posterior medians from fitted model", cex.main = 0.9)
  abline(0,1,lwd=.5, col="gray")

new_dogs_mu_logit_ab <- c(2.4, 1.3)
new_dogs_sigma_ab <- c(0.32, 0.40)
new_dogs_rho_ab <- 0
new_dogs_Sigma_ab <- 
  diag(new_dogs_sigma_ab) %*% 
  rbind(c(1, new_dogs_rho_ab), c(new_dogs_rho_ab, 1)) %*% 
  diag(new_dogs_sigma_ab)

J <- 30
new_dogs_ab <- invlogit(mvrnorm(J, new_dogs_mu_logit_ab, new_dogs_Sigma_ab))
a <- new_dogs_ab[, 1]
b <- new_dogs_ab[, 2]
T <- 25
new_dogs <- array(NA, c(J, T))
for (j in 1:J) {
  prev_shock <- 0
  prev_avoid <-  0
  new_dogs[j, 1] <- 1
  for (t in 2:T) {
    prev_shock = prev_shock + new_dogs[j, t - 1]
    prev_avoid = prev_avoid + 1 - new_dogs[j, t - 1]
    p = a[j] ^ prev_shock * b[j] ^ prev_avoid
    new_dogs[j, t] <- rbinom(1, 1, p)
  }
}
new_dogs_data <- list(y = new_dogs, J = J, T = T)
#| results: hide
new_fit_5 <- dogs_5$sample(data = new_dogs_data, refresh = 0)
#'
print(new_fit_5, variables = c("mu_logit_ab", "sigma_logit_ab", "Omega_logit_ab[1,2]",
                               "a[1]", "a[2]", "b[1]", "b[2]"))

#| label: fig-dogs_parameters
#| fig-height: 4
#| fig-width: 4
par(pty = "s", mar = c(3, 3.5, 2, 1), mgp = c(2, 0.5, 0), tck = -0.01)
plot(a, b, xlim = c(0.55, 1), ylim = c(0.55, 1), 
     xaxs = "i", yaxs = "i",  xlab = "a", ylab = "b", pch = 20, cex = 0.6, 
     main = "Simulated parameters", cex.main = 0.9)
abline(0, 1, lwd = 0.5, col = "gray")


#| label: fig-dogs_data
#| fig-height: 4
#| fig-width: 3
par(pty = "m", mar = c(1, 2, 2, 1))
plot_dogs(new_dogs, main = "Simulated data", cex.main = 0.9)

#| label: fig-dogs_calibration
#| fig-height: 4
#| fig-width: 4
post <- as_draws_rvars(new_fit_5$draws())
par(pty = "s", mar = c(3, 3.5, 2, 1), mgp = c(2, 0.5, 0), tck = -0.01)
plot(0, 0, xlim = c(0.55, 1), ylim = c(0.55, 1),
     xlab = "Posterior inference", ylab = "True parameter value", 
     xaxs = "i", yaxs = "i", pch = 20, cex = 0.6, 
     main = "Calibration check of posterior intervals", cex.main = 0.9)
abline(0, 1, lwd = 0.5, col = "gray")
for (j in 1:J){
  points(median(post$a[j]), a[j], pch = 20, cex = 0.6, col = "blue")
  lines(quantile(post$a[j], c(0.25, 0.75)), rep(a[j], 2), lwd = 0.5, col = "blue")
  points(median(post$b[j]), b[j], pch = 20, cex = 0.6, col = "red")
  lines(quantile(post$b[j], c(0.25, 0.75)), rep(b[j], 2), lwd = 0.5, col = "red")
}

J <- 300
new_dogs_ab <- invlogit(mvrnorm(J, new_dogs_mu_logit_ab, new_dogs_Sigma_ab))
a <- new_dogs_ab[, 1]
b <- new_dogs_ab[, 2]
T <- 25
new_dogs <- array(NA, c(J, T))
for (j in 1:J) {
  prev_shock <- 0
  prev_avoid <-  0
  new_dogs[j, 1] <- 1
  for (t in 2:T) {
    prev_shock = prev_shock + new_dogs[j, t - 1]
    prev_avoid = prev_avoid + 1 - new_dogs[j, t - 1]
    p = a[j] ^ prev_shock * b[j] ^ prev_avoid
    new_dogs[j, t] <- rbinom(1, 1, p)
  }
}
new_dogs_data <- list(y = new_dogs, J = J, T = T)
#| results: hide
new_fit_5 <- dogs_5$sample(data = new_dogs_data, refresh = 0)
#'
print(new_fit_5, variables = c("mu_logit_ab", "sigma_logit_ab", "Omega_logit_ab[1,2]",
                               "a[1]", "a[2]", "b[1]", "b[2]"))


T <- 50
new_dogs <- array(NA, c(J, T))
for (j in 1:J){
  prev_shock <- 0
  prev_avoid <-  0
  new_dogs[j,1] <- 1
  for (t in 2:T){
    prev_shock = prev_shock + new_dogs[j,t-1]
    prev_avoid = prev_avoid + 1 - new_dogs[j,t-1]
    p = a[j]^prev_shock * b[j]^prev_avoid
    new_dogs[j,t] <- rbinom(1, 1, p)
  }
}
new_dogs_data <- list(y = new_dogs, J = J, T = T)
#| results: hide
new_fit_5 <- dogs_5$sample(data = new_dogs_data, refresh = 0)
#'
print(new_fit_5, variables = c("mu_logit_ab", "sigma_logit_ab", "Omega_logit_ab[1,2]",
                               "a[1]", "a[2]", "b[1]", "b[2]"))

#| label: fig-dogs_data_50
#| fig-height: 4
#| fig-width: 6
par(pty = "m", mar = c(1, 2, 2, 1))
plot_dogs(new_dogs, main = "Simulated data:  50 trials", cex.main = 0.9)

#| label: fig-dogs_calibration_50
#| fig-height: 4
#| fig-width: 4
post <- as_draws_rvars(new_fit_5$draws())
par(pty = "s", mar = c(3, 3.5, 2, 1), mgp = c(2, 0.5, 0), tck = -0.01)
plot(0, 0, xlim = c(0.55, 1), ylim = c(0.55, 1), 
     xlab = "Posterior inference", ylab = "True parameter value", 
     xaxs = "i", yaxs = "i", pch = 20, cex = 0.6, 
     main = "Calibration check based on 50 trials", cex.main = 0.9)
abline(0, 1, lwd = 0.5, col = "gray")
for (j in 1:J){
  points(median(post$a[j]), a[j], pch = 20, cex = 0.6, col = "blue")
  lines(quantile(post$a[j], c(0.25, 0.75)), rep(a[j], 2), lwd = 0.5, col = "blue")
  points(median(post$b[j]), b[j], pch = 20, cex = 0.6, col = "red")
  lines(quantile(post$b[j], c(0.25, 0.75)), rep(b[j], 2), lwd = 0.5, col = "red")
}

#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' # Licenses {.unnumbered}
#'
#' * Code &copy; 2023--2025, Andrew Gelman, licensed under BSD-3.
#' * Text &copy; 2023--2025, Andrew Gelman, licensed under CC-BY-NC 4.0.

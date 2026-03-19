library("cmdstanr")
library("posterior")

normal_model <- cmdstan_model("normal.stan")
cauchy_model <- cmdstan_model("cauchy.stan")
bias_model <- cmdstan_model("bias.stan")

set.seed(123)

N <- 100
theta_ <- 2
sigma_ <- 10
y <- rnorm(N, theta_, sigma_)
data <- list(N = N, y = y)
fit_1 <- normal_model$sample(data, refresh = 0)
print(fit_1)

pdf("conflict_1a.pdf", height = 3, width = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(c(-4.1, 6.1), c(0, 0.6), xlab = expression(theta), 
     ylab = "", xaxs = "i", yaxs = "i", yaxt = "n", bty = "l", type = "n")
axis(2, lwd = 4, at = c(-1, 10), col = "white")
curve(dnorm(x, 0, 1), col = "black", n = 1000, add = TRUE)
curve(dnorm(x, 2, 1), col = "black", n = 1000, add = TRUE)
curve(dnorm(x, 1, 1 / sqrt(2)), col = "red", n = 1000, add = TRUE)
text(-1, dnorm(-1, 0, 1), "prior", col = "black")
text(3, dnorm(3, 2, 1), "likelihood", col = "black")
text(1, 0.5, "posterior", col = "red")
dev.off()

pdf("conflict_1b.pdf", height = 3, width = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(c(-4.1, 14.1), c(0, 0.6), xlab = expression(theta), 
     ylab = "", xaxs = "i", yaxs = "i", yaxt = "n", bty = "l", type = "n")
axis(2, lwd = 4, at = c(-1, 10), col = "white")
curve(dnorm(x, 0, 1), col = "black", n = 1000, add = TRUE)
curve(dnorm(x, 10, 1), col = "black", n = 1000, add = TRUE)
curve(dnorm(x, 5, 1 / sqrt(2)), col = "red", n = 1000, add = TRUE)
text(0, 0.25, "prior", col = "black")
text(10, 0.25, "likelihood", col = "black")
text(5, 0.45, "posterior", col = "red")
dev.off()


y <- (y - mean(y)) / sd(y) * 10 + 2
data <- list(N = N, y = y)
fit_2 <- cauchy_model$sample(data, refresh = 0)
print(fit_2)

y <- y - mean(y) + 10
data <- list(N = N, y = y)
fit_2 <- cauchy_model$sample(data, refresh = 0)
print(fit_2)

## Analytic calculation

post_mean_and_sd <- function(ybar){
  theta_grid <- seq(-5, ybar+5, 0.1)
  n_grid <- length(theta_grid)
  post_unnorm <- dcauchy(theta_grid, 0, 1) * dnorm(ybar, theta_grid, 1)
  post_mean <- sum(post_unnorm*theta_grid)/sum(post_unnorm)
  post_sd <- sqrt(sum(post_unnorm*(theta_grid-post_mean)^2)/sum(post_unnorm))
  c(post_mean, post_sd)
}

ybar_grid <- c(seq(0.05, 1, .05), seq(1.1, 10, .1))
post <- sapply(ybar_grid, post_mean_and_sd)
post_mean <- post[1,]
post_sd <- post[2,]
shrink <- 1 - post_mean/ybar_grid
print(post)
print(shrink)

pdf("conflict_2a.pdf", height = 3, width = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0), tck = -.01)
plot(c(-4.1, 6.1), c(0, 0.5), xlab = expression(theta), 
     ylab = "", xaxs = "i", yaxs = "i", yaxt = "n", bty = "l", type = "n")
axis(2, lwd = 4, at = c(-1, 10), col = "white")
curve(dcauchy(x, 0, 1), from = -10, to = 15, col = "black", n = 1000, add = TRUE)
curve(dnorm(x, 2, 1), col = "black", n = 1000, add = TRUE)
post_unnorm <- function(x) {
  dcauchy(x, 0, 1) * dnorm(x, 2, 1)
}
norm_const <- sum(post_unnorm(seq(-10, 10, 0.1))) * 0.1
curve(post_unnorm(x) / norm_const, col = "red", n = 1000, add = TRUE)
text(-1.1, 0.25, "prior", col = "black")
text(3.6, 0.3, "likelihood", col = "black")
text(1, 0.46, "posterior", col = "red")
dev.off()


pdf("conflict_2b.pdf", height = 3, width = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0), tck = -.01)
plot(c(-4.1, 14.1), c(0, 0.5), xlab = expression(theta), 
     ylab = "", xaxs = "i", yaxs = "i", yaxt = "n", bty = "l", type = "n")
axis(2, lwd = 4, at = c(-1, 10), col = "white")
curve(dcauchy(x, 0, 1), from = -10, to = 15, col = "black", n = 1000, add = TRUE)
curve(dnorm(x, 10, 1), col = "black", n = 1000, add = TRUE)
post_unnorm <- function(x) {
  dcauchy(x, 0, 1) * dnorm(x, 10, 1)
}
norm_const <- sum(post_unnorm(seq(-10, 25, 0.1))) * 0.1
curve(post_unnorm(x) / norm_const, col = "red", n = 1000, add = TRUE)
text(0, 0.35, "prior", col = "black")
text(12.3, 0.27, "likelihood", col = "black")
text(7.9, 0.32, "posterior", col = "red")
dev.off()

pdf("conflict_2c.pdf", height = 3, width = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(c(0, ybar_grid), c(shrink[1], shrink), 
     xlim = range(0, ybar_grid), ylim = c(0, 0.51), xaxs = "i", yaxs = "i", 
     type ="l", xlab = expression(bar(y)), ylab = "Shrinkage factor", bty = "l")
abline(0.5, 0)
text(8.2, 0.45, "Shrinkage factor for\nnormal-normal model")
text(8.2, 0.09, "Shrinkage factor for\nnormal-Cauchy model")
dev.off()

pdf("conflict_2d.pdf", height = 3, width = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(c(-rev(ybar_grid), ybar_grid), c(rev(shrink), shrink), 
     xlim = range(-ybar_grid, ybar_grid), ylim = c(0, 0.51), 
     xaxs = "i", yaxs = "i", type = "l", xlab = expression(bar(y)), 
     ylab = "Shrinkage factor", bty = "l")
abline(0.5, 0)
text(6, 0.45, "Shrinkage factor for\nnormal-normal model")
text(6, 0.06, "Shrinkage factor for\nnormal-Cauchy model")
dev.off()

y <- y - mean(y) + 2
data <- list(N = N, y = y)
fit_3a <- bias_model$sample(data, refresh = 0)
print(fit_3a)

y <- y - mean(y) + 10
data <- list(N = N, y = y)
fit_3b <- bias_model$sample(data, refresh = 0)
print(fit_3b)


ybar_grid <- seq(.5, 10, .25)
n_grid <- length(ybar_grid)
theta_mean <- rep(NA, n_grid)
bias_mean <- rep(NA, n_grid)
for (i in 1:n_grid) {
  y <- y - mean(y) + ybar_grid[i]
  data <- list(N = N, y = y)
  fit_3 <- bias_model$sample(data, iter_sampling = 1e4, refresh = 0)
  theta <- fit_3$draws(format = "df")$theta
  bias <- fit_3$draws(format = "df")$bias
  theta_mean[i] <- mean(theta)
  bias_mean[i] <- mean(bias)
}

pdf("conflict_3.pdf", height = 4, width = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(range(0, ybar_grid), range(0, theta_mean, bias_mean), xaxs = "i", yaxs = "i", 
     xlab = expression(bar(y)), ylab = "", bty = "l", type = "n")
lines(c(0, ybar_grid), c(0, theta_mean))
lines(c(0, ybar_grid), c(0, bias_mean))
text(8.2, 1, expression(paste("Posterior mean of ", theta)))
text(8.8, 6, "Posterior\nmean of bias", adj = 1)
dev.off()

library("cmdstanr")
options(mc.cores = 4)
library("posterior")
set.seed(123)

a <- 50; b <- 2; sigma <- 10
N <- 100
x <- runif(N, 0, 10)
y <- rnorm(N, a + b * x, sigma)
fake <- list(N = N, x = x, y = y)
linear <- cmdstan_model("linear.stan")
fit <- linear$sample(data = fake)
sims <- as_draws_rvars(fit$draws())
print(quantile(sims$a, c(0.1, 0.9)))
print(quantile(sims$b, c(0.1, 0.9)))
print(quantile(sims$a / sims$b, c(0.1, 0.9)))

pdf("linear_case_study_1a.pdf", width = 5, height = 4)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .3, 0), tck = -.01)
intervals <- quantile(sims$y_tilde, c(0.1, 0.9))
plot(x, y, pch = 20, cex = .5, xlab = "x", ylab = "y", bty = "l", 
     main = "Data and posterior draws of fitted regression line", 
     cex.main = .9)
for (k in sample(ndraws(sims), 50)) {
  curve(extract_variable(sims, "a")[k] + extract_variable(sims, "b")[k] * x, 
        add = TRUE, lwd = .5, col = "blue")
}
points(x, y, pch = 20, cex = .5)
dev.off()

x_tilde <- 20
sims$y_tilde <- rvar_rng(rnorm, 1, mean = sims$a + sims$b * x_tilde, sd = sims$sigma)
print(quantile(sims$y_tilde, c(0.1, 0.9)))

linear_with_pred <- cmdstan_model("linear_with_pred.stan")
fake_2 <- list(N = N, x = x, y = y, N_tilde = 1, x_tilde = 20)
fit_2 <- linear_with_pred$sample(data = fake_2)
sims_2 <- as_draws_rvars(fit_2$draws())
print(quantile(sims_2$y_tilde, c(0.1, 0.9)))

N_tilde <- 100
x_tilde <- seq(-20, 20, length = N_tilde)
fake_3 <- list(N = N, x = x, y = y, N_tilde = N_tilde, x_tilde = x_tilde)
fit_3 <- linear_with_pred$sample(data = fake_3)
sims_3 <- as_draws_rvars(fit_3$draws())
print(quantile(sims_3$y_tilde, c(0.1, 0.9)))

pdf("linear_case_study_1b.pdf", width = 5, height = 4)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .3, 0), tck = -.01)
intervals <- quantile(sims_3$y_tilde, c(0.1, 0.9))
plot(range(x_tilde), range(intervals), 
     xlab = expression(tilde(x)), 
     ylab = expression(paste("Simulated 80% predictive intervals for ", tilde(y))), 
     bty = "l", type = "n",
     main = "Data and 80% predictive intervals", cex.main = .9)
points(x, y, pch = 20, cex = .5)
for (k in 1:2){
  lines(x_tilde, intervals[k, ], col = "red")
}
dev.off()


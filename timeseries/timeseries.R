library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(scales)
library(arm)
library(posterior)
library(cmdstanr)
options(mc.cores = 4)
set.seed(123)

series <- matrix(
  scan(root("timeseries", "data", "Series1000.txt")),
  nrow = 1000,
  ncol = 135,
  byrow = TRUE
)
T <- 135
N <- 1000

pdf(root("timeseries", "series_1.pdf"), height = 3.4, width = 5.5)
par(mar = c(3, 3, 2, 0), tck = -.01, mgp = c(1.5, .5, 0))
plot(c(1, T), range(series), 
     xaxs = "i", bty = "l", type = "n", 
     xlab = "Time", ylab = "y")
for (n in 1:N){
  lines(1:T, series[n,], lwd = .5, col = alpha("black", 0.2))
}
dev.off()

slope <- rep(NA, N)
se <- rep(NA, N)
for (n in 1:N){
  data <- series[n,]
  time <- 1:T
  fit <- lm(data ~ time)
  slope[n] <- 100*coef(fit)["time"]
  se[n] <- 100*se.coef(fit)["time"]
}

pdf(root("timeseries", "series_2.pdf"), height = 3.6, width = 5.5)
par(mar = c(3, 3, 2, 0), tck = -.01, mgp = c(1.5, .5, 0))
plot(slope, se, 
     ylim = c(0, 1.05 * max(se)), yaxs = "i", bty = "l", 
     xlab = "Estmated slope", ylab = "SE", pch = 20, cex = .5)
dev.off()

pdf(root("timeseries", "series_3.pdf"), height = 3.6, width = 5.5)
par(mar = c(3, 3, 2, 0), tck = -.01, mgp = c(1.5, .5, 0))
hist(
  slope,
  xlab = "Estimated slope",
  breaks = seq(floor(10 * min(slope)), ceiling(10 * max(slope))) / 10,
  main = ""
)
dev.off()


y <- slope
K <- 3
mu <- c(0, -1, 1)
data <- list(y = y, K = K, N = N, mu = mu)
mod <- cmdstan_model(root("timeseries", "mixture.stan"))
fit_mix <- mod$sample(data = data, refresh = 0)
print(fit_mix)

## New model
mod <- cmdstan_model(root("timeseries", "mixture_2.stan"))
fit_mix <- mod$sample(data = data, refresh = 0)
prob_sims <- as_draws_rvars(fit_mix$draws())
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


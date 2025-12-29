library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
library(arm)
library(rstanarm)
library(posterior)

model <- cmdstan_model(root("misc", "chapter_07", "section_07_02", "mrp_demo.stan"))
model_quadratic <- cmdstan_model(root("misc", "chapter_07", "section_07_02", "mrp_demo_quadratic.stan"))
model_logit <- cmdstan_model(root("misc", "chapter_07", "section_07_02", "mrp_demo_logit.stan"))

set.seed(123)
N <- 200
x <- runif(N, 0, 10)
z <- rbinom(N, 1, invlogit(x - 5))

b <- c(0.1, 0.2, -3, 0.4)
sigma <- 0.5
y <- rnorm(N, b[1] + b[2] * x + b[3] * z + b[4] * x * z, sigma)

pdf(root("misc", "chapter_07", "section_07_02", "causal_poststrat_1.pdf"), 
    height = 4, width = 5)
par(mar = c(3, 3, 1, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(x, y, 
     xlab = "Pre-treatment predictor, x", ylab = "Outcome, y", 
     cex = .8, pch = ifelse(z == 0, 20, 1), bty = "l")
dev.off()

# Population info for poststratification
N_pop <- 5000
x_pop <- rnorm(N_pop, 6, 2)  # purposely choosing a population that looks different from the data
n_pop <- rep(1, N_pop)

# Fit model
data <- list(N = N, K = 4, x = x, y = y, z = z, 
             N_pop = N_pop, x_pop = x_pop, n_pop = n_pop)
fit <- model$sample(data = data, parallel_chains = 4, refresh = 0)
print(fit)

# Quadratic model
data_quadratic <- data
data_quadratic$K <- 6
fit_quadratic <- model_quadratic$sample(data = data_quadratic, parallel_chains = 4, refresh = 0)
print(fit_quadratic)

draws <- as_draws_df(fit_quadratic)

pdf(root("misc", "chapter_07", "section_07_02", "causal_poststrat_2.pdf"), 
    height = 3, width = 7.5)
par(mfrow = c(1, 2))
par(mar = c(3, 5, 1, 3), mgp = c(1.7, .5, 0), tck = -.01)
plot(draws$"b[6]", draws$SATE, 
     ylim = c(-1.6, -0.5), 
     xlab = expression(beta[6]), ylab = "SATE", 
     xaxt = "n", yaxt = "n", 
     pch = 20, cex = .2, bty = "l", col = "blue")
axis(1, c(-.05, 0, .05))
axis(2, c(-1.5, -1, -.5))
plot(draws$"b[6]", draws$PATE, 
     ylim = c(-1.1, 0), 
     xlab = expression(beta[6]), ylab = "PATE", 
     xaxt = "n", yaxt = "n", 
     pch = 20, cex = .2, bty = "l", col = "red")
axis(1, c(-.05, 0, .05))
axis(2, c(-1, -.5, 0))
dev.off()

# logit model
y_binary <- ifelse(y > 0, 1, 0)
data_logit <- data
data_logit$y <- y_binary
fit_logit <- model_logit$sample(data = data_logit, parallel_chains = 4, refresh = 0)
print(fit_logit)

# Fit and display a series of estimates

df <- data.frame(x, y_binary, z)
df_pop <- data.frame(x_pop)
df_0 <- data.frame(x, y_binary, z = 0)
df_1 <- data.frame(x, y_binary, z = 1)
df_pop_0 <- data.frame(x = x_pop, z = 0)
df_pop_1 <- data.frame(x = x_pop, z = 1)

fit <- as.list(rep(NA, 4))
fit[[1]] <- stan_glm(y_binary ~ z, family = binomial(link = "logit"), data = df, refresh = 0)
fit[[2]] <- stan_glm(y_binary ~ x + z, family = binomial(link = "logit"), data = df, refresh = 0)
fit[[3]] <- stan_glm(y_binary ~ x * z, family = binomial(link = "logit"), data = df, refresh = 0)
fit[[4]] <- stan_glm(y_binary ~ (x + I(x^2)) * z, family = binomial(link = "logit"), data = df, refresh = 0)

ate <- array(NA, c(4, 2, 2), dimnames = list(paste("Model", 1:4), c("SATE", "PATE"), c("est", "se")))
for (i in 1:length(fit)){
  cat("Model", i, "\n")
  print(fit[[i]])
  sate <- rowMeans(posterior_epred(fit[[i]], newdata = df_1) - posterior_epred(fit[[i]], newdata = df_0))
  ate[i, 1, 1] <- mean(sate)
  ate[i, 1, 2] <- sd(sate)
  pate <- rowMeans(posterior_epred(fit[[i]], newdata = df_pop_1) - posterior_epred(fit[[i]], newdata = df_pop_0))
  ate[i, 2, 1] <- mean(pate)
  ate[i, 2, 2] <- sd(pate)
}

print(ate)

pdf(root("misc", "chapter_07", "section_07_02", "causal_extrap_1.pdf"), 
    height = 4, width = 6)
par(mar = c(3, 3, 1, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(c(1, 4), range(ate[, , "est"]), 
     ylim = c(-0.5, 0.5), 
     xlab = "Model", ylab = "Estimated treatment effects", 
     xaxt = "n", bty = "l", type = "n")
axis(1, 1:4)
colors <- c("blue", "red")
x_shift <- c(-.01, .01)
label_shift <- c(-.05, .05)
for (j in 1:2){
  points(1:4 + x_shift[j], ate[, j, 1], pch = 20, col = colors[j])
  lines(1:4 + x_shift[j], ate[, j, 1], pch = 20, col = colors[j])
  for (i in 1:4){
    lines(rep(i, 2) + x_shift[j], ate[i, j, 1] + c(-1, 1) * ate[i, j, 2], lwd = 2, col = colors[j])
    lines(rep(i, 2) + x_shift[j], ate[i, j, 1] + c(-2, 2) * ate[i, j, 2], lwd = .5, col = colors[j])
  }
  text(3.5, ate[4, j, 1] + label_shift[j], dimnames(ate)[[2]][j], col = colors[j])
}
dev.off()


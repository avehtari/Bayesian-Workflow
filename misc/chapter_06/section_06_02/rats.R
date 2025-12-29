library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
options(mc.cores = 4)

logistic_1 <- cmdstan_model(root("misc", "chapter_06", "section_06_02", "logistic_1.stan"),
                            pedantic = TRUE)
logistic_2 <- cmdstan_model(root("misc", "chapter_06", "section_06_02", "logistic_2.stan"),
                            pedantic = TRUE)
logistic_3 <- cmdstan_model(root("misc", "chapter_06", "section_06_02", "logistic_3.stan"),
                            pedantic = TRUE)

set.seed(123)

rats <- read.table(root("misc", "chapter_06", "section_06_02", "rats.txt"), header = TRUE)
x <- rats$dose
n <- rats$n
y <- rats$y
rats_data <- list(J = nrow(rats), x = x, n = n, y = y)

fit_1 <- logistic_1$sample(data = rats_data, refresh = 0)
print(fit_1)

optimize_1 <- logistic_1$optimize(data = rats_data)
print(optimize_1)
mle_1 <- optimize_1$mle()

sims_1 <- fit_1$draws(format = "df")
logit <- qlogis
invlogit <- plogis

pdf(root("misc", "chapter_06", "section_06_02", "rats_1.pdf"),
    height = 3, width = 5)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0), tck = -.01)
plot(range(x), c(0, 1), 
     xlab = "dose (log g/ml)", ylab = "Pr (death)", 
     bty = "l", type = "n")
for (s in sample(nrow(sims_1), 20)) {
  curve(invlogit(sims_1$a[s] + sims_1$b[s] * x), 
        from = -1, to = 1, 
        col = "red", lwd = 0.5, 
        add = TRUE)
}
curve(invlogit(mle_1["a"] + mle_1["b"] * x), 
      from = -1, to = 1, 
      col = "blue", lwd = 2, 
      add = TRUE)
points(x, y / n, pch = 20, cex = 2)
dev.off()


fit_2 <- logistic_2$sample(
  data = c(rats_data, mu_a = 0, mu_b = 0, sigma_a = 5, sigma_b = 5), 
  refresh = 0
)
print(fit_2)
sims_2 <- fit_2$draws(format = "df")

pdf(root("misc", "chapter_06", "section_06_02", "rats_2.pdf"), 
    height = 3, width = 5)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0), tck = -.01)
plot(range(x), c(0, 1), 
     xlab = "dose (log g/ml)", ylab = "Pr (death)", 
     bty = "l", type = "n")
for (s in sample(nrow(sims_2), 20)) {
  curve(invlogit(sims_2$a[s] + sims_2$b[s] * x), 
        from = -1, to = 1, 
        col = "red", lwd = 0.5, 
        add = TRUE)
}
a_post_mean <- mean(sims_2$a)
b_post_mean <- mean(sims_2$b)
curve(invlogit(a_post_mean + b_post_mean * x), 
      from = -1, to = 1, 
      col = "blue", lwd = 2, 
      add = TRUE)
points(x, y / n, pch = 20, cex = 2)
dev.off()

fit_3 <- logistic_3$sample(
  data = c(rats_data, mu_a = 0, mu_b = 0, sigma_a = 5, sigma_b = 5), 
  refresh = 0
)
print(fit_3)
sims_3 <- fit_3$draws(format = "df")

pdf(root("misc", "chapter_06", "section_06_02", "rats_3.pdf"), 
    height = 3, width = 5)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0), tck = -.01)
plot(range(x), c(0, 1), 
     xlab = "dose (log g/ml)", ylab = "Pr (death)", 
     bty = "l", type = "n")
for (s in sample(nrow(sims_2), 20)) {
  curve(invlogit(sims_3$a[s] + sims_3$b[s] * x), 
        from = -1, to = 1, 
        col = "red", lwd = 0.5, 
        add = TRUE)
}
a_post_mean <- mean(sims_3$a)
b_post_mean <- mean(sims_3$b)
curve(invlogit(a_post_mean + b_post_mean * x), 
      from = -1, to = 1, 
      col = "blue", lwd = 2, 
      add = TRUE)
points(x, y / n, pch = 20, cex = 2)
dev.off()

x_grid <- seq(-1, 1, .01)
n_grid <-  length(x_grid)
Ey_grid <- rep(NA, n_grid)
for (i in 1:n_grid){
  Ey_grid[i] <- mean(invlogit(sims_3$a + sims_3$b * x_grid[i]))
}

pdf(root("misc", "chapter_06", "section_06_02", "rats_3_scatterplot.pdf"), 
    height = 3, width = 3)
par(pty = "s", mar = c(3, 3, 1, 1), mgp = c(1.7, .5, 0), tck = -.01, bty = "l")
plot(sims_3$a, sims_3$b, 
     xlab = "a", ylab = "b", 
     xaxt = "n", yaxt = "n", 
     pch = 20, col = "red", cex = .1)
axis(1, seq(-2, 4, 2))
axis(2, seq(0, 20, 10))
points(mean(sims_3$a), mean(sims_3$b), pch = 20, col = "blue", cex = 1)
dev.off()


pdf(root("misc", "chapter_06", "section_06_02", "rats_3_curves.pdf"), 
    height = 3, width = 5)
par(mar = c(3, 3, 1, 1), mgp = c(1.5, .5, 0), tck = -.01)
plot(range(x), c(0, 1), 
     xlab = "dose (log g/ml)", ylab = "Pr (death)", 
     bty = "l", type = "n")
b_post_mean <- mean(sims_3$b)
curve(invlogit(a_post_mean + b_post_mean * x), 
      from = -1, to = 1, 
      col = "blue", add = TRUE)
lines(x_grid, Ey_grid, col = "red")
points(x, y / n, pch = 20, cex = 2)
text(0.1, 0.85, expression(paste("logit"^{-1}, (hat(a) + hat(b) * x))), 
     adj = 1, cex = .9, col = "blue")
text(0.16, 0.72, "E(y|x,data)", adj = 0, cex = .9, col = "red")
dev.off()

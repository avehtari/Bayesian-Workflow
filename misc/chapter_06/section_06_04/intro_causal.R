library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(rstanarm)
SEED <- 123
set.seed(SEED)

## Simulate student abilities, midterm, and final exam scores
n <- 500
true_ability <- rnorm(n, 50, 20)
noise_1 <- rnorm(n, 0, 5)
noise_2 <- rnorm(n, 0, 5)
midterm <- true_ability + noise_1
final <- true_ability + noise_2
exams <- data.frame(midterm, final)
b <- 20^2 / (20^2 + 5^2)
a <- 50 - 50 * b

fit_1 <- stan_glm(final ~ midterm, data = exams, refresh = 0, seed = SEED)
print(fit_1)


pdf(root("misc", "chapter_06", "section_06_04", "students_1a.pdf"), 
    width = 4, height = 4)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, 0.5, 0), tck = -.01)
par(pty = "s")
plot(midterm, final, xlab = "Midterm exam score", ylab = "Final exam score", 
     main = "If nobody gets the treatment", cex.main = 0.9,
     xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i", 
     pch = 20, cex = 0.5)
grid(5, 5)
abline(a, b)
dev.off()

## Simulate a hypothetical experiment

z_randomized <- rbinom(n, 1, 0.5)
z <- z_randomized
theta <- 10
y <- final + theta * z
exams <- data.frame(midterm, z, y)

fit_2 <- stan_glm(y ~ midterm + z, data = exams, refresh = 0, seed = SEED)
print(fit_2)

pdf(root("misc", "chapter_06", "section_06_04", "students_1b.pdf"), 
    width = 4, height = 4)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
par(pty = "s")
plot(midterm, y, xlab = "Midterm exam score", ylab = "Final exam score",
     xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i", type = "n", 
     main = "Balanced treatment assignment", cex.main = 0.9)
points(midterm[z == 0], y[z == 0], col = "blue", pch = 20, cex = 0.5)
points(midterm[z == 1], y[z == 1], col = "red", pch = 20, cex = 0.5)
grid(5, 5)
abline(a, b, col = "blue")
abline(a + 10, b, col = "red")
dev.off()

# Difference and standard error
diff <- mean(y[z == 1]) - mean(y[z == 0])
n_0 <- sum(z == 0)
n_1 <- sum(z == 1)
se_diff <- sqrt(sd(y[z == 0])^2/n_0 + sd(y[z == 1])^2/n_1)
print(round(c(diff, se_diff), 1))

## Simulate an unbalanced treatment assigment

z_unbalanced <- rbinom(n, 1, invlogit(-(midterm - 50) / 10))
z <- z_unbalanced

pdf(root("misc", "chapter_06", "section_06_04", "students_2a.pdf"), 
    width = 4.5, height = 3.5)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
par(pty = "m")
plot(midterm, ifelse(z == 1, .99, .01), ylim = c(0, 1), 
     xlab = "Pre-test score", ylab = "Pr (z=1)", 
     xaxt = "n", yaxt = "n", yaxs = "i", bty = "l", 
     pch = 20, cex = 0.5)
axis(1, seq(0, 100, 20))
axis(2, seq(0, 1, .5))
text(45, 0.93, "(assigned to treatment group, z=1)", col = "gray30")
text(55, 0.07, "(assigned to control group, z=0)", col = "gray30")
curve(invlogit(-(x - 50) / 20), add = TRUE)
dev.off()

## Plot the data

y <- final + theta * z
exams <- data.frame(midterm, z, y)

fit_3 <- stan_glm(y ~ midterm + z, data = exams, refresh = 0, seed = SEED)
print(fit_3)

pdf(root("misc", "chapter_06", "section_06_04", "students_2b.pdf"), 
    width = 4, height = 4)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
par(pty = "s")
plot(midterm, y, xlab = "Midterm exam score", ylab = "Final exam score", 
     xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i", type = "n", 
     main = "Unalanced treatment assignment", cex.main = 0.9)
points(midterm[z == 0], y[z == 0], col = "blue", pch = 20, cex = 0.5)
points(midterm[z == 1], y[z == 1], col = "red", pch = 20, cex = 0.5)
grid(5, 5)
abline(a, b, col = "blue")
abline(a + 10, b, col = "red")
dev.off()

# Difference and standard error
diff <- mean(y[z == 1]) - mean(y[z == 0])
n_0 <- sum(z == 0)
n_1 <- sum(z == 1)
se_diff <- sqrt(sd(y[z == 0])^2/n_0 + sd(y[z == 1])^2/n_1)
print(round(c(diff, se_diff), 1))

# Regression
fit_3 <- stan_glm(y ~ midterm + z, data = exams, refresh = 0, seed = SEED)
print(fit_3)

## Nonlinear model

final <- 5 + 90 * invlogit((final - 50) / 15)

# Randomized treatment assigment
z <- z_randomized
y <- final + theta * z
exams <- data.frame(midterm, z, y)
fit_4 <- stan_glm(y ~ midterm + z, data = exams, refresh = 0, seed = SEED)
print(fit_4)

pdf(root("misc", "chapter_06", "section_06_04", "students_3a.pdf"), 
    width = 4, height = 4)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
par(pty = "s")
plot(midterm, y, xlab = "Midterm exam score", ylab = "Final exam score", 
     xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i", type = "n", 
     main = "Nonlinear model, balanced assignment", cex.main = 0.9)
points(midterm[z == 0], y[z == 0], col = "blue", pch = 20, cex = 0.5)
points(midterm[z == 1], y[z == 1], col = "red", pch = 20, cex = 0.5)
grid(5, 5)
curve(5 + 90 * invlogit(((a + b * x) - 50) / 15), col = "blue", add = TRUE)
curve(5 + 90 * invlogit(((a + b * x) - 50) / 15) + 10, col = "red", add = TRUE)
dev.off()

# Unbalanced treatment assigment
z <- z_unbalanced
y <- final + theta * z
exams <- data.frame(midterm, z, y)
fit_5 <- stan_glm(y ~ midterm + z, data = exams, refresh = 0, seed = SEED)
print(fit_5)

pdf(root("misc", "chapter_06", "section_06_04", "students_3b.pdf"), 
    width = 4, height = 4)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
par(pty = "s")
plot(midterm, y, xlab = "Midterm exam score", ylab = "Final exam score", 
     xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i", type = "n", 
     main = "Nonlinear model, unbalanced assignment", cex.main = 0.9)
points(midterm[z == 0], y[z == 0], col = "blue", pch = 20, cex = 0.5)
points(midterm[z == 1], y[z == 1], col = "red", pch = 20, cex = 0.5)
grid(5, 5)
curve(5 + 90 * invlogit(((a + b * x) - 50) / 15), col = "blue", add = TRUE)
curve(5 + 90 * invlogit(((a + b * x) - 50) / 15) + 10, col = "red", add = TRUE)
dev.off()

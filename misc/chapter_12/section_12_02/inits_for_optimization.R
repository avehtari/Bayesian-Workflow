library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(cmdstanr)
library(posterior)

# compile Stan model for logistic regression to compare to glm
logit_0 <- cmdstan_model(root("misc", "chapter_12", "section_12_02", "logit_0.stan"))

# logistic regression with just an intercept 
y <- rep(c(1, 0), c(10, 5))
glm(y ~ 1, family = binomial(link = "logit"))

# Check what happens with different starting values 
glm(y ~ 1, family = binomial(link = "logit"), start = 0)
glm(y ~ 1, family = binomial(link = "logit"), start = 2)
glm(y ~ 1, family = binomial(link = "logit"), start = -2)
glm(y ~ 1, family = binomial(link = "logit"), start = 5)

# used later in plots 
clean_estimate <- unlist(coef(glm(y ~ 1, family = binomial(link = "logit"))))

# explore dependence on initial values 
start_grid <- seq(-10, 10, 0.1)
n_grid <- length(start_grid)
est_grid <- array(NA, c(n_grid, 2))
for (i in 1:length(start_grid)) {
  fit_glm <- glm(y ~ 1, family = binomial(link = "logit"), start = start_grid[i])
  fit_stan <- logit_0$optimize(
    data = list(N = length(y), y = y),
    init = draws_df(a = start_grid[i]),
    show_messages = FALSE
  )
  est_grid[i,1] <- coef(fit_glm)
  est_grid[i,2] <- fit_stan$mle()
}

par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(start_grid, est_grid[,1], 
     xlab = "Initial value", ylab = "Parameter estimate", 
     main = "Running glm in R",
     type = "l", bty = "l", log = "y")
dev.off()

par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(start_grid, est_grid[,1], 
     xlim = c(-1.5, 2.2), ylim = clean_estimate + c(-.004, .004), 
     xlab = "Initial value", ylab = "Parameter estimate", 
     main = "Running glm in R",
     type = "l", bty = "l", yaxt = "n")
axis(2, seq(0.67, 0.71, 0.0001))



par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(start_grid, est_grid[,2], 
     ylim = clean_estimate + c(-.004, .004), 
     xlab = "Initial value", ylab = "Parameter estimate", 
     main = "Running Stan's optimizer",
     type = "l", bty = "l", yaxt = "n") 
axis(2, seq(0.67, 0.71, 0.005))
dev.off()


par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.01)
plot(start_grid, est_grid[, 2],
  xlim = c(-1.5, 2.2), ylim = clean_estimate + c(-.00008, .00008),
  xlab = "Initial value", ylab = "Parameter estimate",
  main = "Running Stan's optimizer",
  type = "l", bty = "l", yaxt = "n")
axis(2, seq(0.67, 0.71, 0.0001))
dev.off()



library("cmdstanr")
library("rstanarm")
options(mc.cores = parallel::detectCores())
library("arm")
library("posterior")
library("ellipse")
set.seed(123)

# Compile Stan programs to use throughout the file
logit_0 <- cmdstan_model("logit_0.stan")
logit_prior <- cmdstan_model("logit_prior.stan")
logit_guessing <- cmdstan_model("logit_guessing.stan")
logit_guessing_uncentered <- cmdstan_model("logit_guessing_uncentered.stan")
logit_threshold <- cmdstan_model("logit_threshold.stan")
logit_threshold_prior <- cmdstan_model("logit_threshold_prior.stan")
logit_guessing_multilevel <- cmdstan_model("logit_guessing_multilevel.stan")
logit_guessing_uncentered_multilevel <- cmdstan_model("logit_guessing_uncentered_multilevel.stan")
logit_guessing_multilevel_bivariate <- cmdstan_model("logit_guessing_multilevel_bivariate.stan")
logit_guessing_multilevel_bivariate_cholesky <- cmdstan_model("logit_guessing_multilevel_bivariate_cholesky.stan")
irt_guessing <- cmdstan_model("irt_guessing.stan")
irt_guessing_discrimination <- cmdstan_model("irt_guessing_discrimination.stan")

# Read in data and construct score for each student
responses <- read.csv("final_exam_responses.csv")
answers <- read.csv("final_exam_answers.csv")
J <- nrow(responses)  # number of students
K <- ncol(responses)  # number of items
correct <- array(NA, c(J,K))
for (k in 1:K){
  correct[,k] <- ifelse(responses[,k] == as.character(answers[k]), 1, 0)
}
score <- rowSums(correct)
score_jitt <- score + jitter(rep(0, J), amount=0.3)
score_adj <- (score - mean(score))/sd(score)
score_adj_jitt <- (score_jitt - mean(score))/sd(score)
data <- list(J=J, x=score, y=correct)

item_id_0 <- LETTERS[1:J]  # Only works here because J is no more than 26!

plots_logit <- function(pdf_name, title, xlab, model, data, xvar, guessprob){
  fit <- model$sample(data=data, refresh=0)
  print(fit)
  sims <- fit$draws(format="df")
  a_hat <- median(sims$a)
  b_hat <- median(sims$b)
  
  # Plot fitted parameters
  pdf(paste0(pdf_name, ".pdf"), height=4, width=9)
  par(mfrow=c(1,2), oma=c(0,0,2,0), mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
  plot(range(sims$a), range(sims$b), xlab="a", ylab="b", bty="l", type="n")
  points(sims$a, sims$b, col="red", pch=20, cex=.2)
  points(a_hat, b_hat, pch=20, col="blue", cex=1)
  
  # Plot fitted logistic curves
  plot(range(xvar), c(0, 1), xlab=xlab, ylab="Pr (correct answer)", yaxs="i", bty="l", type="n")
  if (guessprob==0){ 
    for (s in sample(nrow(sims), 100)) {
      curve(invlogit(sims$a[s] + sims$b[s]*x), from = min(xvar) - 1, to = max(xvar) + 1, col="red", lwd=0.5, add=TRUE)
    }
    curve(invlogit(a_hat + b_hat*x), from = min(xvar) - 1, to = max(xvar) + 1, col="blue", lwd=1.5, add=TRUE)
    text(mean(xvar), invlogit(a_hat + b_hat*mean(xvar)) - 0.2,
      paste0("y = invlogit(", pfround(a_hat, 1), " + ", pfround(b_hat, 1), " x)"), cex=.8, adj=0, col="blue")
  }
  else if (is.numeric(guessprob)) {
    for (s in sample(nrow(sims), 100)) {
      curve(guessprob + (1 - guessprob)*invlogit(sims$a[s] + sims$b[s]*x), from = min(xvar) - 1, to = max(x) + 1, col="red", lwd=0.5, add=TRUE)
    }
    curve(guessprob + (1 - guessprob)*invlogit(a_hat + b_hat*x), from = min(xvar) - 1, to = max(xvar) + 1, col="blue", lwd=1.5, add=TRUE)
    text(mean(xvar), guessprob + (1 - guessprob)*invlogit(a_hat + b_hat*mean(xvar)) - 0.2,
      paste0("y = ", guessprob, " + ", 1 - guessprob, " invlogit(", pfround(a_hat, 1), " + ", pfround(b_hat, 1), " x)"), cex=.8, adj=0, col="blue")
  }
  else if (guessprob=="estimated") {
    p0_hat <- median(sims$p0)
    for (s in sample(nrow(sims), 100)) {
      curve(sims$p0[s] + (1 - sims$p0[s])*invlogit(sims$a[s] + sims$b[s]*x), from = min(xvar) - 1, to = max(xvar) + 1, col="red", lwd=0.5, add=TRUE)
    }
    curve(p0_hat + (1 - p0_hat)*invlogit(a_hat + b_hat*x), from = min(xvar) - 1, to = max(xvar) + 1, col="blue", lwd=1.5, add=TRUE)
    text(mean(xvar), p0_hat + (1 - p0_hat)*invlogit(a_hat + b_hat*mean(xvar)) - 0.2,
     paste0("y = ", pfround(p0_hat, 2), " + ", pfround(1 - p0_hat, 2), " invlogit(", pfround(a_hat, 1), " + ", pfround(b_hat, 1), " x)"), cex=.8, adj=0, col="blue")
  }
  points(xvar, 0.5 + 0.98*(data$y - 0.5), cex=.7, pch=20)
  mtext(title, side=3, line=1, outer=TRUE)
  dev.off()
  fit
}

plots_logit_grid <- function(pdf_name, title, xlab, model, data, xvar, item_id, guessprob) {
  pdf(paste0(pdf_name, "_all.pdf"), height=6, width=11)
  par(mfrow=c(4, 6), oma=c(0,0,3.5,0),  mar=c(2.5,2.5,.5,.5), mgp=c(1.5,.5,0), tck=-.01)
  count <- 0
  for (k in order(colSums(correct))){
    print(k)
    count <- count + 1
    data_k <- data
    data_k$y <- data$y[,k]
    fit <- model$sample(data=data_k, refresh=0)
    print(fit)
    sims <- fit$draws(format="df")
    a_hat <- median(sims$a)
    b_hat <- median(sims$b)
    # Plot fitted logistic curves
    plot(range(xvar), c(0, 1), xlab = if (ceiling(count/6)==4) xlab else "", ylab = if (count%%6 == 1) "Pr (correct answer)" else "", yaxs="i", bty="l", xaxt = if (ceiling(count/6)==4) "s" else "n", yaxt="n", type="n")
    if (count%%6 == 1) axis(2, c(0, 0.5, 1))
    for (s in sample(nrow(sims), 100)) {
      curve(guessprob + (1 - guessprob)*invlogit(sims$a[s] + sims$b[s]*x), from = min(xvar) - 1, to = max(xvar) + 1, col="red", lwd=0.5, add=TRUE)
    }
    curve(guessprob + (1 - guessprob)*invlogit(a_hat + b_hat*x), from = min(xvar) - 1, to = max(xvar) + 1, col="blue", lwd=1.5, add=TRUE)
    points(xvar, 0.5 + 0.96*(data_k$y - 0.5), cex=.7, pch=20)
    mtext(paste("item", item_id[k]), side=3, line=0.5, cex=.75)
  }
  mtext(title, side=3, line=1.9, outer=TRUE)
  dev.off()
}

fit_1 <- plots_logit("final_exams_1", "Fit to item A on exam", "Score on exam", logit_0, list(J=data$J, x=data$x, y=data$y[,1]), score_jitt, guessprob=0)

# Add priors
fit_2 <- plots_logit("final_exams_2", "Fit to item A:  rescaled predictor and weakly informative prior", "Standardized exam score", logit_prior, list(J=data$J, x=data$x, y=data$y[,1], mu_a=0, sigma_a=5, mu_b=0, sigma_b=5), score_adj_jitt, guessprob=0)
plots_logit_grid("final_exams_2", "Rescaled predictor and weakly informative prior", "Standardized exam score", logit_prior, c(data, list(mu_a=0, sigma_a=5, mu_b=0, sigma_b=5)), score_adj_jitt, item_id_0, guessprob=0)
fit_2_challenge <- plots_logit("final_exams_2_challenge", "Fit to item G:  rescaled predictor and weakly informative prior", "Standardized exam score", logit_prior, list(J=data$J, x=data$x, y=data$y[,7], mu_a=0, sigma_a=5, mu_b=0, sigma_b=5), score_adj_jitt, guessprob=0)

# Fix the data problems

answers[c(5,14,17)] <- c("d", "d", "c")
correct <- array(NA, c(J,K))
for (k in 1:K){
  correct[,k] <- ifelse(responses[,k] == as.character(answers[k]), 1, 0)
}
score <- rowSums(correct)
score_jitt <- score + jitter(rep(0, J), amount=0.3)
score_adj <- (score - mean(score))/sd(score)
score_adj_jitt <- (score_jitt - mean(score))/sd(score)
data <- list(J=J, x=score, y=correct)
item_id <- rank(colSums(correct), ties="first")

plots_logit_grid("final_exams_3", "After fixing the data problem", "Standardized exam score", logit_prior, c(data, list(mu_a=0, sigma_a=5, mu_b=0, sigma_b=5)), score_adj_jitt, item_id, guessprob=0)

# Allow for guessing
plots_logit_grid("final_exams_4", "Probabilities constrained to range from 0.25 to 1", "Standardized exam score", logit_guessing, c(data, list(mu_a=0, sigma_a=5, mu_b=0, sigma_b=5)), score_adj_jitt, item_id, guessprob=0.25)

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
longdata <- list(N=N, J=J, K=K, student=student, item=item, y=y, x=score)

plots_logit_grid_2 <- function(pdf_name, title, xlab, model, data, xvar, item_id, guessprob) {
  fit <- model$sample(data, refresh=0)
  print(fit)
  a_sims <- as.matrix(fit$draws("a", format="df"))[,1:K]
  b_sims <- as.matrix(fit$draws("b", format="df"))[,1:K]
  a_hat <- apply(a_sims, 2, median)
  b_hat <- apply(b_sims, 2, median)
  se_a <- apply(a_sims, 2, mad)
  se_b <- apply(b_sims, 2, mad)
 
  pdf(paste0(pdf_name, "_all.pdf"), height=6, width=11)
  par(mfrow=c(4,6), oma=c(0,0,3.5,0),  mar=c(2.5,2.5,.5,.5), mgp=c(1.5,.5,0), tck=-.01)

  count <- 0
  for (k in order(colSums(correct))) {
    count <- count + 1
    # Plot fitted logistic curves
    plot(range(xvar), c(0, 1), xlab = if (ceiling(count/6)==4) xlab else "", ylab = if (count%%6 == 1) "Pr (correct answer)" else "", yaxs="i", bty="l", xaxt = if (ceiling(count/6)==4) "s" else "n", yaxt="n", type="n")
    if (count%%6 == 1) axis(2, c(0, 0.5, 1))
    for (s in sample(nrow(a_sims), 100)) {
      curve(guessprob + (1 - guessprob)*invlogit(a_sims[s,k] + b_sims[s,k]*x), from = min(xvar) - 1, to = max(xvar) + 1, col="red", lwd=0.5, add=TRUE)
    }
    curve(guessprob + (1 - guessprob)*invlogit(a_hat[k] + b_hat[k]*x), from = min(xvar) - 1, to = max(xvar) + 1, col="blue", lwd=1.5, add=TRUE)
    points(xvar, 0.5 + 0.96*(data$y[data$item==k] - 0.5), cex=.7, pch=20)
    mtext(paste("item", item_id[k]), side=3, line=0.5, cex=.75)
  }
  mtext(title, side=3, line=1.9, outer=TRUE)
  dev.off()
  
  pdf(paste0(pdf_name, "_scatterplot.pdf"), height=4, width=5)
  par(mar=c(2.5,2.5,.5,.5), mgp=c(1.5,.5,0), tck=-.01)
  x_rng <- range(a_hat - se_a, a_hat + se_a)
  y_rng <- range(b_hat - se_b, b_hat + se_b)
  plot(x_rng, y_rng, xlab="a", ylab="b", bty="l", type="n")
  for (k in 1:K) {
    lines(a_hat[k] + c(-1,1)*0, b_hat[k] + c(-1,1)*se_b[k], col="red", lwd=.5)
    lines(a_hat[k] + c(-1,1)*se_a[k], b_hat[k] + c(-1,1)*0, col="red", lwd=.5)
#    lines(ellipse(cov(cbind(a_sims[,k],b_sims[,k])), centre=c(a_hat[k],b_hat[k]), level=.5), col="red", lwd=.5)
  }
  text(a_hat, b_hat, item_id, col="blue", cex=.9)
  dev.off()
  fit
}

fit_5 <- plots_logit_grid_2("final_exams_5", "Multilevel model, partially pooling across the 24 exam questions", "Standardized exam score", logit_guessing_multilevel, c(longdata, list(mu_mu_a=0, sigma_mu_a=5, mu_mu_b=0, sigma_mu_b=5, mu_sigma_a=5, mu_sigma_b=5)), score_adj_jitt, item_id, guessprob=0.25)
print(fit_5, variables=c("mu_a", "sigma_a", "mu_b", "sigma_b"))

fit_6 <- plots_logit_grid_2("final_exams_6", "Multilevel model with correlation", "Standardized exam score", logit_guessing_multilevel_bivariate, c(longdata, list(mu_mu_ab=c(0,0), sigma_mu_ab=c(5,10), mu_sigma_ab=c(5,10))), score_adj_jitt, item_id, guessprob=0.25)
print(fit_6, variables=c("mu_ab", "sigma_ab", "Omega_ab"))

fit_7 <- plots_logit_grid_2("final_exams_7", "Multilevel model with correlation:  Cholesky parameterization", "Standardized exam score", logit_guessing_multilevel_bivariate_cholesky, c(longdata, list(mu_mu_ab=c(0,0), sigma_mu_ab=c(5,10), mu_sigma_ab=c(5,10))), score_adj_jitt, item_id, guessprob=0.25)
print(fit_7, variables=c("mu_ab", "sigma_ab", "Omega_ab"))

plots_irt <- function(pdf_name, title, model, data, guessprob=0.25, item_id, init=2) {
  pdf(paste0(pdf_name, ".pdf"), height=6, width=11)
  par(mfrow=c(4,6), oma=c(0,0,3,0),  mar=c(2.5,2.5,.5,.5), mgp=c(1.5,.5,0), tck=-.01)
  fit <- model$sample(data, refresh=0, init=init)
  print(fit)
  alpha_sims <- as.matrix(fit$draws("alpha", format="df"))[,1:J]
  beta_sims <- as.matrix(fit$draws("beta", format="df"))[,1:K]
  if ("gamma" %in% names(model$variables()$parameters)) {
    gamma_sims <- as.matrix(fit$draws("gamma", format="df"))[,1:K]
  }
  else {
    gamma_sims <- array(1, dim(beta_sims))
  }
  alpha_hat <- apply(alpha_sims, 2, median)
  beta_hat <- apply(beta_sims, 2, median)
  gamma_hat <- apply(gamma_sims, 2, median)
  x_range <- range(alpha_hat)
  count <- 0
  for (k in order(colSums(correct))){
    count <- count + 1
    # Plot fitted logistic curves
    plot(x_range, c(0, 1), xlab = if (ceiling(count/6)==4) "Estimated student ability" else "", ylab = if (count%%6 == 1) "Pr (correct answer)" else "", yaxs="i", bty="l", xaxt="n", yaxt="n", type="n")
    if (count%%6 == 1) axis(2, c(0, 0.5, 1))
    if (ceiling(count/6)==4) axis(1, c(-1, 0, 1))
    for (s in sample(nrow(alpha_sims), 100)) {
      curve(guessprob + (1 - guessprob)*invlogit(gamma_sims[s,k]*(x - beta_sims[s,k])), from = min(x_range) - 1, to = max(x_range) + 1, col="red", lwd=0.5, add=TRUE)
    }
    curve(guessprob + (1 - guessprob)*invlogit(gamma_hat[k]*(x - beta_hat[k])), from = min(x_range) - 1, to = max(x_range) + 1, col="blue", lwd=1.5, add=TRUE)
    points(alpha_hat, 0.5 + 0.96*(data$y[data$item==k] - 0.5), cex=.7, pch=20)
    mtext(paste("item", item_id[k]), side=3, line=0.5, cex=.75)
  }
  mtext(title, side=3, line=2, outer=TRUE)
  dev.off()
  fit
}

fit_11 <- plots_irt("final_exams_11", "Item-response model", irt_guessing, c(longdata, list(mu_mu_beta=0, sigma_mu_beta=5, mu_sigma_alpha=5, mu_sigma_beta=5)), item_id, guessprob=0.25)
fit_12 <- plots_irt("final_exams_12", "Item-response model with discrimination parameters", irt_guessing_discrimination, c(longdata, list(mu_mu_beta=0, sigma_mu_beta=5, mu_sigma_alpha=5, mu_sigma_beta=5, guessprob=0.25, mu_sigma_gamma=0.5)), item_id, guessprob=0.25)
fit_13 <- plots_irt("final_exams_13", "Item-response model with discrimination parameters", irt_guessing_discrimination, c(longdata, list(mu_mu_beta=0, sigma_mu_beta=5, mu_sigma_alpha=5, mu_sigma_beta=5, mu_sigma_gamma=0.5)), item_id, guessprob=0.25, init=0.1)

alpha_sims <- as.matrix(fit_13$draws("alpha", format="df"))[,1:J]
beta_sims <- as.matrix(fit_13$draws("beta", format="df"))[,1:K]
gamma_sims <- as.matrix(fit_13$draws("gamma", format="df"))[,1:K]
alpha_hat <- apply(alpha_sims, 2, median)
alpha_sd <- apply(alpha_sims, 2, mad)
beta_hat <- apply(beta_sims, 2, median)
beta_sd <- apply(beta_sims, 2, mad)
gamma_hat <- apply(gamma_sims, 2, median)
gamma_sd <- apply(gamma_sims, 2, mad)

pdf("irt_displays_1.pdf", height=4, width=6)
par(mar=c(3,0,0,0), mgp=c(1.5,.2,0), tck=-.01)
rng <- range(alpha_hat - 3*alpha_sd, beta_hat - 3*beta_sd, alpha_hat + 3*alpha_sd, beta_hat + 3*beta_sd)
plot(rng, c(-1,1), ylab="", yaxt="n", bty="n", type="n", xlab="Posterior distributions for student abilitites (above) and item difficulties (below)")
for (j in 1:J){
  curve(dnorm(x, alpha_hat[j], alpha_sd[j]), col="red", add=TRUE)
}
for (k in 1:K){
  curve(-dnorm(x, beta_hat[k], beta_sd[k]), col="red", add=TRUE)
}
dev.off()

pdf("irt_displays_2.pdf", height=3.2, width=4)
par(mar=c(2.5,2.5,.5,.5), mgp=c(1.5,.2,0), tck=-.01)
x_rng <- range(beta_hat - beta_sd, beta_hat + beta_sd)
y_rng <- range(gamma_hat - gamma_sd, gamma_hat + gamma_sd)
plot(x_rng, y_rng, xlab=expression(beta[k]), ylab=expression(gamma[k]), bty="l", type="n")
for (k in 1:K) {
  lines(beta_hat[k] + c(-1,1)*0, gamma_hat[k] + c(-1,1)*gamma_sd[k], col="red", lwd=.5)
  lines(beta_hat[k] + c(-1,1)*beta_sd[k], gamma_hat[k] + c(-1,1)*0, col="red", lwd=.5)
}
text(beta_hat, gamma_hat, item_id, col="blue", cex=.9)
dev.off()

prior_predictive <- function(x, x_jitt, mu_a, sigma_a, mu_b, sigma_b) {
  a <- rnorm(1, mu_a, sigma_a)
  b <- rnorm(1, mu_b, sigma_b)
  y <- rbinom(length(x), 1, invlogit(a + b*x))
  plot(range(x), c(0,1), xlab="x", ylab="y", xaxt ="n", yaxt = "n", yaxs="i", bty="l", type="n")
  axis(1, seq(-2,2,1))
  axis(2, c(0, 1))
  points(x_jitt, 0.5 + 0.96*(y - 0.5), cex=.7, pch=20, col="blue")
}

pdf("multiplechoice_prior_predictive_1.pdf", height=2.5, width=7.5)
par(oma=c(0,0,1.5,0), mfrow=c(2,5), mar=c(3,3,1,1), mgp=c(1.3,.2,0), tck=-.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 0.5, 0, 0.5)
}
mtext("10 prior predictive simulations with a ~ normal(0, 0.5) and b ~ normal(0, 0.5)", side=3, line=.5, outer=TRUE, cex=.7)
dev.off()

pdf("multiplechoice_prior_predictive_2.pdf", height=2.5, width=7.5)
par(oma=c(0,0,1.5,0), mfrow=c(2,5), mar=c(3,3,1,1), mgp=c(1.3,.2,0), tck=-.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 5, 0, 5)
}
mtext("10 prior predictive simulations with a ~ normal(0, 5) and b ~ normal(0, 5)", side=3, line=.5, outer=TRUE, cex=.7)
dev.off()

pdf("multiplechoice_prior_predictive_3.pdf", height=2.5, width=7.5)
par(oma=c(0,0,1.5,0), mfrow=c(2,5), mar=c(3,3,1,1), mgp=c(1.3,.2,0), tck=-.01)
for (loop in 1:10) {
  prior_predictive(score_adj, score_adj_jitt, 0, 50, 0, 50)
}
mtext("10 prior predictive simulations with a ~ normal(0, 50) and b ~ normal(0, 50)", side=3, line=.5, outer=TRUE, cex=.7)
dev.off()

## Breaking the model

set.seed(123)
J <- 32
x <- runif(J, 10, 20)
a_ <- -6
b_ <- 0.4
y <- rbinom(J, 1, 0.25+0.75*invlogit(a_+b_*x))
m_x <- mean(x)
s_x <- sd(x)
x_adj <- (x - m_x)/s_x

break_data <- list(J=J, x=x, y=y, mu_a=0, sigma_a=1000, mu_b=0, sigma_b=1000)
break_1_fit <- logit_guessing_uncentered$sample(data=break_data, refresh=0)
print(break_1_fit)
a <- extract_variable(break_1_fit, "a")
b <- extract_variable(break_1_fit, "b")
n_sims <- length(a)


pdf("break_1.pdf", height=3, width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(x, y, yaxs="i", bty="l", xlab="Exam score", ylab="Pr (correct answer)", type="n")
for (s in sample(n_sims, 20))
  curve(0.25 + 0.75*invlogit(a[s] + b[s]*x), lwd=.5, col="red", add=TRUE)
points(x, 0.5 + 0.985*(y - 0.5), cex=.7, pch=20)
curve(0.25 + 0.75*invlogit(a_ + b_*x), add=TRUE)
dev.off()


break_2_fit <- logit_guessing$sample(data=break_data, refresh=0)
print(break_2_fit)
a <- extract_variable(break_2_fit, "a")
b <- extract_variable(break_2_fit, "b")
n_sims <- length(a)



pdf("break_2.pdf", height=3, width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(x_adj, y, xlim=c(-2,2), yaxs="i", bty="l", xlab="Standardized exam score", ylab="Pr (correct answer)", type="n")
for (s in sample(n_sims, 20))
  curve(0.25 + 0.75*invlogit(a[s] + b[s]*x), lwd=.5, col="red", add=TRUE)
points(x_adj, 0.5 + 0.985*(y - 0.5), cex=.7, pch=20)
curve(0.25 + 0.75*invlogit(a_ + b_*(m_x + s_x*x)), add=TRUE)
dev.off()



fit_5 <- plots_logit_grid_2("final_exams_break_1", "Breaking the model", "Exam score", logit_guessing_multilevel, c(longdata, list(mu_mu_a=0, sigma_mu_a=5, mu_mu_b=0, sigma_mu_b=5, mu_sigma_a=5, mu_sigma_b=5)), score_jitt, item_id, guessprob=0.25)
print(fit_5, variables=c("mu_a", "sigma_a", "mu_b", "sigma_b"))


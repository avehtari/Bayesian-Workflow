plots_logit <- function(pdf_name, title, xlab, model, data, xvar, guessprob){
  fit <- model$sample(data = data, refresh = 0)
  print(fit)
  sims <- fit$draws(format = "df")
  a_hat <- median(sims$a)
  b_hat <- median(sims$b)
  
  # Plot fitted parameters
  pdf(paste0(pdf_name, ".pdf"), height = 4, width = 9)
  par(mfrow = c(1, 2), oma = c(0, 0, 2, 0), mar = c(3, 3, 1, 1), 
      mgp = c(1.5, .5, 0), tck = -.01)
  plot(range(sims$a), range(sims$b), 
       xlab = "a", ylab = "b", 
       bty = "l", type = "n")
  points(sims$a, sims$b, col = "red", pch = 20, cex = .2)
  points(a_hat, b_hat, pch = 20, col = "blue", cex = 1)
  
  # Plot fitted logistic curves
  plot(range(xvar), c(0, 1), 
       xlab = xlab, ylab = "Pr (correct answer)", 
       yaxs = "i", bty = "l", type = "n")
  if (guessprob == 0) {
    for (s in sample(nrow(sims), 100)) {
      curve(invlogit(sims$a[s] + sims$b[s]*x), 
            from = min(xvar) - 1, to = max(xvar) + 1, 
            col = "red", lwd = 0.5, add = TRUE)
    }
    curve(invlogit(a_hat + b_hat*x), 
          from = min(xvar) - 1, to = max(xvar) + 1, 
          col = "blue", lwd = 1.5, add = TRUE)
    text(mean(xvar), invlogit(a_hat + b_hat*mean(xvar)) - 0.2,
         paste0("y = invlogit(", pfround(a_hat, 1), " + ",  pfround(b_hat, 1), " x)"), 
         cex = .8, adj = 0, col = "blue")
  }
  else if (is.numeric(guessprob)) {
    for (s in sample(nrow(sims), 100)) {
      curve(guessprob + (1 - guessprob)*invlogit(sims$a[s] + sims$b[s]*x), 
            from = min(xvar) - 1, to = max(x) + 1, 
            col = "red", lwd = 0.5, add = TRUE)
    }
    curve(guessprob + (1 - guessprob)*invlogit(a_hat + b_hat*x), 
          from = min(xvar) - 1, to = max(xvar) + 1, 
          col = "blue", lwd = 1.5, add = TRUE)
    text(mean(xvar), guessprob + (1 - guessprob)*invlogit(a_hat + b_hat*mean(xvar)) - 0.2,
         paste0("y = ", guessprob, " + ", 1 - guessprob, " invlogit(", pfround(a_hat, 1), " + ", pfround(b_hat, 1), " x)"), 
         cex = .8, adj = 0, col = "blue")
  }
  else if (guessprob == "estimated") {
    p0_hat <- median(sims$p0)
    for (s in sample(nrow(sims), 100)) {
      curve(sims$p0[s] + (1 - sims$p0[s])*invlogit(sims$a[s] + sims$b[s]*x), 
            from = min(xvar) - 1, to = max(xvar) + 1, 
            col = "red", lwd = 0.5, add = TRUE)
    }
    curve(p0_hat + (1 - p0_hat)*invlogit(a_hat + b_hat*x), 
          from = min(xvar) - 1, to = max(xvar) + 1, 
          col = "blue", lwd = 1.5, add = TRUE)
    text(mean(xvar), p0_hat + (1 - p0_hat)*invlogit(a_hat + b_hat*mean(xvar)) - 0.2,
         paste0("y = ", pfround(p0_hat, 2), " + ", pfround(1 - p0_hat, 2), " invlogit(", pfround(a_hat, 1), " + ", pfround(b_hat, 1), " x)"), 
         cex = .8, adj = 0, col = "blue")
  }
  points(xvar, 0.5 + 0.98*(data$y - 0.5), cex = .7, pch = 20)
  mtext(title, side = 3, line = 1, outer = TRUE)
  dev.off()
  fit
}

plots_logit_grid <- function(pdf_name, title, xlab, model, data, xvar, item_id, guessprob) {
  pdf(paste0(pdf_name, "_all.pdf"), height = 6, width = 11)
  par(mfrow = c(4, 6), oma = c(0, 0, 3.5, 0), mar = c(2.5, 2.5, .5, .5), 
      mgp = c(1.5, .5, 0), tck = -.01)
  count <- 0
  for (k in order(colSums(correct))){
    print(k)
    count <- count + 1
    data_k <- data
    data_k$y <- data$y[,k]
    fit <- model$sample(data = data_k, refresh = 0)
    print(fit)
    sims <- fit$draws(format = "df")
    a_hat <- median(sims$a)
    b_hat <- median(sims$b)
    
    # Plot fitted logistic curves
    plot(range(xvar), c(0, 1), 
         xlab = if (ceiling(count/6)==4) xlab else "", 
         ylab = if (count%%6 == 1) "Pr (correct answer)" else "", 
         yaxs = "i", bty = "l",  xaxt = if (ceiling(count/6)==4) "s" else "n", 
         yaxt = "n", type = "n")
    if (count %% 6 == 1) axis(2, c(0, 0.5, 1))
    for (s in sample(nrow(sims), 100)) {
      curve(guessprob + (1 - guessprob)*invlogit(sims$a[s] + sims$b[s]*x), 
            from = min(xvar) - 1, to = max(xvar) + 1, 
            col = "red", lwd = 0.5, add = TRUE)
    }
    curve(guessprob + (1 - guessprob)*invlogit(a_hat + b_hat*x), 
          from = min(xvar) - 1, to = max(xvar) + 1, 
          col = "blue", lwd = 1.5, add = TRUE)
    points(xvar, 0.5 + 0.96*(data_k$y - 0.5), cex = .7, pch = 20)
    mtext(paste("item", item_id[k]), side = 3, line = 0.5, cex = .75)
  }
  mtext(title, side = 3, line = 1.9, outer = TRUE)
  dev.off()
}

plots_logit_grid_2 <- function(pdf_name, title, xlab, model, data, xvar, item_id, guessprob) {
  fit <- model$sample(data, refresh = 0)
  print(fit)
  a_sims <- as.matrix(fit$draws("a", format = "df"))[, 1:K]
  b_sims <- as.matrix(fit$draws("b", format = "df"))[, 1:K]
  a_hat <- apply(a_sims, 2, median)
  b_hat <- apply(b_sims, 2, median)
  se_a <- apply(a_sims, 2, mad)
  se_b <- apply(b_sims, 2, mad)
  
  pdf(paste0(pdf_name, "_all.pdf"), height = 6, width = 11)
  par(mfrow = c(4, 6), oma = c(0, 0, 3.5, 0), mar = c(2.5, 2.5, .5, .5), 
      mgp = c(1.5, .5, 0), tck = -.01)
  
  count <- 0
  for (k in order(colSums(correct))) {
    count <- count + 1
    # Plot fitted logistic curves
    plot(range(xvar), c(0, 1), 
         xlab = if (ceiling(count/6)==4) xlab else "", 
         ylab = if (count%%6 == 1) "Pr (correct answer)" else "", 
         yaxs = "i", bty = "l", xaxt = if (ceiling(count/6)==4) "s" else "n", 
         yaxt = "n", type = "n")
    if (count %% 6 == 1) axis(2, c(0, 0.5, 1))
    for (s in sample(nrow(a_sims), 100)) {
      curve(guessprob + (1 - guessprob)*invlogit(a_sims[s,k] + b_sims[s,k]*x), 
            from = min(xvar) - 1, to = max(xvar) + 1, 
            col = "red", lwd = 0.5, add = TRUE)
    }
    curve(guessprob + (1 - guessprob)*invlogit(a_hat[k] + b_hat[k]*x), 
          from = min(xvar) - 1, to = max(xvar) + 1, 
          col = "blue", lwd = 1.5, add = TRUE)
    points(xvar, 0.5 + 0.96 * (data$y[data$item == k] - 0.5), cex = .7, pch = 20)
    mtext(paste("item", item_id[k]), side = 3, line = 0.5, cex = .75)
  }
  mtext(title, side = 3, line = 1.9, outer = TRUE)
  dev.off()
  
  pdf(paste0(pdf_name, "_scatterplot.pdf"), height = 4, width = 5)
  par(mar = c(2.5, 2.5, .5, .5), mgp = c(1.5, .5, 0), tck = -.01)
  x_rng <- range(a_hat - se_a, a_hat + se_a)
  y_rng <- range(b_hat - se_b, b_hat + se_b)
  plot(x_rng, y_rng, xlab = "a", ylab = "b", bty = "l", type = "n")
  for (k in 1:K) {
    lines(a_hat[k] + c(-1,1)*0, b_hat[k] + c(-1,1)*se_b[k], col = "red", lwd = .5)
    lines(a_hat[k] + c(-1,1)*se_a[k], b_hat[k] + c(-1,1)*0, col = "red", lwd = .5)
  }
  text(a_hat, b_hat, item_id, col = "blue", cex = .9)
  dev.off()
  fit
}

plots_irt <- function(pdf_name, title, model, data, guessprob = 0.25, item_id, init = 2) {
  pdf(paste0(pdf_name, ".pdf"), height = 6, width = 11)
  par(mfrow = c(4, 6), oma = c(0, 0, 3, 0), mar = c(2.5, 2.5, .5, .5), 
      mgp = c(1.5, .5, 0), tck = -.01)
  fit <- model$sample(data, refresh = 0, init = init)
  print(fit)
  alpha_sims <- as.matrix(fit$draws("alpha", format = "df"))[, 1:J]
  beta_sims <- as.matrix(fit$draws("beta", format = "df"))[, 1:K]
  if ("gamma" %in% names(model$variables()$parameters)) {
    gamma_sims <- as.matrix(fit$draws("gamma", format = "df"))[, 1:K]
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
    plot(x_range, c(0, 1), 
         xlab = if (ceiling(count/6)==4) "Estimated student ability" else "", 
         ylab = if (count%%6 == 1) "Pr (correct answer)" else "", 
         yaxs = "i", bty = "l", xaxt = "n", yaxt = "n", type = "n")
    if (count %% 6 == 1) axis(2, c(0, 0.5, 1))
    if (ceiling(count / 6) == 4) axis(1, c(-1, 0, 1))
    for (s in sample(nrow(alpha_sims), 100)) {
      curve(guessprob + (1 - guessprob)*invlogit(gamma_sims[s,k]*(x - beta_sims[s,k])), 
            from = min(x_range) - 1, to = max(x_range) + 1, 
            col = "red", lwd = 0.5, add = TRUE)
    }
    curve(guessprob + (1 - guessprob)*invlogit(gamma_hat[k]*(x - beta_hat[k])), 
          from = min(x_range) - 1, to = max(x_range) + 1, 
          col = "blue", lwd = 1.5, add = TRUE)
    points(alpha_hat, 0.5 + 0.96*(data$y[data$item==k] - 0.5), cex = .7, pch = 20)
    mtext(paste("item", item_id[k]), side = 3, line = 0.5, cex = .75)
  }
  mtext(title, side = 3, line = 2, outer = TRUE)
  dev.off()
  fit
}

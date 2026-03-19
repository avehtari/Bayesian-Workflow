set.seed(123)

invlogit <- plogis
n <- 32
x <- rnorm(n, 0, 1)

prior_predictive <- function(x, mu_a, sigma_a, mu_b, sigma_b) {
  a <- rnorm(1, mu_a, sigma_a)
  b <- rnorm(1, mu_b, sigma_b)
  y <- rbinom(length(x), 1, invlogit(a + b*x))
  plot(range(x), c(0,1), xlab="x", ylab="y", xaxt ="n", yaxt = "n", yaxs="i", bty="l", type="n")
  axis(1, seq(-2,2,1))
  axis(2, c(0, 1))
  points(x, 0.5 + 0.96*(y - 0.5), cex=.7, pch=20, col="blue")
}

pdf("multiplechoice_prior_predictive_1.pdf", height=2.5, width=7.5)
par(oma=c(0,0,1.5,0), mfrow=c(2,5), mar=c(3,3,1,1), mgp=c(1.3,.2,0), tck=-.01)
for (loop in 1:10) {
  prior_predictive(x, 0, 0.5, 0, 0.5)
}
mtext("10 prior predictive simulations with a ~ normal(0, 0.5) and b ~ normal(0, 0.5)", side=3, line=.5, outer=TRUE, cex=.7)
dev.off()

pdf("multiplechoice_prior_predictive_2.pdf", height=2.5, width=7.5)
par(oma=c(0,0,1.5,0), mfrow=c(2,5), mar=c(3,3,1,1), mgp=c(1.3,.2,0), tck=-.01)
for (loop in 1:10) {
  prior_predictive(x, 0, 5, 0, 5)
}
mtext("10 prior predictive simulations with a ~ normal(0, 5) and b ~ normal(0, 5)", side=3, line=.5, outer=TRUE, cex=.7)
dev.off()

pdf("multiplechoice_prior_predictive_3.pdf", height=2.5, width=7.5)
par(oma=c(0,0,1.5,0), mfrow=c(2,5), mar=c(3,3,1,1), mgp=c(1.3,.2,0), tck=-.01)
for (loop in 1:10) {
  prior_predictive(x, 0, 50, 0, 50)
}
mtext("10 prior predictive simulations with a ~ normal(0, 50) and b ~ normal(0, 50)", side=3, line=.5, outer=TRUE, cex=.7)
dev.off()

library("scales")
library("cmdstanr")
set.seed(123)

series <- matrix(scan("Series1000.txt"), nrow=1000, ncol=135, byrow=TRUE)
T <- 135
N <- 1000

pdf("series_1.pdf", height=3.4, width=5.5)
par(mar=c(3,3,2,0), tck=-.01, mgp=c(1.5,.5,0))
plot(c(1,T), range(series), xaxs="i", bty="l", type="n", xlab="Time", ylab="y")
for (n in 1:N){
  lines(1:T, series[n,], lwd=.5, col=alpha("black", 0.2))
}
dev.off()

library("arm")
slope <- rep(NA, N)
se <- rep(NA, N)
for (n in 1:N){
  data <- series[n,]
  time <- 1:T
  fit <- lm(data ~ time)
  slope[n] <- 100*coef(fit)["time"]
  se[n] <- 100*se.coef(fit)["time"]
}

pdf("series_2.pdf", height=3.6, width=5.5)
par(mar=c(3,3,2,0), tck=-.01, mgp=c(1.5,.5,0))
plot(slope, se, ylim=c(0, 1.05*max(se)), yaxs="i", bty="l", xlab="Estmated slope", ylab="SE", pch=20, cex=.5)
dev.off()

pdf("series_3.pdf", height=3.6, width=5.5)
par(mar=c(3,3,2,0), tck=-.01, mgp=c(1.5,.5,0))
hist(slope, xlab="Estimated slope", breaks=seq(floor(10*min(slope)), ceiling(10*max(slope)))/10, main="")
dev.off()

y <- slope
K <- 3
mu <- c(0,-1,1)
mix <- stan("mixture.stan")
print(mix, pars=c("theta", "sigma"))

## New model
mix2 <- stan("mixture_unknown_mu.stan")
print(mix2, pars=c("mu", "theta", "sigma"))

prob_sims <- extract(mix, "p")$p
prob <- array(NA, c(N,K))
for (n in 1:N){
  for (k in 1:K){
    prob[n,k] <- mean(prob_sims[,n,k])
  }
}
max_prob <- apply(prob, 1, max)
choice <- apply(prob, 1, which.max)
expected_correct <- sum(max_prob)
sd_correct <- sqrt(sum(max_prob*(1-max_prob)))
print(table(choice))
pfround(c(expected_correct, sd_correct), 1)

pnorm(expected_correct, 900, sd_correct)
1/pnorm(expected_correct, 900, sd_correct)
1/pnorm(expected_correct, 899.5, sd_correct)

theta_true <- c(.5, .25, .25)
sigma_true <- .4
prob_simple <- array(NA, c(N,K))
for (n in 1:N){
  prob_0 <- theta_true*dnorm(y[n], mu, sigma_true)
  prob_simple[n,] <- prob_0/sum(prob_0)
}
pfround(cbind(prob, prob_simple)[1:10,], 2)

max_prob_simple <- apply(prob_simple, 1, max)
choice_simple <- apply(prob_simple, 1, which.max)
expected_correct_simple <- sum(max_prob_simple)
sd_correct_simple <- sqrt(sum(max_prob_simple*(1 - max_prob_simple)))
print(table(choice_simple))
pfround(c(expected_correct_simple, sd_correct_simple), 1)
1/pnorm(expected_correct_simple, 899.5, sd_correct_simple)

prob_simple_binary <- cbind(prob_simple[,1], prob_simple[,2] + prob_simple[,3])

max_prob_simple_binary <- apply(prob_simple_binary, 1, max)
choice_simple_binary <- apply(prob_simple_binary, 1, which.max)
expected_correct_simple_binary <- sum(max_prob_simple_binary)
sd_correct_simple_binary <- sqrt(sum(max_prob_simple_binary*(1 - max_prob_simple_binary)))
print(table(choice_simple_binary))
pfround(c(expected_correct_simple_binary, sd_correct_simple_binary), 1)

## Check what the individual series look like

noaa <- read.csv("noaa_data.csv", skip=3)

pdf("series_4.pdf", height=8, width=10)
par(mfrow=c(10,10))
par(mar=c(1,1,0,0), tck=-.01, mgp=c(1.5,.5,0))
for (i in 1:100){
  plot(c(1,T), range(series), bty="l", type="n", xlab="", ylab="", xaxt="n", yaxt="n")
  for (n in 10*(i-1) + 1:10){
    lines(1:T, series[n,], lwd=.5)
  }
  lines(1:T, noaa[1:T,2], lwd=.5, col="red")
}
dev.off()


pdf("series_5.pdf", height=4, width=5)
par(mar=c(1,1,0,0), tck=-.01, mgp=c(1.5,.5,0))
for (i in 1:1){
  plot(c(1,T), range(series), bty="l", type="n", xlab="", ylab="", xaxt="n", yaxt="n")
  for (n in 10*(i-1) + 1:10){
    lines(1:T, series[n,], lwd=.5)
  }
  lines(1:T, noaa[1:T,2], lwd=.5, col="red")
}
dev.off()


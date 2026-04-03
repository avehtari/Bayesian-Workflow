# Add Health H3IR1 S35Q1 PHYSICAL ATTRACTIVENESS OF R-W3
# How physically attractive is the respondent?
# Taken from: National Longitudinal Study of Adolescent to Adult Health (Add Health), 1994-2018 [Public Use].
# in DS8

# Raw totals:  1.7% very unattractive, 4.9% unattractive, 45.3% about average, 36.7% attractive, 11.4% very attractive, N = 4877
# Section 24:  Respondent identification number (AID), Was the baby a boy or a girl (H3LB3), 1=boy (687 responses), 2 = girl (644)
# Kanazawa study:  N = 2792 ("Wave III respondents who have had at least one biological child")

library("arm")
library("cmdstanr")
linear_0 <- cmdstan_model("sexratio_linear_0.stan", pedantic=TRUE)
linear_prior <- cmdstan_model("sexratio_linear_prior.stan", pedantic=TRUE)
set.seed(123)

x <- seq(-2,2,1)
y <- c(50, 44, 50, 47, 56)
sexratio_data <- list(N = length(x), x = x, y = y)

display(lm(y~x))

fit_0 <- linear_0$sample(data=sexratio_data, refresh=0)
print(fit_0)
sims_0 <- fit_0$draws(format="df")

optimize_0 <- linear_0$optimize(data=sexratio_data)
print(optimize_0)
mle_0 <- optimize_0$mle()

fit_1 <-  linear_prior$sample(
  data = c(sexratio_data, mu_a = 45.8, sigma_a = 0.5, mu_b = 0, sigma_b = 0.2),
  refresh = 0
)
print(fit_1)
sims_1 <- fit_1$draws(format="df")

pdf("sexratio_bayes.pdf", height=4, width=10)
par(mfrow=c(1,2), mar=c(3,3,3,2), mgp=c(1.7,.5,0), tck=-.01)
plot(x, y, ylim=c(43, 57), xlab="Attractiveness of parent", ylab="Percentage of girl babies", 
     bty="l", yaxt="n", main="Least-squares estimate",  pch=19, cex=1)
axis(2, c(45,50,55), paste(c(45,50,55), "%", sep=""))
#for (s in sample(nrow(sims_0))) {
#  curve(sims_0$a[s] + sims_0$b[s]*x, col="red", lwd=.5, add=TRUE)
#}
curve(mle_0["a"] + mle_0["b"]*x, col="blue", lwd=2, add=TRUE)
text(1, 49.2, paste("y = ", fround(mle_0["a"], 2), " + ", fround(mle_0["b"], 2), " x", sep=""), col="blue")

plot(x, y, ylim=c(43, 57), xlab="Attractiveness of parent", ylab="Percentage of girl babies", 
     bty="l", yaxt="n", main="Bayes estimate with informative prior",  pch=19, cex=1)
axis(2, c(45,50,55), paste(c(45,50,55), "%", sep=""))
#for (s in sample(nrow(sims_1))) {
#  curve(sims_1$a[s] + sims_1$b[s]*x, col="red", lwd=.5, add=TRUE)
#}
curve(median(sims_1$a) + median(sims_1$b)*x, col="blue", lwd=2, add=TRUE)
text(1, 45, paste("y = ", fround(median(sims_1$a), 2), " + ", fround(median(sims_1$b), 2), " x", sep=""), col="blue")
dev.off()

# Simulation

n <- round(2792*c(.02,.05,.45,.37,.11))
p <- 0.488
n_sims <- 1000
y_sim <- array(NA, c(n_sims, length(n)))
for (k in 1:n_sims){
  y_sim[k,] <- rbinom(length(n), n, p)
}

pdf("sexratio_rep_1.pdf", height=5, width=8)
par(mfrow=c(4,5))
par(mar=c(2.5,3,.5,.5), mgp=c(1.5,.3,0), tck=-.01)
for (k in 1:20){
  plot((-2):2, 100*y_sim[k,]/n, ylim=c(43, 57), pch=20, xaxt="n", 
       xlab=if (k>15) "Attractiveness of parent" else "", 
       ylab=if (k%%5==1) "% girls" else "", yaxt="n", bty="l")
  if (k%%5==1)
    axis(2, c(45, 50, 55), c("45%", "50%", "55%"))
  else
    axis(2, c(45, 50, 55), rep("", 3))
  if (k > 15)
    axis(1, (-2):2)
  else
    axis(1, (-2):2, rep("", 5))
  abline(lm(I(100*y_sim[k,]/n) ~ I((-2):2)), col="blue")
}
dev.off()




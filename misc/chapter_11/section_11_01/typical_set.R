# approximate typical set via simulations 
library(mvtnorm)
set.seed(123)
n_sims <- 1e4
d <- 100
theta <- rmvnorm(n_sims, rep(0, d), diag(d))
lp <- dmvnorm(theta, rep(0, d), diag(d), log=TRUE)
H <- - mean(lp)

print(quantile(lp, c(0.005, 0.995)))
typical <- lp > quantile(lp, 0.005) & lp < quantile(lp, 0.995)
theta_typical <- theta[typical,]


# plot of 99% typical set for 100-dimensional normal 
# using direct calculation on grid 
par(mar = c(3, 0, 2, 1), mgp = c(2, .5, 0), tck = -.01)
df <- 100
d <- 100
lp <- function(p, d = 100) {
  -0.5*(d*log(2*pi) + qchisq(p, d))
}
d_lp <- function(x, d = 100) {
  dchisq(-2*x - d*log(2*pi), d)
}
p_grid <- c(
  seq(0.00001, 0.0001, length = 100),
  seq(0.0001, 0.001, length = 100),
  seq(0.001, 0.9999, length = 1000)
)
lp_grid <- c(-0.5*d*log(2*pi), lp(p_grid))
d_lp_grid <- c(0, d_lp(lp(p_grid)))
plot(lp_grid, d_lp_grid, 
     ylab="", xlab=expression(paste("log p(", theta, "|y)")), 
     xaxt="n", yaxt="n", yaxs="i", 
     type="l", bty="n")
axis(1, c(-180, -162.0, -141.9, -125.6, -91.9))
p_typical <- seq(0.005, 0.995, length = 1000)
polygon(lp(c(min(p_typical), p_typical, max(p_typical))), 
        c(0, d_lp(lp(p_typical)), 0), col = "gray")
polygon(lp(c(min(p_typical), p_typical, max(p_typical))), 
        c(0, d_lp(lp(p_typical)), 0), density = 0, col = "black")
text(-141, 0.4*d_lp(-141), "99%\ntypical set")

plot_text <- function(text, ...) {
  plot(0, 0, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", type="n")
  text(0, 0, text, ...)
}

pdf("logit_priors.pdf", height=6, width=6)
par(oma=c(0,0,1,0), mfrow=c(5,3), mar=c(3,0,2,4), mgp=c(1.5,.3,0), tck=-.01)
theta_grid <- seq(-20,20,.01)
p_grid <- invlogit(theta_grid)
for (sigma in c(0.25, 0.5, 1, 2, 4)) {
  if (sigma==0.25) plot_text(expression(paste(sigma, " = ", 0.25)), adj=0, cex=1.1)
  else if (sigma==0.5) plot_text(expression(paste(sigma, " = ", 0.5)), adj=0, cex=1.1)
  else if (sigma==1) plot_text(expression(paste(sigma, " = ", 1)), adj=0, cex=1.1)
  else if (sigma==2) plot_text(expression(paste(sigma, " = ", 2)), adj=0, cex=1.1)
  else if (sigma==4) plot_text(expression(paste(sigma, " = ", 4)), adj=0, cex=1.1)
  p_p <- dlogis(theta_grid, 0, sigma)/dlogis(theta_grid, 0, 1)
  p_theta <- dlogis(theta_grid, 0, sigma)
  plot(theta_grid, p_theta, xlim=range(theta_grid), ylim=c(0, 1.05*max(p_theta)), type="l", bty="n", yaxs="i", xlab=expression(phi), yaxt="n", ylab="")
  subset <- p_grid > 0.01 & p_grid < 0.99
  plot(p_grid[subset], p_p[subset], xlim=c(0,1), ylim=c(0, 1.05*max(p_p[subset])), type="l", bty="n", yaxs="i", xlab="p", yaxt="n", ylab="")
}
mtext("                                          Logit scale                              Probability scale", 3, 0, outer=TRUE, cex=.8)         
dev.off()

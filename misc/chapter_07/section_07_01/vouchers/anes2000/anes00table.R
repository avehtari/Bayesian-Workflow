
setwd("d:/working/voucher/anes2000")
load("anes00MRP.Rdata")
K <- length(unique(releth))
tmp00 <- table(y, releth)
N00 <- sum(tmp00)
n00 <- apply(tmp00, 2, sum)
theta00 <- rep(NA, K)
for(k in 1:K){
  theta00[k] <- mean(y[releth==k], na.rm=TRUE)
}
se.theta00 <- sqrt(theta00*(1-theta00)/n00)

setwd("d:/working/voucher/anes2004")
load("anes04MRP.Rdata")
K <- length(unique(releth))
tmp04 <- table(y, releth)
N04 <- sum(tmp04)
n04 <- apply(tmp04, 2, sum)
theta04 <- rep(NA, K)
for(k in 1:K){
  theta04[k] <- mean(y[releth==k], na.rm=TRUE)
}
se.theta04 <- sqrt(theta04*(1-theta04)/n04)

tag <- c("White Catholics","White Evangelicals","White Non-Evang.\n  Protestants", "White Other/\n  No Religion", "Blacks", "Hispanics", "Other Races")


setwd("d:/working/voucher/anes2000")
pdf("rawplot.pdf", width=7, height=5)
par(mar=c(1, 8, 3, 1), mgp=c(2, 0.2, 0), tcl=-0.2, oma=c(0,0,2,0))
idx <- 1:7
plot(x=theta00, y=idx, ylim=c(7,1), xlim=c(0.34, 0.6), type="n",
  xlab="", ylab="", axes=FALSE)
axis(3, seq(0.35,0.6, 0.05), c("35%", "40%", "45%", "50%", "55%", "60%"))
text(0.25, 1:7, tag, las=2, adj=0, xpd=TRUE)
segments(theta00 + se.theta00, 1:7, theta00 - se.theta00, 1:7)
segments(theta04 + se.theta04, (idx+0.1), theta04 - se.theta04, (idx+0.1), col="gray60")
points(theta00, idx, pch=19)
points(theta04, idx+0.1, col="gray60", pch=19)
abline(v=0.5, lty=2)
mtext("Percentage who support voucher in 2000 & 2004, estimated from\nAnnenberg National Election Surveys",
  side=3, outer=TRUE, cex=1.3, font=2, line=-1)
dev.off()

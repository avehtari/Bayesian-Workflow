setwd("d:/working/voucher/anes2000")
library(arm)
library(car)
library(colorspace)

load("anes00MRP.Rdata")
est <- ybar.weighted.2way
stnames <- rownames(est)
inclab <- colnames(est)

blankplot <- function (words, xpos=.5, ypos=.5,cex=1){
  plot (c(0,1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n",type="n")
  text (xpos, ypos, words, cex=cex, xpd=TRUE)
}
  


pdf("incomeall.pdf", width=7, height=7)
a <- layout(mat=matrix(c(0, 26, 27, 28, 29, 30,
            31, 1, 2, 3, 4, 5,
            32, 6, 7, 8, 9, 10,
            33, 11, 12, 13, 14, 15,
            34, 16, 17, 18, 19, 20,
            35, 21, 22, 23, 24, 25), 6, 6, byrow=TRUE), 
            width=c(4, 5, 5, 5, 5, 5),
            height=c(2, 5, 5, 5, 5, 5))
#layout.show(a)
par(mgp=c(1.5,0.25,0), mar=c(0,0,0,0), tcl=-0.2,pty="m", oma=c(1.5,0,1,1))
for(i in 1:5){
  for(j in 1:5){
    x <- est[,i]
    y <- est[,j]
    plot(x, y, xlim=c(0,1), ylim=c(0,1), type="n",
      axes=FALSE, frame.plot=TRUE, xaxs="i", yaxs="i",
      xlab="", ylab="")
    #abline(a=0, b=1, lty=2)
    if(i==5 & j %in% c(1,2,3,4,5)){
      axis(1, c(0.25, .5, .75), c(0.25, "", .75))#, labels=FALSE)
      axis(2, c(0.25, .5, .75), labels=FALSE)
    }
    if(j==1 & i%in% c(1,2,3,4,5)){
      axis(1, c(0.25, .5, .75), labels=FALSE)
      axis(2, c(0.25, .5, .75), c(0.25, "", .75))#, labels=FALSE)
    } else {
      axis(1, c(0.25, .5, .75), labels=FALSE)
      axis(2, c(0.25, .5, .75), labels=FALSE)    
    }
    abline(lm(y~x), col=2)
    text(x, y, stnames, cex=0.7)
  }
}
par(mgp=c(1.5,0.25,0), mar=c(1,1,1,1), tcl=-0.2, pty="m")
for(i in 1:5){
  blankplot(inclab[i], ypos=0, cex=1.1)
}
for(i in 1:5){
  blankplot(inclab[i],xpos=.3, cex=1.1)
}
mtext("Estimated School Voucher Opinions by State between Income Groups", side=3, outer=TRUE, line=-1, cex=1.2, font=2)
dev.off()


pdf("incomeLowHigh.pdf", width=5, height=5)
par(mar=c(3,3,3,1), mgp=c(1.7,.25,0), tcl=-0.2)
low <- rowMeans(est[,1:2])
high <- rowMeans(est[,4:5])
plot(high, low, xlim=c(.2,.8), ylim=c(.2,.8), type="n",
      axes=FALSE, frame.plot=TRUE, xaxs="i", yaxs="i",
      xlab="Rich Support for Vouchers (>75K)", ylab="Poor Support for Vouchers (0-40K)",
      main="Support for Vouchers for Rich and Poor")
axis(1, c(0.2, .5, .8))
axis(2, c(0.2, .5, .8))
abline(lm(low ~ high), col=2)
text(high, low, stnames, cex=0.8)
dev.off()

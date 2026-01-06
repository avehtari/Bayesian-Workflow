
setwd("D:/working/voucher/anes2000/")

load("anes00reg.Rdata")
feff.00 <- fixef(fit)
feff.p00 <- fixef(fit1)
feff.n00 <- fixef(fit2)
se.feff.00 <- se.fixef(fit)
se.feff.p00 <- se.fixef(fit1)
se.feff.n00 <- se.fixef(fit2)

reff.00 <- ranef(fit)
reff.p00 <- ranef(fit1)
reff.n00 <- ranef(fit2)
se.reff.00 <- se.ranef(fit)
se.reff.p00 <- se.ranef(fit1)
se.reff.n00 <- se.ranef(fit2)


feff.m00 <- fixef(M)
feff.mp00 <- fixef(M1)
feff.mn00 <- fixef(M2)
se.feff.m00 <- se.fixef(M)
se.feff.mp00 <- se.fixef(M1)
se.feff.mn00 <- se.fixef(M2)

reff.m00 <- ranef(M)
reff.mp00 <- ranef(M1)
reff.mn00 <- ranef(M2)
se.reff.m00 <- se.ranef(M)
se.reff.mp00 <- se.ranef(M1)
se.reff.mn00 <- se.ranef(M2)


setwd("D:/working/voucher/anes2004/")
load("anes04reg.Rdata")
feff.04 <- fixef(fit)
feff.p04 <- fixef(fit1)
feff.n04 <- fixef(fit2)
se.feff.04 <- se.fixef(fit)
se.feff.p04 <- se.fixef(fit1)
se.feff.n04 <- se.fixef(fit2)

reff.04 <- ranef(fit)
reff.p04 <- ranef(fit1)
reff.n04 <- ranef(fit2)
se.reff.04 <- se.ranef(fit)
se.reff.p04 <- se.ranef(fit1)
se.reff.n04 <- se.ranef(fit2)



feff.m04 <- fixef(M)
feff.mp04 <- fixef(M1)
feff.mn04 <- fixef(M2)
se.feff.m04 <- se.fixef(M)
se.feff.mp04 <- se.fixef(M1)
se.feff.mn04 <- se.fixef(M2)

reff.m04 <- ranef(M)
reff.mp04 <- ranef(M1)
reff.mn04 <- ranef(M2)
se.reff.m04 <- se.ranef(M)
se.reff.mp04 <- se.ranef(M1)
se.reff.mn04 <- se.ranef(M2)



coefnames <- c("(Intercept)", "z.Education", "z.Income", "Repub. Vote", "Parents",
  "Pub. Employer", "Rural", "Male", "z.Democrats", "z.Age", "Income:Repub.Vote")




pdf("wholefit.pdf", width=5, height=4.5)
par(mgp=c(2,.2,0), tcl=-0.2, oma=c(0,0,2,1))
coefplot(feff.00, se.feff.00, varnames=coefnames, mar=c(0,7,2,0), main="", ylim=c(11,1), xlim=c(-1,1), col.pts=1)
coefplot(feff.04, se.feff.04, add=TRUE, offset=0.2, col.pts=2)
mtext("Individual Level Estimates with the Whole Dataset", side=3, outer=TRUE)
dev.off()

pdf("parentfit.pdf", width=7.5, height=4)
layout(mat=matrix(c(1,2,3), 1, 3, byrow=TRUE), width=c(1,2,2))
par(mgp=c(2,.2,0), tcl=-0.2, oma=c(0,0,1,1))
coefplot(feff.p00, se.feff.p00, plot=FALSE, xlim=c(0,1), ylim=c(10,1))
text(-0.5, 1:length(coefnames[-5]), coefnames[-5], adj=0, xpd=TRUE, cex=1.1)
coefplot(feff.p00, se.feff.p00, varnames="", col.pts=1, main="Parents", mar=c(0,1,4,0), xlim=c(-1,1), ylim=c(10,1))
coefplot(feff.p04, se.feff.p04, add=TRUE, col.pts=2)
coefplot(feff.n00, se.feff.n00, varnames="", col.pts=1, main="Nonparents", mar=c(0,1,4,0), xlim=c(-1,1), ylim=c(10,1))
coefplot(feff.n04, se.feff.n04, add=TRUE, col.pts=2)
mtext("Individual Level Estimates subset by parent", side=3, line=-1, outer=TRUE)
dev.off()


pdf("wholefitM.pdf", width=5, height=4)
par(mgp=c(2,.2,0), tcl=-0.2, oma=c(0,0,2,1))
coefplot(feff.m00, se.feff.m00, varnames=coefnames[-2], mar=c(0,7,2,0), main="", ylim=c(10,1), xlim=c(-1,1), col.pts=1)
coefplot(feff.m04, se.feff.m04, add=TRUE, offset=0.2, col.pts=2)
mtext("Individual Level Estimates with the Whole Dataset", side=3, outer=TRUE)
dev.off()

pdf("parentfitM.pdf", width=7.5, height=3.5)
layout(mat=matrix(c(1,2,3), 1, 3, byrow=TRUE), width=c(1,2,2))
par(mgp=c(2,.2,0), tcl=-0.2, oma=c(0,0,1,1))
coefplot(feff.mp00, se.feff.mp00, plot=FALSE, xlim=c(0,1), ylim=c(9,1))
text(-0.5, 1:length(coefnames[-c(2,5)]), coefnames[-c(2,5)], adj=0, xpd=TRUE, cex=1.1)
coefplot(feff.mp00, se.feff.mp00, varnames="", col.pts=1, main="Parents", mar=c(0,1,4,0), xlim=c(-1,1), ylim=c(9,1))
coefplot(feff.mp04, se.feff.mp04, add=TRUE, col.pts=2)
coefplot(feff.mn00, se.feff.mn00, varnames="", col.pts=1, main="Nonparents", mar=c(0,1,4,0), xlim=c(-1,1), ylim=c(9,1))
coefplot(feff.mn04, se.feff.mn04, add=TRUE, col.pts=2)
mtext("Individual Level Estimates subset by parent", side=3, line=-1, outer=TRUE)
dev.off()

#> names(reff.00)
#[1] "state"           "region * releth" "democrats"       "releth"         
#[5] "age"             "region"          "inc"             "edu"  

coefplot(reff.00[["state"]][,1], se.reff.00[["state"]][,1])

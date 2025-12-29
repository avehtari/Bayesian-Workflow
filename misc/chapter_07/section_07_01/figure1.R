# graphs for simple example in ape paper

par (mar=c(5,5,4,2)+.1)
x <- seq(-7,8,.01)
y1 <- invlogit (x-1)
y2 <- invlogit (x-2)
plot (range(x), c(0,1), xaxs="i", yaxs="i",
      xlab=expression (v),
      ylab=expression (paste ("E(y", group("|", list(u,v), ")"))),
      cex.axis=2, cex.lab=2, type="n")
lines (x, y1)
lines (x, y2)
text (2.3, .9, expression(u==1), cex=2)
text (5.5, .93, expression(u==0), cex=2)
data <- c(rnorm(20,-2,.7), rnorm(20,5,.7))
points (data, rep(.01,40), pch="|", cex=2)
ape <-mean (invlogit(data-1)-invlogit(data-2))
pe.central <- invlogit(mean(data)-1)-invlogit(mean(data)-2)
text (2.5, .5, paste ("avg pred effect =", round(ape,2)),
      adj=0, cex=2)
text (1.5, .3, paste ("pred effect at E(v) =", round(pe.central,2)),
      adj=0, cex=2)
dev.off()


par (mar=c(5,5,4,2)+.1)
x <- seq(-7,8,.01)
y1 <- invlogit (x-1)
y2 <- invlogit (x-2)
plot (range(x), c(0,1), xaxs="i", yaxs="i",
      xlab=expression (v),
      ylab=expression (paste ("E(y", group("|", list(u,v), ")"))),
      cex.axis=2, cex.lab=2, type="n")
lines (x, y1)
lines (x, y2)
text (2.3, .9, expression(u==1), cex=2)
text (5.5, .93, expression(u==0), cex=2)
data <- c(rnorm(20,1.5,.7), rnorm(20,-6,.7))
points (data, rep(.01,40), pch="|-", cex=2)
ape <-mean (invlogit(data-1)-invlogit(data-2))
pe.central <- invlogit(mean(data)-1)-invlogit(mean(data)-2)
text (2.5, .5, paste ("avg pred effect =", round(ape,2)),
      adj=0, cex=2)
text (1.5, .3, paste ("pred effect at E(v) =", round(pe.central,2)),
      adj=0, cex=2)

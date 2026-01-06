setwd("d:/working/voucher/anes2004")
library(arm)
library(car)
library(colorspace)

#=========================================================================
pal <- function(col, border = "light gray"){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = NA)
}

statemaps <- function (a, grayscale=FALSE, ...){
 library (maps)
 if (length(a)==51){
   no.dc <- c(1:8,10:51)
   a <- a[no.dc]
 }
 if (length(a)==50){
   lower48 <- state.abb!="AK" & state.abb!="HI"
   a <- a[lower48]
 }
 else if (length(a)!=48){
  stop ("wrong number of states")
 }
 mapping <- list (1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20:22,23:24,25,26,27,28,29,30,31,32,33,34:37,38:40,41,42,43,44,45,46,47,48,49,50,51,52,53:55,56:60,61,62,63)
 # for (i in 1:length(mapping)) print(regions[mapping[[i]]])
 a.long <- rep (NA, 63)
 projection <- "bonne"
 for (i in 1:48){
   a.long[mapping[[i]]] <- a[i]
 }
  if (grayscale){
    a.long.scaled <- .95*(a.long-min(a,na.rm=TRUE))/(max(a,na.rm=TRUE)-min(a,na.rm=TRUE))
    shades <- a.long.scaled
    not.dc <- !is.na(a.long.scaled)
    shades[not.dc] <- gray (shades[not.dc])
    map('state', proj=projection, param=25, lty=0, ...)
    m <- map('state', proj=projection, param=25, fill=TRUE, plot=FALSE)
    polygon(m$x,m$y, col=shades, lwd=0.5, border="gray85")
  }
  else {
    map('state', proj=projection, param=25, lty=0, ...)
    m <- map('state', proj=projection, param=25, fill=TRUE, plot=FALSE)
    polygon(m$x,m$y, col=a.long, lwd=0.5, border="gray85")
  }
} 

colorblind <- function (a, lo=min(a,na.rm=TRUE), hi=max(a,na.rm=TRUE),
                     n.shades=1001){
  a <- pmin (hi, pmax (lo, a))
  cols <-  diverge_hcl(n.shades, h = c(60, 190), c = 200, l = c(10, 90))
  #cols <- colorScale ("dark golden rod4", "light grey", "turquoise4", n.shades)
  return (cols[floor (1 + (n.shades-1)*(a-lo)/(hi-lo))])
}
#=========================================================================


dat <- read.dta("naes04.dta")
favvoucher <- dat$favvoucher#ifelse(dat$favvoucher < 3, 0,  ifelse(dat$favvoucher > 3, 1, NA))
inc <- dat$inc
eth <- dat$ethn
st <- dat$st

complete <- !is.na(st)&!is.na(inc)&!is.na(eth)&!is.na(favvoucher)
st <- st[complete]
inc <- inc[complete]
eth <- eth[complete]
y <- favvoucher[complete]
y <- ifelse(y > 3, 1, ifelse(y < 3, 0, NA))

 
vote2004 <- read.table ("2004data.dat", header=TRUE)
vote2004.long <- rbind (vote2004[1:8,], c(184,19), vote2004[9:50,])
rvote04 <- 1 - as.vector (vote2004.long[,1]/rowSums(vote2004.long))


votes08 <- read.csv ("2008ElectionResult.csv")
obama08 <- votes08[,"vote_Obama"]
mccain08 <- votes08[,"vote_McCain"]
rvote08 <- mccain08/(obama08+mccain08)

votes <- read.dta ("state vote and income, 68-00.dta")
income2000 <- votes[votes[,"st_year"]==2000, "st_income"]
state.income <- c (income2000[1:8],NA,income2000[9:50])
z.state.income <- rescale (state.income)

n.state <- 51
n.inc <- 5
n.eth <- 4
region <- c(3,4,4,3,4,4,1,1,5,3,3,4,4,2,2,2,2,3,3,1,1,1,2,2,3,2,4,2,4,1,1,4,1,3,2,2,3,4,1,1,3,2,3,3,4,1,3,4,1,2,4)
state.abb.long <- c (state.abb[1:8],"DC",state.abb[9:50])
state.name.long <- c (state.name[1:8],"Washington, D.C.",state.name[9:50])
  
namelist <- list (state.abb.long,
                  c("0-20K","20-40K","40-75K","75-150K",">150K"),
                  c("White_Catholics","White_Evangelicals","White_Non-Evang. Protestants","White_Other_or_No_Religion", "Blacks", "Hispanics", "Other_Races"),
                  c("unweighted","weighted"))
namelist2 <- list (namelist[[1]], namelist[[2]])
namelist3 <- list (namelist[[1]], namelist[[2]], namelist[[3]])

state <- rep (NA, length(y))
for (i in 1:51){
  state[st==state.abb.long[i]] <- i
}

relig <- dat$relig[complete]
relig[is.na(relig)] <- 7
bagain <- dat$bagain[complete]
bagain[is.na(bagain)] <- 0
releth <- ifelse (relig==2 & eth==1, 1,
            ifelse ((relig==1 & eth==1 & bagain==1) | (relig==4 & eth==1), 2,
              ifelse (relig==1 & eth==1 & bagain==0, 3,
                ifelse (eth==1, 4,
                  ifelse (eth==2, 5,
                    ifelse (eth==3, 6, 7))))))
n.releth <- max (releth)


n.releth <- max (releth)
pop.weight <- rep(1, length(y))

wmean <- function (a, w=rep(1,length(a)), subset=rep(TRUE,length(a))){
  keep <- !is.na(a) & !is.na(w) & !is.na(subset) & subset
  sum ((w*a)[keep])/sum(w[keep])
}


ybar.weighted <- array (NA, c(n.state,n.inc,n.releth), dimnames=namelist3)
n <- array (NA, c(n.state,n.inc,n.releth), dimnames=namelist3)
design.effect.by.cell <- array (NA, c(n.state,n.inc,n.releth), dimnames=namelist3)
for (i in (1:n.state)){
  if (!(state.abb.long[i] %in% c("AK","DC","HI"))){
    for (j in 1:n.inc){
      for (k in 1:n.releth){
        ok <- state==i & inc==j & releth==k
        ybar.weighted[i,j,k] <- wmean (y, pop.weight, ok)
        n[i,j,k] <- sum (ok)
        design.effect.by.cell[i,j,k] <- 1 + var(pop.weight[ok]/mean(pop.weight[ok]))
      }
    }
  }
}
ybar.weighted[is.nan(ybar.weighted)] <- NA
design.effect <- wmean (design.effect.by.cell, n, n>1)
n.eff <- n/design.effect

load ("prop.voters.Rdata")
prop.voters0 <- prop.voters
prop.voters <- NA + ybar.weighted
prop.voters[,,5:7] <- prop.voters0[,,2:4]
for (i in 1:n.state){
  for (j in 1:n.inc){
    ok <- state==i & inc==j
    counts <- rep (NA, 4)
    for (k in 1:4){
      counts[k] <- sum (ok & releth==k)
    }
    prop.voters[i,j,1:4] <- prop.voters0[i,j,"Whites"]*counts/sum(counts)
  }
}


cross <- function (...){
  a <- list (...)
  n.factors <- length (a)
  multiplier <- 1
  b <- 0
  for (j in 1:n.factors){
    b <- b + multiplier*a[[j]]
    multiplier <- multiplier * round (10^(ceiling(log10(max(unique(a[[j]]))))))
  }
  return (b)
}

# Put the data vectors together
count <- 0
state.v <- rep (NA, n.state*n.inc*n.releth)
inc.v <- rep (NA, n.state*n.inc*n.releth)
releth.v <- rep (NA, n.state*n.inc*n.releth)
for (k in 1:n.releth){
  for (j in 1:n.inc){
    for (i in 1:n.state){
      count <- count + 1
      state.v[count] <- i
      inc.v[count] <- j
      releth.v[count] <- k
    }
  }
}
z.inc.v <- rescale (inc.v)
z.state.inc.v <- z.state.income[state.v]
rvote04.v <- rvote04[state.v]
region.v <- region[state.v]
inc.region.v <- cross (inc.v, region.v)
inc.state.v <- cross (inc.v, state.v)
inc.releth.v <- cross (inc.v, releth.v)
releth.region.v <- cross (releth.v, region.v)
releth.state.v <- cross (releth.v, state.v)

n.eff.v <- as.vector (n.eff)
ybarw.v <- as.vector (ybar.weighted)
ybarw.v[n.eff.v==0] <- 0.5
data.v <- cbind (ybarw.v*n.eff.v, (1-ybarw.v)*n.eff.v)

threeway.fit <- glmer (data.v ~ z.inc.v*z.state.inc.v + 
  z.inc.v*rvote04.v + 
  (1 + z.inc.v | region.v) + 
  (1 + z.inc.v | state.v) + 
  (1 | inc.region.v) + 
  (1 + z.inc.v| releth.v) + 
  (1 | inc.releth.v) + 
  (1 | inc.state.v) + 
  (1 | releth.region.v) + 
  (1 | releth.state.v) + 
  (1 | inc.v) , family=quasibinomial(link="logit"))
display (threeway.fit)

fit <- rep (NA, n.state*n.inc*n.releth)
fit[!is.na(ybarw.v)] <- fitted (threeway.fit)

theta.hat.v <- fit

# Put the estimated cell means back as a three-way array
theta.hat.unshifted <- array (theta.hat.v, c(n.state,n.inc,n.releth),
                              dimnames=namelist3)

# (5) For estimating opinion on a survey question, no renormalization needed

theta.hat <- theta.hat.unshifted
  

# Create the totals, summing over all ethnic groups
ybar.weighted.2way <- array (NA, c(n.state,n.inc), dimnames=namelist2)
theta.hat.2way <- array (NA, c(n.state,n.inc), dimnames=namelist2)
n.eff.2way <- array (NA, c(n.state,n.inc), dimnames=namelist2)
for (i in 1:n.state){
  for (j in 1:n.inc){
    ybar.weighted.2way[i,j] <- wmean (ybar.weighted[i,j,], prop.voters[i,j,])
    theta.hat.2way[i,j] <- wmean (theta.hat[i,j,], prop.voters[i,j,])
    n.eff.2way[i,j] <- sum (n.eff[i,j,])
  }
}
  

# Make the maps for all voters and whites only
inclabel0 <- c ("under $30,000",
                "between $20,000 and $40,000",
                "between $40,000 and $60,000",
                "between $60,000 and $100,000",
                "over $100,000")
inclabel1 <- c("Income under $20,000", "$20-40,000", "$40-75,000", "$75-150,000", "Over $150,000" )

# Function for a blank plot
blankplot <- function (words, cex=1){
  plot (c(0,1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n",type="n")
  text (.5, .5, words, cex=cex)
}
  
prop.voters.sum <- array (NA, c(n.state, n.releth))
for (i in 1:n.state){
  for (k in 1:n.releth){
    prop.voters.sum[i,k] <- sum (prop.voters[i,,k])
  }
}



# Make the map  
#png (paste ("vouchermaps", year, ".png", sep=""), height=1000, width=1000)
pdf (paste ("vouchermaps", 2004, ".pdf", sep=""), height=14.5, width=14.5)
#====================
v.idx <- seq(1, 8)
h.idx <- c(0, seq(9, 13))
main.idx <- matrix(seq(14, 53), 8, 5)
pal.idx <- c(0,rep(54, 5))
foot.idx <- c(0, rep(55, 5))
total.mat <- rbind(c(h.idx), cbind(v.idx, main.idx), pal.idx, foot.idx)
par (mar=c(0,0,0,0), oma=c(0,0,3,0))
layout(total.mat, width=c(3, rep(4, 5)), height=c(1,rep(4,8),1.5))
#=====================
national.avg <- wmean(theta.hat,prop.voters)
for (k in 0:7){
  blankplot (c("All voters", "White\nCatholics", "White evangelical\nProtestants", "White non-evang.\nProtestants", "White other/\nno religion", "Blacks", "Hispanics", "Other races")[k+1], cex=2)
}
for(j in 1:n.inc){
  blankplot(inclabel1[j], cex=2)
}
for (j in 1:n.inc){
  for (k in 0:7){
    if (k==0){
      statemaps (colorblind (-(theta.hat.2way[,j]-national.avg), lo=-.25, hi=.25))
    }
    else {
      statemaps (ifelse (prop.voters[,j,k]>.01, colorblind (-(theta.hat[,j,k] - national.avg), lo=-.25, hi=.25), "white"))
    }
  }
}
par(mar=c(2,0,1,0))
pal( diverge_hcl(1001, h = c(190, 60), c = 200, l = c(10, 90)))
hi <- max (0, 0.25, na.rm=TRUE)
lo <- min (0, -0.25, na.rm=TRUE)
text(1,   -1, paste (round ((national.avg + hi)*100), "%", sep=""), xpd=TRUE, cex=2)
text(0.5, -1, paste (round(national.avg*100), "%", sep=""), xpd=TRUE, cex=2)
text(0,   -1, paste (round ((national.avg + lo)*100), "%", sep=""), xpd=TRUE, cex=2)
 
#text(1, -1, "Yes", xpd=TRUE, cex=2)
#text(0, -1, "No", xpd=TRUE, cex=2)
par(mar=c(0,0,0,0))
blankplot ("The state is left blank where a category represents less than 1% of the voters of a state", cex=1.8)
mtext("2004: Do you support school vouchers?", side=3, line=0.5, xpd=TRUE, cex=2, outer=TRUE)
#mtext (paste (year, ":  State-level support (orange) or opposition (green) on school vouchers, relative to the national average of ", round (100*national.avg), "% support", sep=""), side=3, outer=TRUE, line=1, cex=1.3)
dev.off()




# (7) Make graphs for the 50 states

for (k in 1:7){
  pdf (paste ("voucherplots", 2004, " ", namelist[[3]][k], ".pdf", sep=""), height=15,width=15)
  par (mar=c(2,2,1,0), tck=0, mgp=c(1.5,.5,0), oma=c(0,0,4,.5))
  graph.dims <- c(7,7)
  par (mfrow=graph.dims)
  sort.state <- rev (order (rvote04))
  count <- 0
  for (i in sort.state){
    if (!(state.abb.long[i] %in% c("AK","HI","DC"))){
      count <- count + 1
      if (!is.na(prop.voters.sum[i,k])&prop.voters.sum[i,k] > .1){
        plot (c(1,5), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n",
              yaxs="i")
        if (count%%graph.dims[2]==1)
          axis (2, c(.02,.5,.95), c("0","50%","100%"), cex.axis=1.1)
        if (count > (48-graph.dims[2]))
          axis (1, c(1.2,3,4.8), c("poor","mid","rich"), cex.axis=1.2)
        abline (.5, 0, col="gray", lwd=.5)
        fit <- theta.hat[i,,k]
        pts <- ybar.weighted[i,,k]
        neff <- n.eff[i,,k]
        se <- sqrt (fit*(1-fit)/neff)
        text (3, .9, state.name.long[i], cex=1.3)
        for (j in 1:n.inc){
          lines (rep (j, 2), pts[j] + c(-1,1)*se[j], lwd=.5)
          points (j, min (.98, max (.02, pts[j])), pch=20, cex=1.5)
        }
        lines (1:n.inc, fit, lwd=2)
      }
      else {
        plot (c(1,5), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n", yaxs="i", type="n")
        text (3, .9, state.name.long[i], cex=1.3)
      }
    }
  }
  mtext (paste (2004, ":  Raw and estimated support for school vouchers within each state among", namelist[[3]][k]), outer=TRUE, cex=1, side=3, line=1)
  dev.off()
}

## Reverse-engineering the Fivethirtyeight forecast
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(rjson)

sims_538 <- fromJSON(file= root("misc", "chapter_08", "section_08_01", "simmed-maps.json"))
states <- sims_538$states
n_sims <- length(sims_538$maps)
sims <- array(NA, c(n_sims, 59), dimnames=list(NULL, c("", "Trump", "Biden", states)))
for (i in 1:n_sims){
  sims[i,] <- sims_538$maps[[i]]
}
state_sims <- sims[,4:59]
trump_share <- (state_sims/100 + 1)/2
biden_wins <- trump_loses <- state_sims < 0
trump_wins <- biden_loses <- state_sims > 0

weird_outcome <- (sims[,"AL"] < 0) | (sims[,"AR"] < 0) | (sims[,"KY"] < 0) | (sims[,"LA"] < 0) | (sims[,"IN"] < 0) | (sims[,"KS"] < 0) | (sims [,"ND"] < 0)

round(apply(biden_wins, 2, mean), 2)

condition <- biden_loses[,"NJ"]
round(apply(biden_wins[condition,], 2, mean), 2)

## New Jersey and Alaska

round(mean(biden_wins[,"AK"][biden_loses[,"NJ"]]), 2)
round(mean(biden_wins[,"AK"][biden_wins[,"NJ"]]), 2)

par(mar=c(3,3,1,1), mgp=c(1.7, .5, 0), tck=-.01)
par(pty="s")
rng <- range(trump_share[,c("NJ", "AK")])
plot(rng, rng, xlab="Trump vote share in New Jersey", ylab="Trump vote share in Alaska", main="40,000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[,"NJ"], trump_share[,"AK"], pch=20, cex=0.1)
text(0.65, 0.25, "Trump wins NJ", col="darkred", cex=0.8)
text(0.35, 0.25, "Trump loses NJ", col="black", cex=0.8)

subset <- sample(n_sims, 1000)
rng <- range(trump_share[,c("NJ", "AK")])
plot(rng, rng, xlab="Trump vote share in New Jersey", ylab="Trump vote share in Alaska", main="Only 1000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[subset,"NJ"], trump_share[subset,"AK"], pch=20, cex=0.1)
text(0.65, 0.25, "Trump wins NJ", col="darkred", cex=0.8)
text(0.35, 0.25, "Trump loses NJ", col="black", cex=0.8)

round(cor(trump_share[,"AK"], trump_share[,"NJ"]), 2)

## New Jersey and Pennsylvania

round(mean(biden_wins[,"PA"][biden_loses[,"NJ"]]), 2)
round(mean(biden_wins[,"PA"][biden_wins[,"NJ"]]), 2)

par(mar=c(3,3,1,1), mgp=c(1.7, .5, 0), tck=-.01)
par(pty="s")
rng <- range(trump_share[,c("NJ", "PA")])
plot(rng, rng, xlab="Trump vote share in New Jersey", ylab="Trump vote share in Pennsylvania", main="40,000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[,"NJ"], trump_share[,"PA"], pch=20, cex=0.1)
text(0.6, 0.25, "Trump wins NJ", col="darkred", cex=0.8)
text(0.4, 0.25, "Trump loses NJ", col="black", cex=0.8)

subset <- sample(n_sims, 1000)
rng <- range(trump_share[,c("NJ", "PA")])
plot(rng, rng, xlab="Trump vote share in New Jersey", ylab="Trump vote share in Pennsylvania", main="Only 1000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[subset,"NJ"], trump_share[subset,"PA"], pch=20, cex=0.1)
text(0.6, 0.25, "Trump wins NJ", col="darkred", cex=0.8)
text(0.4, 0.25, "Trump loses NJ", col="black", cex=0.8)

round(cor(trump_share[,"PA"], trump_share[,"NJ"]), 2)


## Alabama and Mississippi

round(mean(biden_wins[,"MS"][biden_loses[,"AL"]]), 2)
round(mean(biden_wins[,"MS"][biden_wins[,"AL"]]), 2)

par(mar=c(3,3,1,1), mgp=c(1.7, .5, 0), tck=-.01)
par(pty="s")
rng <- range(trump_share[,c("AL", "MS")])
plot(rng, rng, xlab="Trump vote share in Alabama", ylab="Trump vote share in Mississippi", main="40,000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[,"AL"], trump_share[,"MS"], pch=20, cex=0.1)
text(0.6, 0.3, "Trump wins AL", col="darkred", cex=0.8)
text(0.4, 0.3, "Trump loses AL", col="black", cex=0.8)

subset <- sample(n_sims, 1000)
rng <- range(trump_share[,c("AL", "MS")])
plot(rng, rng, xlab="Trump vote share in Alabama", ylab="Trump vote share in Mississippi", main="Only 1000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[subset,"AL"], trump_share[subset,"MS"], pch=20, cex=0.1)
text(0.6, 0.3, "Trump wins AL", col="darkred", cex=0.8)
text(0.4, 0.3, "Trump loses AL", col="black", cex=0.8)

round(cor(trump_share[,"MS"], trump_share[,"AL"]), 2)


## Wisconsin and Pennsylvania

round(mean(trump_wins[,"PA"][trump_wins[,"WI"]]), 2)
round(mean(trump_wins[,"PA"][biden_wins[,"WI"]]), 2)

par(mar=c(3,3,1,1), mgp=c(1.7, .5, 0), tck=-.01)
par(pty="s")
rng <- range(trump_share[,c("WI", "PA")])
plot(rng, rng, xlab="Trump vote share in Wisconsin", ylab="Trump vote share in Pennsylvania", main="40,000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[,"WI"], trump_share[,"PA"], pch=20, cex=0.1)
text(0.55, 0.35, "Trump wins WI", col="darkred", cex=0.8)
text(0.45, 0.35, "Trump loses WI", col="black", cex=0.8)

subset <- sample(n_sims, 1000)
rng <- range(trump_share[,c("WI", "PA")])
plot(rng, rng, xlab="Trump vote share in Wisconsin", ylab="Trump vote share in Pennsylvania", main="Only 1000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[subset,"WI"], trump_share[subset,"PA"], pch=20, cex=0.1)
text(0.6, 0.35, "Trump wins WI", col="darkred", cex=0.8)
text(0.4, 0.35, "Trump loses WI", col="black", cex=0.8)

round(cor(trump_share[,"PA"], trump_share[,"WI"]), 2)

## Look at some correlations

some_states <- c("AK","WA","WI","OH","PA","NJ","VA","GA","FL","AL","MS")
cor_mat <- cor(trump_share[,some_states])
image(cor_mat[,rev(1:nrow(cor_mat))], xaxt="n", yaxt="n")
axis(1, seq(0, 1, length=length(some_states)), some_states, tck=0, cex.axis=0.8)
axis(2, seq(0, 1, length=length(some_states)), rev(some_states), tck=0, cex.axis=0.8, las=1)

round(cor(trump_share[,"WA"], trump_share[,"MS"]), 2)


## Washington and Mississippi
round(mean(trump_wins[,"MS"] [trump_wins[,"WA"]]), 2)
round(mean(trump_wins[,"MS"] [biden_wins[,"WA"]]), 2)

round(mean(biden_wins[,"MS"][biden_loses[,"WA"]]), 2)
round(mean(biden_wins[,"MS"][biden_wins[,"WA"]]), 2)

par(mar=c(3,3,1,1), mgp=c(1.7, .5, 0), tck=-.01)
par(pty="s")
rng <- range(trump_share[,c("WA", "MS")])
plot(rng, rng, xlab="Trump vote share in Washington", ylab="Trump vote share in Mississippi", main="40,000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[,"WA"], trump_share[,"MS"], pch=20, cex=0.1)
text(0.65, 0.3, "Trump wins WA", col="darkred", cex=0.8)
text(0.35, 0.3, "Trump loses WA", col="black", cex=0.8)

subset <- sample(n_sims, 1000)
rng <- range(trump_share[,c("WA", "MS")])
plot(rng, rng, xlab="Trump vote share in Washington", ylab="Trump vote share in Mississippi", main="Only 1000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[subset,"WA"], trump_share[subset,"MS"], pch=20, cex=0.1)
text(0.65, 0.3, "Trump wins WA", col="darkred", cex=0.8)
text(0.35, 0.3, "Trump loses WA", col="black", cex=0.8)

round(cor(trump_share[,"MS"], trump_share[,"WA"]), 2)

round(trump_share[1:10,c("MS","WA")], 2)

## Check first 5 lines of data

check <- rbind(sims_538$maps[[1]], sims_538$maps[[2]], sims_538$maps[[3]], sims_538$maps[[4]], sims_538$maps[[5]])
colnames(check) <- c(rep("", 3), sims_538$states)
round(check)

round(mean(cor(trump_share)), 2)

## New York and Pennsylvania

round(mean(trump_wins[,"PA"][trump_wins[,"NY"]]), 2)
round(mean(trump_wins[,"PA"][biden_wins[,"NY"]]), 2)
round(mean(trump_wins[,"PA"][trump_share[,"NY"] > 0.45]), 2)


par(mar=c(3,3,1,1), mgp=c(1.7, .5, 0), tck=-.01)
par(pty="s")
rng <- range(trump_share[,c("NY", "PA")])
plot(rng, rng, xlab="Trump vote share in New York", ylab="Trump vote share in Pennsylvania", main="40,000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[,"NY"], trump_share[,"PA"], pch=20, cex=0.1)
text(0.56, 0.25, "Trump wins NY", col="darkred", cex=0.8)
text(0.35, 0.25, "Trump loses NY", col="black", cex=0.8)

subset <- sample(n_sims, 1000)
rng <- range(trump_share[,c("NJ", "AK")])
plot(rng, rng, xlab="Trump vote share in New Jersey", ylab="Trump vote share in Alaska", main="Only 1000 simulation draws", cex.main=0.9, bty="l", type="n")
polygon(c(0.5,0.5,1,1), c(0,1,1,0), border=NA, col="pink")
points(trump_share[subset,"NJ"], trump_share[subset,"AK"], pch=20, cex=0.1)
text(0.65, 0.25, "Trump wins NJ", col="darkred", cex=0.8)
text(0.35, 0.25, "Trump loses NJ", col="black", cex=0.8)

round(cor(trump_share[,"AK"], trump_share[,"NJ"]), 2)


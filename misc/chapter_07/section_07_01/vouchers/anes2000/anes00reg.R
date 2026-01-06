setwd("D:/working/voucher/anes2000/")
library(arm)
library(foreign)
library(car)




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

naes00 <- read.dta("naes00new.dta")
vote2004 <- read.table ("2004data.dat", header=TRUE)
vote2004.long <- rbind (vote2004[1:8,], c(184,19), vote2004[9:50,])
rvote04 <- 1 - as.vector (vote2004.long[,1]/rowSums(vote2004.long))
state.abb.long <- c (state.abb[1:8],"DC",state.abb[9:50])
state.name.long <- c (state.name[1:8],"Washington, D.C.",state.name[9:50])
n <- length(naes00[,1])

votes <- read.dta ("state vote and income, 68-00.dta")
rvote00 <- votes[votes[,"st_year"]==2000, "st_repshare"]
state.rvote <- c (rvote[1:8],NA,rvote[9:50])




region <- c(3,4,4,3,4,4,1,1,5,3,3,4,4,2,2,2,2,3,3,1,1,1,2,2,3,2,4,2,4,1,1,4,1,3,2,2,3,4,1,1,3,2,3,3,4,1,3,4,1,2,4)

temp <- rep(NA, n)
for(r in 1:5){
  idx <- grep(r, region)
  temp[naes00$stnum %in% idx] <- r
}

region <- temp
inc <- naes00$inc
z.inc <- rescale(naes00$inc)
y <- naes00$favvoucher
state <- naes00$stnum
relig <- naes00$relig
eth <- naes00$ethn
relig[is.na(relig)] <- 7
bagain <- naes00$bagain
bagain[is.na(bagain)] <- 0
releth <- ifelse (relig==2 & eth==1, 1,
              ifelse ((relig==1 & eth==1 & bagain==1) | (relig==4 & eth==1), 2,
                ifelse (relig==1 & eth==1 & bagain==0, 3,
                  ifelse (eth==1, 4,
                    ifelse (eth==2, 5,
                      ifelse (eth==3, 6, 7))))))
rvote04 <- rvote04[naes00$stnum]
rvote00 <- rvote00[naes00$stnum]
parents <- naes00$parents
pubempl <- naes00$pubempl
rural <- naes00$rural
male <- naes00$male
democrats <- naes00$pid
z.dem <- rescale(democrats)#recode(democrats, "1=-3;2=-2;3=-1;4=0;5=1;6=2;7=3")
age <- naes00$age
z.age <- rescale(age)
edu <- naes00$edu
z.edu <- rescale(edu)
married <- naes00$married
z.releth <- rescale(releth)
regionreleth <- cross(region, releth)
region2 <- ifelse(rvote00 > 0.55, "rep", 
            ifelse(rvote00 < 0.45, "dem",
            ifelse(rvote00 >=0.45 & rvote00<=0.55, "mid",
              NA)))

fit <- glmer(y ~ z.edu + z.inc + rvote + parents + pubempl + rural + male + 
  z.dem + z.age + z.inc*rvote +  
  (1 + z.inc | releth) + (1|inc) + (1|state) + (1|region) + 
  (1|democrats) + (1|age) + (1|edu) + 
  (1|regionreleth), family=binomial)

  
fit1 <- glmer(y ~ z.edu + z.inc + rvote + pubempl + rural + male + 
  z.dem + z.age + z.inc*rvote +  
  (1 + z.inc | releth) + (1|inc) + (1|state) + (1|region) + 
  (1|democrats) + (1|age) + (1|edu) +
  (1|regionreleth), family=binomial, subset=parents==1)

fit2 <- glmer(y ~ z.edu + z.inc + rvote + pubempl + rural + male + 
  z.dem + z.age + z.inc*rvote +  
  (1 + z.inc | releth) + (1|inc) + (1|state) + (1|region) + 
  (1|democrats) + (1|age) + (1|edu) +
  (1|regionreleth), family=binomial, subset=parents==0)

M <- glmer(y ~ z.inc + rvote + parents + pubempl + rural + male + 
  z.dem + z.age + z.inc*rvote +  
  (1 + z.inc | releth) + (1|inc) + (1|state) + (1|region) + 
  (1|democrats) + (1|age) + 
  (1|regionreleth), family=binomial)

M1 <- glmer(y ~ z.inc + rvote + pubempl + rural + male + 
  z.dem + z.age + z.inc*rvote +  
  (1 + z.inc | releth) + (1|inc) + (1|state) + (1|region) + 
  (1|democrats) + (1|age) + 
  (1|regionreleth), family=binomial, subset=parents==1)

M2 <- glmer(y ~ z.inc + rvote + pubempl + rural + male + 
  z.dem + z.age + z.inc*rvote +  
  (1 + z.inc | releth) + (1|inc) + (1|state) + (1|region) + 
  (1|democrats) + (1|age) + 
  (1|regionreleth), family=binomial, subset=parents==0)

MM <- glmer(y ~ z.edu + z.inc + rvote04 + pubempl + rural + male + z.dem + z.age + parents +
  (1|state) + (1+z.inc|releth), family=binomial)



MM <- glmer(y ~ pubempl + rural + male + z.age + z.dem*parents + (z.dem*parents|inc) + (z.dem*parents|releth), family=binomial)


save.image("anes00reg.Rdata")

load("anes00reg.Rdata")

setwd("D:/working/voucher/anes2004/")
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


naes04 <- read.dta("naes04new.dta")
votes08 <- read.csv ("2008ElectionResult.csv")
obama08 <- votes08[,"vote_Obama"]
mccain08 <- votes08[,"vote_McCain"]
rvote08 <- mccain08/(obama08+mccain08)
state.abb.long <- c (state.abb[1:8],"DC",state.abb[9:50])
state.name.long <- c (state.name[1:8],"Washington, D.C.",state.name[9:50])
n <- length(naes04[,1])

region <- c(3,4,4,3,4,4,1,1,5,3,3,4,4,2,2,2,2,3,3,1,1,1,2,2,3,2,4,2,4,1,1,4,1,3,2,2,3,4,1,1,3,2,3,3,4,1,3,4,1,2,4)

temp <- rep(NA, n)
for(r in 1:5){
  idx <- grep(r, region)
  temp[naes04$stnum %in% idx] <- r
}

region <- temp
inc <- naes04$inc
z.inc <- rescale(naes04$inc)
y <- naes04$favvoucher
y <- ifelse(y > 3, 1, ifelse(y < 3, 0, NA))
state <- naes04$stnum
relig <- naes04$relig
eth <- naes04$ethn
relig[is.na(relig)] <- 7
bagain <- naes04$bagain
bagain[is.na(bagain)] <- 0
releth <- ifelse (relig==2 & eth==1, 1,
              ifelse ((relig==1 & eth==1 & bagain==1) | (relig==4 & eth==1), 2,
                ifelse (relig==1 & eth==1 & bagain==0, 3,
                  ifelse (eth==1, 4,
                    ifelse (eth==2, 5,
                      ifelse (eth==3, 6, 7))))))
rvote <- rvote08[naes04$stnum]
parents <- naes04$parents
pubempl <- naes04$pubempl
rural <- naes04$rural
male <- naes04$male
democrats <- naes04$pid
z.dem <- rescale(democrats)#recode(democrats, "1=-3;2=-2;3=-1;4=0;5=1;6=2;7=3")
age <- naes04$age
z.age <- rescale(age)
edu <- naes04$educ
z.edu <- rescale(edu)



regionreleth <- cross(region, releth)

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




save.image("anes04reg.Rdata")

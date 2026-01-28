births <- read.table("births.txt", header=TRUE)
mean_age_45_54 <- function(yr){
  ages <- 45:54
  ok <- births$year %in% (yr - ages)
  return(sum(births$births[ok]*rev(ages))/sum(births$births[ok]))
}
for (yr in 1989:2015) print(mean_age_45_54(yr))

## Calculation
print((.5/10)* (.006423 - .003064)/.003064)

## from life table

deathpr_by_age <- c(.003064, .003322, .003589, .003863, .004148, .004458, .004800, .005165, .005554, .005971)

deathpr_male <- c(.003244, .003571, .003926, .004309, .004719, .005156, .005622, .006121, .006656, .007222)
deathpr_female <- c(.002069, .002270, .002486, .002716, .002960, .003226, .003505, .003779, .004040, .004301)
                    
## sum it up

pop <- read.csv("US-EST00INT-ALLDATA.csv")
years <- 1989:2013
deathpr_1 <- rep(NA, length(years))
deathpr_2 <- rep(NA, length(years))
for (i in 1:length(years)){
  ages_in_2000 <- (2000 - years[i]) + (45:54)
  ok <- pop[,"AGE"] %in% ages_in_2000 & pop[,"MONTH"]==4 & pop[,"YEAR"]==2000
  pop_male <- pop[ok,"NHWA_MALE"]
  pop_female <- pop[ok,"NHWA_FEMALE"]
  print(c(weighted.mean(45:54, pop_male), weighted.mean(45:54, pop_female)))
  deathpr_1[i] <- weighted.mean(deathpr_by_age, pop_male + pop_female)
  deathpr_2[i] <- sum(deathpr_male* pop_male + deathpr_female*pop_female)/sum(pop_male + pop_female)
}

pdf("births.pdf", height=5.5, width=7)
par(mar=c(4,4,3,0), mgp=c(2.2,.5,0), tck=-.01)
plot(years, deathpr_2/deathpr_2[1], type="l", bty="l", xlab="Year", ylab="Death rate (compared to rate in 1989) ", main="Approx increase in death rate among 45-54-year-old whites,\n expected just from the changing age composition of this group", col="red", lwd=2, cex.axis=1.1, cex.lab=1.1)
grid(col="gray")
dev.off()

deaton <- read.table("deaton.txt", header=TRUE)

ages_all <- 35:64
ages_decade <- list(35:44, 45:54, 55:64)
years_1 <- 1999:2013
mort_data <- as.list(rep(NA,3))
group_names <- c("Non-Hispanic white", "Hispanic white", "African American")
mort_data[[1]] <- read.table("white_nonhisp_death_rates_from_1999_to_2013_by_sex.txt", header=TRUE)
mort_data[[2]] <- read.table("white_hisp_death_rates_from_1999_to_2013_by_sex.txt", header=TRUE)
mort_data[[3]] <- read.table("black_death_rates_from_1999_to_2013_by_sex.txt", header=TRUE)


raw_death_rate <- array(NA, c(length(years_1), 3, 3))
male_raw_death_rate <- array(NA, c(length(years_1), 3, 3))
female_raw_death_rate <- array(NA, c(length(years_1), 3, 3))
avg_death_rate <- array(NA, c(length(years_1), 3, 3))
male_avg_death_rate <- array(NA, c(length(years_1), 3, 3))
female_avg_death_rate <- array(NA, c(length(years_1), 3, 3))
for (k in 1:3){
  data <- mort_data[[k]]
  male <- data[,"Male"]==1
  for (j in 1:3){
    for (i in 1:length(years_1)){
      ok <- data[,"Year"]==years_1[i] & data[,"Age"] %in% ages_decade[[j]]
      raw_death_rate[i,j,k] <- 1e5*sum(data[ok,"Deaths"])/sum(data[ok,"Population"])
      male_raw_death_rate[i,j,k] <- 1e5*sum(data[ok&male,"Deaths"])/sum(data[ok&male,"Population"])
      female_raw_death_rate[i,j,k] <- 1e5*sum(data[ok&!male,"Deaths"])/sum(data[ok&!male,"Population"])
      avg_death_rate[i,j,k] <- mean(data[ok,"Rate"])
      male_avg_death_rate[i,j,k] <- mean(data[ok&male,"Rate"])
      female_avg_death_rate[i,j,k] <- mean(data[ok&!male,"Rate"])
    }
  }
}

for (k in 1:3){
  data <- mort_data[[k]]
  pdf(paste("death_rates_by_age_and_eth_", k, ".pdf", sep=""), height=11, width=8)
  par(mfrow=c(7,5), mar=c(2.5, 2.5, 2, .2), mgp=c(1.5,.3,0), tck=-.01, oma=c(0,0,3,0))
  years_1 <- 1999:2013
  for (i in 1:length(ages_all)){
    ok <- data[,"Age"]==ages_all[i]
    male <- data["Male"]==1
    male_deaths <- data[ok&male,"Deaths"]
    female_deaths <- data[ok&!male,"Deaths"]
    male_population <- data[ok&male,"Population"]
    female_population <- data[ok&!male,"Population"]
    male_mort <- male_deaths/male_population
    female_mort <- female_deaths/female_population
    total_mort <- (male_deaths + female_deaths)/(male_population + female_population)
    plot(years_1, total_mort/total_mort[1], xaxt="n", yaxt="n", ylim=range(.65,1.25), type="n", bty="n", xaxs="i", yaxs="i", xlab="", ylab=if (i%%5==1) "Relative death rate" else "", main=paste("age", ages_all[i]))
    lines(years_1, male_mort/male_mort[1], col="blue")
    lines(years_1, female_mort/female_mort[1], col="red")
    axis(1, seq(1990,2020,5))
    axis(2, seq(.6,1.2,.2))
    abline(h=1)
    grid(col="gray")
  }
  for (j in 1:3){
    plot(years_1, avg_death_rate[,j,k]/avg_death_rate[1,j,k], xaxt="n", yaxt="n", ylim=range(.65,1.25), type="n", bty="n", xaxs="i", yaxs="i", xlab="", ylab=if (j==1) "Relative death rate" else "", main=paste("Age-adj, ", min(ages_decade[[j]]), "-", max(ages_decade[[j]]), sep=""))
    lines(years_1, male_avg_death_rate[,j,k]/male_avg_death_rate[1,j,k], col="blue")
    lines(years_1, female_avg_death_rate[,j,k]/female_avg_death_rate[1,j,k], col="red")
    axis(1, seq(1990,2020,5))
    axis(2, seq(.6,1.2,.2))
    abline(h=1, col="gray")
  }
  mtext(paste(group_names[k], "women and men: trends in death rates since 1999"), side=3, outer=TRUE, line=1)
  par(mar=c(0,0,0,0))
  plot(c(-1,1), c(-1,1), xaxt="n", xlab="", yaxt="n", ylab="", bty="n", type="n")
  plot(c(-1,1), c(-1,1), xaxt="n", xlab="", yaxt="n", ylab="", bty="n", type="n")
  text(0, .5, paste("Red lines show\ntrends for women."), col="red")
  text(0, -.2, paste("Blue lines show\ntrends for men."), col="blue")
  dev.off()
}

pdf("effect_of_age_adj.pdf", height=6, width=7)
par(mfrow=c(3,3), mar=c(2.5, 2.5, 2, .2), mgp=c(1.5,.3,0), tck=-.01, oma=c(0,0,4,0))
text_pos <- array(NA, c(2,2,3,3))
text_pos[1,1,,] <- cbind(c(2008, 2003, 2010), c(2005, 2011, 2005), c(2005, 2007, 2005))
text_pos[2,1,,] <- cbind(c(2005, 2004, 2007), c(2004, 2008, 2006), c(2004, 2006, 2006))
text_pos[1,2,,] <- cbind(c(1.04, 1.06, .88), c(.91, .85, .88), c(.90, .90, .84))
text_pos[2,2,,] <- cbind(c(1.02, 1.03, .86), c(.86, .85, .93), c(.82, .80, .90))
for (k in 1:3){
  for (j in 1:3){
    rng <- range(avg_death_rate[,j,k]/avg_death_rate[1,j,k], raw_death_rate[,j,k]/raw_death_rate[1,j,k])
    plot(years_1, avg_death_rate[,j,k]/avg_death_rate[1,j,k], ylim=rng, xaxt="n", type="l", bty="l", xaxs="i", xlab="", ylab=if (j==1) "Death rate relative to 1999" else "", main=paste(group_names[k], " age ", min(ages_decade[[j]]), "-", max(ages_decade[[j]]), sep=""))
    lines(years_1, raw_death_rate[,j,k]/raw_death_rate[1,j,k], lty=2)
    abline(h=1, col="gray")
    axis(1, seq(1990,2020,5))
    text(text_pos[1,1,j,k], text_pos[1,2,j,k], "Raw", cex=.9)
    text(text_pos[2,1,j,k], text_pos[2,2,j,k], "Adjusted", cex=.9)
  }
}
mtext("Effects of age adjustment on trends in death rates by decade of age\n(Note:  these graphs are on different scales)", side=3, line=1, outer=TRUE)
dev.off()


pdf("decades.pdf", height=6, width=7)
par(mfrow=c(3,3), mar=c(2.5, 2.5, 2, .2), mgp=c(1.5,.3,0), tck=-.01, oma=c(0,0,4,0))
for (k in 1:3){
  for (j in 1:3){
    plot(years_1, avg_death_rate[,j,k]/avg_death_rate[1,j,k], xaxt="n", yaxt="n", ylim=range(.7, 1.1), type="n", bty="n", xaxs="i", yaxs="i", xlab="", ylab=if (j==1) "Relative death rate" else "", main=paste(group_names[k], ", ", min(ages_decade[[j]]), "-", max(ages_decade[[j]]), sep=""))
    lines(years_1, male_avg_death_rate[,j,k]/male_avg_death_rate[1,j,k], col="blue")
    lines(years_1, female_avg_death_rate[,j,k]/female_avg_death_rate[1,j,k], col="red")
    axis(1, seq(1990,2020,5))
    axis(2, seq(.7, 1.2, .1))
    abline(h=1, col="gray")
  }
}
mtext("Age-adjusted trends in death rate for 10-year bins", side=3, line=1, outer=TRUE)
dev.off()

pdf("focus_group.pdf", height=7, width=7)
par(mar=c(2.5, 3, 3, .2), mgp=c(2,.5,0), tck=-.01)
plot(years_1, avg_death_rate[,2,1],  ylim=c(382, 416), xaxt="n", yaxt="n", type="l", bty="l", xaxs="i", xlab="", ylab="Death rate per 100,000", main="AGE-ADJUSTED death rates for non-Hispanic whites aged 45-54")
axis(1, seq(1990,2020,5))
axis(2, seq(390, 420, 10))
grid(col="gray")
dev.off()


pdf("focus_group_2.pdf", height=6, width=7)
par(mar=c(2.5, 3, 3, .2), mgp=c(2,.5,0), tck=-.01)
plot(years_1, raw_death_rate[,2,1],  ylim=c(382,  416), xaxt="n", yaxt="n", type="l", bty="l", xaxs="i", xlab="", ylab="Death rate per 100,000", main="RAW death rates for non-Hispanic whites aged 45-54")
axis(1, seq(1990,2020,5))
axis(2, seq(390, 420, 10))
grid(col="gray")
dev.off()


pdf("focus_group_3.pdf", height=6, width=7)
par(mar=c(2.5, 3, 3, .2), mgp=c(2,.5,0), tck=-.01)
plot(range(years_1), c(1, 1.1), xaxt="n", yaxt="n", type="n", bty="l", xaxs="i", xlab="", ylab="Death rate relative to 1999", main="Age-adjusted death rates for non-Hispanic whites aged 45-54:\nTrends for women and men")
lines(years_1, male_avg_death_rate[,2,1]/male_avg_death_rate[1,2,1], col="blue")
lines(years_1, female_avg_death_rate[,2,1]/female_avg_death_rate[1,2,1], col="red")
axis(1, seq(1990,2020,5))
axis(2, seq(1, 1.1, .05))
text(2011.5, 1.075, "Women", col="red")
text(2010.5, 1.02, "Men", col="blue")
grid(col="gray")
dev.off()


pdf("focus_group_3.pdf", height=6, width=7)
par(mar=c(2.5, 3, 3, .2), mgp=c(2,.5,0), tck=-.01)
plot(range(years_1), c(1, 1.15), xaxt="n", yaxt="n", type="n", bty="l", xaxs="i", xlab="", ylab="Death rate relative to 1999", main="RAW death rates for non-Hispanic whites aged 45-54:\nTrends for women and men")
lines(years_1, male_raw_death_rate[,2,1]/male_raw_death_rate[1,2,1], col="blue")
lines(years_1, female_raw_death_rate[,2,1]/female_raw_death_rate[1,2,1], col="red")
axis(1, seq(1990,2020,5))
axis(2, seq(1, 1.2, .05))
text(2011.5, 1.11, "Women", col="red")
text(2010.5, 1.045, "Men", col="blue")
grid(col="gray")
dev.off()

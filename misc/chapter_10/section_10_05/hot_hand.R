expt <- function(n_players=26, n_shots=100, p=0.5) {
  p_after_hits <- rep(NA, n_players)
  p_after_misses <- rep(NA, n_players)
  for (j in 1:n_players) {
    shots <- rbinom(n_shots, 1, p)
    shots_before <- shots[1:(n_shots - 3)] + shots[2:(n_shots - 2)] + shots[3:(n_shots - 1)]
    shots_after <- shots[4:n_shots]
    p_after_hits[j] <- mean(shots_after[shots_before == 3])
    p_after_misses[j] <- mean(shots_after[shots_before == 0])
  }
  ok <- !is.na(p_after_hits + p_after_misses)
  c(mean(p_after_hits[ok]), mean(p_after_misses[ok]))
}

n_reps <- 1000
result <- array(NA, c(n_reps, 2))

for (k in 1:n_reps){
  result[k,] <- expt()
}
print(colMeans(result))
hist(result[,1] - result[,2])

for (k in 1:n_reps){
  result[k,] <- expt(p = 0.5 + 0.2*sin((1:100)*2*pi/100))
}
print(colMeans(result))
hist(result[,1] - result[,2])

for (k in 1:n_reps){
  result[k,] <- result[k,] <- expt(p = 0.5 + 0.1*sin((1:100)*2*pi/100))
}
print(colMeans(result))
hist(result[,1] - result[,2])

for (k in 1:n_reps){
  result[k,] <- result[k,] <- expt(p = 0.5 + 0.115*sin((1:100)*2*pi/100))
}
print(colMeans(result))
hist(result[,1] - result[,2])



expt2 <- function(n_players=26, n_shots=100, p=0.5) {
  p_after_hits <- rep(NA, n_players)
  p_after_misses <- rep(NA, n_players)
  for (j in 1:n_players) {
    shots <- rbinom(n_shots, 1, p)
    shots_before <- shots[1:(n_shots - 1)]
    shots_after <- shots[2:n_shots]
    p_after_hits[j] <- mean(shots_after[shots_before == 1])
    p_after_misses[j] <- mean(shots_after[shots_before == 0])
  }
  ok <- !is.na(p_after_hits + p_after_misses)
  c(mean(p_after_hits[ok]), mean(p_after_misses[ok]))
}

n_reps <- 1000
result <- array(NA, c(n_reps, 2))

for (k in 1:n_reps){
  result[k,] <- expt2()
}
print(colMeans(result))
hist(result[,1] - result[,2])

for (k in 1:n_reps){
  result[k,] <- result[k,] <- expt2(p = 0.5 + 0.115*sin((1:100)*2*pi/100))
}
print(colMeans(result))
hist(result[,1] - result[,2])


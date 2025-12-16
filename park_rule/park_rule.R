library("arm")
library("posterior")
library("rstanarm")
library("cmdstanr")

park <- read.csv("data/park.csv")
N <- nrow(park)
respondents <- sort(unique(park$submission_id))
J <- length(respondents)
respondent <- rep(NA, N)
for (j in 1:J){
  respondent[park$submission_id == respondents[j]] <- j
}
items <- sort(unique(park$question_id))
K <- length(items)
item <- park$question_id
y <- park$answer
n_responses <- rep(NA, J)
for (j in 1:J){
  n_responses[j] <- sum(respondent==j)
}
male_name <- park$male
white_name <- park$white
n_responses_full <- n_responses[respondent]


data_matrix <- data.frame(y, respondent, item, male_name, white_name, n_responses_full)
fit_lme4 <- glmer(y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_responses_full, family=binomial(link="logit"), data=data_matrix)
display(fit_lme4)

n_skipped <- K - n_responses

n_skipped_full <- n_skipped[respondent]
data_matrix <- data.frame(y, respondent, item, male_name, white_name, n_skipped_full)
fit_lme4 <- glmer(y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_skipped_full, family=binomial(link="logit"), data=data_matrix)
display(fit_lme4)

mean(y)
mean(y[male_name==0 & white_name==0 & n_skipped_full==0])
sum(male_name==0 & white_name==0 & n_skipped_full==0)
invlogit(-2.42 + 0.07*mean(male_name) + 0.08*mean(white_name) + 0.03*mean(n_skipped_full))
mean(predict(fit_lme4, type="response"))

set.seed(123)
a_respondent_sim <- rnorm(J, 0, sqrt(VarCorr(fit_lme4)$respondent))
a_item_sim <- rnorm(K, 0, sqrt(VarCorr(fit_lme4)$item))
b_sim <- fixef(fit_lme4)

X <- cbind(1, male_name, white_name, n_skipped_full)
p_sim <- invlogit(a_respondent_sim[respondent] + a_item_sim[item] + X %*% b_sim)
y_sim <- rbinom(N, 1, p_sim)

data_sim <- data.frame(data_matrix, y_sim)
fit_lme4_sim <- glmer(y_sim ~ (1 | item) + (1 | respondent) + male_name + white_name + n_responses_full, family=binomial(link="logit"), data=data_sim)
display(fit_lme4_sim)

a_item_hat <- ranef(fit_lme4)$item
print(a_item_hat)

wordings <- read.csv("data/park.txt", header=FALSE)$V2
wordings <- substr(wordings, 2, nchar(wordings)-1)
a_item_hat <- unlist(a_item_hat)
names(a_item_hat) <- wordings

sort(a_item_hat)

item_avg <- rep(NA, K)
for (k in 1:K) {
  item_avg[k] <- mean(y[item==k])
}
plot(item_avg, a_item_hat, pch=20)

plot(logit(item_avg), a_item_hat, type="n")
text(logit(item_avg), a_item_hat, names(a_item_hat), cex=.5)

a_respondent_hat <- unlist(ranef(fit_lme4)$respondent)
respondent_avg <- rep(NA, J)
for (j in 1:J) {
  respondent_avg[j] <- mean(y[respondent==j])
}
plot(respondent_avg, a_respondent_hat, pch=20, cex=.4)

order_of_response <- rep(NA, N)
for (j in 1:J) {
  order_of_response[respondent==j] <- 1:n_responses[j]
}
data_matrix <- data.frame(data_matrix, order_of_response)

fit_lme4 <- glmer(y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_skipped_full + order_of_response, family=binomial(link="logit"), data=data_matrix)
display(fit_lme4)

# rstanarm
options(mc.cores = 4)

fit_rstanarm <- stan_glmer(y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_skipped_full + order_of_response, family=binomial(link="logit"), data=data_matrix)
print(fit_rstanarm)

fit_rstanarm <- stan_glmer(y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_skipped_full + order_of_response, family=binomial(link="logit"), data=data_matrix, iter=200, control=list(max_treedepth=5))
print(fit_rstanarm)

# Stan

X <- cbind(rep(1,N), male_name, white_name, n_skipped_full, order_of_response)
stan_data <- list(N=N, J=J, K=K, L=ncol(X), y=y, respondent=respondent, item=item, X=X)

park_1 <- cmdstan_model("park_1.stan")
fit_1 <- park_1$sample(data=stan_data, parallel_chains=4, iter_warmup=100, iter_sampling=100, max_treedepth=5)
print(fit_1)

fit_1 <- park_1$sample(data=stan_data, parallel_chains=4, max_treedepth=5)
print(fit_1)

park_2 <- cmdstan_model("park_2.stan")
fit_2 <- park_2$sample(data=stan_data, parallel_chains=4, max_treedepth=5)
print(fit_2)

fit_2 <- park_2$sample(data=stan_data, parallel_chains=4)
print(fit_2)

park_3 <- cmdstan_model("park_3.stan")
fit_3 <- park_3$sample(data=stan_data, parallel_chains=4, max_treedepth=5)
print(fit_3)

# Simulated-data experiment
stan_data_sim <- stan_data
stan_data_sim$y <- y_sim
fit_1_sim <- park_1$sample(data=stan_data_sim, parallel_chains=4, max_treedepth=5)
print(fit_1_sim)
fit_2_sim <- park_2$sample(data=stan_data_sim, parallel_chains=4, max_treedepth=5)
print(fit_2_sim)
fit_3_sim <- park_3$sample(data=stan_data_sim, parallel_chains=4, max_treedepth=5)
print(fit_3_sim)

# Fit normal model
park_normal <- list(cmdstan_model("park_1_normal.stan"),
                    cmdstan_model("park_2_normal.stan"),
                    cmdstan_model("park_3_normal.stan"))
fit_normal <- as.list(rep(NA, 3))
for (i in 1:3){
  fit_normal[[i]] <- (park_normal[[i]])$sample(data=stan_data, parallel_chains=4, max_treedepth=5)
  print(fit_normal[[i]])
}

# Simulate from fitted normal model
sims <- as_draws_rvars(fit_normal[[3]]$draws())
a_respondent_sim <- rnorm(J, 0, median(sims$sigma_respondent))
a_item_sim <- rnorm(K, 0, median(sims$sigma_item))
b_sim <- median(sims$b)
sigma_y_sim <- median(sims$sigma_y)
y_sim <- rnorm(N, a_respondent_sim[respondent] + a_item_sim[item] + X %*% b_sim, sigma_y_sim)

stan_data_sim$y <- y_sim
fit_normal_sim <- as.list(rep(NA, 3))
for (i in 1:3){
  fit_normal_sim[[i]] <- (park_normal[[i]])$sample(data=stan_data_sim, parallel_chains=4, max_treedepth=5)
  print(fit_normal_sim[[i]])
}

# ADVI and Pathfinder
advi_3 <- park_3$variational(data=stan_data)
pathfinder_3 <- park_3$pathfinder(data=stan_data)
print(advi_3)
print(pathfinder_3)
print(fit_3)

# Other ADVI's
advi_1 <- park_1$variational(data=stan_data)
advi_2 <- park_2$variational(data=stan_data)
print(advi_1)
print(advi_2)

# ADVI as starting point
advi_1 <- park_1$variational(data=stan_data)
fit_after_advi_1 <- park_1$sample(data=stan_data, init=advi_1, parallel_chains=4, max_treedepth=5)
print(fit_after_advi_1)

advi_2 <- park_2$variational(data=stan_data)
fit_after_advi_2 <- park_2$sample(data=stan_data, init=advi_2, parallel_chains=4, max_treedepth=5)
print(fit_after_advi_2)


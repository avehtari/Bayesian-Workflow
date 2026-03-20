#' ---
#' title: "Sharks and model building with latent variables"
#' author: "Vianey Leos Barajas"
#' date: 2025-12-14
#' date-modified: today
#' date-format: iso
#' format:
#'   html:
#'     number-sections: true
#'     code-copy: true
#'     code-download: true
#'     code-tools: true
#' bibliography: ../casestudies.bib
#' ---

#' This notebook includes the code for Bayesian Workflow book Chapter
#' 26 *Model building with latent variables: Markov models for animal
#' movement*.
#'
#' # Introduction
#'
#' To demonstrate how hidden Markov models (HMMs) are used to model
#' animal movement, we use positional data analyzed by
#' @Towner-Leos-Barajas-Langrock-etal:2016. Here we replicate the
#' model building workflow and demonstrate how this would be done in
#' a Bayesian framework.
#'
#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=TRUE, comment=NA, out.width='95%')

#' **Load packages**
#| cache: FALSE
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(tidyverse)
library(CircStats)
library(patchwork)
library(lubridate)
library(moveHMM)
library(bayesplot)
library(tidyr)
library(cmdstanr)
options(mc.cores = 4)
# CmdStanR output directory makes Quarto cache to work
dir.create(root("sharks", "stan_output"), showWarnings = FALSE)
options(cmdstanr_output_dir = root("sharks", "stan_output"))

#' # Shark movement data
#' The positions, in longitude and latitude, of multiple sharks were taken
#' over time in Gansbaii, South Africa. Some sharks have repeated trackings.
#' Hidden Markov models (HMMs) applied to the positions of animals over time are
#' often first transformed into step lengths and turning angles, and then
#' analyzed. Here we use the R package `moveHMM` to process our tracks.
#' `moveHMM` requires a column named ID to denote individual tracks, with
#' multiple IDs possibly corresponding to the same shark if they are tracked
#' multiple times. Because these are discrete-time HMMs, `moveHMM` also assumes
#' that the observations are taken regularly over time, with missing values
#' filling any times where observations were not available. It will not check
#' this automatically, so it is important for users to prepare their data
#' appropriately. In the following code we import the data, use the function
#' `prepData` to compute step lengths and turning angles. We also compute
#' covariates, such as time of day, sex, and chum that will be used for
#' model extensions.
load(root("sharks/data","whiteshark_trackdata.RData"))
sharks.HMMtracks.df$ID <- sharks.HMMtracks.df$SharksexTrackNo
moveHMM_wsdata <- prepData(sharks.HMMtracks.df[, c("ID", "Long", "Lat")], type = c("LL"), coordNames = c("Long", "Lat"))
sharks.HMMtracks.df$steplength <- moveHMM_wsdata$step
sharks.HMMtracks.df$turnang <- moveHMM_wsdata$angle
ws_HMM_full <- sharks.HMMtracks.df |>
  filter(!SharksexTrackNo %in% c("WSF9 T4", "WSF9 T3 B "))
ws_HMM <- ws_HMM_full[, c("dateTime",
                          "SharkName",
                          "SharksexTrackNo",
                          "steplength",
                          "turnang",
                          "year",
                          "month",
                          "CDB")]

#' For implementation in Stan, we must set the NA's to numeric values, ideally
#' something non-plausible that can be checked in a for loop during likelihood
#' evaluation.
ws_HMM$steplength[is.na(ws_HMM$steplength)] <- -100
ws_HMM$steplength[which(ws_HMM$steplength > 1.5)] <- -100
ws_HMM$turnang[is.na(ws_HMM$turnang)] <- -100

#' We plot two tracks of white shark female 1 that exhibit the distinct
#' movement patterns we are hoping to capture with our HMM:
#| label: fig-sarika_tracks_1_10
#| fig-height: 5
#| fig-width: 4
#| warning: false
ggplot(data = sharks.HMMtracks.df |>
         dplyr::filter(SharksexTrackNo ==
                         c("WSF1 T1", "WSF1 T10")),
       aes(Long, Lat)) +
  geom_path(aes(group = SharksexTrackNo), alpha = 0.5, color = "grey") +
  geom_point(alpha = 0.5) +
  facet_wrap(~SharksexTrackNo) + theme_classic() +
  ylab("Latitude") + xlab("Longitude") + coord_map()

#' # Initial values for an N-state HMM:
#'
#' We create two functions to provide starting values for the mean of the
#' state-dependent distribution for step length. Without specifying starting
#' values, the MCMC may have some trouble initializing at plausible values.
#' There may still be issues with initialization depending on starting
#' values for other parameters but for now we proceed only specifying
#' starting values for the means.
init_fun_mu <- function(no_states, no_chains) {
  mu_init <- list()
  for (n in 1:no_chains) {
    mu_init[[n]] <- list(mu = sort(runif(no_states, min = 0, max = 0.5)))
  }
  mu_init
}
init_fun_logmu <- function(no_states, no_chains) {
  log_mu_init <- list()
  for (n in 1:no_chains) {
    log_mu_init[[n]] <- list(log_mu = log(sort(runif(no_states, min = 0, max = 0.5))))
  }
  log_mu_init
}

#' # 2-state HMM
#' We fit a 2-state HMM to the white shark tracks, assuming that each 
#' track is an independent realization of the same 2-state HMM. Note that in 
#' the corresponding chapter, we modify either the prior distributions for the 
#' state-dependent distribution or specify an ordering on the means for the 
#' step length state-dependent distribution. The Stan code in 
#' `step_turn_hmm.stan` can be modified to mimic different behavior and 
#' reproduce similar results. 
stanHMM_2states <-  list(
  Nstates = 2, 
  Tlen = dim(ws_HMM)[1], 
  track_index = as.numeric(as.factor(ws_HMM$SharksexTrackNo)),
  steplength = ws_HMM$steplength, 
  angle = ws_HMM$turnang,
  Ntracks = length(unique(ws_HMM$SharksexTrackNo))
)

#| label: fit_2stateHMM
#| cache: true
#| results: hide
model_2stateHMM <- cmdstan_model(
  root("sharks", "step_turn_hmm.stan")
)
fit_2stateHMM <- model_2stateHMM$sample(
  data = stanHMM_2states,
  init = init_fun_mu(2,4),
  chains = 4
)

fit_2stateHMM$summary(variables = c("mu", "sigma", "mixp",
                                    "xangle", "yangle",
                                    "tpm", "initial_dist",
                                    "lp__"))


#' Plot MCMC marginal distributions for certain parameters:
#| label: fig-mcmc_hist_by_chain_2state
fit_2stateHMM_draws <- fit_2stateHMM$draws(format = "df",
                                           variables = c("mu", "sigma", "mixp",
                                                         "shape", "rate",
                                                         "xangle", "yangle",
                                                         "kappa", "loc",
                                                         "tpm", "initial_dist",
                                                         "lp__"))
mcmc_hist_by_chain(fit_2stateHMM_draws, regex_pars = "mu", pars = "lp__")

#' Plot state-dependent distributions for step length and turning angle:
#| results: hide
#| warning: false
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
eval_gamma_sdd <- function(post_draws, no_samples) {
  # setup
  xval <- seq(0.001, 1, len = 200)
  state1_dens <- matrix(NA, nrow = 200, ncol = no_samples)
  state2_dens <- matrix(NA, nrow = 200, ncol = no_samples)

  for (j in 1:no_samples) {
    state1_dens[,j] <- dgamma(xval,
                              shape = as.numeric(post_draws[j, "shape[1]"]), 
                              rate = as.numeric(post_draws[j, "rate[1]"]))
    state2_dens[,j] <- dgamma(xval,
                              shape = as.numeric(post_draws[j, "shape[2]"]), 
                              rate = as.numeric(post_draws[j, "rate[2]"]))
  }
  sdd <- data.frame(xval = rep(xval, times = no_samples),
                    state1_dens = c(state1_dens),
                    state2_dens = c(state2_dens),
                    sample = rep(1:no_samples, each= 200))
  sdd
}
eval_vonMises_sdd <- function(post_draws, no_samples) {
  # setup
  xval <- seq(-pi, pi, len = 200)
  state1_dens <- matrix(NA, nrow = 200, ncol = no_samples)
  state2_dens <- matrix(NA, nrow = 200, ncol = no_samples)

  for (j in 1:no_samples) {
    state1_dens[,j] <- dvm(xval,
                           mu = as.numeric(post_draws[j, "loc[1]"]), 
                           kappa = as.numeric(post_draws[j, "kappa[1]"]))
    state2_dens[,j] <- dvm(xval,
                           mu =  as.numeric(post_draws[j, "loc[2]"]), 
                           kappa = as.numeric(post_draws[j, "kappa[2]"]))
  }
  sdd <- data.frame(xval = rep(xval, times = no_samples),
                    state1_dens = c(state1_dens),
                    state2_dens = c(state2_dens),
                    sample = rep(1:no_samples, each= 200))
  sdd
}
sdd_steplength <- eval_gamma_sdd(fit_2stateHMM_draws, no_samples = 1000)
sdd_angle <- eval_vonMises_sdd(fit_2stateHMM_draws, no_samples = 1000)
#| label: fig-hmm_2state_distributions
#| warning: false
sl_sdd <- ggplot(ws_HMM) +
  geom_histogram(aes(steplength, y = after_stat(density)), bins = 100, alpha = 0.2) +
  xlim(-0.01, 1) +
  geom_line(data = sdd_steplength, aes(xval, 0.5 * state1_dens, group = sample), 
            color = cbPalette[2], alpha = 0.2) +
  geom_line(data = sdd_steplength, aes(xval, 0.5 * state2_dens, group = sample), 
            color = cbPalette[3], alpha = 0.2) +
  geom_line(data = sdd_steplength, aes(xval, 0.5 * state1_dens + 0.5 * state2_dens, group = sample), 
            color = "darkgrey", alpha = 0.02) +
  theme_classic() + xlab("") +
  ylab("density") +
  annotate("text", x = 0.3, y = 4, label = "state 1", color = cbPalette[2])  +
  annotate("text", x = 0.45, y = 1.5, label = "state 2", color = cbPalette[3]) +
  ggtitle("step length state-dependent distributions")
angle_sdd <- ggplot(ws_HMM) +
  geom_histogram(aes(turnang, y = after_stat(density)), bins = 100, alpha = 0.2) +
  xlim(-pi, pi) +
  geom_line(data=sdd_angle, aes(xval, 0.5 * state1_dens, group = sample), 
            color = cbPalette[2], alpha = 0.2) +
  geom_line(data = sdd_angle, aes(xval, 0.5 * state2_dens, group = sample), 
            color = cbPalette[3], alpha = 0.2) +
  geom_line(data = sdd_angle, aes(xval, 0.5 * state1_dens + 0.5 * state2_dens, group = sample), 
            color = "darkgrey", alpha = 0.02) +
  theme_classic() + xlab("") +
  ylab("density") +
  annotate("text", x = -2.5, y = 0.25, label = "state 1", color = cbPalette[2])  +
  annotate("text", x = 1.5, y = 0.25, label = "state 2", color = cbPalette[3]) +
  ggtitle("turning angle state-dependent distributions")
sl_sdd + angle_sdd

#' Plot state-decodings onto track using the forward-backward algorithm 
#' for local state-decoding:
state_probs_draws <- fit_2stateHMM$draws(variables =c("state_probs"), format = "draws_matrix")
state_probs_means <- data.frame(state1prob = colMeans(state_probs_draws[1:1000, (4584 + 1):(2 * 4584)]), 
                                state2prob = colMeans(state_probs_draws[1:1000, 1:4584]))
state1_probs_quants <- data.frame(
  state1prob025 = apply(state_probs_draws[1:1000, (4584 + 1):(2 * 4584)], 2, quantile, probs=0.025),
  state1prob975 = apply(state_probs_draws[1:1000, (4584 + 1):(2 * 4584)], 2, quantile, probs=0.975)
)
state2_probs_quants <- data.frame(
  state1prob025 = apply(state_probs_draws[1:1000, 1:(4584)], 2, quantile, probs=0.025),
  state1prob975 = apply(state_probs_draws[1:1000, 1:(4584)], 2, quantile, probs=0.975)
)
ws_HMM_rep <- ws_HMM_full[,c("dateTime",
                             "SharkName",
                             "SharksexTrackNo",
                             "steplength",
                             "turnang",
                             "year",
                             "month",
                             "CDB",
                             "Lat", "Long")]

ws_HMM_rep$state1prob <- state_probs_means[,1]
ws_HMM_rep$state2prob <- state_probs_means[,2]
ws_HMM_rep$state1prob025 <- state1_probs_quants[,1]
ws_HMM_rep$state1prob975 <- state1_probs_quants[,2]
ws_HMM_rep$state2prob025 <- state2_probs_quants[,1]
ws_HMM_rep$state2prob975 <- state2_probs_quants[,2]
#| label: fig-hmm_2state_decodings
#| warning: false
ggplot(data = ws_HMM_rep |>
         dplyr::filter(SharksexTrackNo ==
                         c("WSF1 T1", "WSF1 T10")),
       aes(Long, Lat)) +
  geom_path(aes(group = SharksexTrackNo), alpha = .5, color = "grey") +
  geom_point(aes(color = state1prob)) +
  labs(color = "State 1 \nProbability") +
  scale_color_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~SharksexTrackNo) + theme_classic() +
  ylab("Latitude") + xlab("Longitude") + coord_map()

#| label: fig-hmm_2state_decodings2
#| warning: false
ggplot(data = ws_HMM_rep |>
         dplyr::filter(SharksexTrackNo ==
                         c("WSF1 T1", "WSF1 T10")),
       aes(dateTime, state1prob)) + theme_classic() +
  geom_point() + geom_ribbon(aes(ymin = state1prob025,
                                 ymax = state1prob975), alpha = 0.5) +
  facet_wrap( ~ SharksexTrackNo, nrow = 2, scales = "free_x") +
  ylab("State 1 Probability") + xlab("Time")

#' Plot pseudo-residuals
# plot pseudo residuals
pseudo_residuals <- fit_2stateHMM$draws(variables = "pseudo_residuals",
                                        format = "draws_matrix")
pr_mean <- colMeans(pseudo_residuals[1:1000, ])
pr_025 <- apply(pseudo_residuals[1:1000, ], 2, quantile, probs = 0.025)
pr_975 <- apply(pseudo_residuals[1:1000, ], 2, quantile, probs = 0.975)
pr_missing_index <- which(pr_mean == 0)
df <- data.frame(y = pr_mean[-pr_missing_index])

#| label: fig-hmm_2state_pseudo_residuals
p <- ggplot(df, aes(sample = y))
p + stat_qq() + stat_qq_line() + theme_classic() +
  xlab("Theoretical Quantiles") +
  ylab("Sample Quantiles") +
  ggtitle("Quantile-Quantile Plot")

#' Plot state decoding histograms
#'
#' Simulate data from posterior predictive distribution
#' of fitted 2-state HMM
#'
init_dist <- fit_2stateHMM$draws("initial_dist", format = "draws_array")
tpm <- fit_2stateHMM$draws("tpm", format = "draws_array")
sdd_sl_shape <- fit_2stateHMM$draws("shape", format = "draws_array")
sdd_sl_rate <- fit_2stateHMM$draws("rate", format = "draws_array")
sdd_sl_zeromass <- fit_2stateHMM$draws("mixp", format = "draws_array")
sdd_angle_loc <- fit_2stateHMM$draws("loc", format = "draws_array")
sdd_angle_conc <- fit_2stateHMM$draws("kappa", format = "draws_array")
post_state_samples <- fit_2stateHMM$draws("state_sequence", format = "draws_array")

state_samples <- matrix(NA, nrow = 4584, ncol = 1000)
state_samples[1, 1] <- sample(x = 2:1, size = 1, prob = init_dist[1,1,])
for (j in 1:1000) {
  tpm_iter <- matrix(data = tpm[j, 1, 4:1], nrow = 2)
  state_samples[1,j] <- sample(x = 2:1, size = 1, prob = init_dist[j,1,])
  for (t in 2:4584) {
    if(ws_HMM$SharksexTrackNo[t]==ws_HMM$SharksexTrackNo[t-1]) {
      state_samples[t, j] <- sample(x = 2:1,
                                    size = 1,
                                    prob = tpm_iter[state_samples[t - 1, j], ])
    } else {
      state_samples[t,j] <- sample(x = 2:1, size = 1, prob = init_dist[j, 1, ])
    }
  }
}
sim_state_counts <- apply(state_samples[1:91, ], 2, table)
sim_state_props <- sim_state_counts/91
post_state_samples_counts <- apply(post_state_samples[, 1, 1:91], 1, table)
post_state_samples_props <- post_state_samples_counts/91

ws_HMM_rep$chum <- ifelse(ws_HMM_full$CDB == "x", yes = 1, no = 0) 
ws_HMM_rep$chum[which(is.na(ws_HMM_rep$CDB))] <- 0
chum_ws1_tr1 <- which(ws_HMM_rep$chum[1:91] == 1)
post_state_samples_chumcounts <- apply(post_state_samples[, 1, chum_ws1_tr1], 1, table)
post_state_samples_st2chumprops <- numeric(42)
for (j in 1:1000) {
  post_state_samples_st2chumprops[j] <-
    max(post_state_samples_chumcounts[[j]])/42
}
#| label: fig-hmm_2state_decoding_histograms
ggplot(data = data.frame(y = sim_state_props[1, ]), aes(y)) +
  geom_histogram(bins = 100, fill = "darkgrey") + theme_classic() + xlim(-0.05, 1) +
  annotate("text", label = "Posterior\nPredictive\nSimulations", x = 0.63, y = 200, col = "darkgrey") +
  geom_histogram(data = data.frame(x = post_state_samples_props[1,]),
                 aes(x), bins = 100, fill = "black") +
  annotate("text", label="State\nDecodings", x=0.35, y=200, col="black") +
  geom_histogram(data = data.frame(p = 1 - post_state_samples_st2chumprops), 
                 aes(p), bins = 100, fill = "blue") +
  annotate("text", label = "Chum", x = 0.07, y = 200, col = "blue") +
  xlab("Proportion of State 1 Observations") +
  ylab("Count")


#' # 2-state HMM with covariates
#'
#' Time inhomogeneous HMMs with covariates in transition probability matrix
#' time of day covariates
ws_HMM$tod_cos <- cos((2*pi*(hour(ws_HMM$dateTime)*60 + minute(ws_HMM$dateTime)))/1440)
ws_HMM$tod_sin <- sin((2*pi*(hour(ws_HMM$dateTime)*60 + minute(ws_HMM$dateTime)))/1440)
#' chum covariate
ws_HMM$chum <- ifelse(ws_HMM$CDB == "x", yes = 1, no = 0)
ws_HMM$chum[which(is.na(ws_HMM$CDB))] <- 0
#' sex covariate
ws_HMM$sex_char <- substring(ws_HMM$SharksexTrackNo, 3, 3)
ws_HMM$sex <- ifelse(ws_HMM$sex_char == "F", yes = 0, no = 1)
HMM_covar <- matrix(data = c(rep(1, dim(ws_HMM)[1]),
                             ws_HMM$chum,
                             ws_HMM$sex,
                             ws_HMM$tod_cos,
                             ws_HMM$tod_sin), nrow = dim(ws_HMM)[1], ncol = 5)
stanHMM_2states_covariates <- list(Nstates = 2,
                                   Tlen = dim(ws_HMM)[1],
                                   track_index = as.numeric(as.factor(ws_HMM$SharksexTrackNo)),
                                   steplength = ws_HMM$steplength,
                                   angle = ws_HMM$turnang,
                                   nCovs = 4,
                                   covs = HMM_covar)

#| label: fit_2stateHMM_covariates
#| cache: true
#| results: hide
model_2stateHMM_covariates <- cmdstan_model(
  root("sharks", "step_turn_hmm_covariates.stan")
)
fit_2stateHMM_covariates <- model_2stateHMM_covariates$sample(
  data = stanHMM_2states_covariates,
  init = init_fun_mu(2, 4), 
  chains = 4
)

#' # 2-state HMM covariates in transition probability matrix and individual varying effects

#' HMM with 2 state, covariates in tpm and non-centered parametrization for varying effects
stanHMM_2states_tpmcov_crencp <- list(
  Nstates = 2,
  Tlen = dim(ws_HMM)[1],
  track_index = as.numeric(as.factor(ws_HMM$SharksexTrackNo)),
  steplength = ws_HMM$steplength,
  angle = ws_HMM$turnang,
  nCovs = 4,
  covs = HMM_covar[, -1], 
  Nsharks = length(unique(ws_HMM$SharksexTrackNo)),
  shark_index = as.numeric(as.factor(ws_HMM$SharksexTrackNo))
)
#| label: fit_2stateHMM_tpmcov_crencp
#| cache: true
#| results: hide
model_2stateHMM_tpmcov_crencp <- cmdstan_model(
  root("sharks", "step_turn_hmm_covariates_cre_ncp.stan")
)
fit_2stateHMM_tpmcov_crencp <- model_2stateHMM_tpmcov_crencp$sample(
  data = stanHMM_2states_tpmcov_crencp,
  init = init_fun_mu(2, 1),
  chains = 1
)

#' Plot entries of transition probability matrix with covariates and individual varying effects:
beta <- fit_2stateHMM_tpmcov_crencp$draws("beta", format = "draws_array")
randeff <- fit_2stateHMM_tpmcov_crencp$draws("randeff_tpm", format = "draws_array")
mu_tpm <- fit_2stateHMM_tpmcov_crencp$draws("mu_tpm", format = "draws_array")

# intercept is female, no chum
DM_female_nochum <- cbind(
  1,
  rep(0, 721),
  rep(0, 721),
  cos(2*pi*480:1200/1440),
  sin(2*pi*480:1200/1440)
)
# intercept is female, chum
DM_female_chum <- cbind(
  1,
  rep(1, 721),
  rep(0, 721),
  cos(2*pi*480:1200/1440),
  sin(2*pi*480:1200/1440)
)
DM_male_nochum <- cbind(
  1,
  rep(0, 721),
  rep(1, 721),
  cos(2*pi*480:1200/1440),
  sin(2*pi*480:1200/1440)
)
beta_array <- array(data = NA, dim = c(2, 2, length(480:1200), 100))
grid_tod <- 480:1200
track_factor <- as.factor(ws_HMM$SharksexTrackNo)
track_numeric <- as.numeric(track_factor)
#| label: fig-hmm_2state_tpmcov_crencp_transitions
par(mfrow = c(2, 2))
for (i in 1:2) {
  for (j in 1:2) {
    beta_mat <- rbind(c(mu_tpm[1,1,1],
                        beta[1,1,c(1, 3, 5, 7)]),
                      c(mu_tpm[1,1,2],
                        beta[1,1,c(2, 4, 6, 8)]))

    tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                                   beta = t(beta_mat),
                                   covs = DM_male_nochum)
    plot(grid_tod, tpm[i, j, ], type = "l",
         ylim = c(0, 1), col = "grey",
         lwd = 0.5, xlab = "minute of the day",
         ylab = paste0("Pr(", i, " -> ",j, ")"))
    for (k in 2:100) {

      beta_mat <- rbind(c(mu_tpm[k,1,1],
                          beta[k,1,c(1, 3, 5, 7)]),
                        c(mu_tpm[k,1,2],
                          beta[k,1,c(2, 4, 6, 8)]))
      tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                                     beta = t(beta_mat),
                                     covs = DM_male_nochum)
      lines(grid_tod, tpm[i, j, ], col = "grey", lwd = .1)
    }
  }
}
beta_mat <- rbind(c(mu_tpm[1, 1, 1], 
                    beta[1, 1, c(1, 3, 5, 7)]), 
                  c(mu_tpm[1, 1, 2], 
                    beta[1, 1, c(2, 4, 6, 8)]))
tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                               beta = t(beta_mat),
                               covs = DM_male_nochum)
male_chum_tpm <- data.frame(omega11 = tpm[1, 1, ], 
                            omega12 = tpm[1, 2, ], 
                            omega21 = tpm[2, 1, ], 
                            omega22 = tpm[2, 2, ], 
                            dateTime = 480:1200, 
                            draw = 1)
for (k in 2:500) {
  beta_mat <- rbind(c(mu_tpm[k, 1, 1], 
                      beta[k, 1, c(1, 3, 5, 7)]), 
                    c(mu_tpm[k, 1, 2], 
                      beta[k, 1, c(2, 4, 6, 8)]))
  tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                                 beta = t(beta_mat),
                                 covs = DM_male_nochum)
  male_chum_tpm <- rbind(male_chum_tpm,
                         data.frame(omega11 = tpm[1, 1, ], 
                                    omega12 = tpm[1, 2, ], 
                                    omega21 = tpm[2, 1, ], 
                                    omega22 = tpm[2, 2, ], 
                                    dateTime = 480:1200,
                                    draw = k))
}
colnames(male_chum_tpm)[1] <- "omega[11](t)"
colnames(male_chum_tpm)[2] <- "omega[12](t)"
colnames(male_chum_tpm)[3] <- "omega[21](t)"
colnames(male_chum_tpm)[4] <- "omega[22](t)"
male_chum_tpm_long <- male_chum_tpm |>
  pivot_longer(!c(dateTime, draw), names_to = "omega")
p1 <- ggplot(male_chum_tpm_long, aes(dateTime, value)) +
  geom_line(aes(group = draw), alpha = 0.1, col = "darkgrey") +
  facet_wrap(~omega, nrow = 2, ncol = 2, labeller = label_parsed) +
  theme_classic() + ylab("Probability") +
  xlab("Time") +
  scale_x_continuous(breaks = c(540, 720, 900, 1080),
                     labels = c("09:00", "12:00", "15:00", "18:00")) +
  ggtitle("Male white shark, no chum")

par(mfrow = c(2, 2))
for (i in 1:2) {
  for (j in 1:2) {

    beta_mat <- rbind(c(mu_tpm[1, 1, 1], 
                        beta[1, 1, c(1, 3, 5, 7)]), 
                      c(mu_tpm[1, 1, 1], 
                        beta[1, 1, c(2, 4, 6, 8)]))

    tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                                   beta = t(beta_mat),
                                   covs = DM_female_nochum)

    plot(grid_tod, tpm[i, j, ], type = "l",
         ylim = c(0, 1), col = "grey",
         lwd = 0.5, xlab = "minute of the day",
         ylab = paste0("Pr(", i, " -> ",j, ")"))

    for (k in 2:100) {

      beta_mat <- rbind(c(mu_tpm[k, 1, 1], 
                          beta[k, 1, c(1, 3, 5, 7)]), 
                        c(mu_tpm[k, 1, 2], 
                          beta[k, 1, c(2, 4, 6, 8)]))

      tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                                     beta = t(beta_mat),
                                     covs = DM_female_nochum)
      lines(grid_tod, tpm[i, j, ], col = "grey", lwd = .1)
    }
  }
}

beta_mat <- rbind(c(mu_tpm[1, 1, 1], 
                    beta[1, 1, c(1, 3, 5, 7)]), 
                  c(mu_tpm[1, 1, 2], 
                    beta[1, 1, c(2, 4, 6, 8)]))
tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                               beta = t(beta_mat),
                               covs = DM_male_nochum)
female_chum_tpm <- data.frame(
  omega11 = tpm[1, 1, ], 
  omega12 = tpm[1, 2, ], 
  omega21 = tpm[2, 1, ], 
  omega22 = tpm[2, 2, ], 
  dateTime = 480:1200,
  draw = 1
)
for (k in 2:500) {
  beta_mat <- rbind(c(mu_tpm[1, 1, 1], 
                      beta[k, 1, c(1, 3, 5, 7)]), 
                    c(mu_tpm[1, 1, 2], 
                      beta[k, 1, c(2, 4, 6, 8)]))
  tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                                 beta = t(beta_mat),
                                 covs = DM_female_nochum)
  female_chum_tpm <- rbind(
    female_chum_tpm,
    data.frame(omega11 = tpm[1, 1, ], 
               omega12 = tpm[1, 2, ], 
               omega21 = tpm[2, 1, ], 
               omega22 = tpm[2, 2, ], 
               dateTime = 480:1200,
               draw = k)
  )
}

colnames(female_chum_tpm)[1] <- "omega[11](t)"
colnames(female_chum_tpm)[2] <- "omega[12](t)"
colnames(female_chum_tpm)[3] <- "omega[21](t)"
colnames(female_chum_tpm)[4] <- "omega[22](t)"
female_chum_tpm_long <- female_chum_tpm |>
  pivot_longer(!c(dateTime, draw), names_to = "omega")
p2 <- ggplot(female_chum_tpm_long, aes(dateTime, value)) +
  geom_line(aes(group = draw), alpha = 0.1, col = "darkgrey") +
  facet_wrap(~omega, nrow = 2, ncol = 2, labeller = label_parsed) +
  theme_classic() + ylab("Probability") +
  xlab("Time") +
  scale_x_continuous(breaks = c(540, 720, 900, 1080),
                     labels = c("09:00", "12:00", "15:00", "18:00")) +
  ggtitle("Female white shark, no chum")

wsf1_t1 <- which(ws_HMM$SharksexTrackNo == "WSF1 T1")
wsf1_t1_covar <- HMM_covar[wsf1_t1,]
wsf1_dateTime <- 60*hour(ws_HMM$dateTime[seq_along(wsf1_t1)]) +
  minute(ws_HMM$dateTime[seq_along(wsf1_t1)])
dim(HMM_covar)
beta_mat <- rbind(c(randeff[1, 1, 1], 
                    beta[1, 1, c(1, 3, 5, 7)]), 
                  c(randeff[1, 1, 2], 
                    beta[1, 1, c(2, 4, 6, 8)]))
tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                               beta = t(beta_mat),
                               covs = wsf1_t1_covar)
wsf1_chum_tpm <- data.frame(
  omega11 = tpm[1, 1, ], 
  omega12 = tpm[1, 2, ], 
  omega21 = tpm[2, 1, ], 
  omega22 = tpm[2, 2, ], 
  dateTime = 60*hour(ws_HMM$dateTime[seq_along(wsf1_t1)]) +
    minute(ws_HMM$dateTime[seq_along(wsf1_t1)]),
  draw = 1
)
for (k in 2:500) {
  beta_mat <- rbind(c(mu_tpm[k, 1, 1], 
                      beta[k, 1, c(1, 3, 5, 7)]), 
                    c(mu_tpm[k, 1, 2], 
                      beta[k, 1, c(2, 4, 6, 8)]))
  tpm <- moveHMM:::trMatrix_rcpp(nbStates = 2,
                                 beta = t(beta_mat),
                                 covs = wsf1_t1_covar)
  wsf1_chum_tpm <- rbind(
    wsf1_chum_tpm,
    data.frame(
      omega11 = tpm[1, 1, ], 
      omega12 = tpm[1, 2, ], 
      omega21 = tpm[2, 1, ], 
      omega22 = tpm[2, 2, ], 
      dateTime = 60*hour(ws_HMM$dateTime[seq_along(wsf1_t1)]) +
        minute(ws_HMM$dateTime[seq_along(wsf1_t1)]),
      draw = k
    )
  )
}
colnames(wsf1_chum_tpm)[1] <- "omega[11](t)"
colnames(wsf1_chum_tpm)[2] <- "omega[12](t)"
colnames(wsf1_chum_tpm)[3] <- "omega[21](t)"
colnames(wsf1_chum_tpm)[4] <- "omega[22](t)"
wsf1_chum <- which(HMM_covar[wsf1_t1, 2] ==1)
wsf1_chum_tpm_long <- wsf1_chum_tpm |>
  pivot_longer(!c(dateTime, draw), names_to = "omega")

p3 <- ggplot(wsf1_chum_tpm_long, aes(dateTime, value)) +
  geom_line(aes(group = draw),
            alpha = 0.1, col = "darkgrey") +
  geom_vline(xintercept = wsf1_dateTime[wsf1_chum],
             linetype = 1, col = "lightgrey", alpha = 0.15) +
  facet_wrap(~omega, nrow = 2, ncol = 2, labeller = label_parsed) +
  theme_classic() + ylab("Probability") +
  xlab("Time") +
  scale_x_continuous(breaks = c(540, 720, 900, 1080),
                     labels = c("09:00", "12:00", "15:00", "18:00")) +
  ggtitle("White shark female 1, track 1")

#| label: fig-2state_tpmcov_crencp
(p1 + p2)/p3

#' <br />
#'
#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' # Licenses {.unnumbered}
#'
#' * Code &copy; 2025, Vianey Leos Barajas, licensed under BSD-3.
#' * Text &copy; 2025, Vianey Leos Barajas, licensed under CC-BY-NC 4.0.

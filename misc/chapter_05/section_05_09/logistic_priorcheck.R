#' ---
#' title: "Logistic regression prior predictive checking demo for Bayesian Workflow book"
#' author: "Aki Vehtari"
#' date: 2024-06-27
#' date-modified: today
#' date-format: iso
#' format:
#'   html:
#'     number-sections: true
#'     code-copy: true
#'     code-download: true
#'     code-tools: true
#' bibliography: ../../../casestudies.bib
#' ---
#' 
#' This notebook includes part of the code for the Bayesian Workflow
#' book Section 5.9 about prior predictive checking.

#+ setup, include=FALSE
knitr::opts_chunk$set(
  cache = FALSE,
  message = FALSE,
  error = FALSE,
  warning = FALSE,
  comment = NA,
  out.width = "90%"
)

#' **Load packages**
library(ggplot2)
library(patchwork)

#' # Logistic regression with normal(0, 1) prior and binary predictor
#'
prior_draws <- \(k, S, sd=1) plogis(replicate(
  n=S,
  rnorm(n=1, sd=sd) + t(rnorm(n=k, sd=sd)) %*% rbinom(n=k, size=1, prob=0.5))
  )
S <- 1000*100
prior_draws_1 <- prior_draws(k=2, S=S)
prior_draws_2 <- prior_draws(k=4, S=S)
prior_draws_3 <- prior_draws(k=15, S=S)
data.frame(p = c(prior_draws_1, prior_draws_2, prior_draws_3),
           k = rep(c(2, 4, 15), each = S)) |>
  ggplot(aes(x = p)) +
  geom_histogram(color = "gray60", fill = "gray60", bins = 100) +
  facet_grid(
    cols = vars(k),
    labeller = labeller(
      k = ~ paste(.x, "predictors")
    )
  ) +
  labs(x="Prior predictive probability", y = "") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y   = element_blank()
  )

#ggsave(root("prior_predictive_logistic_new.pdf"), width=8, height=2)

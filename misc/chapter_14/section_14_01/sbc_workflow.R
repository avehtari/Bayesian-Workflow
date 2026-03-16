#' ---
#' title: "Simple simulation based calibration checking demo"
#' author: "Martin Modrák"
#' date: 2020-10-07
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
#' This notebook includes code for the Bayesian Workflow book figures
#' in Chapter 14 *Simulation-based calibration checking*.

#+ setup, include=FALSE
knitr::opts_chunk$set(
  cache = FALSE,
  message = FALSE,
  error = FALSE,
  warning = FALSE,
  comment = NA,
  out.width = "95%"
)

#' 
#' **Load packages**
#| cache: FALSE
library(dplyr)
library(SBC)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))

#' Figure 14.1
#| label: fig-sbc_point_estimate
#| fig-width: 6
#| fig-height: 2.5
set.seed(65123365)
N <- 4000
true_val <- 1.58
df_sbc1 <- data.frame(
  sc = rep(paste0("scenario ", 1:3), each = N),
  mu = c(
    rnorm(N, true_val, 1.3),
    rnorm(N, 0, 1),
    if_else(runif(N) < 0.5, rnorm(N, true_val, 1.3), rnorm(N, 6, 0.75))
  )
)
ggplot(data = df_sbc1, aes(x = mu)) +
  geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = true_val, linewidth = 2, linetype = "dashed", color = "lightblue") +
  facet_wrap(~sc) +
  scale_x_continuous("mu", breaks = c(-3,0,3,6,9)) 
#ggsave("sbc_point_estimate.pdf")

#' Figure 14.3
#| label: fig-sbc_rank_hist_1000
#| fig-width: 7.32
#| fig-height: 4.51
res_1000 <- SBC_example_results("visualizations", n_sims = 1000)
plot_rank_hist(res_1000) +
  scale_x_continuous("rank")
#ggsave(filename = "sbc_rank_hist_1000.pdf", width = 7.32, height = 4.51)

#' Figure 14.4
#| label: fig-sbc_ecdf_diff_1000
#| fig-width: 7.32
#| fig-height: 4.51
plot_ecdf_diff(res_1000) +
  theme(legend.position = "bottom")
#ggsave(filename = "sbc_ecdf_diff_1000.pdf", width = 7.32, height = 4.51)

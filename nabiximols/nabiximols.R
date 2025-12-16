#' ---
#' title: "Nabiximols treatment efficiency"
#' author: "Aki Vehtari"
#' date: 2024-02-23
#' date-modified: today
#' date-format: iso
#' format:
#'   html:
#'     toc: true
#'     toc-location: left
#'     toc-depth: 2
#'     number-sections: true
#'     smooth-scroll: true
#'     theme: readable
#'     code-copy: true
#'     code-download: true
#'     code-tools: true
#'     embed-resources: true
#'     anchor-sections: true
#'     html-math-method: katex
#' bibliography: ../casestudies.bib
#' ---

#' # Introduction
#' 
#' This notebook was inspired by [a question by Llew Mills in Stan
#' Discourse](https://discourse.mc-stan.org/t/comparing-models-with-different-functional-forms-of-the-same-y-variable-using-loo-compare/34161)
#' about comparison of continuous and discrete observation
#' models. This question has been asked several times before, and
#' [an answer is included in CV-FAQ](https://users.aalto.fi/~ave/CV-FAQ.html#13_Can_cross-validation_be_used_to_compare_different_observation_models__response_distributions__likelihoods), too. CV-FAQ answer states
#' 
#' > You can’t compare densities and probabilities directly. Thus you
#' can’t compare model given continuous and discrete observation
#' models, unless you compute probabilities in intervals from the
#' continuous model (also known as discretising continuous model).
#'
#' The answer is complete in theory, but doesn't tell how to do the
#' discretization in practice. The first part of this notebook shows
#' how easy that discretization can be in a special case of counts.
#' 
#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=TRUE, comment=NA, out.width='95%')
#' **Load packages**
#| code-fold: true
#| cache: FALSE
library(tidyr)
library(dplyr)
library(tibble)
library(modelr)
library(loo)
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 1)
library(rstan)
## options(mc.cores = 4)
library(posterior)
options(posterior.num_args=list(digits=2))
library(pillar)
options(pillar.negative = FALSE)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size=14))
library(tidybayes)
library(ggdist)
library(patchwork)
library(tinytable)
options(tinytable_format_num_fmt = "significant_cell", tinytable_format_digits = 2, tinytable_tt_digits=2)
library(matrixStats)
library(reliabilitydiag)
library(priorsense)
set1 <- RColorBrewer::brewer.pal(3, "Set1")

#' # Data
#'
#' Data comes from a study [Nabiximols for the Treatment of Cannabis
#' Dependence: A Randomized Clinical Trial] by @Lintzeris+etal:2019:Nabiximols,
#' and was posted in that Discourse thread by Mills (one of the co-authors).
#' The data includes 128 participants (`id`) in two groups (`group` = `Placebo`
#' or `Nabiximols`). The data are used to illustrate various workflow
#' aspects including comparison of models, but the analysis here do
#' not match exactly the analysis in the paper or in the follow-up
#' paper by @Lintzeris+etal:2020:nabiximols, and further improvements
#' in the analysis could be made by using additional data not included
#' in the data used here.
#'
#' > Participants received 12-week treatment involving weekly clinical
#' reviews, structured counseling, and flexible medication doses—up to
#' 32 sprays daily (tetrahydrocannabinol, 86.4 mg, and cannabidiol, 80
#' mg), dispensed weekly.
#'
#' The number of cannabis used days (`cu`) was for previous $28$ days
#' (`set`) asked after 0, 4, 8, and 12 weeks (`week` as a factor).
#' 

id <- factor(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 12, 12, 13, 14, 15, 16, 16, 17, 18, 18, 18, 18, 19, 20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 23, 24, 24, 24, 24, 25, 25, 25, 25, 26, 27, 27, 28, 28, 28, 28, 29, 30, 30, 30, 30, 31, 31, 32, 32, 32, 32, 33, 33, 33, 34, 34, 34, 35, 35, 36, 36, 37, 37, 37, 37, 38, 39, 39, 39, 39, 40, 40, 40, 41, 42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 45, 45, 46, 46, 46, 46, 47, 47, 47, 47, 48, 48, 49, 49, 49, 50, 50, 50, 50, 51, 51, 51, 52, 52, 52, 52, 53, 53, 53, 53, 54, 54, 55, 55, 55, 55, 56, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 62, 63, 63, 64, 64, 64, 65, 65, 65, 65, 66, 66, 66, 66, 67, 67, 67, 67, 68, 68, 68, 69, 69, 69, 69, 70, 70, 70, 70, 71, 71, 71, 71, 72, 73, 73, 73, 73, 74, 74, 74, 75, 76, 76, 76, 76, 77, 77, 77, 77, 78, 78, 78, 79, 79, 79, 79, 80, 80, 80, 80, 81, 81, 81, 81, 82, 82, 83, 83, 84, 84, 84, 85, 85, 85, 86, 86, 86, 86, 87, 87, 87, 87, 88, 88, 88, 88, 89, 89, 89, 89, 90, 90, 90, 90, 91, 91, 91, 91, 92, 92, 92, 92, 93, 93, 93, 93, 94, 94, 94, 94, 95, 95, 95, 95, 96, 96, 96, 96, 97, 97, 97, 98, 98, 98, 98, 99, 99, 99, 99, 100, 101, 101, 101, 102, 102, 102, 102, 103, 103, 103, 103, 104, 104, 105, 105, 105, 105, 106, 106, 106, 106, 107, 107, 107, 107, 108, 108, 108, 108, 109, 109, 109, 109, 110, 110, 111, 111, 112, 112, 112, 112, 113, 113, 113, 113, 114, 115, 115, 115, 115, 116, 116, 116, 116, 117, 117, 117, 117, 118, 118, 119, 119, 119, 119, 120, 120, 120, 120, 121, 121, 121, 122, 123, 123, 123, 123, 124, 124, 124, 125, 125, 125, 125, 126, 126, 126, 126, 127, 127, 128))
group <- factor(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0),
                                levels = 0:1,
                                labels = c("placebo", "nabiximols"))
week <- factor(c(0, 4, 8, 12, 0, 4, 8, 0, 4, 8, 0, 4, 8, 12, 0, 4, 0, 4, 8, 12, 0, 4, 0, 4, 8, 12, 0, 4, 8, 0, 4, 8, 12, 0, 4, 8, 0, 4, 0, 0, 0, 0, 4, 0, 0, 4, 8, 12, 0, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0, 4, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 0, 4, 0, 4, 8, 12, 0, 0, 4, 8, 12, 0, 4, 0, 4, 8, 12, 0, 4, 8, 0, 4, 12, 0, 8, 0, 4, 0, 4, 8, 12, 0, 0, 4, 8, 12, 0, 4, 8, 0, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0, 4, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0, 4, 12, 0, 4, 8, 12, 0, 4, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0, 4, 8, 12, 0, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 12, 0, 0, 4, 0, 4, 8, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 0, 4, 8, 12, 0, 4, 8, 0, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0, 4, 0, 4, 8, 0, 4, 8, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 0, 4, 8, 12, 0, 4, 8, 12, 0, 0, 4, 8, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0, 4, 0, 4, 8, 12, 0, 4, 8, 12, 0, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 12, 0, 0, 4, 8, 12, 0, 4, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 0))
cu <- c(13, 12, 12, 12, 28, 0, NA, 16, 9, 2, 28, 28, 28, 28, 28, NA, 28, 28, 17, 28, 28, NA, 16, 0, 0, NA, 28, 28, 28, 28, 17, 0, NA, 28, 27, 28, 28, 26, 24, 28, 28, 28, 25, 28, 26, 28, 18, 16, 28, 28, 7, 0, 2, 28, 2, 4, 1, 28, 28, 16, 28, 28, 24, 26, 15, 28, 25, 17, 1, 8, 28, 24, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 20, 28, 28, 28, 28, 12, 28, NA, 17, 15, 14, 28, 0, 28, 28, 28, 0, 0, 0, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 21, 24, 28, 27, 28, 28, 26, NA, 28, NA, 20, 2, 3, 7, 28, 1, 19, 8, 21, 7, 28, 28, 20, 28, 28, 28, 24, 20, 17, 11, 25, 25, 28, 26, 28, 24, 17, 16, 27, 14, 28, 28, 28, 28, 28, 28, 14, 13, 4, 24, 28, 28, 28, 21, 28, 21, 26, 28, 28, 0, 0, 28, 23, 20, 28, 20, 16, 28, 28, 28, 10, 1, 1, 2, 28, 28, 28, 28, 18, 22, 9, 15, 28, 9, 1, 20, 18, 20, 24, 28, 28, 28, 19, 28, 28, 28, 28, 28, 28, 28, 28, 28, 4, 14, 20, 28, 28, 0, 0, 0, 28, 20, 9, 24, 28, 28, 28, 28, 28, 21, 28, 28, 14, 24, 28, 23, 0, 0, 0, 28, NA, 28, NA, 28, 15, NA, 12, 25, NA, 28, 2, 0, 0, 28, 10, 0, 0, 28, 0, 0, 0, 23, 0, 0, 0, 28, 0, 0, 0, 28, 0, 0, 0, 28, 2, 1, 0, 21, 14, 7, 8, 28, 28, 28, 0, 28, 28, 20, 18, 24, 0, 0, 0, 28, 15, NA, 28, 1, 1, 2, 28, 1, 0, 0, 28, 28, 14, 21, 25, 19, 16, 13, 28, 28, 28, 28, 28, 28, 28, 27, 19, 21, 18, 1, 0, 0, 28, 28, 28, 28, 28, 24, 27, 28, 18, 0, 3, 8, 28, 28, 28, 9, 20, 25, 20, 12, 19, 0, 0, 0, 27, 28, 0, 0, 0, 20, 17, 16, 14, 28, 7, 0, 1, 28, 24, 28, 25, 23, 20, 28, 14, 16, 7, 28, 28, 26, 28, 28, 26, 28, 28, 28, 24, 20, 28, 28, 28, 28, 28, 8, 6, 4, 28, 20, 28)
set <- rep(28, length(cu))
cu_df <- data.frame(id, group, week, cu, set)

#' I remove the rows with NA's as they are not used for posteriors anyway
cu_df <- cu_df |>
  drop_na(cu)

#' It's good to visualize the data distribution. There is a clear
#' change in the distribution going from week $0$ to later weeks.
#| label: fig-data
#| fig-height: 4.5
#| fig-width: 8
cu_df |> 
  ggplot(aes(x=cu)) +
  geom_histogram(breaks=seq(-.5,28.5,by=1)) +
  scale_x_continuous(breaks=c(0,10,20,28)) +
  facet_grid(group ~ week, switch="y", axes="all_x",
             labeller = labeller(group = label_value, week = label_both)) +
  labs(y="")

#' # Initial models
#' 
#' Mills provided two `brms` [@Buerkner:2017:brms] models with
#' specified priors. The first one is a normal regression model with
#' varying intercept for each participant (`id`).
#| results: hide
#| cache: true
fit_normal <- brm(formula = cu ~ group*week + (1 | id),
          data = cu_df,
          family = gaussian(),
          prior = c(prior(normal(14, 1.5), class = Intercept),
                    prior(normal(0, 11), class = b),
                    prior(cauchy(1,2), class = sd)),
          save_pars = save_pars(all=TRUE),
          seed = 1234,
          refresh = 0)
fit_normal <- add_criterion(fit_normal, criterion="loo", save_psis=TRUE,
                            moment_match=TRUE)

#' The second provided models is binomial model with the number of
#' trials being $28$ for each outcome (`cu`)
#| results: hide
#| cache: true
fit_binomial <- brm(formula = cu | trials(set)  ~ group*week + (1 | id),
        data = cu_df,
        binomial(link = logit),
        prior = c(prior(normal(0, 1.5), class = Intercept),
                  prior(normal(0, 1), class = b),
                  prior(cauchy(0,2), class = sd)),
        save_pars = save_pars(all=TRUE),
        seed = 1234,
        refresh = 0)
fit_binomial <- add_criterion(fit_binomial, criterion="loo", save_psis=TRUE,
                              moment_match=TRUE)

#' Mills compared the models using PSIS-LOO
#' [@Vehtari+Gelman+Gabry:2017:psisloo;@Vehtari+etal:2024:loo] and asked is this valid:
loo_compare(fit_normal, fit_binomial)

#' # Comparison of discrete and continuous models
#' 
#' At this point I ignore the fact that there many high Pareto-$k$
#' values indicating unreliable PSIS-LOO estimates
#' [@Vehtari+etal:PSIS:2024], and first discuss comparison of
#' predictive probabilities and predictive probability densities.
#'
#' I pick the `id=2` and compute the posterior predictive probability
#' of each possible outcome $(0,1,2,\ldots,28)$. At this point we
#' don't need to consider LOO-CV issues, and can focus on comparison
#' of posterior predictive probabilities. We can use `log_lik()`
#' function to compute the log predictive probabilities. In the
#' following plot the area of each bar is equal to the probability of
#' corresponding count outcome and these probabilities sum to 1. They
#' grey bar illustrates the posterior predictive probability of
#' `cu=12`.
i<-2
cu_df_predi <- cu_df[i,] |> select(!cu) |> expand_grid(cu=0:28)
S<-4000
p <- exp(colLogSumExps(log_lik(fit_binomial,
                               newdata = cu_df_predi))-log(S))
predi <- data.frame(cu=0:28, p=p)
predi |>
  ggplot(aes(x=cu, y=p)) +
  geom_col(width=1,fill="white", color="blue") +
  geom_col(data=predi[cu[i]+1,],width=1,fill="gray", color="blue") +
  labs(y="Probability") +
  scale_x_continuous(breaks=c(0,10,20,28), minor_breaks = 0:28) +
  guides(x = guide_axis(minor.ticks = TRUE))

#' In case of normal model, the observation model is continuous, and
#' we can compute posterior predictive densities. I use `log_lik()` to
#' compute log predictive densities. I compute the predictive density
#' of `cu=12` and predictive density in fine grid for other values
#' (although the observation model is continuous, in practice with
#' computers, we can evaluate the density only in finite number of
#' points). The following plot shows the predictive density and grey
#' line at `cu=12` with end of line corresponding to density at
#' `cu=12`.
p1 <- exp(colLogSumExps(log_lik(fit_normal,
                               newdata = cu_df[i,]))-log(S))
cu_pred <- seq(-15,40,length.out=400)
cu_df_predi_2 <- cu_df[i,] |> select(!cu) |> expand_grid(cu=cu_pred)
p <- exp(colLogSumExps(log_lik(fit_normal,
                               newdata = cu_df_predi_2))-log(S))
predi <- data.frame(cu=cu_pred, p=p)
pn <- predi |>
  ggplot(aes(x=cu, y=p)) +
  geom_line(color="blue") +
  annotate(geom="segment", x=cu_df[i,'cu'], y=0, xend=cu_df[i,'cu'], yend=p1, color="gray") +
  labs(y="Density")+
  scale_x_continuous(breaks=c(0,10,20,28), minor_breaks = 0:28) +
  guides(x = guide_axis(minor.ticks = TRUE))
pn

#' We can't compare probabilities and densities directly, but we can
#' discretize the density to get probabilities. As the outcomes are
#' integers $(0,1,2,\ldots,28)$, we can compute probabilities for
#' intervals $((-0.5, 0.5), (0.5,1.5),
#' (1.5,2.5),\ldots,(27.5,28.5))$. The following plot shows the
#' vertical lines for the edges of these intervals.
p2 <- exp(colLogSumExps(log_lik(fit_normal,
                               newdata = cu_df[i,] |> select(!cu) |> expand_grid(cu=seq(-0.5,28.5,by=1))))-log(S))
pn +
  annotate(geom="segment",x=seq(-0.5,28.5,by=1),y=0,xend=seq(-0.5,28.5,by=1),yend=p2, color="blue") +
  scale_x_continuous(breaks=c(0,10,20,28), minor_breaks = 0:28, lim=c(-0.5,28.5))
  

#' We can then integrate the density over the interval. For example,
#' integrate predictive density from $11.5$ to $12.5$ to get a
#' probability that $11.5 < cu < 12.5$ (it doesn't matter here if the
#' interval ends are open or closed). To simplify the computation, we
#' approximate the density as piecewise constant function shown in the
#' next plot.
p3 <- exp(colLogSumExps(log_lik(fit_normal,
                               newdata = cu_df[i,] |> select(!cu) |> expand_grid(cu=seq(0,28,by=1))))-log(S))
pn <-
  predi |>
  ggplot(aes(x=cu, y=p)) +
  geom_line(color="blue", alpha=0.2) +
  annotate(geom="segment", x=cu_df[i,"cu"], y=0, xend=cu_df[i,"cu"], yend=p1, color="gray") +
  labs(y="Density")+
  scale_x_continuous(breaks=c(0,10,20,28), minor_breaks = 0:28, lim=c(-0.5,28.5)) +
  guides(x = guide_axis(minor.ticks = TRUE))
pn +
  annotate(geom="segment",x=rep(seq(-0.5,28.5,by=1),each=2)[2:59],y=0,xend=rep(seq(-0.5,28.5,by=1),each=2)[2:59],yend=rep(p3, each=2), color="blue") +
  annotate(geom="segment",x=rep(seq(-0.5,28.5,by=1),each=2)[1:58],y=rep(p3, each=2),xend=rep(seq(-0.5,28.5,by=1),each=2)[3:60],yend=rep(p3, each=2), color="blue")

#' Now the probability of each interval is approximated by the height
#' times the width of a bar. The height is the density in the middle
#' of the interval and width of the bar is 1, and thus the probability
#' value is the same as the density value! In this case this is simple
#' as the counts are integers and the distance between counts is
#' 1. The following plot shows the probability of `cu=12`.
#| label: fig-discretized-normal
#| fig-height: 4
#| fig-width: 7
pn +
  annotate(geom="segment",x=rep(seq(-0.5,28.5,by=1),each=2)[2:59],y=0,xend=rep(seq(-0.5,28.5,by=1),each=2)[2:59],yend=rep(p3, each=2), color="blue") +
  annotate(geom="segment",x=rep(seq(-0.5,28.5,by=1),each=2)[1:58],y=rep(p3, each=2),xend=rep(seq(-0.5,28.5,by=1),each=2)[3:60],yend=rep(p3, each=2), color="blue") +
  geom_col(data=data.frame(cu=12,p=p3[12+1]),width=1,fill="gray", color="blue")

#' Thus, in this case, the LOO-ELPD comparison is valid.
#'
#' The normal model predictions are not constrained between $0$ and
#' $28$ (or $-0.5$ and $28.5$), and thus the probabilities for `cu`
#' $\in (0,1,2,\ldots,28)$ do not necessarily sum to 1.
sum(p3)

#' We could switch to truncated normal, but for this data, it would
#' not work well, and as we have better options, I skip illustration of
#' that.
#' 
#' Before continuing with better models, I illustrate the case where
#' the integration interval is not 1, and the comparison would require
#' more work. It is common with continuous targets to normalize the
#' target to have zero mean and unit standard deviation.
cu_scaled_df <- cu_df |>
  mutate(cu_scaled = (cu - mean(cu)) / sd(cu))
#' Fit the normal model with scaled target.
#| results: hide
#| cache: true
fit_normal_scaled <- brm(formula = cu_scaled ~ group*week + (1 | id),
          data = cu_scaled_df,
          family = gaussian(),
          prior = c(prior(normal(0, 1.5), class = Intercept),
                    prior(normal(0, 11), class = b),
                    prior(cauchy(1,2), class = sd)),
          save_pars = save_pars(all=TRUE),
          seed = 1234,
          refresh = 0)
fit_normal_scaled <- add_criterion(fit_normal_scaled, criterion="loo")

#' Now the densities for the scaled target are much higher, and the
#' scaled model seems to be much better.
loo_compare(fit_normal, fit_normal_scaled)

#' However to make a fair comparison, we need to take into account the
#' scaling we used. To get the probabilities for integer counts, the
#' densities need to be multiplied by the discretization binwidth
#' $0.093$. Correspondingly in log scale we add $\log(0.093)$ to each
#' log-density, and the total adjustment is
nrow(cu_df)*log(sd(cu_df$cu))

#' Adding this to `eldp_loo` of the scaled outcome model, makes the
#' two models to have almost the same `elpd_loo` (the small
#' difference comes from not scaling the priors).
#'
#' # Model checking
#' 
#' Let's get back to the first models. The normal model was better than the
#' binomial model in the LOO comparison, but that doesn't mean it is a
#' good model. Before trying to solve the high Pareto-$k$ issues in
#' LOO computation, we can use posterior predictive checking.
#'
#' Histograms of posterior predictive replicates show that the
#' binomial model predicts lower probability for both $0$ and $28$ than
#' what is the observed proportion.
#| label: fig-ppc_hist-binomial
#| fig-height: 4
#| fig-width: 7
pp_check(fit_binomial, type="hist", ndraws=3) +
  scale_x_continuous(breaks = c(0,10,20,28)) +
  theme(axis.line.y=element_blank())
#'
#' We can also look at the marginal distribution of the data and
#' posterior replicates.
#| label: fig-ppc_bars-binomial
#| fig-height: 4
#| fig-width: 7
pp_check(fit_binomial, type="bars", ndraws=4000) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))

#' There is a clear mismatch between data and posterior predictive
#' replicates near edges.
#' 
#' This effect can be seen also in LOO-PIT-ECDF calibration check plot
#' [@Sailynoja+etal:2022:PIT-ECDF] which shows that the pointwise predictive
#' distributions are too narrow.
#| label: fig-ppc_pit_ecdf-binomial
#| fig-height: 4
#| fig-width: 4
pp_check(fit_binomial, type="loo_pit_ecdf")

#' We see that binomial model has too many PIT values near $0$ and $1$,
#' which indicates the posterior predictive intervals are too narrow,
#' and the model is not appropriately handling the overdispersion
#' compared to binomial model, even with included varying intercept
#' term `(1 | id)`.
#'
#' The normal model posterior predictive replicates look even more
#' different than the observed data, including predicting outcomes
#' less than $0$ and larger than $28$,  but the higher probability of $0$ and
#' $28$ makes the difference that the normal model wins in the LOO
#' comparison.
#| label: fig-ppc_hist-normal
#| fig-height: 4
#| fig-width: 7
pp_check(fit_normal, type="hist", ndraws=3) +
  scale_x_continuous(breaks = c(0,10,20,28)) +
  theme(axis.line.y=element_blank())

#' LOO-PIT-ECDF calibration check plot for normal model does indicate
#' some problem with too few PIT values near 0.5, but surprisingly the
#' miscalibration effects at left and right cancel each other out.
#| label: fig-ppc_pit_ecdf-normal
#| fig-height: 4
#| fig-width: 4
pp_check(fit_normal, type="loo_pit_ecdf", moment_match=TRUE)


#' We can also examine PIT values computed by comparing all posterior
#' predictive draws with observations marginally, which in this case
#' show the miscalibration very clearly.
#| fig-height: 4
#| fig-width: 4
marginal_pit <- \(y, x) {
  pit <- vapply(seq_len(length(y)), function(j) {
    sel_min <- x < y[j]
    if (!any(sel_min)) {
      pit <- 0
    } else {
      pit <- mean(sel_min)
    }
    sel_sup <- x == y[j]
    if (any(sel_sup)) {
      # randomized PIT for discrete y (see, e.g., Czado, C., Gneiting, T.,
      # Held, L.: Predictive model assessment for count data.
      # Biometrics 65(4), 1254–1261 (2009).)
      pit_sup <- pit + mean(sel_sup)
      pit <- runif(1, pit, pit_sup)
    }
  pit
  }, FUN.VALUE = 1.0)
}
ppc_pit_ecdf(pit = marginal_pit(cu_df$cu, posterior_predict(fit_normal))) +
  labs(x="Marginal posterior PIT")

#' # Model extension
#' 
#' Binomial model doesn't have overdispersion term, but the normal
#' model does. How about using beta-binomial model which is an
#' overdispersed version of the binomial model.
#' 
#| results: hide
#| cache: true
fit_betabinomial <- brm(formula = cu | trials(set) ~ group*week + (1 | id),
        data = cu_df,
        beta_binomial(),
        prior = c(prior(normal(0, 1.5), class = Intercept),
                  prior(normal(0, 1), class = b),
                  prior(cauchy(0,2), class = sd)),
        save_pars = save_pars(all=TRUE),
        seed = 1234,
        refresh = 0)
fit_betabinomial <- add_criterion(fit_betabinomial, criterion="loo", save_psis=TRUE)

#'
loo_compare(fit_normal, fit_binomial, fit_betabinomial)

#'
#' The beta-binomial model beats the normal model big time. The
#' difference is so big that it is unlikely that fixing Pareto-$k$
#' warnings matter.
#'
#' # Improved LOO computation
#' 
#' Since in case of high Pareto-$k$ warnings, the estimate is usually
#' optimistic, it is sufficient if we use moment matching
#' [@Paananen+etal:2021:implicit] for beta-binomial model.
#| cache: true
fit_betabinomial <- add_criterion(fit_betabinomial, criterion="loo",
                                  save_psis=TRUE,
                                  moment_match=TRUE,
                                  overwrite=TRUE)

#' Moment matching gets all Pareto-$k$ values below the diagnostic
#' threshold. There is negligible change in the comparison results.
loo_compare(fit_normal, fit_binomial, fit_betabinomial)

#' Truncated normal might get closer to beta-binomial than normal, but
#' as I said it had difficulties here, and since beta-binomial is
#' proper discrete observation model, I continue with that.
#' 
#' # Model checking
#'
#' The posterior predictive replicates from beta-binomial look very
#' much like the original data.
#| label: fig-ppc_hist-betabinomial
#| fig-height: 4
#| fig-width: 7
pp_check(fit_betabinomial, type="hist", ndraws=3) +
  scale_x_continuous(breaks = c(0,10,20,28)) +
  theme(axis.line.y=element_blank())

#' We can also look at the marginal distribution of the data and
#' posterior replicates.
#| label: fig-ppc_bars-betabinomial
#| fig-height: 4
#| fig-width: 7
pp_check(fit_betabinomial, type="bars", ndraws=4000) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))

#' This plot looks much better than with binomial model, although data
#' seem to have more variation.
#' 
#' LOO-PIT-ECDF calibration looks very nice and also better than for
#' the normal model.
#| label: fig-loo_pit_ecdf-betabinomial
#| fig-height: 4
#| fig-width: 4
pp_check(fit_betabinomial, type="loo_pit_ecdf", moment_match = TRUE)

#' The hierarchical beta-binomial model can also be used to illustrate
#' the importance for using crossvalidation for computing PIT
#' values. If we use the same observations to fit the model and check
#' the model, the double use of the data leads to too many PIT values
#' near 0.5 which shows an S-shape in the PIT-ECDF plot.
#| label: fig-ppc_pit_ecdf-betabinomial
#| fig-height: 4
#| fig-width: 4
pp_check(fit_betabinomial, type="pit_ecdf")

#' The data included many counts of $0$ and $28$, and we can further
#' check whether we might need to include zero-inflation or
#' $28$-inflation model component. As illustrated in [Roaches case
#' study](https://users.aalto.fi/~ave/modelselection/roaches.html)
#' PIT's may sometimes be weak to detect zero-inflation, but
#' reliability diagram [@Dimitriadis+etal:2021:reliabilitydiag]
#' can focus there in more detail.
#'
#' Calibration check with reliability diagram for $0$ vs others
#| label: fig-calibration-0-betabinomial
#| fig-height: 5
#| fig-width: 5.5
th<-0
rd=reliabilitydiag(EMOS = pmin(E_loo(0+(posterior_predict(fit_betabinomial)>th),loo(fit_betabinomial)$psis_object)$value,1),
                   y = as.numeric(cu_df$cu>th))
autoplot(rd)+
  labs(x="Predicted probability of non-zero",
       y="Conditional event probabilities")+
  bayesplot::theme_default(base_family = "sans", base_size=14)

#' Calibration check with reliability diagram for $28$ vs others
#| label: fig-calibration-28-betabinomial
#| fig-height: 5
#| fig-width: 5.5
th<-27
rd=reliabilitydiag(EMOS = pmin(E_loo(0+(posterior_predict(fit_betabinomial)>th),loo(fit_betabinomial)$psis_object)$value,1),
                   y = as.numeric(cu_df$cu>th))
autoplot(rd)+
  labs(x="Predicted probability of 28",
       y="Conditional event probabilities")+
  bayesplot::theme_default(base_family = "sans", base_size=16)

#' Although the red line did not completely stay within blue envelope,
#' these look good.
#'
#' # Prior sensitivity analysis
#'
#' The priors provided by Mills seems a bit narrow, so I do
#' prior-likelihood sensitivity analysis using powerscaling approach
#' [@Kallioinen+etal:2023:priorsense].
powerscale_sensitivity(fit_betabinomial,
                       variable=variables(as_draws(fit_betabinomial))[1:11]) |>
  tt() |>
  format_tt(num_fmt="decimal")

#' There are prior-data conflicts for all global parameters. It is
#' possible that the priors had been chosen based on substantial prior
#' information and then the posterior should be influenced by the
#' prior, too, and the conflict is not necessarily a problem.
#'
#' # Alternative priors
#' 
#' However, I decided to try slightly wider priors, especially as the
#' data seem to be quite informative
#| results: hide
#| cache: true
fit_betabinomial2 <- brm(formula = cu | trials(set) ~ group*week + (1 | id),
        data = cu_df,
        beta_binomial(),
        prior = c(prior(normal(0, 3), class = Intercept),
                  prior(normal(0, 3), class = b),
                  prior(normal(0, 3), class = sd)),
        save_pars = save_pars(all=TRUE),
        seed = 1234,
        refresh = 0)

#' There are no prior data conflicts.
powerscale_sensitivity(fit_betabinomial2,
                       variable=variables(as_draws(fit_betabinomial2))[1:11]) |>
  tt() |>
  format_tt(num_fmt="decimal")

#' We can also do LOO comparison, which indicates that the wider
#' priors provide a tiny bit better predictive performance.
#| cache: true
fit_betabinomial2 <- add_criterion(fit_betabinomial2, criterion="loo", save_psis=TRUE, moment_match=TRUE, overwrite=TRUE)
#| cache: false
loo_compare(fit_betabinomial, fit_betabinomial2)

#' # Model refinement
#'
#' All the above models included the baseline `week=0` in the
#' interaction term with `group`, which does not make sense as the
#' group should not affect the baseline. I modify the data and models by
#' moving the baseline `week=0` `cu`s to be pre-treatment covariate.
#| results: hide
#| cache: true
cu_df_b <- cu_df |> filter(week != 0) |>
  mutate(week = droplevels(week))
cu_df_b <- left_join(cu_df_b,
                     cu_df |> filter(week == 0) |> select(id, cu),
                     by="id",
                     suffix=c("","_baseline"))
fit_betabinomial2b <- brm(formula = cu | trials(set) ~ group*week + cu_baseline + (1 | id),
        data = cu_df_b,
        beta_binomial(),
        prior = c(prior(normal(0, 3), class = Intercept),
                  prior(normal(0, 3), class = b)),
        save_pars = save_pars(all=TRUE),
        seed = 1234,
        refresh = 0)
fit_betabinomial2b <- add_criterion(fit_betabinomial2b, criterion="loo",
                                    save_psis=TRUE, moment_match=TRUE)

#' To compare the new model against the previous ones, we exclude the
#' pointwise elpds for week 0.
loo2 <- loo(fit_betabinomial2)
loo2$pointwise <- loo2$pointwise[cu_df$week!="0",]
loo_compare(loo2, loo(fit_betabinomial2b))

#' The quick fix to exclude week 0, does not change the `yhash` and we
#' get a warning which we can ignore. The new model is clearly better than
#' the previous best model.
#'
#' By comparing `phi` parameter posteriors, we see that the
#' overdispersion of the beta-binomial in the new model is smaller
#' (smaller `phi` means bigger overdispersion).
as_draws_df(fit_betabinomial2) |>
  subset_draws(variable="phi") |>
  summarise_draws() |>
  tt()
as_draws_df(fit_betabinomial2b) |>
  subset_draws(variable="phi") |>
  summarise_draws() |>
  tt()

#' # Model checking
#' 
#' The marginal distribution of the posterior replicates matches the
#' data better than earlier models.
#| label: fig-ppc_bars-betabinomial2b
#| fig-height: 4
#| fig-width: 7
pp_check(fit_betabinomial2b, type="bars", ndraws=4000) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))

#' 
#' LOO-PIT-ECDF plot looks good:
#| label: fig-ppc_pit_ecdf-betabinomial2b
#| fig-height: 4
#| fig-width: 4
pp_check(fit_betabinomial2b, type="loo_pit_ecdf", moment_match = TRUE)


#' Calibration check with reliability diagrams for $0$ vs others and
#' $28$ vs others look better than for the previous model.
#| label: fig-calibration-0-betabinomial2
#| fig-height: 5
#| fig-width: 5.5
th<-0
rd=reliabilitydiag(EMOS = pmin(E_loo(0+(posterior_predict(fit_betabinomial2b)>th),loo(fit_betabinomial2b)$psis_object)$value,1),
                   y = as.numeric(cu_df_b$cu>th))
autoplot(rd)+
  labs(x="Predicted probability of non-zero",
       y="Conditional event probabilities")+
  bayesplot::theme_default(base_family = "sans", base_size=16)
#'
#| label: fig-calibration-28-betabinomial2
#| fig-height: 5
#| fig-width: 5.5
th<-27
rd=reliabilitydiag(EMOS = pmin(E_loo(0+(posterior_predict(fit_betabinomial2b)>th),loo(fit_betabinomial2b)$psis_object)$value,1),
                   y = as.numeric(cu_df_b$cu>th))
autoplot(rd)+
  labs(x="Predicted probability of 28",
       y="Conditional event probabilities")+
  bayesplot::theme_default(base_family = "sans", base_size=16)

#' # Treatment effect
#'
#' As the treatment is included in interaction term, looking at the
#' univariate posterior marginal is not sufficient for checking
#' whether the treatment helps.
#'
#' We can examine the effect of the treatment by looking at the
#' posterior predictive distribution and the expectation of the
#' posterior predictive distribution for a new person (new `id`).  For
#' each person we can compare the prediction to the counterfactul
#' case, that is, what if that person would have been in the other
#' treatment group. We can compare the difference from being in the
#' placebo to being in the Nabiximols group.
#' 
#' We use `tidybayes` [@Kay:2023:tidybayes] and `ggdist`
#' [@Kay:2024:ggdist] to generate and plot predictions.
#' 
#' We start by examining the posterior predictive distribution. We
#' predict `cu` for a new `id=129` given `cu_baseline=28` (median
#' value). The distribution is wide as it included the aleatoric
#' uncertainty from the unknown `id` specific intercept and
#' beta-binomial distribution.
#| label: fig-posterior_prediction-betabinomial2
cu_df_b |>
  data_grid(group, week, cu_baseline=28, id=129, set=28) |>
  add_predicted_draws(fit_betabinomial2b, allow_new_levels=TRUE) |>
  ggplot(aes(x = .prediction)) +
  facet_grid(group ~ week, switch = "y", axes="all_x",
             labeller = labeller(group = label_value, week = label_both)) +
  stat_dotsinterval(quantiles=100, fill=set1[2], slab_color=set1[2], binwidth=2/3, overflow = "keep")+
  coord_cartesian(expand = FALSE, clip="off") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(x="", y="") +  
  scale_x_continuous(breaks=c(0,10,20,28), lim=c(-0.5,28.5)) +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_blank())

#' The difference between groups is the biggest in week 12. The next
#' plot shows the comparison of predicted `cu` between the groups in
#' different weeks.
#| label: fig-posterior_prediction-diff-betabinomial2
cu_df_b |>
  data_grid(group, week, cu_baseline=28, id=129, set=28) |>
  add_predicted_draws(fit_betabinomial2b, allow_new_levels=TRUE) |>
  compare_levels(.prediction, by=group) |>
  ggplot(aes(x = .prediction)) +
  facet_grid(. ~ week, switch = "y", axes="all_x",
             labeller = labeller(group = label_value, week = label_both)) +
  stat_dotsinterval(quantiles=100, fill=set1[2], slab_color=set1[2], binwidth=2, overflow = "keep")+
  coord_cartesian(expand = FALSE, clip="off") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(x="Difference in cu given placebo vs Nabiximols", y="") +
  geom_vline(xintercept=0, color=set1[1], linewidth=1, alpha=0.3) +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_blank())

#' For a new individual the posterior predictive distribution
#' indicates 90% possibility of smaller `cu` with Nabiximols than with
#' placebo. I have used dot plots with 100 dots, as dot plots are
#' better than kernel density estimates and histograms for showing
#' spikes (as here at 0) and make it easy to estimate tail
#' probabilities by counting the dots.
#'
#' For comparison I plot 1) the ggplot2 default KDE with
#' stat_density, 2) ggdist KDE with stat_slabinterval, and 3) ggplot
#' default histogram with stat_histogram for the week 12
#' difference. The ggplot2 KDE oversmooths a lot, the ggdist KDE and
#' ggplot2 histogram smooth less, but it is not as easy to estimate
#' the tail probabilities as with dots plot.
#| label: fig-kde-kde-hist-comparison-betabinomial2
pp <- cu_df_b |>
  data_grid(group, week=12, cu_baseline=28, id=129, set=28) |>
  add_predicted_draws(fit_betabinomial2b, allow_new_levels=TRUE) |>
  compare_levels(.prediction, by=group) |>
  ggplot(aes(x = .prediction)) +
  coord_cartesian(expand = FALSE, clip="off") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(x="Difference in cu given placebo vs Nabiximols", y="") +
 theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_blank())
(pp + stat_density(fill=set1[2], color=set1[2])+labs(title="KDE ggplot2")+
  geom_vline(xintercept=0, color=set1[1], linewidth=1, alpha=0.3)) +
(pp + stat_slabinterval(fill=set1[2], slab_color=set1[2])+labs(title="KDE ggdist")+
  geom_vline(xintercept=0, color=set1[1], linewidth=1, alpha=0.3)) +
(pp + geom_histogram(fill=set1[2], color=set1[2])+labs(title="Hist. ggplot2")+
   geom_vline(xintercept=0, color=set1[1], linewidth=1, alpha=0.3)) +
  plot_layout(axes = "collect")  

#' 
#' To get more accuracy we can remove the aleatoric uncertainty and
#' focus on expected effect. The following plots shows the expectation
#' of the posterior predictive distribution for `cu` for a new
#' `id=129` given `cu_baseline=28` (expectation as in average over
#' many new individuals). The distribution is much narrower as it
#' includes only the epistemic uncertainty, but not the aleatoric
#' uncertainty (using `re_formula=NA` excludes the uncertainty in the
#' unknown `id` specific intercept and using `_epred` gives just the
#' expectation excluding the uncertainty in beta-binomial
#' distribution).
#| label: fig-epred-betabinomial2
cu_df_b |>
  data_grid(group, week, cu_baseline=28, id=129, set=28) |>
  add_epred_draws(fit_betabinomial2b, re_formula=NA, allow_new_levels=TRUE) |>
  ggplot(aes(x = .epred)) +
  facet_grid(group ~ week, switch = "y", axes="all_x",
             labeller = labeller(group = label_value, week = label_both)) +
  stat_dotsinterval(quantiles=100, fill=set1[2], slab_color=set1[2], layout="swarm", binwidth=1, overflow = "keep")+
  coord_cartesian(expand = FALSE, clip="off") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(x="", y="") +  
  scale_x_continuous(breaks=c(0,10,20,28), lim=c(-0.5,28.5)) +
  ylim(c(0, 1)) +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_blank())

#' The next plot shows the comparison of expected `cu` between the
#' groups in different weeks.
#| label: fig-epred-diff-betabinomial2
cu_df_b |>
  data_grid(group, week, cu_baseline=28, id=129, set=28) |>
  add_epred_draws(fit_betabinomial2b, re_formula=NA, allow_new_levels=TRUE) |>
  compare_levels(.epred, by=group) |>
  ggplot(aes(x = .epred, y=week)) +
  stat_dotsinterval(quantiles=100, fill=set1[2], slab_color=set1[2], layout="swarm")+
  coord_cartesian(expand = FALSE, clip="off") +
  scale_x_continuous(breaks=seq(-20,5,by=5), lim=c(-21,8)) +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(x="Difference in expected cu given placebo vs Nabiximols", y="Week") +  
  geom_vline(xintercept=0, color=set1[1], linewidth=1, alpha=0.3) +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_blank())

#' For new individuals the posterior predictive distribution indicates
#' 99% possibility of smaller expected `cu` with Nabiximols than with
#' placebo. 
cu_df_b |>
  data_grid(group, week=12, cu_baseline=28, id=129, set=28) |>
  add_epred_draws(fit_betabinomial2b, re_formula=NA, allow_new_levels=TRUE) |>
  compare_levels(.epred, by=group) |>
  mutate(p=as.numeric(.epred<0)) |>
  pull("p") |>
  mean()

#' # Sensitivity to model choice
#'
#' Finally we check how much there is difference in the conclusions,
#' if we had used the continuos normal model.  We include also models
#' that are closer to the model reported by
#' @Lintzeris+etal:2019:Nabiximols, which used sum of days of cannabis
#' use across the 12-week trial (the follow-up paper by
#' @Lintzeris+etal:2020:nabiximols considered several 4 week
#' intervals). We build normal and beta-binomial models which match
#' the model in the paper, expect that we don't include site factor as
#' this was not available for us.
#| results: hide
#| cache: true
fit_normal2b <- brm(formula = cu ~ group*week + cu_baseline + (1 | id),
        data = cu_df_b,
        gaussian(),
        prior = c(prior(normal(0, 3), class = Intercept),
                  prior(normal(0, 3), class = b)),
        save_pars = save_pars(all=TRUE),
        seed = 1234,
        refresh = 0)
fit_normal2b <- add_criterion(fit_normal2b, criterion="loo",
                                    save_psis=TRUE, moment_match=TRUE)

#'
#' LOO-PIT-ECDF hints slight miscalibration
#| label: fig-ppc_pit_ecdf-normal2b
#| fig-height: 4
#| fig-width: 4
pp_check(fit_normal2b, type="loo_pit_ecdf")

#'
#' Marginal PIT-ECDF shows clear miscalibration
#| fig-height: 4
#| fig-width: 4
ppc_pit_ecdf(pit = marginal_pit(cu_df$cu, posterior_predict(fit_normal2b))) +
  labs(x="Marginal posterior PIT")

#' The beta-binomial model beats the normal model big time. 
loo_compare(fit_normal2b, fit_betabinomial2b)

#' Next we build models predicting total cu for all weeks (except 0), by dropping
#' week covariate.
#' 
#| results: hide
#| cache: true
cu_df_c <- cu_df_b |>
  group_by(id,group,cu_baseline) |>
  summarise(cu_total=sum(cu), set_total=sum(set)) |>
  as.data.frame()
fit_normal2c <- brm(formula = cu_total ~ group + cu_baseline,
        data = cu_df_c,
        family = gaussian(),
        prior = c(prior(normal(0, 3), class = Intercept),
                  prior(normal(0, 3), class = b)),
        save_pars = save_pars(all=TRUE),
        seed = 1234,
        refresh = 0)
fit_normal2c <- add_criterion(fit_normal2c, criterion="loo",
                                    save_psis=TRUE, moment_match=TRUE)
#| results: hide
#| cache: true
fit_betabinomial2c <- brm(formula = cu_total | trials(set_total) ~ group + cu_baseline,
        data = cu_df_c,
        beta_binomial(),
        prior = c(prior(normal(0, 3), class = Intercept),
                  prior(normal(0, 3), class = b)),
        save_pars = save_pars(all=TRUE),
        seed = 1234,
        refresh = 0)
fit_betabinomial2c <- add_criterion(fit_betabinomial2c, criterion="loo",
                                    save_psis=TRUE, moment_match=TRUE)

#' We're not able to compare the models predicting total count in 12
#' weeks and models predicting counts in three 4-week period using
#' `elpd_loo` as the targets are different.
#'

#' We plot here the difference either for the last 4 week period or
#' for the all weeks depending on the model.
#' 

dat_bb2b <- cu_df_b |>
  data_grid(group, week=12, cu_baseline=28, id=129, set=28) |>
  add_epred_draws(fit_betabinomial2b, re_formula=NA, allow_new_levels=TRUE) |>
  compare_levels(.epred, by=group) |>
  mutate(model="beta-binomial model\nweeks 9-12")

dat_n2b <- cu_df_b |>
  data_grid(group, week=12, cu_baseline=28, id=129, set=28) |>
  add_epred_draws(fit_normal2b, re_formula=NA, allow_new_levels=TRUE) |>
  compare_levels(.epred, by=group) |>
  mutate(model="normal model\nweeks 9-12")

dat_n2c <- cu_df_c |>
  data_grid(group, cu_baseline=28, id=129, set_total=84) |>
  add_epred_draws(fit_normal2c, re_formula=NA, allow_new_levels=TRUE) |>
  compare_levels(.epred, by=group) |>
  mutate(model="normal model\nweeks 1-12")

dat_bb2c <- cu_df_c |>
  data_grid(group, cu_baseline=28, id=129, set_total=84) |>
  add_epred_draws(fit_betabinomial2c, re_formula=NA, allow_new_levels=TRUE) |>
  compare_levels(.epred, by=group) |>
  mutate(model="beta-binomial model\nweeks 1-12")

#| label: fig-epred-diff-4models
rbind(dat_bb2b, dat_n2b, dat_bb2c, dat_n2c) |>
  ggplot(aes(x = .epred, y=model)) +
  stat_dotsinterval(quantiles=100, layout="swarm", fill=set1[2], slab_color=set1[2])+
  coord_cartesian(expand = FALSE, clip="off") +
  theme(strip.background = element_blank(), strip.placement = "outside",
        legend.position="none") +
  scale_x_continuous(lim=c(-27,7), breaks=seq(-25,5,by=5)) +
  labs(x="Difference in expected cu given placebo vs Nabiximols", y="") +  
  geom_vline(xintercept=0, color=set1[1], linewidth=1, alpha=0.3) +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_blank())

#' The normal models underestimate the magnitude of the change and are
#' overconfident having much narrower posteriors. When we compare
#' posteriors for effect in weeks 9-12 given the models which included
#' week as a covariate, the conclusion about benefit of Nabiximols is
#' the same with normal and beta-binomial model despite that the
#' normal model underestimates the effect. When we compare posteriors
#' for effect in weeks 1-12 given the models which did not include
#' week as a covariate, the conclusion on benefit of Nabiximols is
#' clearly different. In this case, we know based on LOO-CV
#' comparisons and posterior predictive checking that the
#' beta-binomial has much better performance and matches the data much
#' better, and we can safely drop the normal models from
#' consideration. The conclusions using the models with week as a
#' covariate or all counts summed together are similar, although
#' focusing in the last four week period gives sharper posterior which
#' is plausible as the effect seems to be increasing in time.
#'
#' The results here are illustrating that the choice of data model can
#' matter for the conclusions, but also the best models shown here
#' could be even further improved by using the additional data used by
#' @Lintzeris+etal:2020:nabiximols.
#' 

#' # Predictive performance comparison
#'
#' Above, posterior predictive checking, LOO predictive checking, and
#' LOO-ELPD comparisons were important for model checking and
#' comparison, and only after we trusted that there is no significant
#' discrepancy between the model predictions and the data, we did look
#' at the treatment effect. We used LOO model comparison to assess the
#' the more elaborate models did provide better predictions, but we did
#' not use it above for assessing the predictive performance of using
#' the treatment effect in the model.
#'
#' Let's now build the mode without treatment group variable and check
#' whether there is significant difference in predictive performance.
fit_betabinomial3b <- brm(cu | trials(set) ~ week + cu_baseline + (1 | id),
        data = cu_df_b,
        beta_binomial(),
        prior = c(prior(normal(0, 9), class = Intercept),
                  prior(normal(0, 3), class = b)),
        save_pars = save_pars(all=TRUE),
        seed = 1234,
        refresh = 0)
fit_betabinomial3b <- add_criterion(fit_betabinomial3b, criterion="loo", save_psis=TRUE, moment_match=TRUE)

#' We compare the models with and without group variable:
loo_compare(fit_betabinomial2b, fit_betabinomial3b)

#' There are a few high Pareto k values for PSIS-LOO even after moment
#' matching, and for final results we would use also `reloo=TRUE`
#' argument, but the conclusions are likely to be the same.
#' 
#' There is a negligible difference in predictive performance for
#' predicting cannabis use in a 4-week period for a new person. This
#' is probably due to
#' 1. the actual effect not being very far from 0,
#' 2. the aleatoric uncertainty (modelled by the beta-binomial)
#'    for a new 4-week period is big, that is the predictive
#'    distribution is very wide and due to the constrained range
#'    has also thick tails (actually U shape), which makes the log
#'    score not to be sensitive in tails.
#'
#' As the predictive distribution is wide with thick tails, we can
#' also focus on comparing absolute error of using means of the
#' predictive distributions as point estimates. This approach has the
#' benefit that we can improve accuracy of the finite MCMC sample by
#' dropping the sampling from the varying intercept population
#' distribution and from the beta-binomial distribution (using
#' `posterior_epred(, re_formula=NA)`. `E_loo` is used to go from
#' posterior predictive draws to the mean of leave-one-out predictive
#' distribution.
ae2 <- abs(cu_df_b$cu - E_loo(posterior_epred(fit_betabinomial2b, re_formula=NA),
                              loo(fit_betabinomial2b)$psis_object, type="mean")$value)
ae3 <- abs(cu_df_b$cu - E_loo(posterior_epred(fit_betabinomial3b, re_formula=NA),
                              loo(fit_betabinomial3b)$psis_object, type="mean")$value)
#' Probability that the leave-one-out predictive absolute error is
#' smaller with model 2b (with group variable) than with model 3b
#' (without group variable) using normal approximation.
pnorm(0, mean(ae2-ae3), sd(ae2-ae3)/sqrt(257)) |> round(2)

#' By dropping out the aleatoric part from the LOO-CV predictions we
#' do see clear difference in the predictive performance if the
#' treatment group variable is dropped, but still the randomness in
#' the actual observations makes the probability of the difference to
#' be further away from 1 than when we focused in the expected
#' counterfactual predictions above.
#' 

#' 
#' 
#' We don't usually use Bayes factors for many reasons, but for
#' completness we used `bridgesampling`
#' [@Gronau:2020:bridgesampling] to get estimated Bayes factor.
(br2b <- bridgesampling::bridge_sampler(fit_betabinomial2b, silent=TRUE))
(br3b <- bridgesampling::bridge_sampler(fit_betabinomial3b, silent=TRUE))
bridgesampling::bf(br2b, br3b)

#' Like LOO-CV, Bayes factor does not see difference between models
#' with or without treatment group variable. Bayes factor is a ratio
#' of marginal likelihoods, and marginal likelihoods can be formulated
#' as sequential log score predictive performance quantity, which
#' explains why BF is also weak to rank models in case of high
#' aleatoric uncertainty.

#' <BR />
#' 
#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' <br />
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2024, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2024, Aki Vehtari, licensed under CC-BY-NC 4.0.

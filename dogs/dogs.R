#' ---
#' title: "Stochastic learning in dogs"
#' author: "Aki Vehtari and Andrew Gelman"
#' date: 2024-04-20
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
#' 
#' This notebook includes the `brms` code code for the Bayesian Workflow book
#' Chapter 18 *Posterior predictive checking: Stochastic learning in dogs*.
#' 
#' # Introduction
#'
#' This notebook is a remake of the Andrew Gelman's analysis of
#' stochastic learning in dogs data by @Bush+Mosteller:1955. Andrew
#' wrote his models in [Stan language](https://mc-stan.org/), and
#' here we use [`brms`](https://paul-buerkner.github.io/brms/)
#' and add some further diagnostics.
#'
#+ setup, include=FALSE
knitr::opts_chunk$set(
  cache=FALSE,
  message=FALSE,
  error=FALSE,
  warning=TRUE,
  comment=NA,
  out.width='95%'
)

#' **Load packages**
#| cache: FALSE
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 16))
library(patchwork)
library(tidybayes)
library(RColorBrewer)
library(priorsense)
library(loo)
library(posterior)
library(brms)
options(brms.backend = "cmdstanr")
options(mc.cores = 4)
library(dplyr)
library(tibble)
library(matrixStats)
library(tinytable)
options(
  tinytable_format_num_fmt = "significant_cell",
  tinytable_format_digits = 2,
  tinytable_tt_digits = 2
)
library(reliabilitydiag)
library(Iso)
library(assertthat)

#' # Data
#'
#' Data comes originally from a book by @Bush+Mosteller:1955. The data
#' come from 30 dogs in a stochastic learning experiment.  Each dog
#' was put in a cage where it would be shocked if it did not jump out
#' in time, a few seconds after a light goes on (we do not advocate
#' giving shocks to dogs).  After 25 tries, all of the dogs learned to
#' jump and avoid the shock.  Bush and Mosteller then posited a
#' two-parameter model which allowed different amounts of learning
#' from shocks and avoidances:
#' $$
#' \Pr(\mathrm{shock}) = a^\mathrm{(\# \ of \ previous \ shocks)}\,
#'                       b^\mathrm{(\# \ of \ previous \ avoidances)}
#' $$
#' Under this model, the probability of being shocked starts at 1,
#' which is appropriate, as there is no reason the dogs should know at
#' the start that the light would precede a shock, and indeed any dogs
#' that jumped before the first trial were excluded from the
#' experiment.  From then on, as long as the parameters $a$ and $b$
#' are between 0 and 1, the probability of shock gradually declines
#' over time.
dogs <- read.table(root("dogs", "data", "dogs.dat"), skip = 2)
shock <- ifelse(as.matrix(dogs[, 2:26]) == "S", 1, 0)
#' We create a data frame
dogs_df <- data.frame(
        shock = as.numeric(shock),
        dog = rep(1:nrow(shock), times = ncol(shock)),
        time = rep(1:ncol(shock), each = nrow(shock)))
#' The original data includes the first event which is always shock,
#' and as there is no uncertainty about that event, it can be removed
#' from the data without any loss of information.
dogs_df <- dogs_df |> dplyr::filter(time > 1)
#' Add the number of previous shocks and avoidances as covariates
dogs_df$prev_shock <- as.numeric(rowCumsums(shock)[, 1:(ncol(shock) - 1)])
dogs_df$prev_avoid <- as.numeric(rowCumsums(1 - shock)[, 1:(ncol(shock) - 1)])

#' # Model 0: Logistic regression
#'
#' We start with a simple logistic regression model
#' $$
#' \Pr(\mathrm{shock}) = {\mathrm{logit}^{-1}(\alpha + \beta t)}
#' $$
#' where $t$ denotes time. By default, `brms` assigns
#' $\mathrm{Student}_3(0, 2.5)$ prior on intercept $\alpha$, and we
#' add $\mathrm{normal}(0, 1)$ prior on coefficent $\beta$,
#' $$
#' \begin{aligned}
#' \alpha & \sim \mathrm{Student}_3(0, 2.5)\\
#' \beta & \sim \mathrm{normal}(0, 1).
#' \end{aligned}
#' $$
#| label: bfit_0
#| results: "hide"
bfit_0 <- brm(shock ~ time, family = bernoulli(),
              prior = prior(normal(0,1)),
              data = dogs_df, refresh = 0)
bfit_0 <- add_criterion(bfit_0, criterion = "loo")

#' # Model 0h: Hierarchical logistic regression
#'
#' Instead of doing model checking for the simple logistic regression, we build a
#' hierarchical model so that each dog has their own parameters. We number this as
#' 0h, so that the rest of model numbers follow Andrew's numbering.
#'
#' $$
#' \Pr(\mathrm{shock}) = {\mathrm{logit}^{-1}(\alpha_j + \beta_j t)},
#' $$
#' where $\alpha_j$ and $\beta_j$ are the parameters for dog $j$. `brms` uses following default priors, except we added the weakly informative normal prior for $\beta_0$.
#' $$
#' \begin{aligned}
#' \left(\begin{array}{c}\alpha_j \\ \beta_j\end{array} \right) & \sim \mathrm{MVN}(\mu_{\alpha,\beta}, \Sigma_{\alpha,\beta})\\
#' \mu_{\alpha} & \sim \mathrm{Student}_3(0, 2.5)\\
#' \mu_{\beta} & \sim \mathrm{normal}(0, 1)\\
#' \Sigma_{\alpha,\beta} & = \left(\begin{array}{cc}\sigma_\alpha & 0 \\ 0 & \sigma_\beta\end{array}\right) Q_{\alpha,\beta} \left(\begin{array}{cc}\sigma_\alpha & 0 \\ 0 & \sigma_\beta\end{array}\right) \\
#' \sigma_\alpha,\sigma_\beta & \sim \mathrm{Student}^{+}_3(0, 2.5) \\
#' Q_{\alpha,\beta} & \sim \mathrm{LKJ}(1).
#' \end{aligned}
#' $$
#| label: bfit_0h
#| results: "hide"
bfit_0h <- brm(shock ~ time + (time | dog), family = bernoulli(),
               prior = prior(normal(0,1)),
               data = dogs_df, refresh = 0, adapt_delta = 0.95)
bfit_0h <- add_criterion(bfit_0h, criterion = "loo", save_psis = TRUE)

#' ## Model comparison
#'
#' We compare the simple and hierarchical logistic regression using
#' PSIS-LOO-CV [@Vehtari+Gelman+Gabry:2017:psisloo], and as the latter
#' is clearly better, we can skip model checking for the simpler
#' model.
loo_compare(bfit_0, bfit_0h) |>
        as.data.frame() |>
        tibble::rownames_to_column("model") |>
        select(model, elpd_diff, se_diff) |>
        tt()

#' ## Visualize predictions
#'
#' Visualize the model fit for 9 first dogs with some help from
#' `tidybayes` [@Kay:2023:tidybayes].
#| label: fig-pred9-0h
plot_pred9_0h <- dogs_df |>
  dplyr::filter(dog <= 9) |>
        add_linpred_draws(bfit_0h, transform = TRUE) |>
        ggplot(aes(x = time, y = .linpred)) +
        stat_lineribbon(.width = c(.95), alpha = 1 / 2, color = brewer.pal(5, "Blues")[[5]]) +
  scale_fill_brewer() +
  geom_point(data = dplyr::filter(dogs_df, dog<=9), aes(x = time, y = shock, group = dog)) +
  facet_wrap( ~ dog, labeller = label_both) +
  scale_y_continuous(breaks = c(0, 1)) +
  theme(legend.position = "none") +
  labs(y = "Shocks and predicted probability of shock", title="Model 0h")
plot_pred9_0h

#' ## Predictive calibration check
#'
#' Examine how well the leave-one-out predictive probabilities
#' (computed with `loo_epred()`) are
#' calibrated using PAV-adjusted calibration plot
#' [@Dimitriadis+etal:2021:reliabilitydiag] implemented in
#' `reliabilitydiag`. Looks quite good.
#| label: fig-calib-0h
rd <- reliabilitydiag(EMOS = loo_epred(bfit_0h), y = dogs_df$shock)
plot_calib_0h <- autoplot(rd) +
        labs(
                x = "Predicted (LOO)",
                y = "Conditional event probabilities",
                title = "Model 0h"
        ) +
  bayesplot::theme_default(base_family = "sans", base_size = 16)
plot_calib_0h

#' ## Residual plots
#'
#' When can use PAV-adjustment to make also residual plots with
#' respect to covariates (details in forthcoming paper).
#| code-fold: true
ppc_pava_residual <-
  function(y,
           epred,
           x = NULL,
           prob = .9,
           n.boot = 1000,
           interval_geom = "ribbon",
           alpha = .2,
           ...) {
    require("Iso")
    require("ggplot2")
    assertthat::assert_that(is.numeric(y))
    assertthat::assert_that(is.numeric(epred))
    assertthat::are_equal(length(epred), length(y))
    # Compute predictive means and their ordering.
    yrep_bar_order <- order(epred)

    # If no x value supplied, plot residuals against predictive means.
      if (is.null(x)) {
        x <- epred
      }

      cep_df <- seq_len(n.boot) |>
        lapply(\(i) data.frame(cep = Iso::pava(rbinom(
          length(y), 1, epred
        )[yrep_bar_order]),
        id_ = yrep_bar_order)) |>
        dplyr::bind_rows() |>
        dplyr::group_by(id_) |>
        dplyr::summarise(upper = quantile(cep, .5 + .5 * prob),
                         lower = quantile(cep, .5 * (1 - prob)))
      cep_df$yrep_bar <- epred
      cep_df[yrep_bar_order, "cep"] <- Iso::pava(y[yrep_bar_order])
    cep_df$x <- x
    cep_df$ymax <- cep_df$upper - cep_df$yrep_bar
    cep_df$ymin <- cep_df$lower - cep_df$yrep_bar
    bw <- .5 * bw.SJ(cep_df$x)
    w <- sapply(cep_df$x, \(x_i) dnorm(cep_df$x, x_i, bw))
    cep_df$ymaxs <- (t(w) %*% cep_df$ymax) / colSums(w)
    cep_df$ymins <- (t(w) %*% cep_df$ymin) / colSums(w)

    ggplot(cep_df,
           aes(
             y = cep - yrep_bar,
             ymax = ymaxs,
             ymin = ymins,
             x = x
           )) +
      geom_hline(yintercept=0, alpha=0.3) +
      stat_identity(aes(
        colour = TRUE,
        fill = TRUE
      ),
      alpha = alpha,
      geom = interval_geom,
      ...) +
      geom_point(aes(colour = (cep >= lower) & (cep <= upper))) +
      scale_colour_discrete(aesthetics = c("fill" , "colour")) +
      theme(legend.position = "none") +
      labs(y = "PAVA Residual")
  }

#' PAV-adjusted residual plot looks reasonable.
#| label: fig-ppc_pava_residual-0h
ppc_pava_residual(dogs_df$shock,
                  loo_epred(bfit_0h),
                  jitter(dogs_df$time,0.3)) +
  labs(x = "Time")

#' ## Posterior predictive checking
#'
#' Visual posterior predictive checking plotting predicted shocks and
#' avoidances by ordering the dogs with last observed shock.
#' 
#| code-fold: true
pred_logit <- function(fit) {
  matrix(c(rep(1, 30),
           posterior_predict(fit, ndraws = 1) |>
             as.numeric()), nrow = 30, ncol = 25)
}
ppc_shocks <- function(shock, title) {
  expand.grid(dog = rev(1:30), time = 1:25) |>
    mutate(shock = as.numeric(shock[order(apply(shock, 1, \(x) max(which(x == 1)))), ])) |>
  ggplot(aes(time, dog, fill = shock)) +
  geom_tile() +
  ## coord_fixed() +
  scale_fill_gradient(low = "#ffffc8", high = "#7c0025") +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.6)) +
    labs(y = title)
}

#| label: fig-ppc_shocks-0h
#| fig-height: 3
#| fig-width: 8
ppc_shocks(shock, "Real dogs") +
  ppc_shocks(pred_logit(bfit_0h), "Model 0h: hier. logit")

#' Posterior predictive checking using mean number of switches between
#' shocks and avoidances as the test statistic.
mean_switches <- function(shock) {
  shock |> rowDiffs() |> abs() |> rowSums() |> mean()
}
yrep <- replicate(100, mean_switches(pred_logit(bfit_0h)))
#| label: fig-ppc-mean_switches-0h
ppc_stat(mean_switches(shock), matrix(yrep, nrow=length(yrep)), stat="identity")

#' ## Prior-likelihood sensitivity analysis
#'
#' Using `priorsense` package for powerscaling prior-likelihood
#' sensitivity analysis [@Kallioinen+etal:2023:priorsense] shows that
#' the data are informative and there is no need to think more about
#' priors unless we do happen to have easily available strong prior
#' information.
powerscale_sensitivity(bfit_0h, variable = variables(as_draws(bfit_0h))[1:6]) |>
        tt() |>
        format_tt(num_fmt = "decimal")

#' # Model 1: 1-parameter log model
#'
#' Instead of going straight to the 2-parameter log model by @Bush+Mosteller:1955,
#' we test one parameter model which by construction gives probability 1 at time $t=1$.
#' We assign a uniform prior ($\matrhm{beta}(1,1) is uniform from $0$ to $1$) on $a$.
#' $$
#' \begin{aligned}
#' \Pr(\mathrm{shock}) & = a^{(t - 1)}\\
#' a & \sim \mathrm{uniform}(0,1).
#' \end{aligned}
#' $$
#| label: bfit_1
#| results: "hide"
bfit_1 <- brm(bf(shock ~ a^(time - 1), a ~ 1, nl = TRUE),
              family = bernoulli(link = "identity"),
              prior = prior(beta(1, 1), nlpar = "a", lb = 0, ub = 1),
              data = dogs_df, refresh = 0)
bfit_1 <- add_criterion(bfit_1, criterion = "loo")

#' ## Model comparison
#'
#' PSIS-LOO-CV comparison shows that the 1-parameter log model is worse
#' than either logistic regression model. Thus we don't examine it
#' further and move to more elaborate log models.
loo_compare(bfit_0, bfit_0h, bfit_1) |>
        as.data.frame() |>
        tibble::rownames_to_column("model") |>
        select(model, elpd_diff, se_diff) |>
        tt()

#' # Model 2: 2-parameter log model
#'
#' This is the original 2-parameter model proposed by @Bush+Mosteller:1955:
#' $$
#' \begin{aligned}
#' \Pr(\mathrm{shock}) & = a^{x_{1jt}}\,b^{x_{2jt}}\\
#' a,b & \sim \mathrm{uniform}(0,1),
#' \end{aligned}
#' $$
#' where $x_{1jt}$ and $x_{1jt}$ are the number of previous shocks and
#' avoidances, respectively, in trials $1,\ldots,t-1$ for dog $j$.
#| label: bfit_2
#| results: "hide"
bfit_2 <- brm(bf(shock ~ a^prev_shock * b^prev_avoid,
                a ~ 1, b ~ 1, nl = TRUE),
              family = bernoulli(link = "identity"),
              prior = c(prior(beta(1, 1), nlpar = "a", lb = 0, ub = 1),
                        prior(beta(1, 1), nlpar = "b", lb = 0, ub = 1)),
              data = dogs_df, refresh=0)
bfit_2 <- add_criterion(bfit_2, criterion = "loo")

#' ## Model comparison
#'
#' PSIS-LOO-CV comparison shows that the 2-parameter log model is worse
#' than the hierarchical logistic regression model, but not
#' significantly. We don't examine this model further, as we can make
#' the comparison more fair by using hierarchical 2-parameter log
#' model.
loo_compare(bfit_0h, bfit_2) |>
        as.data.frame() |>
        tibble::rownames_to_column("model") |>
        select(model, elpd_diff, se_diff) |>
        tt()

#' # Model 3: hierarchical 1-parameter log model
#'
#' We could go directly to hierarchical 2-parameter log model, but for
#' completeness compared to Andrew's analysis, we include hierarchical
#' 1-parameter log model, too. Each dog has now it's own parameter $a_j$
#' with a normal hierarchical prior on $\mathrm{logit}(a_j)$.
#' $$
#' \begin{aligned}
#' \Pr(\mathrm{shock}) & = a_j^{(t - 1)}\\
#' \mathrm{logit}(a_j) & \sim \mathrm{normal}(\mu_{\mathrm{logit}(a)}, \sigma_{\mathrm{logit}(a)})\\
#' \mu_{\mathrm{logit}(a)} & \sim \mathrm{Student}_3(0, 2.5)\\
#' \sigma_{\mathrm{logit}(a)} & \sim \mathrm{Student}^{+}_3(0, 2.5)
#' \end{aligned}
#' $$
#' where $a_j$ is parameter for dog $j$.
#| label: bfit_3
#| results: "hide"
inv_logit <- function(x) 1 / (1 + exp(-x))
bfit_3 <- brm(bf(shock ~ inv_logit(etaa)^(time - 1),
                 etaa ~ (1 | dog), nl = TRUE),
              family = bernoulli(link = "identity"),
              prior = prior(student_t(3, 0, 2.5), nlpar = "etaa"),
              data = dogs_df, refresh = 0)
bfit_3 <- add_criterion(bfit_3, criterion = "loo")

#' ## Model comparison
#' 
#' PSIS-LOO-CV comparison shows that the hierarchical 1-parameter log model
#' is worse than the hierarchcial logistic regression model and
#' non-hierarchcial 2-parameter log model, and thus we don't examine
#' this model further.
loo_compare(bfit_0h, bfit_2, bfit_3) |>
        as.data.frame() |>
        tibble::rownames_to_column("model") |>
        select(model, elpd_diff, se_diff) |>
        tt()

#' # Model 4: hierarchical 2-parameter log model
#'
#' Each dog has now it's own parameters $a_j$ and $b_j$ with a
#' multivariate normal hierarchical prior.
#' $$
#' \begin{aligned}
#' \Pr(\mathrm{shock}) & = a_j^{x_{1jt}}\,b_j^{x_{2jt}}\\
#' \left(\begin{array}{c}\mathrm{logit}(a)_j \\ \mathrm{logit}(b)_j\end{array} \right) & \sim \mathrm{MVN}(\mu_{\mathrm{logit}(a),\mathrm{logit}(b)}, \Sigma_{\mathrm{logit}(a),\mathrm{logit}(b)})\\
#' \mu_{\mathrm{logit}(a)}, \mu_{\mathrm{logit}(b)} & \sim \mathrm{Student}_3(0, 2.5)\\
#' \Sigma_{\mathrm{logit}(a),\mathrm{logit}(b)} & = \left(\begin{array}{cc}\sigma_\mathrm{logit}(a) & 0 \\ 0 & \sigma_\mathrm{logit}(b)\end{array}\right) Q_{\mathrm{logit}(a),\mathrm{logit}(b)} \left(\begin{array}{cc}\sigma_\mathrm{logit}(a) & 0 \\ 0 & \sigma_\mathrm{logit}(b)\end{array}\right) \\
#' \sigma_{\mathrm{logit}(a)}, \sigma_{\mathrm{logit}(b)} & \sim \mathrm{Student}^{+}_3(0, 2.5)\\
#' Q_{\mathrm{logit}(a),\mathrm{logit}(b)} & \sim \mathrm{LKJ}(1).
#' \end{aligned}
#' $$
#| label: bfit_4
#| results: "hide"
bfit_4 <- brm(bf(shock ~ inv_logit(etaa)^prev_shock * inv_logit(etab)^prev_avoid,
                 mvbind(etaa, etab) ~ (1 |p| dog), nl=TRUE),
              family = bernoulli(link = "identity"),
              prior = c(prior(student_t(3, 0, 2.5), nlpar = "etaa"),
                        prior(student_t(3, 0, 2.5), nlpar = "etab")),
              data = dogs_df, refresh = 0)
bfit_4 <- add_criterion(bfit_4, criterion = "loo", save_psis = TRUE)

#' ## Model comparison
#' 
#' PSIS-LOO-CV comparison shows that hierarchical 2-parameter log model is
#' worse than the hierarchical logistic regression, although not
#' significantly. While adding hierarchy to logistic regression
#' improved the predictive performance significantly, adding hierarchy
#' to 2-parameter log model has a very small effect. This is probably
#' due to the fact that the 2-parameter log model was already able to
#' take into account the variation in the shock and avoidances by
#' using them as covariates.
loo_compare(bfit_0h, bfit_2, bfit_4) |>
        as.data.frame() |>
        tibble::rownames_to_column("model") |>
        select(model, elpd_diff, se_diff) |>
        tt()

#' ## Visualize predictions
#' 
#' Visualize the model fit for 9 first dogs. These look different from
#' the logistic regression. Specifically we tend to see a sharper drop
#' after the first avoidance, which makes sense as the magnitude of
#' drop in probability after repeated shocks diminishes, but the first
#' avoidance provides another big drop.
#| label: fig-pred9-4
plot_pred9_4 <-dogs_df |>
  dplyr::filter(dog <= 9) |>
        add_linpred_draws(bfit_4, transform = TRUE) |>
        ggplot(aes(x = time, y = .linpred)) +
        stat_lineribbon(.width = c(.95), alpha = 1 / 2, 
                        color = brewer.pal(5, "Blues")[[5]]) +
  scale_fill_brewer() +
  geom_point(data = dplyr::filter(dogs_df, dog <= 9), 
             aes(x = time, y = shock, group = dog)) +
  facet_wrap( ~ dog, labeller = label_both) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c(0, 1)) +
  labs(y = "Shocks and predicted probability of shock", title = "Model 4")
plot_pred9_4
  
#' We can compare these two the predictions from the non-hierarchial
#' 2-parameter log model, and we can see that even if the model is
#' non-hierarchical, dog-specific number of shocks and avoidances do
#' make the model fit to vary by dog. This can explain why adding
#' hierarchy to the model does not improve the predictive performance.
#| label: fig-pred9-2
dogs_df |>
  dplyr::filter(dog <= 9) |>
        add_linpred_draws(bfit_2, transform = TRUE) |>
        ggplot(aes(x = time, y = .linpred)) +
        stat_lineribbon(.width = c(.95), alpha = 1 / 2, 
                        color = brewer.pal(5, "Blues")[[5]]) +
  scale_fill_brewer() +
  geom_point(data = dplyr::filter(dogs_df, dog<=9), 
             aes(x = time, y = shock, group = dog)) +
  facet_wrap( ~ dog, labeller = label_both) +
  theme(legend.position = "none") +
  labs(y = "Shocks and predicted probability of shock")


#' ## Predictive calibration check
#' 
#' Examine how well the leave-one-out predictive probabilities from
#' hierarchical 2-parameter log model are calibrated. Looks quite good.
#| label: fig-calib-4
rd <- reliabilitydiag(EMOS = loo_epred(bfit_4), y = dogs_df$shock)
plot_calib_4 <- autoplot(rd) +
  labs(x = "Predicted (LOO)", 
       y = "Conditional event probabilities", 
       title = "Model 4")+
  bayesplot::theme_default(base_family = "sans", base_size = 16)
plot_calib_4

#' ## Residual plots
#' 
#' PAV-adjusted residual plot looks reasonable.
#| label: fig-ppc_pava_residual-4
ppc_pava_residual(dogs_df$shock,
                  loo_epred(bfit_4),
                  jitter(dogs_df$time,0.3)) +
  labs(x="Time")

#' ## Posterior predictive checking
#'
#' Visual posterior predictive checking
#'
#| code-fold: true
pred_log <- function(fit) {
  pred_shock <- matrix(rep(1, 30), nrow = 30)
  for (t in 2:25) {
    dogs_df_pred <- dplyr::filter(dogs_df, time == t)
    dogs_df_pred$prev_pred_shock <- as.numeric(rowSums(pred_shock))
    dogs_df_pred$prev_avoid <- as.numeric(rowSums(1 - pred_shock))
    pred <- posterior_predict(fit, ndraws = 1, newdata = dogs_df_pred)
    pred_shock <- cbind(pred_shock, as.numeric(pred))
  }
  pred_shock
}
#| label: fig-ppc_shocks-4
#| fig-height: 4
#| fig-width: 5
ppc_shocks(pred_log(bfit_4), "PPsims from M4:\nhier logit model") 

#' Posterior predictive checking using mean number of switches between
#' shocks and avoidances as the test statistic.
#| label: fig-ppc_mean_switches-4
yrep <- replicate(100, mean_switches(pred_log(bfit_4)))
ppc_stat(mean_switches(shock), matrix(yrep, nrow=length(yrep)), stat="identity")

#' ## Prior-likelihood sensitivity analysis
#'
#' Using `priorsense` package for prior-likelihood sensitivity
#' analysis shows that the data are informative and there is no need
#' to think more about priors unless we do happen to have easily
#' available strong prior information.
powerscale_sensitivity(bfit_4, variable = variables(as_draws(bfit_4))[1:5]) |>
  tt() |>
  format_tt(num_fmt="decimal")


#' # Comparing posterior predictions
#'
#' We can compare the posterior predictions by overlaying them in the
#' same plot. We see the biggest difference is visible in time $t=2$,
#' but the differences are that small that given the small size of the
#' data, there is no clear preference for either model.
#| label: fig-pred-0h-4
dogs_df |>
  mutate(epred_0h = colMeans(posterior_epred(bfit_0h)),
         epred_4 = colMeans(posterior_epred(bfit_4))) |>
  ggplot(aes(x = time, group = dog)) +
  geom_line(aes(y = epred_0h, color = "Model 0h"), alpha = 0.5) +
  geom_line(aes(y = epred_4, color = "Model 4"), alpha = 0.5) +
  scale_y_continuous(lim = c(0, 1)) +
  labs(x="Time",y="Posterior predictive probability")

#' # Leave-future-out cross-validation
#'
#' Above we used LOO-CV which leaves out only one observation, but
#' maybe we should examine the predictive performance only for the
#' future observations. As the model fits are quite fast, we can do
#' brute force leave.future-out cross-validation. We start with
#' fitting the model with $t<5$ and predict the outcomes at $t=5$, and
#' then fit with one more time point and repeat.
#'
#| results: "hide"
ll_0h <- matrix(nrow=4000,ncol=0)
for (t in 5:25) {
  ll_0h <- cbind(ll_0h, rowSums(log_lik(update(bfit_0h, newdata = dplyr::filter(dogs_df, time < t)),
                                        newdata = dplyr::filter(dogs_df, time == t))))
}
ll_4 <- matrix(nrow=4000,ncol=0)
for (t in 5:25) {
  ll_4 <- cbind(ll_4, rowSums(log_lik(update(bfit_4, newdata = dplyr::filter(dogs_df, time < t)),
                                      newdata = dplyr::filter(dogs_df, time == t))))
}

#' There is no difference in predictive performance between the
#' hierarchical logistic regression and hierarchical log model.
loo_compare(list(bfit_0h=elpd(ll_0h),bfit_4=elpd(ll_4))) |>
        as.data.frame() |>
        tibble::rownames_to_column("model") |>
        select(model, elpd_diff, se_diff) |>
        tt()

#' # Are the data informative on two parameters of 2-parameter log model?
#'
#' One assumption about the 2-parameter log model is that if the
#' posteriors of $a$ and $b$ are different, then the dogs learn a
#' different amount from shocks and avoidances.
#| label: fig-mcmc_areas-4
as_draws_df(bfit_4) |>
  mutate(a = inv_logit(b_etaa_Intercept),
         b = inv_logit(b_etab_Intercept)) |>
  mcmc_areas(pars = c("a", "b"))

#' The posteriors are clearly different, but is this due to different
#' learning rate? We can test this by generating simulated data from
#' the logistic regression model, which is noty making assumption
#' about different learning rates from shocks and avoidances.
#'
#' Generate shocks and avoidances using the simple logistic regression
#| label: bfit_4s
#| results: "hide"
dogs_df_pred_0 <- dogs_df
pred_shock_0 <- matrix(c(rep(1, 30), posterior_predict(bfit_0, ndraws = 1) |> as.numeric()),
                       nrow = 30, ncol = 25)
dogs_df_pred_0$prev_shock_0 <- as.numeric(rowCumsums(pred_shock_0)[, 1:(ncol(pred_shock_0) - 1)])
dogs_df_pred_0$prev_avoid <- as.numeric(rowCumsums(1-pred_shock_0)[, 1:(ncol(pred_shock_0) - 1)])
bfit_4s <- brm(bf(shock ~ inv_logit(etaa)^prev_shock * inv_logit(etab)^prev_avoid,
                 mvbind(etaa, etab) ~ (1 |p| dog), nl = TRUE),
              family = bernoulli(link = "identity"),
              data = dogs_df_pred_0, refresh=0)
#' The posterior for a and b are different!
#| label: fig-mcmc_areas-4s
as_draws_df(bfit_4s) |>
  mutate(a = inv_logit(b_etaa_Intercept),
         b = inv_logit(b_etab_Intercept)) |>
  mcmc_areas(pars = c("a", "b"))

#' Generate shocks and avoidances using the hierarchical logistic regression
#| label: bfit_4sh
#| results: "hide"
dogs_df_pred_0h <- dogs_df
pred_shock_0h <- matrix(c(rep(1, 30), posterior_predict(bfit_0h, ndraws = 1) |> as.numeric()),
                       nrow = 30, ncol = 25)
dogs_df_pred_0h$prev_shock_0h <- as.numeric(rowCumsums(pred_shock_0h)[, 1:(ncol(pred_shock_0h) - 1)])
dogs_df_pred_0h$prev_avoid <- as.numeric(rowCumsums(1-pred_shock_0h)[, 1:(ncol(pred_shock_0h) - 1)])
bfit_4sh <- brm(bf(shock ~ inv_logit(etaa)^prev_shock * inv_logit(etab)^prev_avoid,
                mvbind(etaa, etab) ~ (1 |p| dog), nl=TRUE),
                family = bernoulli(link = "identity"),
                data = dogs_df_pred_0h, refresh = 0)

#' The posterior for a and b are different!
#| label: fig-mcmc_areas-4sh
as_draws_df(bfit_4sh) |>
  mutate(a = inv_logit(b_etaa_Intercept),
         b = inv_logit(b_etab_Intercept)) |>
  mcmc_areas(pars = c("a", "b")) +
  labs(title = "Model 4: hier. logit, simulated data from Model 0h")

#' # Conclusion
#'
#' Based on model comparison and various model checking diagnostics,
#' data are not informative to make difference between hierarchical
#' logistic regression model, 2-parameter log model, and hierarchical
#' 2-parameter log model. Although the visual predictive checking in
#' dogs example is historically interesting, it seems that it is not
#' sufficient for making difference between some plausible models.
#'
#'
#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' # Licenses {.unnumbered}
#'
#' * Code &copy; 2023--2025, Aki Vehtari and Andrew Gelman, licensed under BSD-3.
#' * Text &copy; 2023--2025, Aki Vehtari and Andrew Gelman, licensed under CC-BY-NC 4.0.

#' ---
#' title: "Student grades and variable selection with projpred"
#' author: "Aki Vehtari"
#' date: 2023-12-14
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

#| label: setup
#| include: FALSE
knitr::opts_chunk$set(
  message = FALSE,
  error = FALSE,
  warning = FALSE,
  comment = NA,
  out.width = "90%",
  cache = TRUE
)

#' ## Portugal students data
#'
#' We work with an example of predicting mathematics and Portuguese
#' exam grades for a sample of high school students in Portugal. The
#' same data was used in Chapter 12 of Regression and Other Stories
#' book to illustrate different models for regression coefficients .
#'
#' We predict the students' final-year median exam grade in
#' mathematics (n=407) and Portuguese (n=657) given a large number of
#' potentially relevant predictors: student's school, student's sex,
#' student's age, student's home address type, family size, parents'
#' cohabitation status, mother's education, father's education,
#' home-to-school travel time, weekly study time, number of past class
#' failures, extra educational support, extra paid classes within the
#' course subject, extra-curricular activities, whether the student
#' attended nursery school, whether the student wants to take higher
#' education, internet access at home, whether the student has a
#' romantic relationship, quality of family relationships, free time
#' after school, going out with friends, weekday alcohol consumption,
#' weekend alcohol consumption, current health status, and number of
#' school absences.
#'
#' ## Variable selection
#'
#' If we would care only about the predictive performance, we would
#' not need to do variable selection, but we would use all the
#' variables and a sensible joint prior. Here we are interested in
#' finding the smallest set of variables that provide similar
#' predictive performance as using all the variables (and sensible
#' prior). This helps to improve explainability and to design further
#' studies that could include also interventions. We are not
#' considering causal structure, and the selected variables are
#' unlikely to have direct causal effect, but the selected variables
#' that have high predictive relevance are such that their role in
#' causal graph should be eventually considered.
#'
#' We first build models with all predictors, and then use projection
#' predictive variable selection [@Piironen+etal:projpred:2020;
#' @McLatchie+etal:2025:projpred_workflow] implemented in R package
#' [`projpred`](https://cran.r-project.org/package=projpred). We also
#' demonstrate use of subsampling LOO with difference estimator
#' [@Magnusson+etal:2020:bigloocomparison] to speed-up model size
#' selection.
#'

#' **Load packages**
#| cache: FALSE
library(brms)
library(cmdstanr)
options(brms.backend = "cmdstanr")
options(mc.cores = 4)
library(posterior)
options(digits = 2, posterior.digits = 2, 
        pillar.neg = FALSE, pillar.subtle = FALSE, pillar.sigfig = 2)
library(loo)
library(projpred)
if (interactive()) {
  library(progressr)
  options(projpred.use_progressr = TRUE)
  handlers(global = TRUE)
}
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 16))
library(RColorBrewer)
set1 <- brewer.pal(7, "Set1")
library(khroma)
library(latex2exp)
library(tinytable)
options(tinytable_format_num_fmt = "significant_cell", tinytable_format_digits = 2,
        tinytable_tt_digits = 2)
library(dplyr)
library(matrixStats)
library(patchwork)
library(ggdist)
library(doFuture)
library(doRNG)
library(tictoc)


#' ## Data preparation and visualization
#'
#' Get the data from Regression and Other Stories R package.
student <- read.csv(url('https://raw.githubusercontent.com/avehtari/ROS-Examples/master/Student/data/student-merged-all.csv'))
#' List the predictors to be used.
predictors <- c("school","sex","age","address","famsize","Pstatus", "Medu","Fedu",
                "traveltime","studytime","failures","schoolsup", "famsup","paid",
                "activities", "nursery", "higher", "internet", "romantic","famrel",
                "freetime","goout","Dalc","Walc","health", "absences")
p <- length(predictors)

#' The data includes 3 grades for both mathematics and Portuguese.  To
#' reduce the variability in the outcome we use median grades based on
#' those three exams for each topic. We select only students with
#' non-zero grades.
grades <- c("G1mat", "G2mat", "G3mat", "G1por", "G2por", "G3por")
student <- student %>%
  mutate(across(matches("G[1-3]..."), ~na_if(.,0))) %>%
  mutate(Gmat = rowMedians(as.matrix(select(., matches("G.mat"))), na.rm = TRUE), 
         Gpor = rowMedians(as.matrix(select(., matches("G.por"))), na.rm = TRUE))
student_Gmat <- subset(student, is.finite(Gmat), select = c("Gmat", predictors))
student_Gmat <- student_Gmat[is.finite(rowMeans(student_Gmat)),]
student_Gpor <- subset(student, is.finite(Gpor), select = c("Gpor", predictors))
(nmat <- nrow(student_Gmat))
head(student_Gmat) |> tt()
(npor <- nrow(student_Gpor))
head(student_Gpor) |> tt()

#' The following plot shows the distributions of median math and
#' Portuguese exam scores for each student.
#| label: fig-data-histograms
#| fig-height: 4
#| fig-width: 8
p1 <- ggplot(student_Gmat, aes(x = Gmat)) +
  geom_dots() + 
  labs(x = "Median math exam score")
p2 <- ggplot(student_Gpor, aes(x = Gpor)) +
  geom_dots() + 
  labs(x = "Median Portuguese exam score")
(p1 + p2) * scale_x_continuous(lim = c(0, 20)) * 
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

#' We standardize all to have standard deviation 1, to make the
#' comparison of relevances easier and also to make the definition of
#' priors easier.
studentstd_Gmat <- student_Gmat
studentstd_Gmat[, predictors] <- scale(student_Gmat[, predictors])
studentstd_Gpor <- student_Gpor
studentstd_Gpor[, predictors] <- scale(student_Gpor[, predictors])

#' ## Default uniform prior on coefficients
#'
#' Before variable selection, we want to build a good model with all
#' covariates. We first illustrate that default priors may be
#' bad when we have many predictors.
#'
#' By default `brms` uses uniform ("flat") priors on regression coefficients.
#| results: hide
#| cache: true
fitm_u <- brm(Gmat ~ ., data = studentstd_Gmat)

#' If we compare posterior-$R^2$ (`bayes_R2()`) and LOO-$R^2$
#' (`loo_R2()`) [@Gelman+etal:2019:BayesR2], we see that the
#' posterior-$R^2$ is much higher, which means that the posterior
#' estimate for the residual variance is strongly underestimated and
#' the model has overfitted the data.
#'
#' We first re-define `bayes_R2()` to match the paper.
bayes_R2.brmsfit <- function(fit, summary=TRUE, probs = c(0.025, 0.975)) {
  mupred <- brms::posterior_epred(fit)
  var_mupred <- apply(mupred, 1, var)
  sigma2 <- as.matrix(fit, variable = c("sigma"))^2
  R2 <- var_mupred / (var_mupred + sigma2)
  colnames(R2) <- "R2"
  if (summary) {
    R2 <- posterior_summary(R2, probs = probs)
  }
  R2
}
bayes_R2(fitm_u) |> as.data.frame() |> tt()
loo_R2(fitm_u) |> as.data.frame() |> tt()

#' We plot the marginal posteriors for coefficients. For many
#' coefficients the posterior is quite wide.
#| label: fig-fitmu-mcmc_areas
#| fig-height: 6
#| fig-width: 8
#| cache: false
drawsmu <- as_draws_df(fitm_u, variable = paste0("b_",predictors)) |>
  set_variables(predictors)
p <- mcmc_areas(drawsmu, prob_outer = 0.98, area_method = "scaled height") +
  xlim(c(-1.5,1.5))
p <- p + scale_y_discrete(limits = rev(levels(p$data$parameter)))
p

#' ## Piranha theorem
#' 
#' The common proper prior choice for coefficients is independent wide
#' normal prior. Piranha theorem states that it is not possible that
#' all predictors would be independent and have large coefficients at
#' the same time [@Tosh+etal:2024:piranha]. Thus when we have many
#' predictors we should include prior information that not all
#' coefficients can be large. One option is simply to divide the
#' normal prior variance with the number of predictors which keeps the
#' total variance constant (assuming the predictors are
#' normalized). Other option is to use more elaborate joint priors on
#' coefficients.
#' 

#' ## Implied priors on $R^2$ and R2D2 prior
#'
#' In regression analysis cases we might assume that data is noisy and
#' it is unlikely that we would get almost perfect fit. We can measure
#' the proportion of variance explained by the model with $R^2$.  To
#' understand what is the implied prior on $R^2$ given different
#' priors on coefficients, we can simply sample from the prior and
#' compute Bayesian $R^2$ using the prior draws.  The Bayesian-$R^2$
#' depends only on the model parameters, and thus can be used without
#' computing residuals that depend on data.
#'
#' If we have some prior information about $R^2$ we can use R2D2 prior
#' [@Zhang+etal:2022:R2D2] to first define a prior directly on $R^2$,
#' and then the prior is propagated to coefficients so that despite
#' the number of predictors the prior on $R^2$ stays constant. As
#' $R^2$ depends also on the residual scale, R2D2 prior is a joint
#' prior on coefficients and residual sigma.
#'
#' Although we can fit models with uniform prior on coefficients and
#' get proper posterior, we can't sample from improper unbounded
#' uniform prior. Thus, in the following we consider only models with
#' proper priors. For all the following models we use normal+(0, 3)
#' prior (where normal+ indicates normal distribution constrained to
#' be positive) for residual scale sigma, which is very weak prior as
#' the standard deviation of the whole data is 3.3. We use four
#' different priors for coefficients:
#'
#' 1. Independent normal(0, 2.5) prior which is a proper prior used as default by
#'    `rstanarm` and considered to be weakly informative for a single
#'    coefficient.
#' 2. Independent scaled normal prior. If we assume that many
#'    predictors may each have small relevance, we can scale the
#'    independent priors so that the sum of the prior variance stays
#'    reasonable. In this case we have 26 predictors and could have a
#'    prior guess that proportion of the explained variance is near
#'    $0.3$. Then a simple approach would be to assign independent
#'    priors to regression coefficients with mean $0$ and standard
#'    deviation $\sqrt{0.3/26}\operatorname{sd}(y)$.
#' 3. Regularized horseshoe prior [@Piironen+Vehtari:2017:rhs] which
#'    is a joint prior for the coefficients depending also on the
#'    residual scale, and it can be used to define sparsity assuming
#'    prior with the expected number of relevant coefficients. Here we
#'    guess that maybe 6 coefficients are relevant and set the global
#'    scale according to that. The result is not sensitive to the
#'    exact value for the prior guess of the number of relevant
#'    coefficients as it states just the mean for the prior. Regularized
#'    horseshoe has good prior predictive behavior when more variables
#'    are added.
#' 4. R2D2 prior, which has the benefit that it first defines the
#'    prior directly on $R^2$, and then the prior is propagated to the
#'    coefficients. The R2D2 prior is predictively consistent, so that
#'    the prior on $R^2$ stays constant as the number of predictors
#'    increases. We assign the R2D2 prior with mean 1/3 and precision
#'    3, which corresponds to the $\Beta(1,2)$ distribution on $R^2$
#'    implying that higher $R^2$ values are less likely. We set the
#'    concentration parameter to 1/2, which implies we assume some of
#'    the coefficients can be big and some small. R2D2 prior implementation
#'    in `brms` assumes the predictors have been standardized to have
#'    unit variance.
#'     
#| label: fitm
#| results: hide
#| cache: true
# we sample from both posterior and prior
# normal(0, 2.5)
fitm_n1 <- brm(Gmat ~ ., data = studentstd_Gmat,
               prior=c(prior(normal(0, 2.5), class = b)),
               warmup = 1000, iter = 5000,
               refresh = 0)
fitm_n1p <- update(fitm_n1, sample_prior = "only")
# we also sample from a truncated prior to get more draws in the
# interesting region of R^2<0.5 needed for zoomed plot
fitm_n1pt <- update(fitm_n1, sample_prior = "only",
                    prior = c(prior(normal(0, 2.5), class = b),
                              prior(student_t(3, 0, 3), lb=5, class = sigma)),
                    refresh = 0)

# normal(0, sqrt(0.3/26)*sd(y))
scale_b <- sqrt(0.3/26) * sd(studentstd_Gmat$Gmat)
fitm_n2 <- brm(Gmat ~ ., data = studentstd_Gmat,
               prior=c(prior(normal(0, scale_b), class = b)),
               warmup = 1000, iter = 5000,
               stanvars = stanvar(scale_b, name = "scale_b"), 
               refresh = 0)
fitm_n2p <- update(fitm_n2, sample_prior = "only")

# horseshoe
p <- length(predictors)
p0 <- 6
scale_slab <- sd(studentstd_Gmat$Gmat)/sqrt(p0)*sqrt(0.3)
scale_global <- p0 / (p - p0) / sqrt(nrow(studentstd_Gmat))
fitm_hs <- brm(Gmat ~ ., data = studentstd_Gmat,
               prior = c(prior(horseshoe(scale_global = scale_global,
                                         scale_slab = scale_slab), class = b)),
               warmup = 1000, iter = 5000,
               refresh = 0)
fitm_hsp <- update(fitm_hs, sample_prior = "only")

# R2D2
fitm <- brm(Gmat ~ ., data = studentstd_Gmat,
            prior = c(prior(R2D2(mean_R2 = 1/3, prec_R2 = 3, cons_D2 = 1/2), class = b)), 
            warmup = 1000, iter = 5000,
            refresh = 0)
fitmp <- update(fitm, sample_prior = "only")

#| label: prior-plots
#| cache: true
# plot prior on R^2 in (0,1)
types <- factor(c("Prior", "Posterior"), levels = c("Prior", "Posterior"))
types <- rep(types, times = 4)
priornames <- factor(c("Wide normal","Scaled normal","RHS","R2D2"), 
                     levels = c("Wide normal", "Scaled normal", "RHS", "R2D2"))
priornames <- rep(priornames, each = 2)
clr <- colour("bright", names = FALSE)(7)
fits <- list(fitm_n1p, fitm_n1, fitm_n2p, fitm_n2, fitm_hsp, fitm_hs, fitmp, fitm)
pp1 <- lapply(1:8, \(i) data.frame(
                R2 = as.numeric(bayes_R2(fits[[i]], summary = FALSE)),
                type = types[i],
                priorname = priornames[i]
              )) |> 
  bind_rows() |>
  filter(type=="Prior") |>
  ggplot(aes(x=R2, color=priorname)) +
  stat_slab(expand = TRUE, trim = FALSE, alpha = .6, fill = NA, adjust = 2) +
  xlim(c(0,1)) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 16), 
        legend.position = "none") +
  scale_color_bright() +
  labs(x = TeX("$R^2$"), title = TeX("a) Implied prior on $R^2$")) +
  annotate(geom = "text", x = 0.94, y = 0.33, label = "Wide normal", 
           hjust = 1, color = clr[1], size = 5) + 
  annotate(geom = "text", x = 0.85, y = 0.12, label = "Scaled normal", 
           hjust = 1, color = clr[2], size = 5) +
  annotate(geom = "text", x = 0.06, y = 0.33, label = "Horseshoe", 
           hjust = 0, color = clr[3], size = 5) +
  annotate(geom = "text", x = 0.15, y = 0.12, label = "R2D2", 
           hjust = 0, color = clr[4], size = 5) 

# plot prior on R^2 in narrower range
fits <- list(fitm_n1pt, fitm_n1, fitm_n2p, fitm_n2, fitm_hsp, fitm_hs, fitmp, fitm)
pp2 <- lapply(1:8, \(i) data.frame(
                R2 = as.numeric(bayes_R2(fits[[i]], summary = FALSE)),
                type = types[i],
                priorname = priornames[i]
              )) |> 
  bind_rows() |>
  filter(type == "Prior") |> 
  ggplot(aes(x = R2, color = priorname)) +
  stat_slab(expand = TRUE, trim = FALSE, alpha = 0.5, fill = NA, adjust = 2) +
  xlim(c(0.04,.419)) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 16), 
        legend.position = "none") +
  scale_color_bright() +
  labs(x = TeX("$R^2$"), title = TeX("b) Implied prior on $R^2$")) +
  annotate(geom = "text", x = 0.415, y = 0.61, label = "Wide normal", 
           hjust = 1, color = clr[1], size = 5) + 
  annotate(geom = "text", x = 0.415, y = 0.34, label = "Scaled normal", 
           hjust = 1, color = clr[2], size = 5) +
  annotate(geom = "text", x = 0.06, y = 0.87, label = "Horseshoe", 
           hjust = 0, color = clr[3], size = 5) +
  annotate(geom = "text", x = 0.06, y = 0.28, label = "R2D2", 
           hjust = 0, color = clr[4], size = 5) 

# plot R^2 posterior in narrower range
pp3 <- lapply(1:8, \(i) data.frame(
                R2 = as.numeric(bayes_R2(fits[[i]], summary = FALSE)),
                type = types[i],
                priorname = priornames[i]
              )) |> 
  bind_rows() |>
  filter(type == "Posterior") |> 
  ggplot(aes(x = R2, color = priorname)) +
  stat_slab(expand = TRUE, trim = FALSE, alpha = 0.5, fill = NA, adjust = 2) +
  xlim(c(0.04,.42)) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 16), 
        legend.position = "none") +
  scale_color_bright() +
  labs(x = TeX("$R^2$"), title = TeX("c) Posterior-$R^2$")) +
  annotate(geom = "text", x = 0.28, y = 0.95, label = "Wide normal", 
           hjust = 0, color = clr[1], size = 5) + 
  annotate(geom = "text", x = 0.26, y = 0.95, label = "Scaled normal", 
           hjust = 1, color = clr[2], size = 5) +
  annotate(geom = "text", x = 0.19, y = 0.73, label = "Horseshoe", 
           hjust = 1, color = clr[3], size = 5) +
  annotate(geom = "text", x = 0.21, y = 0.84, label = "R2D2", 
           hjust = 1, color = clr[4], size = 5) 

# plot LOO-R^2 in narrower range
looR2 <- lapply(1:8, \(i) data.frame(
                  R2 = as.numeric(loo_R2(fits[[i]], summary = FALSE)),
                  type = types[i],
                  priorname = priornames[i]
                )) |> bind_rows()
pp4 <- looR2 |>
  filter(type == "Posterior") |> 
  ggplot(aes(x = R2, color = priorname)) +
  stat_slab(expand = TRUE, trim = FALSE, alpha = 0.5, fill = NA, adjust = 2) +
  xlim(c(0.04,.42)) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 16), 
        legend.position = "none") +
  scale_color_bright() +
  labs(x = TeX("$R^2$"), title = TeX("d) LOO-$R^2$")) +
  annotate(geom = "text", x = 0.175, y = 0.72, label = "Wide normal", 
           hjust = 1, color = clr[1], size = 5) + 
  annotate(geom = "text", x = 0.235, y = 0.72, label = "Scaled normal", 
           hjust = 0, color = clr[2], size = 5) +
  annotate(geom = "text", x = 0.175, y = 0.9, label = "Horseshoe", 
           hjust = 1, color = clr[3], size = 5) +
  annotate(geom = "text", x = 0.23, y = 0.9, label = "R2D2", 
           hjust = 0, color = clr[4], size = 5) 

#' ::: {.content-visible when-format="html"}
#| label: fig-implied-R2-prior-posterior-loo-html
#| cache: false
pp1
pp2
pp3
pp4
#' :::
#'
#' ::: {.content-visible unless-format="html"}
#| label: fig-implied-R2-prior-posterior-loo
#| fig-height: 4
#| fig-width: 13
#| cache: false
pp1 + pp2 + pp3 + pp4 + plot_layout(ncol = 4)
#' :::

#' Subplot a) shows the implied priors on $R^2$. We see that
#' independent wide normal prior on coefficients leads to a prior
#' strongly favoring values near 1. Independent scaled normal and R2D2
#' priors imply relatively flat prior on $R^2$. Regularized horseshoe
#' prior implies prior favoring values near 0. Looking at the implied
#' priors on the whole range can be misleading, and what is more
#' important is what is the behavior of the priors where the
#' likelihood is not very close to 0. Subplot b) shows the implied
#' priors on $R^2$ in a shorter range, and subplots c) and d) show
#' posterior-$R^2$ and LOO-$R^2$, respectively in the same range. Wide
#' normal prior favors larger $R^2$ values which then pushes the
#' posterior towards higher $R^2$ values and makes the model ovefit
#' which causes LOO-$R^2$ to be lower than with priors that do not
#' favor higher $R^2$ values. Scaled normal, regularized horseshe and
#' R2D2 priors all slightly favor smaller $R^2$ values, which is a
#' sensible prior assumption, and the posterior have not been pushed
#' towards higher $R^2$ values. LOO-$R^2$ results are quite similar
#' with scaled normal, regularized horseshe and R2D2 priors, although
#' with regularized horseshe and R2D2 priors there is slightly less
#' uncertainty. Regularized horseshoe is favoring smaller $R^2$ values
#' more than scaled normal and R2D2, but as the likelihood has thinner
#' tail towards smaller $R^2$ values (is more informative in that
#' direction) there is not much difference in the posteriors.
#'
#' Although in this case scaled normal, regularized horseshoe, and
#' R2D2 all produce quite similar $R^2$ results, in general we favor
#' R2D2 prior as it is easiest to define the prior on
#' $R^2$. Regularized horseshoe can be easier to define when our prior
#' is about the sparsity in the coefficients. We now continue with the
#' R2D2 prior.
#'
#' In addition of LOO-$R^2$ shown above, we can compare LOO estimated
#' expected log predictive scores. With R2D2 prior the predictive
#' performance is better than with wide normal prior with probability
#' 0.94. The difference to scaled normal and horseshoe prior is
#' smaller.
loo_compare(list(
  `Wide normal` = loo(fitm_n1),
  `Scaled normal` = loo(fitm_n2),
  `Horseshoe` = loo(fitm_hs),
  `R2D2` = loo(fitm)
))

#' ## Marginal posteriors of coefficients
#' 
#' We plot the marginal posteriors for coefficients. For many
#' coefficients the posterior has been shrunk close to 0. Some
#' marginal posteriors are wide.
#| label: fitm-fewer-draws
#| results: hide
#| cache: true
# refit with fewer posterior draws
fitm <- update(fitm, iter = 2000)

#| label: fig-fitm-mcmc_areas
#| fig-height: 6
#| fig-width: 8
drawsm <- as_draws_df(fitm, variable = paste0('b_', predictors)) |> 
  set_variables(predictors)
p <- mcmc_areas(drawsm, prob_outer=0.98, area_method = "scaled height") +
  xlim(c(-1.5,1.5))
p <- p + scale_y_discrete(limits = rev(levels(p$data$parameter)))
p

#' We check the bivariate marginal for Fedu and Medu coefficients, and
#' see that while the univariate marginals overlap with 0, jointly
#' there is not much posterior mass near 0. This is due to Fedu and
#' Medu being collinear. Collinearity of predictors, make it difficult
#' to infer the predictor relevance from the marginal posteriors.
#| label: fig-fitm-mcmc_scatter-Fedu-Medu
#| fig-height: 3
#| fig-width: 3
mcmc_scatter(drawsm, pars = c("Fedu","Medu"), size = 1, alpha = 0.1) +
  vline_0(linetype = "dashed") +
  hline_0(linetype = "dashed")

#' ## Model checking
#'
#' We're using a normal observation model, although we know that the
#' exam scores are in a bounded range. The posterior predictive checking
#' shows that we sometimes predict exam scores higher than 20, but the
#' discrepancy is minor.
#| label: fig-fitm-ppc-hist
#| fig-height: 3
#| fig-width: 6
pp_check(fitm, type = "hist", ndraws = 5)

#' LOO-PIT-ECDF plots shows that otherwise the normal model is quite well
#' calibrated.
#| label: fig-fitm-ppc-pitecdf
#| fig-height: 4
#| fig-width: 5
pp_check(fitm, type = "loo_pit_ecdf")

#' We could use truncated normal as more accurate model, but for
#' example beta-binomial model cannot be used for median exam scores
#' as some of the median scores are not integers. A fancier model
#' could model the three exams hierarchically, but as the normal model
#' is not that bad, we now continue with it.
#' 
#' ## Projection predictive variable selection
#'
#' We use projection predictive variable selection implemented in
#' `projpred` R package to find the minimal set of predictors that can
#' provide similar predictive performance as all predictor jointly.
#' By default `projpred` starts from the intercept only model and uses
#' forward search to find in which order to add predictors to minimize
#' the divergence from the full model predictive distribution. 
#'
#' ### Math exam scores
#' 
#' We start with doing fast PSIS-LOO-CV only for the full data search path.
#'
#| label: vselm_fast
#| results: hide
#| cache: true
vselm_fast <- cv_varsel(fitm, nterms_max = 27, validate_search = FALSE)

#' The following plot shows the relevance order of the predictors and
#' estimated predictive performance given those variables. As the
#' search can overfit and we didn't cross-validate the search, the
#' performance estimates can go above the reference model
#' performance. However, this plot helps as to see that 10 or fewer
#' predictors would be sufficient.
#'
#| label: fig-vselm_fast
#| fig-width:  8
#| fig-height: 6
#| out-width: "100%"
#| cache: false
plot(vselm_fast, stats = c("elpd", "R2"), deltas = "mixed", 
     text_angle = 45, alpha = 0.1,  size_position = "primary_x_top", 
     show_cv_proportions = FALSE) +
  geom_vline(xintercept = seq(0, 25, by = 5), colour = "black", alpha = 0.1)


#' Next we repeat the search, but now cross-validate the search,
#' too. We repeat the search with PSIS-LOO-CV criterion only for
#' `nloo=50` folds, and combine the result with the fast PSIS-LOO
#' result using difference estimator
#' [@Magnusson+etal:2020:bigloocomparison]. The use of sub-sampling LOO
#' affects only the model size selection and given that is stable,
#' the projection for the selected model is as good as with computationally
#' more expensive search validation. Based on the previous
#' quick result, we search only up to models of size 10. 
#' With my laptop and 8 parallel workers, this takes less than 5min.
#| label: vselm-subsampling-loo
#| results: hide
#| cache: true
registerDoFuture()
plan(multisession, workers = 8)
tic()
vselm <- cv_varsel(fitm, nterms_max = 10, validate_search = TRUE,
                   refit_prj = TRUE, nloo = 50,
                   parallel = TRUE, verbose = TRUE)
toc()
plan(sequential)

#' The following plot shows the relevance order of the predictors and
#' estimated predictive performance given those variables. The order
#' is the same as in the previous plot, but now the predictive
#' performance estimates are taking into account search and have
#' smaller bias. It seems using just four predictors can provide the
#' similar predictive performance as using all the predictors.
#| label: fig-vselm
#| fig-width:  8
#| fig-height: 5
#| out-width: "100%"
#| cache: false
plot(vselm, stats = c("elpd", "R2"), deltas = "mixed", 
     text_angle = 45, alpha = 0.1, size_position = "primary_x_top", 
     show_cv_proportions = FALSE) +
  geom_vline(xintercept = seq(0, 10, by = 5), colour = "black", alpha = 0.1)

#' `projpred` can also provide suggestion for the sufficient model size.
(nselm <- suggest_size(vselm))

#' Form the projected posterior for the selected model.
#' 
#| results: hide
#| cache: true
rankm <- ranking(vselm, nterms = nselm)
projm <- project(vselm, nterms = nselm)
drawsm_proj <- as_draws_df(projm) |>
  subset_draws(variable = paste0('b_', rankm$fulldata[1:nselm])) |> 
  set_variables(variable = rankm$fulldata[1:nselm])

#' The marginals of the projected posterior are all clearly away from 0.
#| label: fig-fitm-proj-mcmc_areas
#| fig-height: 2
#| fig-width: 8
mcmc_areas(drawsm_proj, prob_outer = 0.98, area_method = "scaled height")

#' The following plot shows the stability of the search over the
#' different LOO-CV folds. The numbers indicate the proportion of
#' folds, the specific predictor was included at latest on the given
#' model size.
#| label: fig-vselm-cv_proportions
#| fig-width:  8
#| fig-height: 5
#| out-width: "95%"
plot(cv_proportions(rankm, cumulate = TRUE))

#' ### Portuguese exam scores
#'
#' We repeat the same, but predicting grade for Portuguese instead of mathematics
#'
#' Fit a model with R2D2 prior with mean 1/3 and precision 3.
#| results: hide
#| cache: true
fitp <- brm(Gpor ~ ., data = studentstd_Gpor,
            prior = c(prior(R2D2(mean_R2 = 1/3, prec_R2 = 3, cons_D2 = 1/2), class = b)),
            refresh = 0)

#' Compare posterior-$R^2$ and LOO-$R^2$. We see that Portuguese grade
#' is easier to predict given the predictors (but there is still a lot
#' of unexplained variance).
#' 
#| cache: TRUE
fitp <- add_criterion(fitp, criterion = "loo")
#+
bayes_R2(fitp) |> round(2)
loo_R2(fitp) |> round(2)

#' Plot marginal posteriors of coefficients
#| label: fig-fitp-mcmc_areas
#| fig-height: 6
#| fig-width: 8
drawsp <- as_draws_df(fitp, variable = paste0("b_",predictors)) |>
  set_variables(predictors)
p <- mcmc_areas(drawsp, prob_outer = 0.98, area_method = "scaled height") +
  xlim(c(-1.5,1.5))
p <- p + scale_y_discrete(limits = rev(levels(p$data$parameter)))
p

#' We use projection predictive variable selection with fast LOO-CV
#' only for the full data search path.
#| label: vselp_fast
#| results: hide
#| cache: true
vselp_fast <- cv_varsel(fitp, nterms_max = 27, validate_search = FALSE)

#' The following plot shows the relevance order of the predictors and
#' estimated predictive performance given those variables. As there is
#' some overfitting in the search and we didn't cross-validate the
#' search, the performance estimates scan go above the reference model
#' performance. However, this plot helps as to see that 10 or fewer
#' predictors would be sufficient.
#| label: fig-vselp_fast
#| fig-width:  8
#| fig-height: 6
#| out-width: "100%"
#| cache: false
plot(vselp_fast, stats = c("elpd","R2"), deltas = "mixed", 
     text_angle = 45, alpha = 0.1,  size_position = "primary_x_top", 
     show_cv_proportions = FALSE) +
  geom_vline(xintercept = seq(0, 25, by = 5), colour = "black", alpha = 0.1)

#' Next we repeat the search, but now cross-validate the search,
#' too. We repeat the search with PSIS-LOO-CV criterion only for
#' `nloo=50` folds, and combine the result with the fast PSIS-LOO
#' result using difference estimator
#' [@Magnusson+etal:2020:bigloocomparison].   The use of sub-sampling LOO
#' affects only the model size selection and given that is stable,
#' the projection for the selected model is as good as with computationally
#' more expensive search validation. Based on the previous
#' quick result, we search only up to models of size 10. With my
#' laptop and 8 parallel workers, this takes less than 5min.
#| label: vselp-subsampling-loo
#| results: hide
#| cache: true
registerDoFuture()
plan(multisession, workers = 8)
tic()
vselp <- cv_varsel(fitp, nterms_max = 10, validate_search = TRUE, 
                   refit_prj = TRUE, nloo = 50, 
                   parallel = TRUE)
toc()
plan(sequential)

#' The following plot shows the relevance order of the predictors and
#' estimated predictive performance given those variables. The order
#' is the same as in the previous plot, but now the predictive
#' performance estimates are taking into account search and have
#' smaller bias. It seems using just seven predictors can provide the
#' similar predictive performance as using all the predictors.
#| label: fig-vselp
#| fig-width:  8
#| fig-height: 5
#| out-width: "100%"
#| cache: false
plot(vselp, stats = c("elpd", "R2"), deltas = "mixed", 
     text_angle = 45, alpha = 0.1, size_position = "primary_x_top", 
     show_cv_proportions = FALSE) +
  geom_vline(xintercept = seq(0, 10, by = 5), colour = "black", alpha = 0.1)

#' `projpred` can also provide suggestion for the sufficient model size.
(nselp <- suggest_size(vselp))

#' Form the projected posterior for the selected model.
#' 
#| results: hide
#| cache: true
rankp <- ranking(vselp, nterms = nselp)
projp <- project(vselp, nterms = nselp)
drawsp_proj <- as_draws_df(projp) |>
  subset_draws(variable = paste0("b_", rankp$fulldata[1:nselp])) |>
  set_variables(variable = rankp$fulldata[1:nselp])

#' The marginals of the projected posterior are all clearly away from 0.
#| label: fig-fitp-proj-mcmc_areas
#| fig-height: 2.5
#| fig-width: 8
mcmc_areas(drawsp_proj, prob_outer = 0.98, area_method = "scaled height")

#' The following plot shows the stability of the search over the
#' different LOO-CV folds. The numbers indicate the proportion of
#' folds, the specific predictor was included at latest on the given
#' model size.
#| label: fig-vselp-cv_proportions
#| fig-width: 8
#| fig-height: 5
#| out-width: "95%"
plot(cv_proportions(rankp, cumulate = TRUE))

#' ### Using the selected model
#' 
#' For further predictions we can use the projected draws. Due to how
#' different packages work, sometimes it can be easier to rerun MCMC
#' conditionally on the selected variables. This gives a slightly
#' different result, but when the reference model has been good the
#' difference tends to be small, and the main benefit form using
#' `projpred` is still that the selection process itself has not
#' caused overfitting and selection of spurious covariates.
#' 



#'
#' ## References {.unnumbered}
#'
#' <div id="refs"></div>
#'
#' ## Licenses {.unnumbered}
#' 
#' * Code &copy; 2023-2024, Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2023-2024, Aki Vehtari, licensed under CC-BY-NC 4.0.
#' 
#' ## Original Computing Environment {.unnumbered}
#' 
sessionInfo()

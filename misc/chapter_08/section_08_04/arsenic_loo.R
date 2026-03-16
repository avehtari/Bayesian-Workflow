#' ---
#' title: "Pointwise LOO-CV comparison demo"
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
#' book Section 8.4 about the influence of individual data points.

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
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(rstanarm)
library(loo)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans", base_size = 14))
library(khroma)
library(latex2exp)
library(patchwork)

wells <- read.csv(root("misc/chapter_08/section_08_04","wells.csv"))

fit_3 <- stan_glm(switch ~ dist100 + arsenic,
                  family = binomial(link = "logit"),
                  data = wells,
                  refresh = 0)
(loo3 <- loo(fit_3, save_psis = TRUE))

fit_3l <- stan_glm(switch ~ dist100 + log10(arsenic),
                   family = binomial(link = "logit"),
                   data = wells,
                   refresh = 0)
(loo3l <- loo(fit_3l, save_psis = TRUE))

wells$elpd_loo3=pointwise(loo3, "elpd_loo")
wells$elpd_loo3l=pointwise(loo3l, "elpd_loo")

p1 <- wells |> ggplot(aes(x=elpd_loo3, y=elpd_loo3l, color=factor(switch), shape=factor(switch))) +
  coord_fixed() +
  geom_abline(color="gray") +
  xlim(range(cbind(wells$elpd_loo3,wells$elpd_loo3l)))+
  ylim(range(cbind(wells$elpd_loo3,wells$elpd_loo3l)))+
  scale_color_bright(labels = c("Didn't switch","Switched")) +
  scale_shape_manual(values=c(1,3), labels = c("Didn't switch","Switched"))+
  geom_point(alpha=1) +
  theme(legend.position = "inside", legend.position.inside = c(0.7, 0.2),
        legend.title = element_blank(), legend.box.background = element_rect(color="grey")) +
  labs(x=TeX("$LOO_1$"), y=TeX("$LOO_2$"))
p1

p2 <- wells |> ggplot(aes(x=log10(arsenic), y=elpd_loo3-elpd_loo3l, color=factor(switch), shape=factor(switch))) +
  geom_hline(yintercept=0,  color="gray") +
  scale_color_bright(labels = c("Didn't switch","Switched")) +
  scale_shape_manual(values=c(1,3), labels = c("Didn't switch","Switched"))+
  geom_point(alpha=1) +
  theme(legend.position = "inside", legend.position.inside = c(0.3, 0.2),
        legend.title = element_blank(), legend.box.background = element_rect(color="grey")) +
  labs(x=TeX("$log_{10}(arsenic)$"), y=TeX("$LOO_1-LOO_2$"))
p2

p1 + p2

#ggsave("arsenic_loo_comparison_new.pdf", width=8.2, height=4)

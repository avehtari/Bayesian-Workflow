#' ---
#' title: "Kilpisjärvi prior predictive checking demo for Bayesian Workflow book"
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
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(ggplot2)
library(patchwork)

#' # Linear Gaussian model
#' 
#' Use the Kilpisjärvi summer month temperatures 1952--2022 data
#' recorded by Finnish Meteorological Institute (CC BY 4.0).
data_kilpis <- read.delim(root("digits/data", "kilpisjarvi-summer-temp-2022.csv"), sep = ";")
data_lin <-list(N = nrow(data_kilpis),
             x = data_kilpis$year,
             xpred = 2023,
             y = data_kilpis[,5])


data_lin_priors <- c(list(
    pmualpha_c = 0,     # prior mean for average temperature
    psalpha = 100,        # weakly informative
    pmubeta = 0,         # a priori incr. and decr. as likely
    psbeta = 100,   # avg temp prob does does not incr. more than a degree per 10 years:  setting this to +/-3 sd's
    pssigma = 1),        # setting sd of total variation in summer average temperatures to 1 degree implies that +/- 3 sd's is +/-3 degrees: 
  data_lin)


ns=1000
xc=data_kilpis$year-mean(data_kilpis$year)
n=length(xc)
mupp=rnorm(ns, mean=0, sd=100) + matrix(rnorm(ns*n, mean=0, sd=100),nrow=ns)%*%xc

theme_set(bayesplot::theme_default(base_family = "sans", base_size=15))

p1<-ggplot() +
    geom_abline(
        intercept = rnorm(100, mean=0, sd=100),
        slope = rnorm(100, mean=0, sd=100),
        linewidth = 0.1,
        color = "black",
        alpha = 0.2
    ) +
  scale_x_continuous(limits=c(-31,31),breaks=c(-22.5,-12.5,-2.5,7.5,17.7,27.5),labels=c(1960,1970,1980,1990,2000,2010))+
  ylim(y=c(-10000,10000))
## +
##   ggtitle("Models drawn from a wide prior")

p2<-ggplot() +
    geom_abline(
        intercept = rnorm(100, mean=0, sd=10)+10,
        slope = rnorm(100, mean=0, sd=.1/3),
        linewidth = 0.1,
        color = "black",
        alpha = 0.2
    ) +
  scale_x_continuous(limits=c(-31,31),breaks=c(-22.5,-12.5,-2.5,7.5,17.7,27.5),labels=c(1960,1970,1980,1990,2000,2010))+
  ylim(y=c(-20,40)) ## +
  ## ggtitle("Models drawn from an informative prior")

#| label: fig-kilpisjarvi_priorcheck
#| width: 9
#| height: 3
p1+p2

#ggsave(root("Digits","kilpisjarvi_priorcheck.pdf"), width=9, height=3)

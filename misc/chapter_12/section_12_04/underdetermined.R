#' 
#' **Load packages**
#| cache: FALSE
library(rprojroot)
root <- has_file(".Bayesian-Workflow-root")$make_fix_file()
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(patchwork)
library(cmdstanr)
options(mc.cores = 4)
library(posterior)

#' Simple analytic 2D posterior
y = c(-0.6906332,  1.2873201, -1.6422984,  2.0089285, 0.1347772 -0.8905993 -1.5295488 -0.8171846)
underdetermined_theory <- crossing(
  mu = seq(-3, 3, length.out = 200),
  log_sigma = seq(-2, 2, length.out = 200),
  data.frame(y_id = 1:length(y), y = y), n_data = c(1,2,8)
) |>
  filter(y_id <= n_data) |>
  mutate(log_density_point = dnorm(y, mu, exp(log_sigma), log = TRUE)) |>
  group_by(mu, log_sigma, n_data) |>
  summarise(log_density = sum(log_density_point), .groups = "drop") |>
  group_by(n_data) |>
  mutate(rel_density = exp(log_density - max(log_density))) |>
  ggplot(aes(x = mu, y = log_sigma, z = rel_density)) +
  geom_contour() + 
  facet_wrap(~n_data, nrow = 1, labeller = label_bquote(cols = paste(N, " = ", .(n_data)) )) +
  scale_y_continuous("log(sigma)")
#| label: fig-underdetermined_theory
underdetermined_theory

#' Lotka-Volterra model with Hudson bay lynx-hare data
lynx_hare_df <- read.csv(
  root("misc/chapter_12/section_12_04", "hudson-bay-lynx-hare.csv"),
  comment.char="#"
)
lv_model <- cmdstan_model(
  root("misc/chapter_12/section_12_04","lotka-volterra.stan")
)

#' Full data
lynx_hare_data_full <- with(
  lynx_hare_df,
  list(
    N = length(lynx_hare_df$Year) - 1,
    ts = 1:N,
    y_init = c(Hare[1], Lynx[1]),
    y = as.matrix(lynx_hare_df[2:(N + 1), c(3,2)]) # hare, lynx order
  )
)

#| label: fit_full
#| results: hide
fit_full <- lv_model$sample(data = lynx_hare_data_full,
                            init = 0.1,
                            refresh = 0)

drws_full <- bind_draws(fit_full$draws(variables = c("sigma[1]","theta[1]"), format="df"),
                        fit_full$sampler_diagnostics(format="df")) |>
  mutate_variables(n_data = lynx_hare_data_full$N + 1)


#' Mid sized data
set.seed(589762569)
N_mid <- 8
ts <- 1:N_mid
y_init <- c(lynx_hare_df$Hare[1], lynx_hare_df$Lynx[1])
y <- as.matrix(lynx_hare_df[2:(N_mid + 1), 2:3])
y <- cbind(y[ , 2], y[ , 1]); # hare, lynx order
lynx_hare_data_mid <- list(N = N_mid, ts = ts, y_init = y_init, y = y)

N <- 8
lynx_hare_data_mid <- with(
  lynx_hare_df,
  list(
    N = N,
    ts = 1:N,
    y_init = c(Hare[1], Lynx[1]),
    y = as.matrix(lynx_hare_df[2:(N + 1), c(3,2)]) # hare, lynx order
  )
)

#| label: fit_mid
#| results: hide
fit_mid <- lv_model$sample(data = lynx_hare_data_mid,
                           init = 0.1,
                           refresh = 0)

drws_mid <- bind_draws(fit_mid$draws(variables = c("sigma[1]","theta[1]"), format="df"),
                        fit_mid$sampler_diagnostics(format="df")) |>
  mutate_variables(n_data = lynx_hare_data_mid$N + 1)

#' Small data
set.seed(3248856)
N <- 5
lynx_hare_data_small <- with(
  lynx_hare_df,
  list(
    N = N,
    ts = 1:N,
    y_init = c(Hare[1], Lynx[1]),
    y = as.matrix(lynx_hare_df[2:(N + 1), c(3,2)]) # hare, lynx order
  )
)

#| label: fit_small
#| results: hid
fit_small <- lv_model$sample(data = lynx_hare_data_small,
                               init = 0.1,
                               refresh = 0)

drws_small <- bind_draws(fit_small$draws(variables = c("sigma[1]","theta[1]"), format="df"),
                        fit_small$sampler_diagnostics(format="df")) |>
  mutate_variables(n_data = lynx_hare_data_small$N + 1)


underdetermined_lv <- bind_draws(drws_full, drws_mid, drws_small, along = "draw") |>
  filter(divergent__==0) |>
  ggplot(aes(x = `theta[1]`, y = log(`sigma[1]`))) +
  geom_point(shape = 20, color = bayesplot:::get_color("m"), fill = bayesplot:::get_color("lh"), alpha = 0.1, size = 2) +
  geom_point(data = drws_small |> filter(divergent__==1),
             shape = 23, fill = "red", color = "white", size = 2) +
  facet_wrap(~n_data, nrow = 1, labeller = label_bquote(cols = paste(N, " = ", .(n_data)) )) +
  scale_y_continuous("log(sigma[1])") +
  scale_x_continuous("log(theta[1])")
underdetermined_lv

underdetermined_all <-
  (underdetermined_theory + labs(tag = "A")) |
  (undetermined_lv + labs(tag = "B"))
underdetermined_all

#ggsave("underdetermined.pdf", underdetermined_all, width = 8, height = 3)


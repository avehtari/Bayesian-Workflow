# heavily inspired by # downloaded from https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
# Sleepstudy example
library(brms)
library(tidyverse)

#-------------  Sleepstudy
# Reaction times (ms) to a series of tests after sleep deprivation.
# subjects sleep max 3 h for a number of days
data(sleepstudy, package = "lme4")
# str(sleepstudy)
# summary(sleepstudy)
sleepstudy <- sleepstudy %>% 
  as_tibble() %>% 
  mutate(Subject = as.character(Subject))

xlab <- "Days of sleep deprivation"
ylab <- "Average reaction time (ms)"

# No pooling
prior_nop <- prior(normal(5, 0.55), class = Intercept) +
  prior(normal(0, 0.2), class = b) +
  prior(exponential(3), class = sigma)

fit_nop_empty <- brm(Reaction ~ Days, data = sleepstudy,
                     family = lognormal(), prior = prior_nop, 
                     chains = 0)

Subjects <- unique(sleepstudy$Subject)
coefs_nop <- lapply(Subjects, function(S) {
  fit <- update(fit_nop_empty, newdata = filter(sleepstudy, Subject == S))
  out <- matrix(fixef(fit)[, 1], nrow = 1)
  out <- as.data.frame(out) 
  names(out) <- c("Intercept", "Slope")
  out$Subject <- S
  return(out)
})
coefs_nop <- Reduce(rbind, coefs_nop) %>%
  add_column(Model = "Independent models")

# Partial Pooling
prior_pp <- prior(normal(5, 0.55), class = Intercept) +
  prior(normal(0, 0.2), class = b) +
  prior(exponential(6), class = sigma) +
  prior(exponential(6), class = sd, group = Subject, coef = Intercept) +
  prior(exponential(10), class = sd, group = Subject, coef = Days) +
  prior(lkj(1), class = cor)

fit_pp <- brm(Reaction ~ 1 + Days + (1 + Days | Subject), 
              data = sleepstudy, family = lognormal(),
              prior = prior_pp)

coefs_pp <- coef(fit_pp)[["Subject"]][, 1, ] %>% 
  as.data.frame() %>% 
  rownames_to_column("Subject") %>% 
  rename(Slope = Days) %>% 
  add_column(Model = "Joint multilevel model")

coefs_all <- bind_rows(coefs_nop, coefs_pp)

gg <- ggplot(coefs_all) + 
  aes(x = Intercept, y = Slope, color = Model, shape = Model) + 
  geom_point(size = 2) + 
  # geom_point(
  #   data = df_shrinkage, 
  #   size = 3,
  #   show.legend = FALSE
  # ) + 
  geom_path(
    aes(group = Subject, color = NULL), 
    arrow = arrow(length = unit(.02, "npc")),
    show.legend = FALSE
  ) + 
  ggrepel::geom_text_repel(
    aes(label = Subject, color = NULL), 
    data = coefs_nop,
    show.legend = FALSE
  ) + 
  theme_bw() +
  theme(
    legend.position = "bottom"
    # legend.justification = "right"
  ) + 
  # ggtitle("Pooling of regression parameters") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_shape_manual(values = c(15:18)) +
  scale_color_viridis_d()

gg

ggsave("plots/sleep_multilevel_shrinkage_ln.pdf", height = 4, width = 7)

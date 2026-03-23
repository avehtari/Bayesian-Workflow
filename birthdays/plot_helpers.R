#' # Plot helpers for birthdays.R
#'
#' Shared color aliases, layer lists, holiday annotations, subplot
#' factory functions, and patchwork layout helpers.  Source this file
#' once near the top of birthdays.R after `set1` is defined.

# ── Step 1: Named color aliases ───────────────────────────────────────────────
#' col_fit  — model fit line / points  (Set1 red,  set1[1])
#' col_data — raw data scatter points  (Set1 blue, set1[2])
col_fit  <- set1[1]
col_data <- set1[2]

# ── Step 2: Shared layer lists ────────────────────────────────────────────────
#' Each variable is a list() of ggplot2 layer / scale objects.
#' Add any of them to a ggplot with `+`, e.g.:
#'   p + layers_hline100 + layers_labs_date

#' Raw data scatter background (alpha = 0.2, used in every subplot)
layers_data_scatter <- list(
  geom_point(color = col_data, alpha = 0.2)
)

#' Horizontal reference line at y = 100 (relative births scale)
layers_hline100 <- list(
  geom_hline(yintercept = 100, color = "gray")
)

#' Horizontal reference line at y = 1 (ratio scale, used in pf2c)
layers_hline1 <- list(
  geom_hline(yintercept = 1, color = "gray")
)

#' Standard axis labels for date-vs-births panels
layers_labs_date <- list(
  labs(x = "Date", y = "Relative number of births")
)

#' Month-tick x-axis for day-of-year panels (pf2, pf2b, pf2c)
layers_scale_doy <- list(
  scale_x_date(date_breaks = "1 month", date_labels = "%b")
)

#' Weekday x-axis labels (pf3)
layers_scale_weekday <- list(
  scale_x_continuous(
    breaks = 1:7,
    labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  )
)

# ── Step 3: Holiday annotation lists ─────────────────────────────────────────
#' holiday_labels_fixed(Ef4vec)
#'   Returns 7 annotate("text") calls for fixed US holidays.
#'   Ef4vec  — 366-element numeric vector (one value per day of 1988,
#'              a leap year) used to position labels on the y axis.
#'   The offsets (+/-) keep labels clear of the fit line.
holiday_labels_fixed <- function(Ef4vec) {
  list(
    annotate("text", x = as.Date("1988-01-01"),      y = Ef4vec[1]   - 1,   label = "New year"),
    annotate("text", x = as.Date("1988-02-14"),      y = Ef4vec[45]  + 1.5, label = "Valentine's day"),
    annotate("text", x = as.Date("1988-02-29"),      y = Ef4vec[60]  - 2.5, label = "Leap day"),
    annotate("text", x = as.Date("1988-04-01"),      y = Ef4vec[92]  - 1.5, label = "April 1st"),
    annotate("text", x = as.Date("1988-07-04")+5,      y = Ef4vec[186] - 1.5, label = "Independence day"),
    annotate("text", x = as.Date("1988-10-31"),      y = Ef4vec[305] - 1.5, label = "Halloween"),
    annotate("text", x = as.Date("1988-12-24"),      y = Ef4vec[360] - 2,   label = "Christmas")
  )
}

#' holiday_labels_floating(Ef4vec)
#'   Fixed 7 + 3 floating US holidays (Memorial day, Labor day,
#'   Thanksgiving) as annotate("text") calls.
holiday_labels_floating <- function(Ef4vec) {
  c(
    holiday_labels_fixed(Ef4vec),
    list(
      annotate("text", x = as.Date("1988-05-30"), y = Ef4vec[151] - 2, label = "Memorial day"),
      annotate("text", x = as.Date("1988-09-05"), y = Ef4vec[249] - 1.5, label = "Labor day"),
      annotate("text", x = as.Date("1988-11-24"), y = Ef4vec[329] - 1,   label = "Thanksgiving")
    )
  )
}

#' holiday_labels_fixed_ratio(Ef4vec)
#'   Same 7 fixed holidays for pf2c (ratio-scale y, offsets in ~0.02
#'   units rather than 1-2 birth units).
holiday_labels_fixed_ratio <- function(Ef4vec) {
  list(
    annotate("text", x = as.Date("1988-01-01") + 2,  y = Ef4vec[1]   - .025,  label = "New year"),
    annotate("text", x = as.Date("1988-02-14"),       y = Ef4vec[45]  + .025,  label = "Valentine's day"),
    annotate("text", x = as.Date("1988-02-29"),       y = Ef4vec[60]  - .03,  label = "Leap day"),
    annotate("text", x = as.Date("1988-04-01"),       y = Ef4vec[92]  - .025, label = "April 1st"),
    annotate("text", x = as.Date("1988-07-04") + 6,  y = Ef4vec[186] - .02,  label = "Independence day"),
    annotate("text", x = as.Date("1988-10-31"),       y = Ef4vec[305] - .025,  label = "Halloween"),
    annotate("text", x = as.Date("1988-12-24"),       y = Ef4vec[360] - .025, label = "Christmas")
  )
}

#' holiday_labels_floating_ratio(Ef4vec)
#'   Fixed 7 + 3 floating for pf2c (ratio scale).
holiday_labels_floating_ratio <- function(Ef4vec) {
  c(
    holiday_labels_fixed_ratio(Ef4vec),
    list(
      annotate("text", x = as.Date("1988-05-30") - 5, y = Ef4vec[151] - .03, label = "Memorial day"),
      annotate("text", x = as.Date("1988-09-05"),      y = Ef4vec[249] - .02, label = "Labor day"),
      annotate("text", x = as.Date("1988-11-24"),      y = Ef4vec[329] - .02, label = "Thanksgiving")
    )
  )
}

# ── Step 4: Plot factory functions ────────────────────────────────────────────

#' make_pf(data, Ef_vec, fit_geom = c("line", "point"))
#'   Overall fit vs raw data (full time series).
#'   data     — birthdays data frame (must contain columns date,
#'               births_relative100)
#'   Ef_vec   — numeric vector of posterior median fit values
#'               (same length as nrow(data))
#'   fit_geom — "line" (models 1-7) or "point" (model 8+, where the
#'               fit is day-tagged and a line would be misleading)
make_pf <- function(data, Ef_vec, fit_geom = c("line", "point")) {
  fit_geom <- match.arg(fit_geom)
  p <- data |>
    mutate(.Ef = Ef_vec) |>
    ggplot(aes(x = date, y = births_relative100)) +
    layers_data_scatter +
    layers_labs_date
  if (fit_geom == "line") {
    p <- p + geom_line(aes(y = .Ef), color = col_fit, alpha = 0.75)
  } else {
    p <- p + geom_point(aes(y = .Ef), color = col_fit, alpha = 0.2)
  }
  p
}

#' make_pf1(data, Ef1_vec)
#'   Slow-trend component vs raw data (full time series).
#'   Arguments as in make_pf(); always uses geom_line for the fit.
make_pf1 <- function(data, Ef1_vec) {
  data |>
    mutate(.Ef1 = Ef1_vec) |>
    ggplot(aes(x = date, y = births_relative100)) +
    layers_data_scatter +
    geom_line(aes(y = .Ef1), color = col_fit) +
    layers_hline100 +
    layers_labs_date
}

#' make_pf2(data, Ef2_vec, date_breaks = "1 month")
#'   Day-of-year seasonal component: data aggregated by day_of_year2,
#'   fit overlaid as a line.
#'   data        — birthdays data frame (columns: births_relative100,
#'                  day_of_year2)
#'   Ef2_vec     — numeric vector of fit values (same length as
#'                  nrow(data)); will be averaged over day_of_year2
#'   date_breaks — passed to scale_x_date; use "2 month" for the
#'                  compact 4-panel layout (model 4b fit)
make_pf2 <- function(data, Ef2_vec, date_breaks = "1 month") {
  data |>
    mutate(.Ef2 = Ef2_vec) |>
    group_by(day_of_year2) |>
    summarise(
      meanbirths = mean(births_relative100),
      meanEf2    = mean(.Ef2)
    ) |>
    ggplot(aes(
      x = as.Date("1987-12-31") + day_of_year2,
      y = meanbirths
    )) +
    layers_data_scatter +
    scale_x_date(date_breaks = date_breaks, date_labels = "%b") +
    geom_line(aes(y = meanEf2), color = col_fit) +
    layers_hline100 +
    layers_labs_date
}

#' make_pf3(data, Ef_dow_vec)
#'   Day-of-week component: data scatter + fit line on integer x axis.
#'   data        — birthdays data frame (columns: day_of_week,
#'                  births_relative100)
#'   Ef_dow_vec  — numeric vector of length 7 (one value per weekday)
make_pf3 <- function(data, Ef_dow_vec) {
  ggplot(data = data, aes(x = day_of_week, y = births_relative100)) +
    layers_data_scatter +
    layers_scale_weekday +
    geom_line(
      data = data.frame(x = 1:7, y = Ef_dow_vec),
      aes(x = x, y = y),
      color = col_fit
    ) +
    layers_hline100 +
    layers_labs_date
}

#' make_pf3b(data, Ef3_vec, Ef1_vec, weekday_labels = FALSE)
#'   Time-varying weekday magnitude component.
#'   data           — birthdays data frame (columns: date,
#'                     births_relative100, id)
#'   Ef3_vec        — time-varying weekday effect (same length as
#'                     nrow(data)); interpretation differs by model:
#'                     * models 4–7: raw weekday residuals, y axis is
#'                       births_relative100 / Ef1 / Ef2 * 100 * 100
#'                     * model 8+:   Ef3 * Ef1 / 100, y axis is
#'                       births_relative100
#'   Ef1_vec        — slow-trend vector (needed for model 8+ formula)
#'   weekday_labels — if TRUE, annotate the final four days with
#'                     Mon/Tue/Sat/Sun text labels (model 8+)
make_pf3b <- function(data, Ef3_vec, Ef1_vec = NULL,
                      weekday_labels = FALSE) {
  if (is.null(Ef1_vec)) {
    # Models 4–7: y = residuals after slow trend and day-of-year
    p <- data |>
      mutate(.Ef3 = Ef3_vec) |>
      ggplot(aes(x = date, y = births_relative100)) +
      layers_data_scatter +
      geom_point(aes(y = .Ef3), color = col_fit, size = 0.1) +
      layers_hline100 +
      layers_labs_date
  } else {
    # Model 8+: time-varying weekday × slow trend, plotted vs raw
    combined <- Ef3_vec * Ef1_vec / 100
    N <- nrow(data)
    p <- data |>
      mutate(.Ef3_combined = combined) |>
      ggplot(aes(x = date, y = births_relative100)) +
      layers_data_scatter +
      geom_point(aes(y = .Ef3_combined), color = col_fit, size = 0.1) +
      layers_hline100 +
      layers_labs_date
    if (weekday_labels) {
      p <- p +
        annotate(
          "text",
          x     = as.Date("1989-08-01"),
          y     = combined[c((N - 5):(N - 4), N, N - 6)],
          label = c("Mon", "Tue", "Sat", "Sun")
        )
    }
  }
  p
}

#' make_pf2b(Ef4_vec, f13, holidays = c("fixed", "floating"))
#'   Annotated day-of-year smooth line for year 1988.
#'   Ef4_vec  — 366-element numeric vector of day-of-year effect
#'               (on the births_relative100 / 100 scale, i.e. centred
#'               at 100, not 1)
#'   f13      — data frame with columns date & y for the 13th-of-month
#'               highlight points
#'   holidays — which holiday set to annotate
make_pf2b <- function(Ef4_vec, f13,
                      holidays = c("fixed", "floating")) {
  holidays <- match.arg(holidays)
  annots <- if (holidays == "fixed") {
    holiday_labels_fixed(Ef4_vec)
  } else {
    holiday_labels_floating(Ef4_vec)
  }
  data.frame(x = as.Date("1988-01-01") + 0:365, y = Ef4_vec) |>
    ggplot(aes(x = x, y = y)) +
    geom_line(color = col_fit) +
    layers_scale_doy +
    layers_hline100 +
    layers_labs_date +
    annots +
    geom_point(data = f13, aes(x = date, y = y), size = 3, shape = 1)
}

#' make_pf2c(Ef4_vec, Ef4r_rvar, f13,
#'            holidays = c("fixed", "floating"))
#'   Day-of-year effect with 90% posterior intervals (ggdist).
#'   Ef4_vec   — 366-element posterior median vector (ratio scale, ~1)
#'   Ef4r_rvar — rvar of the same quantity (from as_draws_rvars)
#'   f13       — data frame with columns date & y (ratio scale) for
#'                the 13th-of-month highlight points
#'   holidays  — which holiday set to annotate
make_pf2c <- function(Ef4_vec, Ef4r_rvar, f13,
                      holidays = c("fixed", "floating")) {
  holidays <- match.arg(holidays)
  annots <- if (holidays == "fixed") {
    holiday_labels_fixed_ratio(Ef4_vec)
  } else {
    holiday_labels_floating_ratio(Ef4_vec)
  }
  data.frame(
    x     = as.Date("1988-01-01") + 0:365,
    y     = Ef4_vec,
    ydist = Ef4r_rvar
  ) |>
    ggplot(aes(x = x, y = y, ydist = ydist)) +
    stat_pointinterval(.width = 0.9, color = col_fit,
                       alpha = 0.3, size = 0.5) +
    layers_scale_doy +
    layers_hline1 +
    layers_labs_date +
    annots +
    geom_point(data = f13, aes(x = date, y = y),
               inherit.aes = FALSE, size = 3, shape = 1)
}

#' make_pth_vs_fit(sm, sp, max.overlaps = 10)
#'   Scatter plot comparing Pathfinder and MCMC mean ± sd for each
#'   parameter.  Replaces the 11+ verbatim copies in birthdays.R.
#'   sm           — summarise_draws() result for MCMC draws
#'   sp           — summarise_draws() result for Pathfinder draws
#'   max.overlaps — passed to geom_text_repel(); use 20 for models
#'                   with many parameters (e.g. 8rhs second comparison)
make_pth_vs_fit <- function(sm, sp, max.overlaps = 10) {
  ggplot(
    data = NULL,
    aes(
      x = sm$mean, xmin = sm$mean - sm$sd, xmax = sm$mean + sm$sd,
      y = sp$mean, ymin = sp$mean - sp$sd, ymax = sp$mean + sp$sd,
      label = sm$variable
    )
  ) +
    geom_point(color = 4) +
    geom_errorbar(width = 0, color = 4) +
    geom_errorbarh(height = 0, color = 4) +
    geom_text_repel(max.overlaps = max.overlaps) +
    geom_abline(linetype = "dotted") +
    labs(x = "MCMC mean and sd", y = "Pathfinder mean and sd")
}

#' make_sd_ratio(sm, sp)
#'   Diagnostic scatter: sd(Pathfinder) / sd(MCMC) per parameter.
#'   sm, sp — summarise_draws() results (MCMC and Pathfinder)
make_sd_ratio <- function(sm, sp) {
  data.frame(varid = seq_len(nrow(sp)), sd_ratio = sp$sd / sm$sd) |>
    ggplot(aes(x = varid, y = sd_ratio)) +
    geom_point() +
    scale_y_log10() +
    labs(x = "Variable number", y = "sd(pathfinder) / sd(MCMC)") +
    geom_hline(yintercept = 1, color = "gray")
}

# ── Step 6: Patchwork layout helpers ─────────────────────────────────────────
#' compose_3panel(pf, pf1, pf2)
#'   Layout for models 1–2: full fit on top, components below.
#'     pf / (pf1 + pf2)
compose_3panel <- function(pf, pf1, pf2) {
  pf / (pf1 + pf2)
}

#' compose_4panel(pf, pf1, pf2, pf3_or_pf3b)
#'   Layout for models 3–4: 2×2 grid.
#'     (pf + pf1) / (pf2 + pf3_or_pf3b)
compose_4panel <- function(pf, pf1, pf2, pf3_or_pf3b) {
  (pf + pf1) / (pf2 + pf3_or_pf3b)
}

#' compose_6panel(pf, pf1, pf2, pf3_or_pf3b, pf2b)
#'   Layout for models 5–8: 3-row layout.
#'     (pf + pf1) / (pf2 + pf3_or_pf3b) / pf2b
compose_6panel <- function(pf, pf1, pf2, pf3_or_pf3b, pf2b) {
  (pf + pf1) / (pf2 + pf3_or_pf3b) / pf2b
}

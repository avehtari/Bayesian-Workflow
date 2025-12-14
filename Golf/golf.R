#' ---
#' title: "Model building and expansion for golf putting"
#' author: "Andrew Gelman and Aki Vehtari"
#' date: 2019-09-24
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
#' 
#' # Setup  {.unnumbered}
#' 
#+ setup, include=FALSE
knitr::opts_chunk$set(cache=FALSE, message=FALSE, error=FALSE, warning=FALSE, comment=NA, out.width='95%')

#' 
#' **Load packages**
#| code-fold: true
#| cache: FALSE
library("cmdstanr")
library("posterior")
print_file <- function(file) {
  cat(paste(readLines(file), "\n", sep=""), sep="")
}
logit <- qlogis
invlogit <- plogis
fround <- function (x, digits) format(round(x, digits), nsmall = digits)

#' # Introduction
#' 
#' The following graph shows data from professional golfers on the proportion of successful
#' putts as a function of distance from the hole [from @Berry:1995].  Unsurprisingly, the
#' probability of making the shot declines as a function of distance:
golf <- read.table("data/golf_data.txt", header = TRUE, skip = 2)
x <- golf$x
y <- golf$y
n <- golf$n
J <- length(y)
r <- (1.68 / 2) / 12
R <- (4.25 / 2) / 12
se <- sqrt((y / n) * (1 - y / n) / n)
#| label: fig-golf-data
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(0, 0, xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.02), 
     xaxs = "i", yaxs = "i", bty = "l", 
     xlab = "Distance from hole (feet)", 
     ylab = "Probability of success", 
     main = "Data on putts in pro golf", 
     type = "n")
points(x, y / n, pch = 20, col = "blue")
segments(x, y / n + se, x, y / n - se, lwd = .5, col = "blue")
text(x + .4,
     y / n + se + .02,
     paste(y, "/", n, sep = ""),
     cex = .6,
     col = "gray40")

#' The error bars associated with each point $j$ in the above graph are
#' simple classical standard deviations,
#' $\sqrt{\hat{p}_j(1-\hat{p}_j)/n_j}$, where $\hat{p}_j=y_j/n_j$ is the
#' success rate for putts taken at distance $x_j$.
#'
#' # Logistic regression
#'
#' Can we model the probability of success in golf putting as a function
#' of distance from the hole?  Given usual statistical practice, the
#' natural starting point would be logistic regression:
#'
#' $$
#' y_j\sim\mathrm{binomial}(n_j, \mathrm{logit}^{-1}(a + bx_j)),
#' \text{ for } j=1,\dots, J.
#' $$
#' In Stan, this is:
print_file("golf_logistic.stan")

#' The code in the above model block is (implicitly) vectorized, so
#' that it is mathematically equivalent to modeling each data point,
#' `y[i] ~ binomial_logit(n[i], a + b*x[i])`.  The vectorized code is
#' more compact (no need to write a loop, or to include the
#' subscripts) and faster (because of more efficient gradient
#' evaluations).
#'
#' We fit the model to the data:
#| label: golf_logistic.stan
#| results: hide
golf_data <- list(x = x, y = y, n = n, J = J)
model_1 <- cmdstan_model("golf_logistic.stan")
fit_1 <- model_1$sample(data = golf_data, refresh = 0)

#' And here is the result:
draws_1 <- fit_1$draws(format = "df")
a_sim <- draws_1$a
b_sim <- draws_1$b
a_hat <- median(a_sim)
b_hat <- median(b_sim)
n_sims <- nrow(draws_1)
print(fit_1)

#' Going through the columns of the above table: Stan has computed the
#' posterior means $\pm$ standard deviations of $a$ and $b$ to be
#' `{r} sprintf("%.2f", mean(a_sim))` $\pm$ `{r} sprintf("%.2f", sd(a_sim))`
#' and `{r} sprintf("%.2f", mean(b_sim))` $\pm$ `{r} sprintf("%.2f", sd(b_sim))`,
#' respectively. The Monte Carlo standard error of the mean of each of
#' these parameters is 0 (to two decimal places), indicating that the
#' simulations have run long enough to estimate the posterior means
#' precisely.  The posterior quantiles give a sense of the uncertainty
#' in the parameters, with 50% posterior intervals of
#' [`{r} sprintf("%.2f", quantile(a_sim, 0.25))`,
#' `{r} sprintf("%.2f", quantile(a_sim, 0.75))`] and
#' [`{r} sprintf("%.2f", quantile(b_sim, 0.25))`,
#' `{r} sprintf("%.2f", quantile(b_sim, 0.75))`] for $a$ and $b$, respectively.
#' Finally, the values of $\widehat{R}$ near 1 tell us that the
#' simulations from Stan's four simulated chains have mixed well.
#'
#' The following graph shows the fit plotted along with the data:
#| label: fig-golf-fit1
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(0, 0, xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.02), 
     xaxs = "i", yaxs = "i", bty = "l", 
     xlab = "Distance from hole (feet)", 
     ylab = "Probability of success", 
     main = "Fitted logistic regression", 
     type = "n")
for (i in sample(n_sims, 10)) {
  curve(invlogit(a_sim[i] + b_sim[i]*x),
        from = 0, to = 1.1 * max(x), 
        wd = 0.5, add = TRUE, col = "green")
}
curve(invlogit(a_hat + b_hat*x), from=0, to=1.1*max(x), add=TRUE)
points(x, y/n, pch=20, col="blue")
segments(x, y/n + se, x, y/n-se, lwd=.5, col="blue")
text(11, .57, paste("Logistic regression,\n    a = ",
                      fround(a_hat, 2), ", b = ", fround(b_hat, 2), sep=""))

#' The black line shows the fit corresponding to the posterior median
#' estimates of the parameters $a$ and $b$; the green lines show 10 draws
#' from the posterior distribution.
#'
#' In this example, posterior uncertainties in the fits are small, and
#' for simplicity we will just plot point estimates based on posterior
#' median parameter estimates for the remaining models.  Our focus
#' here is on the sequence of models that we fit, not so much on
#' uncertainty in particular model fits.
#'
#' # Modeling from first principles
#'
#' As an alternative to logistic regression, we shall build a model from
#' first principles and fit it to the data.
#'
#' The graph below shows a simplified sketch of a golf shot.  The
#' dotted line represents the angle within which the ball of radius
#' $r$ must be hit so that it falls within the hole of radius $R$.
#' This threshold angle is $\sin^{-1}((R-r)/x)$.  The graph, which is
#' not to scale, is intended to illustrate the geometry of the ball
#' needing to go into the hole.
#| label: fig-golf-sketch-1
#| fig-height:  2
#| fig-width:  7
par(mar=c(0,0,0,0))
dist <- 2
r_plot <- r
R_plot <- R
plot(0, 0, xlim = c(-R_plot, dist + 3 * R_plot), ylim = c(-2 * R_plot, 2 * R_plot), 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n", 
     xlab = "", ylab = "", type = "n", asp = 1)
symbols(0, 0, circles=r_plot, inches=FALSE, add=TRUE)
symbols(dist, 0, circles = R_plot - r_plot, inches = FALSE, lty = 2, add = TRUE)
symbols(dist, 0, circles = R_plot, inches = FALSE, add = TRUE)
curve(0 * x, from = 0, to = dist, add = TRUE)
curve(((R_plot - r_plot) / dist) * x, from = 0, to = dist, lty = 2, add = TRUE)
curve(-((R_plot - r_plot) / dist) * x, from = 0, to = dist, lty = 2, add = TRUE)
text(0.5 * dist, -1.5 * R_plot, "x")
arrows(0.5 * dist + 0.05, -1.5 * R_plot, dist, -1.5 * R_plot, 2, length = .1)
arrows(0.5 * dist - 0.05, -1.5 * R_plot, 0, -1.5 * R_plot, 2, length = .1)
text(dist + 1.2 * R_plot, .5 * R_plot, "R")
arrows(dist + 1.2 * R_plot, .7 * R_plot, dist + 1.2 * R_plot, R_plot, length = .05)
arrows(dist + 1.2 * R_plot, .3 * R_plot, dist + 1.2 * R_plot, 0, length = .05)
text(0, r_plot / 2, "r")

#' The next step is to model human error.  We assume that the golfer is
#' attempting to hit the ball completely straight but that many small
#' factors interfere with this goal, so that the actual angle follows a
#' normal distribution centered at 0 with some standard deviation
#' $\sigma$.
#| label: fig-golf-sketch-2
#| fig-height:  3
#| fig-width:  7
par(mar=c(3, 3, 0, 0), mgp = c(1.7, .5, 0), tck = -.02)
plot(0, 0, xlim = c(-4, 4), ylim = c(0, 1.02), 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n", 
     xlab = "Angle of shot", ylab = "", type = "n")
axis(1, seq(-4,4), 
     c("", "", expression(-2*sigma), "", 
       0, "", expression(2*sigma),"", ""))
curve(dnorm(x) / dnorm(0), add = TRUE)

#' The probability the ball goes in the hole is then the probability
#' that the angle is less than the threshold; that is,
#' $\mathrm{Pr}\left(|\mathrm{angle}| < \sin^{-1}((R-r)/x)\right) =
#' 2\Phi\left(\frac{\sin^{-1}((R-r)/x)}{\sigma}\right) - 1$, where
#' $\Phi$ is the cumulative normal distribution function.  The only
#' unknown parameter in this model is $\sigma$, the standard deviation
#' of the distribution of shot angles. Stan (and, for that matter, R)
#' computes trigonometry using angles in radians, so at the end of our
#' calculations we will need to multiply by $180/\pi$ to convert to
#' degrees, which are more interpretable by humans.
#' 
#' Our model then has two parts:
#' \begin{align}
#' y_j &\sim \mathrm{binomial}(n_j, p_j)\\
#' p_j &= 2\Phi\left(\frac{\sin^{-1}((R-r)/x_j)}{\sigma}\right) - 1 , \text{ for } j=1,\dots, J.
#' \end{align}
#' Here is a graph showing the curve for some potential values of the parameter $\sigma$.
#| label: fig-golf-angle-curves
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(0, 0, xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.02), 
     xaxs = "i", yaxs = "i", bty = "l", 
     xlab = "Distance from hole (feet)", ylab = "Probability of success", 
     main = expression(paste("Modeled Pr(success) for different values of ", sigma)), 
     type = "n")
sigma_degrees_plot <- c(0.5, 1, 2, 5, 20)
x_text <- c(15, 10, 6, 4, 2)
for (i in 1:length(sigma_degrees_plot)){
  sigma <- (pi / 180) * sigma_degrees_plot[i]
  x_grid <- seq(R - r, 1.1 * max(x), .01)
  p_grid <- 2 * pnorm(asin((R - r) / x_grid) / sigma) - 1
  lines(c(0, R - r, x_grid), c(1, 1, p_grid))
  text(x_text[i] + 0.7,
       2 * pnorm(asin((R - r) / x_text[i]) / sigma) - 1,
       bquote(sigma == .(sigma_degrees_plot[i]) * degree),
       adj = 0)
}

#' The highest curve on the graph corresponds to $\sigma=0.5^\circ$:
#' if golfers could control the angles of their putts to an accuracy
#' of approximately half a degree, they would have a very high
#' probability of success, making over 80\% of their ten-foot putts,
#' over 50\% of their fifteen-foot putts, and so on.  At the other
#' extreme, the lowest plotted curve corresponds to $\sigma=20^\circ$:
#' if your putts could be off as high as 20 degrees, then you would be
#' highly inaccurate, missing more than half of your two-foot
#' putts. When fitting the model in Stan, the program moves around the
#' space of $\sigma$, sampling from the posterior distribution.
#'
#' We now write the Stan model in preparation to estimating $\sigma$:
print_file("golf_angle_binomial.stan")

#' In the transformed data block above, the `./` in the calculation of
#' p corresponds to componentwise division in this vectorized
#' computation.
#' 
#' The data $J,n,x,y$ have already been set up; we just need to define
#' $r$ and $R$ (the golf ball and hole have diameters 1.68 and 4.25
#' inches, respectively), and run the Stan model:
#| label: golf_angle_binomial.stan
#| results: hide
r <- (1.68/2)/12
R <- (4.25/2)/12
golf_data <- c(golf_data, r=r, R=R)
model_2 <- cmdstan_model("golf_angle_binomial.stan")
fit_2 <- model_2$sample(data = golf_data, refresh = 0)

#' Here is the result:
draws_2 <- fit_2$draws(format = "df")
sigma_sim <- draws_2$sigma
sigma_degrees_sim <- draws_2$sigma_degrees
sigma_hat <- median(sigma_sim)
print(fit_2)

#' The model has a single parameter, $\sigma$.  From the output, we
#' find that Stan has computed the posterior mean of $\sigma$ to be 
#' `{r} sprintf("%.2f", mean(sigma_sim))`.  Multiplying this by $180/\pi$,
#' this comes to `{r} sprintf("%.2f", mean(sigma_degrees_sim))` degrees.
#' The Monte Carlo standard error of the mean is 0 (to two decimal
#' places), indicating that the simulations have run long enough to
#' estimate the posterior mean precisely.  The posterior standard
#' deviation is calculated at `{r} sprintf("%.2f",
#' sd(sigma_degrees_sim))` degrees, indicating that $\sigma$ itself
#' has been estimated with high precision, which makes sense given the
#' large number of data points and the simplicity of the model.  The
#' precise posterior distribution of $\sigma$ can also be seen from
#' the narrow range of the posterior quantiles.  Finally,
#' $\widehat{R}$ is near 1, telling us that the simulations from
#' Stan's four simulated chains have mixed well.
#'
#' We next plot the data and the fitted model (here using the
#' posterior median of $\sigma$ but in this case the uncertainty is so
#' narrow that any reasonable posterior summary would give essentially
#' the same result), along with the logistic regression fitted
#' earlier:
#| label: fig-golf-fit-1-2
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(0, 0, xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.02), 
     xaxs = "i", yaxs = "i", bty = "l", 
     xlab = "Distance from hole (feet)", 
     ylab = "Probability of success", 
     main = "Two models fit to the golf putting data", type = "n")
segments(x, y / n + se, x, y / n - se, lwd = .5)
curve(invlogit(a_hat + b_hat*x), 
      from = 0, to = 1.1 * max(x), 
      add = TRUE)
x_grid <- seq(R - r, 1.1 * max(x), .01)
p_grid <- 2 * pnorm(asin((R - r) / x_grid) / sigma_hat) - 1
lines(c(0, R - r, x_grid), c(1, 1, p_grid), col = "blue")
points(x, y / n, pch = 20, col = "blue")
text(10.3, .58, "Logistic regression")
text(18.5, .24, "Geometry-based model", col="blue")

#' The custom nonlinear model fits the data much better.  This is not to
#' say that the model is perfect---any experience of golf will reveal
#' that the angle is not the only factor determining whether the ball
#' goes in the hole---but it seems like a useful start.
#'
#' # Testing the fitted model on new data
#'
#' Several years after fitting the above model, we were presented with
#' a newer and more comprehensive dataset on professional golf putting
#' [@Broadie:2018].  For simplicity we'll just look here at the
#' summary data, probabilities of the ball going into the hole for
#' shots up to 75 feet from the hole.  The graph below shows these new
#' data (in red), along with our earlier dataset (in blue) and the
#' already-fit geometry-based model from before, extending to the
#' range of the new data.
golf_new <- read.table("data/golf_data_new.txt", header = TRUE, skip = 2)
#| label: fig-golf-fit-2-new
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02), 
     xaxs = "i", yaxs = "i", bty = "l", 
     xlab = "Distance from hole (feet)", 
     ylab = "Probability of success", 
     main = "Checking already-fit model to new data")
x_grid <- seq(R - r, 1.1 * max(golf_new$x), .01)
p_grid <- 2 * pnorm(asin((R - r) / x_grid) / sigma_hat) - 1
lines(c(0, R - r, x_grid), c(1, 1, p_grid), col = "blue")
points(golf$x, golf$y / golf$n, pch = 20, col = "blue")
points(golf_new$x, golf_new$y/golf_new$n, pch = 20, col = "red")
legend(60, 0.4, legend = c("Old data", "New data"), col = c("blue", "red"), pch = 20) 

#' Comparing the two datasets in the range 0-20 feet, the success rate
#' is similar for longer putts but is much higher than before for the
#' short putts. This could be a measurement issue, if the distances to
#' the hole are only approximate for the old data, and it could also
#' be that golfers are better than they used to be.
#'
#' Beyond 20 feet, the empirical success rates become lower than would
#' be predicted by the old model. These are much more difficult
#' attempts, even after accounting for the increased angular precision
#' required as distance goes up.  In addition, the new data look
#' smoother, which perhaps is a reflection of more comprehensive data
#' collection.
#'
#' # A new model accounting for how hard the ball is hit
#'
#' To get the ball in the hole, the angle isn’t the only thing you
#' need to control; you also need to hit the ball just hard enough.
#'
#' @Broadie:2018 added this feature to the geometric model by
#' introducing another parameter corresponding to the golfer's control
#' over distance. Supposing $u$ is the distance that golfer's shot
#' would travel if there were no hole, the assumption is that the putt
#' will go in if (a) the angle allows the ball to go over the hole,
#' and (b) $u$ is in the range $(x,x+3)$. That is the ball must be hit
#' hard enough to reach the whole but not go too far. Factor (a) is
#' what we have considered earlier; we must now add factor (b).
#'
#' The following sketch, which is not to scale, illustrates the need
#' for the distance as angle as well as the angle of the shot to be in
#' some range, in this case the gray zone which represents the
#' trajectories for which the ball would reach the hole and stay in
#' it.
#| label: fig-golf-sketch-3
#| fig-height:  2
#| fig-width: 7
par(mar = c(0, 0, 0, 0))
dist <- 2
r_plot <- r
R_plot <- R
distance_tolerance <- 0.6
plot(0, 0, 
     xlim = c(-R_plot, dist + 3 * R_plot + 1.5 * distance_tolerance), 
     ylim = c(-2 * R_plot, 2 * R_plot), 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n", 
     xlab = "", ylab = "", type = "n", asp = 1)
polygon(
  c(dist, dist, dist + distance_tolerance, dist + distance_tolerance),
  c(
    R_plot - r_plot,
    -(R_plot - r_plot),
    -(R_plot - r_plot) * (dist + distance_tolerance) / dist,
    (R_plot - r_plot) * (dist + distance_tolerance) / dist
  ),
  border = NA,
  col = "gray"
)
symbols(0, 0, circles = r_plot, inches = FALSE, add = TRUE)
symbols(dist, 0, circles = R_plot, inches = FALSE, add = TRUE)
symbols(dist, 0, circles = R_plot - r_plot, inches = FALSE, 
        lty = 2, bg = "gray", add = TRUE)
curve(((R_plot - r_plot) / dist) * x, from = 0, to = dist + 1.5 * distance_tolerance, 
      lty = 2, add = TRUE)
curve(-((R_plot - r_plot) / dist) * x, from = 0, to = dist + 1.5 * distance_tolerance, 
      lty = 2, add = TRUE)
text(0.5 * dist, -1.5 * R_plot, "x")
arrows(0.5 * dist + 0.05, -1.5 * R_plot, dist, -1.5 * R_plot, 2, length = .1)
arrows(0.5 * dist - 0.05, -1.5 * R_plot, 0, -1.5 * R_plot, 2, length = .1)

#' Suppose that a golfer will aim to hit the ball one foot past the
#' hole but with a multiplicative error in the shot's potential
#' distance, so that $u=(x+1)(1+\epsilon)$, where the error $\epsilon$
#' has a normal distribution with mean 0 and standard deviation
#' $\sigma_{\rm distance}$. In statistics notation, this model is,
#' $u\sim\normal(x+1,(x+1)\sigma_{\rm distance})$, and the distance is
#' acceptable if $u\in [x,x+3]$, an event that has probability
#' $\Phi\left(\frac{2}{(x+1)\sigma_{\rm distance}}\right)-\Phi\left(\frac{-1}{(x+1)\sigma_{\rm distance}}\right).$
#'
#' Putting these together, the probability a shot goes in becomes,
#' $\left(2\Phi\left(\frac{\sin^{-1}((R-r)/x)}{\sigma_{\rm
#' angle}}\right) -
#' 1\right)\left(\Phi\left(\frac{2}{(x+1)\,\sigma_{\rm
#' distance}}\right) - \Phi\left(\frac{-1}{(x+1)\,\sigma_{\rm
#' distance}}\right)\right)$, where we have renamed the parameter
#' $\sigma$ from our earlier model to $\sigma_{\rm angle}$ to
#' distinguish it from the new $\sigma_{\rm distance}$ parameter.  We
#' write the new model in Stan, giving it the name
#' `golf_angle_distance_binomial.stan` to convey that it accounts both
#' for angle and distance:
print_file("golf_angle_distance_binomial.stan")

#' The result is a model with two parameters, $\sigma_{\rm angle}$ and
#' $\sigma_{\rm distance}$. Even this improved geometry-based model is
#' a gross oversimplification of putting, and the average distances in
#' the binned data are not the exact distances for each shot.  But it
#' should be an advance on the earlier one-parameter model; the next
#' step is to see how it fits the data.
#'
#' Here we have defined `overshot` and `distance_tolerance` as data,
#' which @Broadie:2018 has specified as 1 and 3 feet, respectively.
#' We might wonder why if the distance range is 3 feet, the overshot
#' is not 1.5 feet. One reason could be that it is riskier to hit the
#' ball too hard than too soft.  In addition we assigned weakly
#' informative half-normal(0,1) priors on the scale parameters,
#' $\sigma_{\rm angle}$ and $\sigma_{\rm distance}$, which are
#' required in this case to keep the computations stable.
#'
#' We fit the model to the new dataset.
#| label: golf_angle_distance_binomial.stan
#| results: hide
overshot <- 1
distance_tolerance <- 3
golf_new_data <- list(
  x = golf_new$x,
  y = golf_new$y,
  n = golf_new$n,
  J = nrow(golf_new),
  r = r,
  R = R,
  overshot = overshot,
  distance_tolerance = distance_tolerance
)
model_3 <- cmdstan_model("golf_angle_distance_binomial.stan")
fit_3 <- model_3$sample(data = golf_new_data, refresh = 0)

#' Here is the result:
print(fit_3)

#' There is poor convergence. Very high Rhat and very low ESSs
#' indicate multimodality. We try initializing MCMC with Pathfinder
#| results: hide
pth_3 <- model_3$pathfinder(data = golf_new_data, refresh = 0,
                            num_paths=20, max_lbfgs_iters=100)
fit_3 <- model_3$sample(data = golf_new_data, refresh = 0, init = pth_3)

#' Here is the result:
print(fit_3)
#' Now the convergence looks fine.

#' We graph the new data and the fitted model:
#| label: fig-golf-fit-3
draws_3 <- fit_3$draws(format = "df")
sigma_angle_hat <- median(draws_3$sigma_angle)
sigma_distance_hat <- median(draws_3$sigma_distance)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02), 
     xaxs = "i", yaxs = "i", bty = "l", 
     xlab = "Distance from hole (feet)", 
     ylab = "Probability of success", 
     main = "Checking model fit", type = "n")
x_grid <- seq(R - r, 1.1 * max(golf_new$x), .01)
p_angle_grid <- (2 * pnorm(asin((R - r) / x_grid) / sigma_angle_hat) - 1)
p_distance_grid <- pnorm((distance_tolerance - overshot) / ((x_grid + overshot) * sigma_distance_hat)) -
  pnorm((-overshot) / ((x_grid + overshot) * sigma_distance_hat))
lines(c(0, R-r, x_grid), c(1, 1, p_angle_grid*p_distance_grid), col="red")
points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")

#' There are problems with the fit in the middle of the range of $x$.
#' We suspect this is a problem with the binomial error model, as it
#' tries harder to fit points where the counts are higher.  Look at
#' how closely the fitted curve hugs the data at the very lowest
#' values of $x$.
#'
#' Here are the first few rows of the data:
print(golf_new[1:5,])

#' With such large values of $n_j$, the binomial likelihood enforces
#' an extremely close fit at these first few points, and that drives
#' the entire fit of the model.
#'
#' # Expanding the model by including a fudge factor
#' 
#' To fix this problem we took the data model, $y_j \sim
#' \mathrm{binomial}(n_j, p_j)$, and added an independent error term to
#' each observation.  There is no easy way to add error directly to
#' the binomial distribution---we could replace it with its
#' overdispersed generalization, the beta-binomial, but this would not
#' be appropriate here because the variance for each data point $i$
#' would still be roughly proportional to the sample size $n_j$, and
#' our whole point here is to get away from this assumption and allow
#' for model misspecification---so instead we first approximate the
#' binomial data distribution by a normal and then add independent
#' variance; thus: $$y_j/n_j \sim \mathrm{normal}\left(p_j,
#' \sqrt{p_j(1-p_j)/n_j + \sigma_y^2}\right),$$ To write this in Stan
#' there are some complications:
#'
#' * $y$ and $n$ are integer variables, which we convert to vectors so
#'   that we can multiply and divide them.
#'
#' * To perform componentwise multiplication or division using
#'   vectors, you need to use `.*` or `./` so that San knows not to
#'   try to perform vector/matrix multiplication and division.  Stan
#'   is opposite from R in this way: Stan defaults to vector/matrix
#'   operations and has to be told otherwise, whereas R defaults to
#'   componentwise operations, and vector/matrix multiplication in R
#'   is indicated using the `%*%` operator.
#'
#' We implement these via the following new code in the transformed data block:
#' ```
#'   vector[J] raw_proportions = to_vector(y) ./ to_vector(n);
#' ```
#' And then in the model block:
#' ```
#'   raw_proportions ~ normal(p, sqrt(p .* (1-p) ./ to_vector(n) + sigma_y^2));
#' ```
#'
#' To complete the model we add $\sigma_y$ to the parameters block and
#' assign it a weakly informative half-normal(0,1) prior
#' distribution. Here's the new Stan program:
print_file("golf_angle_distance_normal.stan")

#' We now fit this model to the data:
#| label: golf_angle_distance_normal.stan
#| results: hide
model_4 <- cmdstan_model("golf_angle_distance_normal.stan")
fit_4 <- model_4$sample(data = golf_new_data, refresh = 0)

#' Here is the result
draws_4 <- fit_4$draws(format = "df")
print(fit_4)

#' The new parameter estimates are:
#' 
#' * $\sigma_{\rm angle}$ is estimated at `{r} sprintf("%.2f",
#'   mean(draws_4$sigma_angle))`, which when corresponds to
#'   $\sigma_{\rm degrees}=$ `{r} sprintf("%.1f",
#'   mean(draws_4$sigma_degrees))`.  According to the fitted
#'   model, there is a standard deviation of `{r} sprintf("%.1f",
#'   mean(draws_4$sigma_degrees))` degree in the angles of
#'   putts taken by pro golfers.  The estimate of $\sigma_{\rm angle}$
#'   has decreased compared to the earlier model that only had angular
#'   errors.  This makes sense: now that distance errors have been
#'   included in the model, there is no need to explain so many of the
#'   missed shots using errors in angle.
#'
#' * $\sigma_{\rm distance}$ is estimated at `{r} sprintf("%.2f",
#'   mean(draws_4$sigma_distance))`.  According to the fitted
#'   model, there is a standard deviation of 8\% in the errors of
#'   distance.
#'
#' * $\sigma_y$ is estimated at `{r} sprintf("%.3f",
#'   mean(draws_4$sigma_y))`.  The fitted model fits the
#'   aggregate data (success rate as a function of distance) to an
#'   accuracy of `{r} sprintf("%.1f", mean(draws_4$sigma_y)*100)`
#'   percentage points.
#'
#' And now we graph:
#| label: fig-golf-fit-4
draws_4 <- fit_4$draws(format = "df")
sigma_angle_hat <- median(draws_4$sigma_angle)
sigma_distance_hat <- median(draws_4$sigma_distance)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02), 
     xaxs = "i", yaxs = "i", pch = 20, bty = "l", 
     xlab = "Distance from hole (feet)", 
     ylab = "Probability of success", 
     main = "Checking model fit", type = "n")
x_grid <- seq(R - r, 1.1 * max(golf_new$x), .01)
p_angle_grid <- (2 * pnorm(asin((R - r) / x_grid) / sigma_angle_hat) - 1)
p_distance_grid <- pnorm((distance_tolerance - overshot) / ((x_grid + overshot) * sigma_distance_hat)) -
  pnorm((-overshot) / ((x_grid + overshot) * sigma_distance_hat))
lines(c(0, R - r, x_grid), c(1, 1, p_angle_grid * p_distance_grid), col = "red")
points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")

#' We can go further and plot the residuals from this fit.  First we
#' augment the Stan model to compute residuals in the generated
#' quantities block.
#| label: golf_angle_distance_normal_with_resids.stan
#| results: hide
model_4_with_resids <- cmdstan_model("golf_angle_distance_normal_with_resids.stan")
fit_4_with_resids <- model_4_with_resids$sample(data = golf_new_data, refresh = 0)

#' Then we compute the posterior means of the residuals, $y_j/n_j -
#' p_j$, then plot these vs. distance:
#| label: fig-golf-res-4
posterior_mean_residual <- mean(as_draws_rvars(fit_4_with_resids$draws())$residual)
par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
plot(golf_new$x, posterior_mean_residual, xlim = c(0, 1.1 * max(golf_new$x)), 
     xaxs = "i", pch = 20, bty = "l", 
     xlab = "Distance from hole (feet)", 
     ylab = "y/n - fitted E(y/n)", 
     main = "Residuals from fitted model", type = "n")
abline(0, 0, col = "gray", lty = 2)
lines(golf_new$x, posterior_mean_residual)

#' The fit is good, and the residuals show no strong pattern, also
#' they are low in absolute value---the model predicts the success
#' rate to within half a percentage point at most distances,
#' suggesting not that the model is perfect but that there are no
#' clear problems given the current data.
#'
#' The above model fit, but we were bothered by the normal
#' approximation, not so much for these particular data but rather
#' because it was a sort of admission of failure to not be able to
#' directly use the binomial model.
#' 
#' # Binomial with errors in logit scale
#'
#' What we wanted to do was to keep the binomial model and then add
#' the error on the logistic scale or something like that, something
#' to keep the probabilities bounded between 0 and 1. We added an
#' error term on the logistic scale with a scale parameter, `sigma_eta`,
#' estimated from the data.
print_file("golf_angle_distance_binomial_with_logit_errors.stan")

#' We fit the model to the data:
#| label: golf_angle_distance_binomial_with_logit_errors.stan
#| results: hide
model_5 <- cmdstan_model("golf_angle_distance_binomial_with_logit_errors.stan")
fit_5 <- model_5$sample(data = golf_new_data, refresh = 0)

#' Here is the result:
print(fit_5, digits=3)

#' Unfortunately, when we try to fit this new model to our data, the
#' console fills up with warnings and the chains don’t mix. We try
#' initializing with Pathfinder.
#| results: hide
pth_5 <- model_5$pathfinder(data = golf_new_data, refresh = 0,
                            num_paths=20, max_lbfgs_iters=100)
fit_5 <- model_5$sample(data = golf_new_data, refresh = 0, init = pth_5)

#' Here is the result:
print(fit_5)
#' The sampling is slower than for earlier models, but the convergence diagnostic look just fine.

#' We graph the new data and the fitted model:
#| label: fig-golf-fit-5
draws_5 <- fit_5$draws(format = "df")
sigma_angle_hat <- median(draws_5$sigma_angle)
sigma_distance_hat <- median(draws_5$sigma_distance)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(0, 0, xlim=c(0, 1.1*max(golf_new$x)), ylim=c(0, 1.02),
     xaxs="i", yaxs="i", pch=20, bty="l",
     xlab="Distance from hole (feet)",
     ylab="Probability of success",
     main="Checking model fit", type="n")
x_grid <- seq(R-r, 1.1*max(golf_new$x), .01)
p_angle_grid <- (2*pnorm(asin((R-r)/x_grid) / sigma_angle_hat) - 1)
p_distance_grid <- pnorm((distance_tolerance - overshot) / ((x_grid + overshot)*sigma_distance_hat)) -
  pnorm((- overshot) / ((x_grid + overshot)*sigma_distance_hat))
lines(c(0, R-r, x_grid), c(1, 1, p_angle_grid*p_distance_grid), col="red")
points(golf_new$x, golf_new$y/golf_new$n, pch=20, col="red")

#' There are now problems with the fit in the low range of $x$.
#' Clearly the additive error in logit scale is not good.
#' 
#' # Binomial with proportional errors
#'
#' Instead of additive errors in logit scale, we next fit a
#' three-parameter model that scales all the probabilities down from
#' 1. Each observation has its own proportional error term. The key is
#' to make each element of the multiplier vector `(1- epsilon)`
#' positive and less than 1. This eliminates the problem with the
#' boundary and the need for the logit.  The prior distribution for
#' `epsilon` keeps the errors under control.
print_file("golf_angle_distance_binomial_with_proportional_errors.stan")

#' We fit the model to the data:
#| label: golf_angle_distance_binomial_with_proportional_errors.stan
#| results: hide
model_6 <- cmdstan_model("golf_angle_distance_binomial_with_proportional_errors.stan")
fit_6 <- model_6$sample(data = golf_new_data, refresh = 0)

#' Here is the result:
print(fit_6)

#' We graph the new data and the fitted model:
#| label: fig-golf-fit-6
draws_6 <- fit_6$draws(format = "df")
sigma_angle_hat <- median(draws_6$sigma_angle)
sigma_distance_hat <- median(draws_6$sigma_distance)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(0, 0, xlim=c(0, 1.1*max(golf_new$x)), ylim=c(0, 1.02),
     xaxs="i", yaxs="i", pch=20, bty="l",
     xlab="Distance from hole (feet)",
     ylab="Probability of success",
     main="Checking model fit", type="n")
x_grid <- seq(R-r, 1.1*max(golf_new$x), .01)
p_angle_grid <- (2*pnorm(asin((R-r)/x_grid) / sigma_angle_hat) - 1)
p_distance_grid <- pnorm((distance_tolerance - overshot) / ((x_grid + overshot)*sigma_distance_hat)) -
  pnorm((- overshot) / ((x_grid + overshot)*sigma_distance_hat))
lines(c(0, R-r, x_grid), c(1, 1, p_angle_grid*p_distance_grid), col="red")
points(golf_new$x, golf_new$y/golf_new$n, pch=20, col="red")

#| label: fig-golf-res-6
posterior_mean_residual <- mean(as_draws_rvars(fit_6$draws())$residual)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(golf_new$x, posterior_mean_residual, xlim=c(0, 1.1*max(golf_new$x)),
     xaxs="i", pch=20, bty="l",
     xlab="Distance from hole (feet)",
     ylab="y/n - fitted E(y/n)",
     main="Residuals from fitted model", type="n")
abline(0, 0, col="gray", lty=2)
lines(golf_new$x, posterior_mean_residual)

#' The model fit looks good. The residual plots shows a pattern for
#' short distances: for very short distances the model underestimates
#' and then overestimates the probabilities. This pattern could be
#' explained by sensitivity to `distance_tolerance` and `overshot`
#' parameters that were fixed.
#'
#' We change `distance_tolerance` to be a parameter. We assume
#' @Broadie:2018 did use his expertise to choose the value 3 feet. We
#' use log-normal prior with mean log(3) and standard deviation 0.3 to
#' include that expert information but still have the prior to be
#' relatively weak.
#| label: golf_angle_distance_binomial_with_proportional_errors_2.stan
#| results: hide
model_7 <- cmdstan_model("golf_angle_distance_binomial_with_proportional_errors_2.stan")
fit_7 <- model_7$sample(data = golf_new_data, refresh = 0)

#' Here is the result:
print(fit_7)

#' We graph the new data and the fitted model:
#| label: fig-golf-fit-7
draws_7 <- fit_7$draws(format = "df")
sigma_angle_hat <- median(draws_7$sigma_angle)
sigma_distance_hat <- median(draws_7$sigma_distance)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(0, 0, xlim=c(0, 1.1*max(golf_new$x)), ylim=c(0, 1.02),
     xaxs="i", yaxs="i", pch=20, bty="l",
     xlab="Distance from hole (feet)",
     ylab="Probability of success",
     main="Checking model fit", type="n")
x_grid <- seq(R-r, 1.1*max(golf_new$x), .01)
p_angle_grid <- (2*pnorm(asin((R-r)/x_grid) / sigma_angle_hat) - 1)
p_distance_grid <- pnorm((distance_tolerance - overshot) / ((x_grid + overshot)*sigma_distance_hat)) -
  pnorm((- overshot) / ((x_grid + overshot)*sigma_distance_hat))
lines(c(0, R-r, x_grid), c(1, 1, p_angle_grid*p_distance_grid), col="red")
points(golf_new$x, golf_new$y/golf_new$n, pch=20, col="red")

#| label: fig-golf-res-7
posterior_mean_residual <- mean(as_draws_rvars(fit_7$draws())$residual)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(golf_new$x, posterior_mean_residual, xlim=c(0, 1.1*max(golf_new$x)),
     xaxs="i", pch=20, bty="l",
     xlab="Distance from hole (feet)",
     ylab="y/n - fitted E(y/n)",
     main="Residuals from fitted model", type="n")
abline(0, 0, col="gray", lty=2)
lines(golf_new$x, posterior_mean_residual)

#' The residuals are now smaller with standard deviation halved, and
#' there is no obvious pattern anymore. Looking at the posterior of
#' `distance_tolerance` most of the posterior mass is above 3 and the
#' posterior is much narrower than the prior and thus likelihood is
#' informative about it.
fit_7$summary(variables = c("distance_tolerance"))

#' For the final model we change `overshot` to be a parameter, too. We assume
#' @Broadie:2018 did use his expertise to choose the value 1 feet. We
#' use log-normal prior with mean log(3) and standard deviation 0.3 to
#' include that expert information but still have the prior to be
#' relatively weak.
#| label: golf_angle_distance_binomial_with_proportional_errors_3.stan
#| results: hide
model_8 <- cmdstan_model("golf_angle_distance_binomial_with_proportional_errors_3.stan")
fit_8 <- model_8$sample(data = golf_new_data, refresh = 0)

#' Here is the result:
print(fit_8)

#' We graph the new data and the fitted model:
#| label: fig-golf-fit-8
draws_8 <- fit_8$draws(format = "df")
sigma_angle_hat <- median(draws_8$sigma_angle)
sigma_distance_hat <- median(draws_8$sigma_distance)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(0, 0, xlim=c(0, 1.1*max(golf_new$x)), ylim=c(0, 1.02),
     xaxs="i", yaxs="i", pch=20, bty="l",
     xlab="Distance from hole (feet)",
     ylab="Probability of success",
     main="Checking model fit", type="n")
x_grid <- seq(R-r, 1.1*max(golf_new$x), .01)
p_angle_grid <- (2*pnorm(asin((R-r)/x_grid) / sigma_angle_hat) - 1)
p_distance_grid <- pnorm((distance_tolerance - overshot) / ((x_grid + overshot)*sigma_distance_hat)) -
  pnorm((- overshot) / ((x_grid + overshot)*sigma_distance_hat))
lines(c(0, R-r, x_grid), c(1, 1, p_angle_grid*p_distance_grid), col="red")
points(golf_new$x, golf_new$y/golf_new$n, pch=20, col="red")

#| label: fig-golf-res-8
posterior_mean_residual <- mean(as_draws_rvars(fit_8$draws())$residual)
par(mar=c(3,3,2,1), mgp=c(1.7,.5,0), tck=-.02)
plot(golf_new$x, posterior_mean_residual, xlim=c(0, 1.1*max(golf_new$x)),
     xaxs="i", pch=20, bty="l",
     xlab="Distance from hole (feet)",
     ylab="y/n - fitted E(y/n)",
     main="Residuals from fitted model", type="n")
abline(0, 0, col="gray", lty=2)
lines(golf_new$x, posterior_mean_residual)

#' The residuals are very similar with the previous model and with
#' similar standard deviation.  If we examine the posterior marginals
#' of `distance_tolerance` and `overshot` the standard deviations are
#' only half from the prior standard deviations, hinting that
#' likelihood would be weakly informative on these.
fit_8$summary(variables = c("distance_tolerance","overshot"))

#' However, when we examine the bivariate posterior we see strong
#' dependency and their ratio is well informed by the likelihood. The
#' posterior standard deviation of the ratio is one fifth of the prior
#' standard deviation. This also explains why adding `overshot` did
#' not reduce much the residual standard deviation.
#' 
#' As we can expect the residual standard deviation to decrease when
#' we add more parameters, we use also cross-validation to compare the
#' models. As we have one `epsilon` parameter for each observation, we
#' need to use integerated PSIS-LOO. We can compute the integrated
#' `log_lik` with the following stand alone generated quantities code.
gq_ll <- cmdstan_model("golf_log_lik.stan")

#' When calling generated quantities, we pass only the required variables which
#' allowed us to write more compact Stan code for the `log_lik` computation.
#| results: hide
loo_6 <- gq_ll$generate_quantities(fit_6$draws(variables = c("sigma_epsilon",
                                                             "p_angle",
                                                             "p_distance")),
                                   data = golf_new_data)$draws(variables = "log_lik") |>
                                                       loo::loo()
loo_7 <- gq_ll$generate_quantities(fit_7$draws(variables = c("sigma_epsilon",
                                                             "p_angle",
                                                             "p_distance")),
                                   data = golf_new_data)$draws(variables = "log_lik") |>
                                                       loo::loo()
loo_8 <- gq_ll$generate_quantities(fit_8$draws(variables = c("sigma_epsilon",
                                                             "p_angle",
                                                             "p_distance")),
                                   data = golf_new_data)$draws(variables = "log_lik") |>
                                                       loo::loo()

#' As we have integrated out the `epsilon` parameters, the effective
#' number of parameters match the number of remaining parameters,
#' except for the last model where `distance_tolerance` and `overshot`
#' parameters are not well informaed by the likelihood separately.
print(loo_6)
print(loo_7)
print(loo_8)

#' Finally we compare the expected predictive performances
loo::loo_compare(list(`Distance tolerance and overshot fixed`=loo_6,
                      `Distance tolerance parameter and overshot fixed`=loo_7,
                      `Distance tolerance and overshot parameters`=loo_8))

#' Adding `distance_tolerance` parameter significantly improves the
#' performance, but adding `overshot` parameter does not improve the
#' fit. However the model which has both `distance_tolerance` and
#' `overshot` has posterior that is more informative on what can be
#' learned about this aspects of the model.
#' 
#' <br />
#' 
#' # References {.unnumbered}
#'
#' <div id="refs"></div>
#' 
#' # Licenses {.unnumbered}
#' 
#' * Code &copy; 2019-2025, Andrew Gelman and Aki Vehtari, licensed under BSD-3.
#' * Text &copy; 2019-2025, Andrew Gelman and Aki Vehtari, licensed under CC-BY-NC 4.0.

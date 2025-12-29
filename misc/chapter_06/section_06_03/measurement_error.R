library("arm")

set.seed(123)
n <- 1000
x <- runif(n, 0, 10)
a <- 0.2
b <- 0.3
sigma <- 0.5
y <- rnorm(n, a + b * x, sigma)
fake <- data.frame(x, y)

fit_1 <- lm(y ~ x, data = fake)
display(fit_1)


sigma_y <- 1
fake$y_star <- rnorm(n, fake$y, sigma_y)
sigma_x <- 4
fake$x_star <- rnorm(n, fake$x, sigma_x)


fit_2 <- lm(y_star ~ x, data=fake)
display(fit_2)
fit_3 <- lm(y ~ x_star, data=fake)
display(fit_3)
fit_4 <- lm(y_star ~ x_star, data=fake)
display(fit_4)
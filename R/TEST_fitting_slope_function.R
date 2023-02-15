
# approximating a function for the slope of the relationship
# between distance factor and distance 
# in the potential for apparent interactions functions

# fit a function of the form -beta*x^(-alpha)
# which is equivalent to a rational function -alpha/beta*x

library("tidyverse")

# -------------------------------------------------------------------------

# random points
set.seed(123)
N <- 1000
x <- runif(N,1e-4,.1)
# f_old <- function(x) 50*x^2/(1 + 4*x) # data-simulating function
f <- function(x) -.1*x^(-.0001)
y <- f(x) + rnorm(N, sd=1e-8)

# x <- c(1e-4,1e-1)
# y <- c(-1,-1/100)

point_data <- data.frame(x, y)

ggplot(point_data, aes(x=x, y=y)) + 
  geom_point() + 
  # ylim(-.1,0) +
  ggtitle("data points")

fit.fun <- function(beta, alpha,x) (-beta*x^(-alpha))

# fit nls
nlc <- nls.control(maxiter = 1000)
my.fit <- nls(y ~ fit.fun(mybeta,myalpha,x), control=nlc, data = point_data, start=list(myalpha=1,mybeta=1))
# summary(my.fit)

# -------------------------------------------------------------------------
# plot
ggplot(point_data, aes(x=x, y=y)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method = "nls", 
              formula = y ~ beta * x^(-alpha), 
              se =  FALSE, # this is important 
              method.args = list(start = list(beta = 1, alpha = 1))) + 
  NULL


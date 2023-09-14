library(RTMB)

set.seed(123)
b0 <- 0.3
b1 <- 2
N <- 200L
x <- sort(runif(N))
mu <- exp(b0 + b1 * x)
obs <- rpois(N, lambda = mu)
plot(x, obs)
lines(x, mu)

data <- list(y = obs, x = x)
parameters <- list(b0 = 0, b1 = 0)

f <- function(parms) {
  getAll(data, parms)
  mu <- exp(b0 + b1 * x)
  -sum(RTMB::dpois(y, lambda = mu, log = TRUE))
}

obj <- MakeADFun(f, parameters)
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$convergence
sdrep <- sdreport(obj)
summary(sdrep)


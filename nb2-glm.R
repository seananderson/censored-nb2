library(RTMB)

set.seed(123)
b0 <- 0.3
b1 <- 2
ln_phi <- 0.5
N <- 200L
x <- sort(runif(N))
mu <- exp(b0 + b1 * x)
obs <- rnbinom(N, size = exp(ln_phi), mu = mu)
plot(x, obs)
lines(x, mu)

data <- list(y = obs, x = x)
parameters <- list(b0 = 0, b1 = 0, ln_phi = 0)

f <- function(parms) {
  getAll(data, parms)
  mu <- exp(b0 + b1 * x)
  s1 <- log(mu)
  s2 <- 2 * s1 - ln_phi
  -sum(dnbinom_robust(y, s1, s2, log = TRUE))
}

obj <- MakeADFun(f, parameters)
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$convergence
sdrep <- sdreport(obj)
summary(sdrep)

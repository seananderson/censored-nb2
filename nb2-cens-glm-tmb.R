library(TMB)

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

cens <- rbinom(N, 1, 0.1)
cens
data <- list(y = obs, x = x, cens = cens)
# data <- list(y = obs, x = x, cens = rep(0, N)) # test
parameters <- list(b0 = 0.3, b1 = 2, ln_phi = 0.5)

compile("nb2cens.cpp")
dyn.load(dynlib("nb2cens"))
obj <- MakeADFun(data, parameters, DLL = "nb2cens")
obj$fn(obj$par)
obj$gr(obj$par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
sdr

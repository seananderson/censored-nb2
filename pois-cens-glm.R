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

logspace_sub <- function(log_a, log_b) {
  if (log_a < log_b) {
    return(log_b + log1p(-exp(log_a - log_b)))
  } else {
    return(log_a + log1p(-exp(log_b - log_a)))
  }
}

f1 <- function(x) {
  log_b <- x[1]
  log_a <- x[2]
  if (log_a < log_b) {
    return(log_b + log1p(-exp(log_a - log_b)))
  } else {
    return(log_a + log1p(-exp(log_b - log_a)))
  }
}
f1(c(2, 3))

RTMB::TapeConfig(
  comparison = "allow",
)
MakeTape(f1, numeric(2))

f2 <- function(x) {
  logx <- x[1]
  logy <- x[2]
  logx + log1p(-exp(logy - logx))
}

f2(c(3, 2))


MakeTape(f2, numeric(2))


dcenspois_right <- function(x, lambda, log = FALSE) {
  # ll <- RTMB::ppois(x - 1, lambda) # F(lower-1)
  # ll <- stats::ppois(as.numeric(x - 1), as.numeric(lambda))) # F(lower-1)
  print(ll)
  # ll <- 2
  ll <- logspace_sub(0, ll) # 1 - F(lower-1)
  if (log) (ll) else exp(ll)
}

f <- function(parms) {
  getAll(data, parms)
  # b0 <- parms$b0
  # b1 <- parms$b1
  # y <- data$y

  mu <- exp(b0 + b1 * x)
  nll <- 0
  for (i in seq_along(y)) {
    nll <- nll - stats::dpois(y[i], lambda = a, log = TRUE)
    # nll <- nll - dpois(y[i], lambda = mu[i], log = TRUE)
    # nll <- nll - dcenspois_right(y[i], lambda = mu[i], log = TRUE)
  }
  nll
}

RTMB::TapeConfig(comparison = "forbid")
obj <- MakeADFun(f, parameters)
obj$fn(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$convergence
sdrep <- sdreport(obj)
summary(sdrep)


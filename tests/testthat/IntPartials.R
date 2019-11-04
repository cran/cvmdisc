intpartial.exp.a <- function(x, lambda) {
  eval(D(expression(
    lambda * exp(- lambda * x)), "lambda"))
}

intpartial.lnorm.a <- function(x, mu, sigma) {
  eval(D(expression(
    1 / (sqrt(2 * pi) * sigma * x) * exp(-(log(x) - mu) ^ 2 / (2 * sigma ^ 2))),
    "mu"))
}

intpartial.lnorm.b <- function(x, mu, sigma) {
  eval(D(expression(
    1 / (sqrt(2 * pi) * sigma * x) * exp(-(log(x) - mu) ^ 2 / (2 * sigma ^ 2))),
    "sigma"))
}

intpartial.weibull.a <- function(x, a, b) {
  eval(D(expression(
    (a / b) * (x / b) ^ (a - 1) * exp(-(x / b) ^ a)), "a"))
}

intpartial.weibull.b <- function(x, a, b) {
  eval(D(expression(
    (a / b) * (x / b) ^ (a - 1) * exp(-(x / b) ^ a)), "b"))
}

intpartial.gamma.a <- function(x, a, r){
  #a=shape
  #r=rate the default output of pgamma is rate.
  eval(D(expression(
    x ^ (a - 1) * exp(-x * r) * r ^ (a) / gamma(a)), "a"))
}

intpartial.gamma.b <- function(x, a, r) {
  #a=shape
  #r=rate the default output of pgamma is rate.
  eval(D(expression(
    x ^ (a - 1) * exp(-x * r) * r ^ (a) / gamma(a)), "r"))
}

intpartial.norm.a <- function(x, mu, sigma) {
  eval(D(expression(
    1 / (sqrt(2 * pi) * sigma) * exp(-(x - mu) ^ 2 / (2 * sigma ^ 2))),
    "mu"))
}

intpartial.norm.b <- function(x, mu, sigma) {
  eval(D(expression(
    1 / (sqrt(2 * pi) * sigma) * exp(-(x - mu) ^ 2 / (2 * sigma ^ 2))),
    "sigma"))
}

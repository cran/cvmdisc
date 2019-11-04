source("groupFit1.R")

test_that("Asymptotic", {
  #exponential example (Spinelli 2001)
  breaks <- c(0, 2, 6, 10, 14, 18, 22, 26)
  counts <- c(21, 9, 5, 2, 1, 1, 0)

  for (distr in c("unif", "exp", "gamma", "lnorm", "norm", "weibull")) {
    test <- groupFit(breaks, counts, distr = distr)
    expect_true(is.list(test))
  }

  #Poisson Example (Spinelli 1994) p36
  counts <- c(9, 22, 6, 2, 1, 0, 0, 0, 0)
  vals <- 0 : 8
  data <- rep(vals, counts)

  test <- groupFit(data = data, distr = "pois")
  expect_true(is.list(test))

  #Spinelli (1994) P92.
  N=12
  counts= c(185, 1149, 3265, 5475, 6114, 5194, 3067, 1331, 403, 105, 14, 4, 0)
  vals=0:12
  data<- rep(vals, counts)

  test <- groupFit(data = data, distr = "binom", N = N)
  expect_true(is.list(test))
})

test_that("Bootstrap", {
  breaks <- c(0, 2, 6, 10, 14, 18, 22, 26)
  counts <- c(21, 9, 5, 2, 1, 1, 0)

  for (distr in c("unif", "exp", "gamma", "lnorm", "norm", "weibull")) {
    test <- groupFit(breaks, counts, distr = distr, bootstrap = TRUE, numLoops = 100)
    expect_true(is.list(test))
  }

  counts <- c(9, 22, 6, 2, 1, 0, 0, 0, 0)
  vals <- 0 : 8
  data <- rep(vals, counts)

  test <- groupFit(data = data, distr = "pois", bootstrap = TRUE, numLoops = 100)
  expect_true(is.list(test))

  N=12
  counts= c(185, 1149, 3265, 5475, 6114, 5194, 3067, 1331, 403, 105, 14, 4, 0)
  vals=0:12
  data<- rep(vals, counts)

  test <- groupFit(data = data, distr = "binom", N = N, bootstrap = TRUE, numLoops = 100)
  expect_true(is.list(test))
})

test_that("Errors", {
  counts <- c(9, 22, 6, 2, 1, 0, 0, 0, 0)
  vals <- 0 : 8
  data <- rep(vals, counts)

  expect_error(groupFit(data = data, discrete = 18, distr = "pois"))
  expect_error(groupFit(data = data, discrete = TRUE, distr = not_pois))
  expect_error(groupFit(data = data, distr = "pois", known = TRUE))
  expect_error(groupFit(data = data, distr = "binom"))
  expect_error(groupFit(data = data, distr = "binom", N = NULL))
  expect_error(groupFit(data = data, distr = "binom", N = 0.5))
  expect_error(groupFit(data = data, distr = "binom", N = 2))

  counts <- c(21, 9, 5, 2, 1, 1, 0)

  breaks <- c(8, 7, 6, 5, 4, 3, 2, 1)
  expect_error(groupFit(breaks, counts, distr = "gamma"))

  breaks <- c(1, 2, 3, 4, 4, 5, 6, 7)
  expect_error(groupFit(breaks, counts, distr = "gamma"))

  breaks <- c(0, 2, 6, 10, 14, 18, 22, 26)
  counts <- c(1, 2, 3)
  expect_error(groupFit(breaks, counts, distr = "gamma"))

  counts <- c(-1, 2, 3, 4, -5, 6, 7)
  expect_error(groupFit(breaks, counts, distr = "gamma"))

  counts <- c(0, 0, 0, 0, 0, 0, 1)
  expect_error(groupFit(breaks, counts, distr = "gamma"))
})

test_that("Other", {
  breaks <- c(0, 2, 6, 10, 14, 18, 22, 26)
  counts <- c(21, 9, 5, 2, 1, 1, 0)

  test <- groupFit(breaks, counts, distr = "exp",
                   initials = 0.26)
  expect_true(is.list(test))
  test <- groupFit(breaks, counts, distr = "exp",
                   params = 0.26, known = TRUE)
  expect_true(is.list(test))
  test <- groupFit(breaks, counts, distr = "exp",
                   pave = TRUE)
  expect_true(is.list(test))
  test <- groupFit(breaks, counts, distr = "exp",
                   params = 0.26, known = TRUE, pave = TRUE)
  expect_true(is.list(test))
  test <- groupFit(breaks, counts, distr = "exp",
                   imhof = TRUE)
  expect_true(is.list(test))
  test <- groupFit(breaks, counts, distr = "exp",
                   params = 0.26, known = TRUE, imhof = TRUE)
  expect_true(is.list(test))
  test <- groupFit(breaks, counts, distr = "exp",
                   pave = TRUE, imhof = TRUE)
  expect_true(is.list(test))
  test <- groupFit(breaks, counts, distr = "exp",
                   params = 0.26, known = TRUE, pave = TRUE, imhof = TRUE)
  expect_true(is.list(test))
})

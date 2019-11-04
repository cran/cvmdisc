groupFit <- function(breaks, counts, data, discrete, distr,
                     N, params, initials, bootstrap = FALSE,  numLoops = 5000,
                     known = FALSE,  T = FALSE, imhof = TRUE) {
  names <- c("A-squared", "W-squared", "U-squared",  "Chi-squared")
  distr <- tolower(distr)
  if (distr == "binomial") {
    distr <- "binom"
  } else if (distr == "poisson") {
    distr <- "pois"
  } else if (distr == "exponential") {
    distr <- "exp"
  } else if (is.element(distr, c("l norm", "l-norm", "lognormal", "log normal", "log-normal"))) {
    distr <- "lnorm"
  } else if (distr == "normal") {
    distr <- "norm"
  } else if (distr == "uniform") {
    distr <- "unif"
  }
  if (distr == "unif") {
    known <- TRUE
    discrete <- FALSE
  }
  if (missing(discrete)) {
    if (is.element(distr, c("binom", "pois")))
      discrete <- TRUE
    else discrete <- FALSE
  } else {
    if (!is.logical(discrete))
      stop("discrete must be logical. (TRUE or FALSE)")
  }

  if (discrete == TRUE) {
    if (!(is.element(distr, c("binom", "pois")))) {
      stop(paste("Discrete function ", distr, " is not currently supported."))
    }
    n <- length(data)
    plim <- 1 / (n * 10 ^ 3)
    if (distr == "pois") {
      if (known == TRUE) {
        if (missing(params)) {
          stop("Parameters are known but not provided")
        }
        hats <- params
      } else {
        hats <- mean(data)
      }
      phat <- dpois(0:(2*max(data)), hats)
      k <- max(which(phat > plim))
      phat <- phat[1 : k]
    } else if (distr == "binom") {
      if (!(is.numeric(N))) {
        stop(print("Parameter 'N' must be numeric."))
      }
      if (!(N == round(N)) | N < 1) {
        stop(print("Paraneter 'N' must be a positive integer."))
      }
      if (N < max(data)) {
        stop("The number of trials should be greater than or equal
             to the greatest number of successes.")
      }
      if (known == TRUE) {
        if (missing(params)) {
          stop("Parameters are known but not provided")
        }
        hats <- params
      } else {
        hats <- sum(data) / (N * n)
      }
      k <- N + 1
      phat <- dbinom(0:N, N, hats)
    }
    breaks <- seq(-0.5, (k - 0.5), 1)
    counts <- table(cut(data, breaks))
  } else {
    if (is.unsorted(breaks, strictly = TRUE)) {
      if (is.unsorted(breaks)) {
        stop(paste("Vector 'breaks' is not in increasing order."))
      } else {
        stop(paste("Vector 'breaks' has duplicate values."))
      }
    }
    if (!(length(counts) == length(breaks)-1)) {
      stop(paste("The length of vector 'counts' must be one less
                 than the length of vector 'breaks'."))
    }
    if (length(counts) != length(counts[counts >= 0])) {
      stop(paste("Vector 'counts' must contain non-negative values only."))
    }
    k <- length(counts)
    n <- sum(counts)
    plim <- 1 / (n * 10 ^ 3)
    phat <- rep(0, k)
    if (distr == "unif") {
      phat <- diff(breaks) / (breaks[ (k + 1) ] - breaks[1])
    } else {
      breaks[length(breaks)] <- Inf
      if (is.element(distr, c("exp", "gamma", "lnorm", "weibull"))) {
        breaks[1] <- 1e-7
      } else if (distr == "norm") {
        breaks[1] <- -Inf
      }
      pdistr <- paste("p", distr, sep = "")
      if (!exists(pdistr)) {
        stop(paste("Function '", pdistr, "' is not defined.", sep = ""))
      }
      if ((length(counts[counts > 0]) < 1 && distr == "exp")
          || length(counts[counts > 0]) < 2 && !(distr == "exp")) {
        stop("Length of non-zero counts must be greater
             than the degrees of freedom of the distribution")
      }
      if (known == TRUE) {
        if (missing(params)) {
          stop("Parameters are known but not provided")
        }
        hats <- params
      } else {
        if (!missing(initials)) {
          hats <- distrFit(breaks, counts, distr, initials)
        } else {
          hats <- distrFit(breaks, counts, distr)
        }
      }
      phat <- diff(do.call(pdistr, c(list(breaks), as.list(as.numeric(hats)))))
    }
  }
  k_max <- max(which(phat > plim))
  phat <- phat[phat > plim]
  k_trunc <- length(phat)
  k_min <- k_max - k_trunc
  counts <- counts[(k_min + 1) : k_max]
  stats <- cvmTest(counts, phat, T)
  if (bootstrap == FALSE) {
    if (any(n*phat < 5)) {
      warning("Chi-squared approximation may be incorrect.")
    }
    Hj <- stats$Hj
    k1 <- length(Hj)
    I <- diag(k_trunc) 
    A <- matrix(rep(1, k_trunc ^ 2), k_trunc)
    for (i in 1 : k_trunc) {
      for (j in 1 : k_trunc) {
        if (j > i) {
          A[i, j] <- 0
        }
      }
    }

    invA <- diag(k_trunc) - rbind(rep(0, k_trunc),
                                  cbind(diag( (k_trunc-1) ),
                                        rep(0, (k_trunc - 1))))
    D <- diag(phat)
    invD <- diag(1 / phat)
    one <- rep(1, k_trunc)
    tt <- I - D %*% one %*% t(one)
    Kmat <- diag(1 / (Hj * (1 - Hj)))

    if (known == FALSE) {
      if (discrete == FALSE){
        derp.alpha <- rep(0, k_trunc)
        derp.beta <- rep(0, k_trunc)
        intpartial.a <- paste("intpartial.", distr, ".a", sep = "")
        if (!exists(intpartial.a)) {
          stop(paste("Function '", intpartial.a,
                     "' is not defined.", sep = ""))
        }
        if (distr == "exp") {
          for (i in 1 : k_trunc) {
            derp.alpha[i] <- integrate(get(intpartial.a),
                                      lower = breaks[k_min + i],
                                      upper = breaks[k_min + i + 1],
                                      hats[1])$value
          }
          bhat <- derp.alpha
        } else {
          intpartial.b <- paste("intpartial.", distr, ".b", sep = "")
          if (!exists(intpartial.b)) {
            stop(paste("Function '", intpartial.b,
                       "' is not defined.", sep = ""))
          }
          for (i in 1 : k_trunc) {
            derp.alpha[i] <- integrate(get(intpartial.a),
                                       lower = breaks[k_min + i],
                                       upper = breaks[k_min + i + 1],
                                       hats[1], hats[2])$value
            derp.beta[i] <- integrate(get(intpartial.b),
                                      lower = breaks[k_min + i],
                                      upper = breaks[k_min + i + 1],
                                      hats[1], hats[2])$value
          }
          bhat <- matrix(c(derp.alpha, derp.beta), k_trunc, 2)
        }
        vhat <- solve(t(bhat) %*% invD %*% bhat)
      }
    }
    SigO <- D-phat %*% t(phat)

    if (known == TRUE) {
      SigU = A %*% SigO %*% t(A)
    } else {
      if (distr == "pois") {
        gj <- phat / hats * (0 : (k_trunc - 1) - hats)
        SigD <- SigO - hats * (gj %*% t(gj))
      } else if (distr == "binom") {
        gj <- ((0 : (k_trunc - 1)) - N * hats) * phat / (hats * (1 - hats))
        SigD <- SigO - (hats * (1 - hats)) / N * (gj %*% t(gj))
      } else {
        SigD <- SigO - bhat %*% vhat %*% t(bhat)
      }
      SigU <- A %*% SigD %*% t(A)
    }
    if (T == TRUE) {
      phat1 <- (phat + c(phat[-1], phat[1])) / 2
      D <- diag(phat1)
    }
    MsigW2 <- D %*% SigU
    MsigU2 <- (I - D %*% one %*% t(one)) %*% D %*%
              (I - one %*% t(one) %*% D) %*% SigU
    MsigA2 <- t((D[1 : k1, 1 : k1] %*% Kmat) ^ 0.5) %*%
              SigU[1 : k1, 1 : k1] %*% (D[1 : k1, 1 : k1] %*% Kmat) ^ 0.5

    A2pval <- cvmPval(stats$Asq, MsigA2, imhof)
    W2pval <- cvmPval(stats$Wsq, MsigW2, imhof)
    U2pval <- cvmPval(stats$Usq, MsigU2, imhof)
    if (known == TRUE) {
      paraEst <- 0
    } else {
      paraEst <- length(hats)
    }
    X2pval <- 1 - pchisq(stats$Chisq, k - paraEst - 1)
    pvals <- c(A2pval, W2pval, U2pval, X2pval)
    pvals <- as.data.frame(pvals, row.names = names)
  } else if (bootstrap == TRUE) {
    if (discrete == TRUE) {
      if (distr == "pois") {
        bootdata <- matrix(rpois(n * numLoops, hats), ncol = numLoops)
        if (known == FALSE) {
          boothats <- colMeans(bootdata)
          N <- qpois(1 - plim, max(boothats))
          bootphats <- outer(0 : N, boothats, FUN = dpois)
        } else {
          N <- qpois(1 - plim, hats)
        }
      } else if (distr == "binom") {
        bootdata <- matrix(rbinom(n * numLoops, N, hats), ncol = numLoops)
        if (known == FALSE) {
          boothats <- colSums(bootdata) / (N * n)
          bootphats <- outer(0 : N, boothats,
                             FUN = function(x, y)dbinom(x, prob = y,
                                                        size = N))
        }
      }
      bootbreaks <- seq(-0.5, (N + 0.5), 1)
      bootcounts <- apply(bootdata, 2,
                          FUN = function(x)table(cut(x, bootbreaks)))
      numEff <- numLoops
    } else {
      bootcounts <- rmultinom(numLoops, n, prob = phat)
      bootbreaks <- breaks[(k_min + 1) : (k_max + 1)]
      if (distr != "unif") {
        smallcounts <- which((apply(bootcounts, 2,
                                    function(x)length(x[x != 0])))
                             < (length(hats) + 1))
        numS <- length(smallcounts)
        numEff <- numLoops- numS

        if (numS != 0) {
          bootcounts <- bootcounts[, -smallcounts]
        }
      } else {
        numEff <- numLoops
      }
      if (known == FALSE) {
        boothats <- matrix(apply(bootcounts, 2,
                                 function(y)distrFit(bootbreaks,
                                                     y, distr,
                                                     initials = hats)),
                           ncol = numEff)

        bootphats<- apply(boothats, 2,
                          function(x)diff(do.call(pdistr,
                                                  c(list(bootbreaks),
                                                    as.list(x)))))
      }
    }
    if (known == TRUE){
      bootstats <- matrix(unlist(apply(bootcounts, 2,
                                       function(x)cvmTest(x, phat, T))),
                          ncol = numEff)
    } else if (known == FALSE) {
      bootstats <- mapply(function(x, y) {cvmTest(x, y, T)},
                          as.data.frame(bootcounts), as.data.frame(bootphats))
    }
 
    pvals <- mapply(function(x, y){length(x[x > y]) / numEff},
                    as.data.frame(t(bootstats[1 : 4, ])),
                    as.data.frame(stats[1 : 4]))
    pvals <- as.data.frame(pvals, row.names = names)
  }
  test_statistics <- unlist(stats[1 : 4])
  test_statistics <- as.data.frame(test_statistics, row.names = names)

  if (distr == "unif") {
    reslist <- list(stats = test_statistics, pvals = pvals)
  } else if (known == FALSE) {
    reslist <- list(estimates = hats, stats = test_statistics,
                    pvals = pvals)
  } else if (known == TRUE) {
    reslist <- list("Known Parameter" = hats, stats = test_statistics,
                    pvals = pvals)
  }
  return(reslist)
}
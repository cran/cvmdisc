#' groupFit
#'
#' Fits grouped continuous data or discrete data to a distribution and computes Cramer-Von Mises Goodness of Fit statistics and p-values.
#'
#' @param breaks Vector defining the breaks in each group.
#' @param counts Vector containing the frequency of counts in each group.
#' @param data Vector containing values of discrete random variables.
#' @param discrete Logical. Is the distribution discrete?
#' @param distr Character; the name of the distribution users want to fit the data to.
#'
#'    Included continuous distributions are: "exp", "gamma", "lnorm", "norm", "unif" and "weibull".
#'
#'    Included discrete distributions are: "binom" and "pois".
#'
#'    User defined distributions are supported as "user".
#'
#'    Short-hand or full spelling of these distributions will be recognized, case insensitive.
#'
#' @param N Number of trials, used only for the binomial distribution.
#' @param params Vector of distribution parameters. This is only required when known == TRUE. groupFit will estimate the parameters if known == FALSE.
#' @param initials Vector of distribution parameters to use as starting points for calculating MLEs.
#' @param pfixed Vector of known probabilities for corresponding counts vector. pfixed must be provided when distr = "user".
#' @param bootstrap Logical. Should p-values be calculated via bootstrapping?
#' @param numLoops Number of Bootstrap iterations. Set to be 5000 by default.
#' @param pave Logical. Set to be FALSE if the probabilities in groups are used; set to TRUE if the average of probabilities of groups j and j + 1 are used. See page 2 of Spinelli (2001) for more details.
#' @param known Logical. Set to be TRUE if the parameters are known and do not need to be estimated.
#' @param imhof Logical. Set to be TRUE if Imhof's method from the package "CompQuadForm" is used to approximate the null distributions.
#'
#'
#' @details For grouped continuous data: call groupFit with arguments breaks, counts, and distr to fit the data.
#'
#'    For discrete data call groupFit with data and distr to fit the data. If distr = "binom", then be sure to call groupFit with N as well.
#'
#'    Provide initials to suggest starting points for parameter estimation.
#'
#'    Provide params and set known = TRUE to test goodness of fit when parameters are known.
#'
#'    Set bootstrap = TRUE to use bootstrapping to estimate p-values rather than using asymptotic distributions.
#'
#'    Set imhof = FALSE when \eqn{a+bX^2_p} is used to approximate the null distributions. See page 5 of Spinelli (2001) for details.
#'
#'    groupFit can test the fit of user defined discrete distributions. To do so, set distr to "user", and provide the vector pfixed, where each cell contains the probability corresponding to that same cell in counts.
#'
#' @return List containing the components:
#'
#' \item{estimates}{The estimated parameters of the distribution}
#'
#' \item{stats}{Data frame containing the goodness of fit statistics}
#'
#' \item{pvals}{Data frame containing p-values for the goodness of fit statistics}
#'
#' @author Shaun Zheng Sun and Dillon Duncan
#'
#' @seealso \code{\link{distrFit}}: Parameter estimation function
#'
#' @examples
#'
#' #Poisson Example (Spinelli 1994) p36
#'
#' counts <- c(9, 22, 6, 2, 1, 0, 0, 0, 0)
#' vals <- 0:8
#' data <- rep(vals, counts)
#'
#' groupFit(data = data, distr = "pois")
#'
#'
#' # When the parameters are unknown
#' #(Spinelli 1994) p56
#'
#' counts <- c(57, 203, 383, 525, 532, 408, 273, 139, 45, 27, 10, 4, 0, 1, 1)
#' vals <- 0:14
#' data <- rep(vals, counts)
#'
#' (pois_fit <- groupFit(data = data, distr = "pois"))
#'
#' #Binomial example when the parameter is unknown
#' #Spinelli (1994) P92.
#' N=12
#' counts= c(185, 1149, 3265, 5475, 6114, 5194,
#'           3067, 1331, 403, 105, 14, 4, 0)
#' vals <- 0:12
#' data <- rep(vals, counts)
#'
#' (binom_fit <- groupFit(data = data, N = N, distr = "binom"))
#'
#' #When the parameter is assumed known and is equal to 1/3
#'
#' (binom_fit <- groupFit(data = data, N = N, distr = "binom", params = 1/3, known = TRUE))
#'
#' #uniform example (Choulakian, Lockhart and Stephens(1994) Example 2, p8)
#'
#' counts <- c(10, 19, 18, 15, 11, 13, 7, 10, 13, 23, 15, 22)
#'
#' (uni_fit <- groupFit(0:12, counts, distr = "unif"))
#'
#'
#' #uniform example (Choulakian, Lockhart and Stephens(1994) Example 3, p8)
#' counts <- c(1, 4, 11, 4, 0)
#' probability <- c(0.05, 0.3, 0.3, 0.3, 0.05)
#' breaks <- c(0, cumsum(probability))
#'
#' #with bootstrapping
#' (uni_fit1 <- groupFit(breaks, counts, distr = "unif", bootstrap = TRUE, numLoops = 500))
#'
#' #without bootstrapping
#' (uni_fit2 <- groupFit(breaks, counts, distr = "unif", bootstrap = FALSE))
#'
#' #exponential example (Spinelli 2001)
#'
#' breaks <- c(0, 2, 6, 10, 14, 18, 22, 26)
#' counts <- c(21, 9, 5, 2, 1, 1, 0)
#'
#' (exp_fit <- groupFit(breaks, counts, distr = "exp", pave = TRUE))
#'
#'
#' #Example Sun, Stephens & Spinelli (2012) set 2.
#' breaks <- c(0, 2, 6, 10, 14, 18, 22, 26)
#' counts <- c(21, 9, 5, 2, 1, 1, 0)
#' breaks[1] <- 1e-6
#' breaks[8] <- Inf
#'
#' (weibull_fit <- groupFit(breaks, counts, distr = "Weibull"))
#'
#'
#' #Example Sun, Stephens & Spinelli (2012) set 3.
#'
#' breaks <- c(0, seq(0.5, 6.5, 1) )
#' counts <- c(32, 12, 3, 6, 0, 0, 1)
#' breaks[1] <- 1e-6
#' breaks[8] <- Inf
#'
#' (weibull_fit <- groupFit(breaks, counts, distr = "Weibull"))
#'
#'
#' #Example Sun, Stephens & Spinelli (2012) set 3.
#' breaks <- c(0, 2, 6, 10, 14, 18, 22, 26)
#' counts <- c(21, 9, 5, 2, 1, 1, 0)
#' breaks[1] <- 1e-6
#' breaks[8] <- Inf
#'
#' (gamma_fit <- groupFit(breaks, counts, distr = "exp"))
#'
#'
#' #Example Sun, Stephens & Spinelli (2012) set 3.
#' breaks <- c(0, seq(0.5, 6.5, 1) )
#' counts <- c(32, 12, 3, 6, 0, 0, 1)
#' breaks[1] <- 1e-6
#' breaks[8] <- Inf
#'
#' (gamma_fit <- groupFit(breaks, counts, distr = "gamma"))
#'
#'
#' #More examples
#'
#' breaks <- c(0, seq(0.5, 6.5, 1) )
#' counts <- table(cut(rgamma(100, 3, 1/3), breaks))
#' breaks[8] <- Inf
#'
#' #setting pave to true
#' (exp_fit <- groupFit(breaks, counts, distr = "exp", initials = 0.2, pave = TRUE))
#'
#' #Setting known to true, with params
#' (gamma_fit <- groupFit(breaks, counts, distr = "gamma",
#'                       params = c(3, 1/3), known = TRUE))
#'
#' #with bootstrapping, specifying the number of loops.
#' (lnorm_fit <- groupFit(breaks, counts, distr = "lnorm",
#'                       bootstrap = TRUE, numLoops = 1000))
#'
#' #fitting with both pave and imhof set to false
#' #by setting imhof to false, we use a+bX^2_p to approximate
#' #the distribution of the goodness-of-fit Statistics
#' (weibull_fit <- groupFit(breaks, counts, distr = "weibull",
#'                         pave = TRUE, imhof = FALSE))
#'
#' #Using the user defined distribution to test for Benford's law
#'
#' #genomic data, Lesperance et al (2016)
#'
#' genomic <- c(48, 14, 12, 6, 18, 5, 7, 8, 9)
#' phat <- log10(1+1/1:9)
#'
#' (fit <- groupFit(counts = genomic, distr = "user", pfixed = phat, imhof = FALSE, pave = TRUE))
#'
#' @references
#' V. Choulakian, R. A. Lockhart & M. A. Stephens (1994). Cramer-vonMises statistics for discrete distributions.
#' The Canadian Journal of Statistics,2 2,125-137.
#'
#' J.J.Spinelli (2001). Testing fit for the grouped exponential distribution.
#' The Canadian Journal of Statistics,29,451-458.
#'
#' J.J. Spinelli (1994). Cramer-vonMises statistics for discrete distributions.
#' Phd thesis. Simon Fraser University, Burnaby, Canada.
#'
#' S. Z. Sun, J.J. Spinelli & M. A. Stephens (2012).  Testing fit for the grouped gamma and weibull distributions.
#' Technical report. Simon Fraser University, Burnaby, Canada.
#'
#' @importFrom stats runif rexp rgamma rlnorm rnorm rweibull rmultinom
#' @importFrom stats pexp pgamma plnorm pnorm pweibull
#' @importFrom stats D dbinom dpois integrate pchisq rbinom rpois sd var qpois
#' @importFrom stats4 coef mle
#' @importFrom CompQuadForm imhof
#'
#' @export

groupFit <- function(breaks, counts, data, discrete, distr,
                     N, params, initials, pfixed, bootstrap = FALSE,  numLoops = 5000,
                     known = FALSE,  pave = FALSE, imhof = TRUE) {
  #tryCatch({
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
  if (distr == "user") { ###########################
    known <- TRUE
    discrete <- TRUE
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
    if (!(is.element(distr, c("binom", "pois", "user")))) {
      stop(paste("Discrete function ", distr, " is not currently supported."))
    }
    if(distr == "user"){
      n <- sum(counts)
    } else {
      n <- length(data)
    }
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
    } else if (distr == "user") { ##########################
      if (!(is.numeric(pfixed))) {
        stop(print("Parameter 'pfixed' must be numeric."))
      }
      if(sum(pfixed) != 1){
        warning("Provided probability 'pfixed' does not sum to 1. They have been automatically normalized.")
        pfixed <- pfixed/sum(pfixed)
      }
      k <- length(pfixed)
      phat <- pfixed
    }
    if (distr != "user") {
      breaks <- seq(-0.5, (k - 0.5), 1)
      counts <- table(cut(data, breaks))
    }
  } else { #for continuous distributions
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
  stats <- cvmTest(counts, phat, pave)
  if (bootstrap == FALSE) {
    if (any(n*phat < 5)) {
      warning("Chi-squared approximation may be incorrect.")
    }
    Hj <- stats$Hj
    k1 <- length(Hj)
    I <- diag(k_trunc)  #############
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
        #gj <- ((0 : (k2 - 1)) - N * hats ) * phat[1 : k2] / (hats * (1 - hats))
        gj <- ((0 : (k_trunc - 1)) - N * hats) * phat / (hats * (1 - hats))
        SigD <- SigO - (hats * (1 - hats)) / N * (gj %*% t(gj))
      } else {
        SigD <- SigO - bhat %*% vhat %*% t(bhat)
      }
      SigU <- A %*% SigD %*% t(A)
    }
    if (pave == TRUE) {
      phat1 <- (phat + c(phat[-1], phat[1])) / 2
      D <- diag(phat1)
    }
    MsigW2 <- D %*% SigU
    MsigU2 <- (I - D %*% one %*% t(one)) %*% D %*%
              (I - one %*% t(one) %*% D) %*% SigU
    MsigA2 <- t((D[1 : k1, 1 : k1] %*% Kmat) ^ 0.5) %*%
              SigU[1 : k1, 1 : k1] %*% (D[1 : k1, 1 : k1] %*% Kmat) ^ 0.5
    #MsigX2<- t(invA)%*%invD%*%invA%*%SigU

    A2pval <- cvmPval(stats$Asq, MsigA2, imhof)
    W2pval <- cvmPval(stats$Wsq, MsigW2, imhof)
    U2pval <- cvmPval(stats$Usq, MsigU2, imhof)
    #X2pval <- cvmPval(stats$Chisq, MsigX2, imhof)
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
          #N=the counts needed to make the cumulative probability less plim  .
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
      } else if (distr == "user") { ##########
        bootcounts <- rmultinom(numLoops, n, prob = phat)
      }
      if (distr != "user") {
        #bootbreaks for both Poisson & Binomial
        bootbreaks <- seq(-0.5, (N + 0.5), 1)
        bootcounts <- apply(bootdata, 2,
                            FUN = function(x)table(cut(x, bootbreaks)))
      }
      numEff <- numLoops
    } else { ##continous distributions
      bootcounts <- rmultinom(numLoops, n, prob = phat)
      bootbreaks <- breaks[(k_min + 1) : (k_max + 1)]
      if (distr != "unif") {
        ##remove bootstrap samples have nonempty cells less
        ##than the number of parameters
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
      if (known == FALSE) { ##continous distributions with unknown parameters
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
    ###computing statistics for all cases, discrete and continuous.
    if (known == TRUE){
      bootstats <- matrix(unlist(apply(bootcounts, 2,
                                       function(x)cvmTest(x, phat, pave))),
                          ncol = numEff)
    } else if (known == FALSE) {
      bootstats <- mapply(function(x, y) {cvmTest(x, y, pave)},
                          as.data.frame(bootcounts), as.data.frame(bootphats))
    }
           #All distributions should generate "bootcounts"

    pvals <- mapply(function(x, y){length(x[x > y]) / numEff},
                    as.data.frame(t(bootstats[1 : 4, ])),
                    as.data.frame(stats[1 : 4]))
    pvals <- as.data.frame(pvals, row.names = names)
  }
  test_statistics <- unlist(stats[1 : 4])
  test_statistics <- as.data.frame(test_statistics, row.names = names)

  if (is.element(distr, c("unif", "user"))) {
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

#' distrFit
#'
#' Finds the Maximum Likelihood Estimates of the parameters in a requested distribution.
#'
#' @param breaks Vector defining the breaks in each group
#' @param counts Vector containing the frequency of counts in each group
#' @param distr Character; the name of the distribution users want to fit the data to
#'
#' distrFit supports all of the continuous distributions supported in \code{\link{groupFit}}:
#' @param initials Vector of initial values for the maximum likelihood estimates.
#'
#' @details
#'
#' distFit uses Maximum Likelihood Estimates to optimize the parameters for a requested distribution.
#'
#' @return distFit returns a vector containing the MLEs.
#'
#' @author Shaun Zheng Sun and Dillon Duncan
#'
#' @seealso \code{\link{groupFit}} for fitting data and providing GoF statistics.
#'
#' @examples
#'
#' #fitting exponential data without initial values (Spinelli 2001)
#'
#' breaks <- c(0, 2, 6, 10, 14, 18, 22, 26)
#' counts <- c(21, 9, 5, 2, 1, 1, 0)
#'
#' (mle1 <- distrFit(breaks, counts, distr = "exp"))
#'
#' #fitting generated data with initial values
#'
#' breaks <- seq(0, 40, 2)
#' counts <- table(cut(rweibull(200, 0.5, 3), breaks))
#'
#' (mle2 <- distrFit(breaks, counts, distr = "weibull", initials = c(0.5, 3)))
#'
#' #fitting generated data to a different distribution
#'
#' breaks <- seq(-100, 100, 5)
#' counts <- table(cut(rcauchy(500, -20, 10), breaks))
#'
#' (mle3 <- distrFit(breaks, counts, distr = "norm"))
#'
#' @export


distrFit <-
function(breaks, counts, distr, initials){
  if (! (length(counts) == length(breaks) - 1)) {
    stop(paste("The length of vector 'counts' must be one less than the
               length of vector 'breaks'."))
  }
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
	}
	pdistr <- paste("p", distr, sep = "")
	if (!exists(pdistr)) {
	  stop(paste("Function '", pdistr, "' is not defined.", sep=""))
	}
	if (sum(counts) == 0) {
		data.est <- counts
	} else {
		mids <- (breaks[-length(breaks)] + breaks[-1]) / 2
		data.est <- rep(mids, counts)
	}
	if (is.na(sd(data.est))) {
	  data.est <- counts
  }
	if (!missing(initials)) {
		if ((!(length(initials) == 2) && !(distr == "exp"))
		    || (distr == "exp" && !(length(initials) == 1))) {
			stop("Vector 'initials' has too many or too few starting values")
		}
		start <- as.numeric(initials)
	}
	min <- 0.01
	if (distr == "exp") {
		ll_exp <- function(a) {
			z <- do.call(pdistr, list(breaks, exp(a)))
			- sum(counts * log(diff(z)))
		}
		if (missing(initials)) {
			start <- 1 / mean(data.est)
		}
		fit <- mle(ll_exp, start = list(a = log(start)))
		hats <- exp(coef(fit))
	} else if (distr == "gamma") {
		ll_gamma <- function(a, b) {
			z <- do.call(pdistr, list(breaks, exp(a), exp(b)))
			#the default output of pgamma is shape and rate.
			- sum(counts * log(diff(z)))
		}
		if (missing(initials)) {
			start <- c(mean(data.est) ^ 2 / var(data.est),
			           mean(data.est) / var(data.est))
		}
		fit <- mle(ll_gamma, start = list(a = log(start[1]), b = log(start[2])))
		hats <- exp(coef(fit))
		} else if (distr == "lnorm") {
		ll_lnorm <- function(a, b) {
			z <- do.call(pdistr, list(breaks, a, exp(b)))
			- sum(counts * log(diff(z)))
		}
		if (missing(initials)) {
		  m <- mean(data.est)
		  v <- var(data.est)
			start <- c(log(m / sqrt((1 + v / m ^ 2))), log(1 + v / m ^ 2))
		}
		fit <- mle(ll_lnorm, start = list(a = start[1], b = log(start[2])))
		hats <- c(coef(fit)[1], exp(coef(fit)[2]))
	} else if (distr == "norm") {
		ll_norm <- function(a, b) {
			z <- do.call(pdistr, list(breaks, a, exp(b)))
			- sum(counts * log(diff(z)))
		}
		if(missing(initials)) {
			start <- c(mean(data.est), sd(data.est))
		}
		fit <- mle(ll_norm, start = list(a = start[1], b = log(start[2])))
		hats <- c(coef(fit)[1], exp(coef(fit)[2]))
	} else if (distr == "weibull") {
		ll_weibull <- function(a, b) {
			z <- do.call(pdistr, list(breaks, exp(a), exp(b)))
			- sum(counts * log(diff(z)))
		}
		if (missing(initials)) {
			start <- c(1, mean(data.est))
		}
		fit <- mle(ll_weibull, start = list(a = log(start[1]), b = log(start[2])))
		hats <- exp(coef(fit))
	} else {
	  stop(paste("Distribution '", distr, "' is not currently supported.", sep = ""))
	}
	return(hats)
}

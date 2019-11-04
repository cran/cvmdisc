#' cvmPval
#'
#' Calculate P values of CVM statistics using their asymptotic distributions
#'
#' @param statistic test statistic
#' @param Msig Matrix used to produce eigenvalues to estimate the asymptotic distributions of the test statistic
#' @param imhof Logical. Set to be FALSE if Imhof's method is NOT used to approximate the null. The package"CompQuadForm" is required if "imhof=TRUE"
#'
#' @details cvmPval is used by groupFit to calculate test statistics for fit distributions.
#'
#' @return P-value for test statistic
#'
#' @author Shaun Zheng Sun and Dillon Duncan
#'
#' @seealso \code{\link{groupFit}}: Data fitting function
#'
#' @examples
#'
#' # A_squared and MsigA derived from
#' #(Choulakian, Lockhart and Stephens(1994) Example 3, p8)
#'
#' A_squared <- 1.172932
#'
#' MsigA <- matrix(c(0.05000, 0.03829, 0.02061, 0.00644,
#'                   0.03829, 0.30000, 0.16153, 0.05050,
#'                   0.02061, 0.16153, 0.30000, 0.09379,
#'                   0.00644, 0.05050, 0.09379, 0.30000),
#'                 nrow = 4, ncol = 4, byrow = TRUE)
#'
#' (U2Pval1 = cvmPval(A_squared, MsigA))
#'
#' U_squared <- 0
#'
#' MsigU <- matrix(c(0.16666667, 0.10540926, 0.0745356, 0.05270463, 0.03333333,
#'                   0.10540926, 0.16666667, 0.1178511, 0.08333333, 0.05270463,
#'                   0.07453560, 0.11785113, 0.1666667, 0.11785113, 0.07453560,
#'                   0.05270463, 0.08333333, 0.1178511, 0.16666667, 0.10540926,
#'                   0.03333333, 0.05270463, 0.0745356, 0.10540926, 0.16666667),
#'                 nrow = 5, ncol = 5, byrow = TRUE)
#'
#' (U2Pval2 = cvmPval(U_squared, MsigU, imhof = FALSE))
#'
#' @export

cvmPval <- function(statistic, Msig, imhof = TRUE){

  eigVals <- eigen(Msig)
  if (imhof == TRUE){
    Pval <- imhof(q = statistic, lambda = eigVals$values)$Qq
  }
  else{
    eigenvalues <- Re(eigVals$values)
    k1w <- sum(eigenvalues)
    k2w <- 2 * sum(eigenvalues ^ 2)
    k3w <- 8 * sum(eigenvalues ^ 3)

    bappw <- k3w / (4 * k2w)
    pappw <- 8 * (k2w ^ 3) / (k3w ^ 2)
    aappw <- k1w - bappw * pappw
    Pval <- 1 - pchisq((statistic - aappw) / bappw, pappw)
  }
  return(Pval)
}

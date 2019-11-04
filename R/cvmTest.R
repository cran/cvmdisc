#' cvmTest
#'
#' Calculate test statistics for grouped data's counts and probabilities
#'
#' @param counts vector containing the frequency of the counts in each group
#' @param p vector of probabilities for counts, often called p-hat
#' @param pave Logical. Set to be FALSE if the probabilities in groups are used; set to TRUE if the average of probabilities of groups j and j+1 are used
#'
#' @details cvmTest is used by groupFit to calculate test statistics for the fitted distributions.
#'
#' @return A list with the components:
#'
#' \item{Asq}{Anderson-Darling test statistic}
#'
#' \item{Wsq}{Cramer-Von-Mises test statistic}
#'
#' \item{Usq}{Watson's test statistic}
#'
#' \item{Chisq}{Pearson's Chi-squared test statistic}
#'
#' \item{Hj}{Estimated cumulative probability}
#'
#' @author Shaun Zheng Sun and Dillon Duncan
#'
#' @seealso \code{\link{groupFit}}: Data fitting function
#'
#' @examples
#'
#' #Choulakian, Lockhart and Stephens (1994)
#' counts <- c(10, 19, 18, 15, 11, 13, 7, 10, 13, 23, 15, 22)
#' phat <- rep(1/12, 12)
#'
#' (stats1 <- cvmTest(counts, phat))
#'
#' #Choulakian, Lockhart and Stephens (1994)
#' counts <- c(1 ,4 ,11 ,4 ,0)
#' phat <- c(0.05, 0.3, 0.3, 0.3, 0.05)
#'
#' (stat2 <- cvmTest(counts, phat))
#'
#' #Utilizing Benford's Law
#'
#' #Setting pave to TRUE
#'
#' #genomic data, Lesperance et al (2016)
#' genomic<- c(48, 14, 12, 6, 18, 5, 7, 8, 9)
#' phat<- log10(1+1/1:9)
#'
#' (stat3 <- cvmTest(genomic, phat, pave = TRUE))
#'
#' @export

cvmTest <- function(counts, p, pave = FALSE){
  n <- sum(counts)
  expcounts <- n * p
  Tj <- cumsum(expcounts)
  Sj <- cumsum(counts)
  Zj <- Sj - Tj
  zbar <- sum(Zj * p)

  if (pave == TRUE) {
    p <- (p + c(p[-1], p[1])) / 2
  }
  Wsq <- sum(Zj ^ 2 * p) / n
  Usq <- sum((Zj - zbar) ^ 2 * p) / n
  Hj <- Tj / n
  Hj <- Hj[Hj < 0.99999999]
  k1 <- length(Hj)

  Asq <- sum(Zj[1:k1] ^ 2 * p[1:k1] / (Hj * (1 - Hj))) / n
  Chisq <- sum((counts - expcounts) ^ 2 / expcounts)
  reslist <- list(Asq = Asq,  Wsq = Wsq, Usq = Usq, Chisq = Chisq, Hj = Hj)
  return(reslist)
}

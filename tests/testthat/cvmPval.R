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
cvmTest <- function(counts, p, T){
  n <- sum(counts)
  expcounts <- n * p
  Tj <- cumsum(expcounts)
  Sj <- cumsum(counts)
  Zj <- Sj - Tj
  zbar <- sum(Zj * p)

  if (T == TRUE) {
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
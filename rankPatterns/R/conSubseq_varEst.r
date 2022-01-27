#' Variance estimation for W1
#'
#' Estimate variance for the W1 measure using Monte Carlo simulations on randomly generated sequence pairs that are independent of each other.
#' @param n length of sequences
#' @param k length of subsequences
#' @param m number of random sequence pairs to generate
#' @param ysign indicating whether only the positive indicators ("pos"), only the negative indicators ("neg"), or both ("both") should be counted
conSubseq_varEst <- function(n, k, m, ysign="both"){
  reps <- replicate(m, conSubseq(rnorm(n), rnorm(n), k, ysign))
  varEst <- var(reps)
  return(varEst)
}
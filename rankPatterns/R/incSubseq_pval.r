#' Asymptotic p-value for W2
#'
#' An overall function that calculates W2 using \code{incSubseq}, estimates the variance using \code{incSubseq_varEst}, computes the normalized statistic T2 and the asymptotic p-value. Returns all four outputs.
#' @param x sequence 1
#' @param y sequence 2
#' @param k length of subsequences
#' @param m number of random sequence pairs to generate
#' @param ysign indicating whether only the positive indicators ("pos"), only the negative indicators ("neg"), or both ("both") should be counted
#' @return \code{coexp_meas} value of W2
#' @return \code{est_var} estimated variance of W2
#' @return \code{coexp_meas_norm} value of T2
#' @return \code{p_val} asymptotic p-value
incSubseq_pval <- function(x, y, k, m, ysign="both"){
  n <- length(x)
  
  # calculate measure
  w <- incSubseq(x, y, k, ysign)
  
  # estimate variance 
  estVar <- incSubseq_varEst(n, k, m, ysign)
  
  # calculate asymptotic p-value 
  if (ysign=="both"){
    mu <- 2*choose(n,k)/gamma(k+1)
  } else {
    mu <- choose(n,k)/gamma(k+1)
  }
  t <- (w-mu)/sqrt(estVar)
  pval <- 1-pnorm(t)
  
  ret <- list(coexp_meas=w, est_var=estVar, coexp_meas_norm=t, p_val=pval)
  return(ret)
}
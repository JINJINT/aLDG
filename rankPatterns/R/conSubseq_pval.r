#' Asymptotic p-value for W1
#'
#' An overall function that calculates W1 using \code{conSubseq}, estimates the variance using \code{conSubseq_varEst}, computes the normalized statistic T1 and the asymptotic p-value. Returns all four outputs.
#' @param x sequence 1
#' @param y sequence 2
#' @param k length of subsequences
#' @param m number of random sequence pairs to generate
#' @param ysign indicating whether only the positive indicators ("pos"), only the negative indicators ("neg"), or both ("both") should be counted
#' @return \code{coexp_meas} value of W1
#' @return \code{est_var} estimated variance of W1
#' @return \code{coexp_meas_norm} value of T1
#' @return \code{p_val} asymptotic p-value
conSubseq_pval <- function(x, y, k, m, ysign="both"){
  n <- length(x)
  
  # calculate measure
  w <- conSubseq(x, y, k, ysign)
  
  # estimate variance 
  estVar <- conSubseq_varEst(n, k, m, ysign)
  
  # calculate asymptotic p-value
  if (ysign=="both"){
    mu <- 2*(n-k+1)/gamma(k+1)
  } else {
    mu <- (n-k+1)/gamma(k+1)
  }
  t <- (w-mu)/sqrt(estVar)
  pval <- 1-pnorm(t)
  
  ret <- list(coexp_meas=w, est_var=estVar, coexp_meas_norm=t, p_val=pval)
  return(ret)
}
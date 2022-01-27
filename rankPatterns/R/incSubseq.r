#' W2
#'
#' A wrapper function that calculates W2 and its variants, calling \code{incSubseq_small_n} for \code{n} less than 5000 and \code{incSubseq_large_n} for larger \code{n}.
#' @param x sequence 1
#' @param y sequence 2
#' @param k length of subsequences
#' @param ysign indicating whether only the positive indicators ("pos"), only the negative indicators ("neg"), or both ("both") should be counted
incSubseq <- function(x,y,k,ysign="both"){
  n <- length(x)
  if (n<=5000){
    ret <- incSubseq_small_n(x,y,k,ysign)
  } else {
    ret <- incSubseq_large_n(x,y,k,ysign)
  }
  return(ret)
}
#' W1
#'
#' A wrapper function that calculates W1 and its variants, calling \code{conSubseq_small_k} for \code{k} less than 100 and \code{conSubseq_large_k} for larger \code{k}.
#' @param x sequence 1
#' @param y sequence 2
#' @param k length of subsequences
#' @param ysign indicating whether only the positive indicators ("pos"), only the negative indicators ("neg"), or both ("both") should be counted
conSubseq <- function(x,y,k,ysign="both"){
  if (k<=100){
    ret <- conSubseq_small_k(x,y,k,ysign)
  } else {
    ret <- conSubseq_large_k(x,y,k,ysign)
  }
  return(ret)
}
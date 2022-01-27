#' W1 for large k
#'
#' Calculate W1 and its variants. The algorithm has a running time of O(klog(k)n), but runs slower than \code{conSubseq_small_k} for small \code{k}; recommended for large \code{k}.
#' @param x sequence 1
#' @param y sequence 2
#' @param k length of subsequences
#' @param ysign indicating whether only the positive indicators ("pos"), only the negative indicators ("neg"), or both ("both") should be counted
conSubseq_large_k <- function(x,y,k,ysign="both"){
  n <- length(x)
  cnt <- 0
  if (ysign=="both"){
    z1 <- .C('conSubseqr', as.double(x), as.double(y), as.integer(k), as.integer(n), cnt = as.double(cnt));
    z2 <- .C('conSubseqr', as.double(x), as.double(-y), as.integer(k), as.integer(n), cnt = as.double(cnt));
    return(z1$cnt + z2$cnt)
  } else if (ysign=="pos"){
    z <- .C('conSubseqr', as.double(x), as.double(y), as.integer(k), as.integer(n), cnt = as.double(cnt));
  } else if (ysign=="neg"){
    z <- .C('conSubseqr', as.double(x), as.double(-y), as.integer(k), as.integer(n), cnt = as.double(cnt));
  }
  return(z$cnt)
}
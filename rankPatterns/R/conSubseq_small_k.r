#' W1 for small k
#'
#' Calculate W1 and its variants. The algorithm has a running time of O(k^2 n) and is recommended for small \code{k}.
#' @param x sequence 1
#' @param y sequence 2
#' @param k length of subsequences
#' @param ysign indicating whether only the positive indicators ("pos"), only the negative indicators ("neg"), or both ("both") should be counted
conSubseq_small_k <- function(x,y,k,ysign="both"){
  n <- length(x)
  pt <- 0
  if (ysign=="both"){
    z1 <- .C('conSubseq_r', as.double(x), as.double(y), as.integer(k), as.integer(n), pt = as.double(pt));
    z2 <- .C('conSubseq_r', as.double(x), as.double(-y), as.integer(k), as.integer(n), pt = as.double(pt));
    return(z1$pt + z2$pt)
  } else if (ysign=="pos"){
    z <- .C('conSubseq_r', as.double(x), as.double(y), as.integer(k), as.integer(n), pt = as.double(pt));
  } else if (ysign=="neg"){
    z <- .C('conSubseq_r', as.double(x), as.double(-y), as.integer(k), as.integer(n), pt = as.double(pt));
  }
  return(z$pt)
}
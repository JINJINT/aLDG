#' W2 for small n
#'
#' Calculate W2 and its variants. The algorithm has a running time of O(k n^2) and is recommended for small \code{n}.
#' @param x sequence 1
#' @param y sequence 2
#' @param k length of subsequences
#' @param ysign indicating whether only the positive indicators ("pos"), only the negative indicators ("neg"), or both ("both") should be counted
incSubseq_small_n <- function(x,y,k,ysign="both"){
  n <- length(x)
  pt <- 0
  if (ysign=="both"){
    pairs1 <- cbind (x,y)
    pairs2 <- cbind (x,-y)
    z1 <- pairs1[order(pairs1[,1]),][,2]
    z2 <- pairs2[order(pairs2[,1]),][,2]
    rst1 <- .C('incSubseq_r', as.double(z1), as.integer(k), as.integer(n), pt = as.double(pt));
    rst2 <- .C('incSubseq_r', as.double(z2), as.integer(k), as.integer(n), pt = as.double(pt));
    return(rst1$pt + rst2$pt)
  } else if (ysign=="pos"){
    pairs <- cbind (x,y)
  } else if (ysign=="neg"){
    pairs <- cbind (x,-y)
  }
  z <- pairs[order(pairs[,1]),][,2]
  rst <- .C('incSubseq_r', as.double(z), as.integer(k), as.integer(n), pt = as.double(pt));
  return(rst$pt)
}
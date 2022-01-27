
#' @export
randcop <-function(Rho, nCells){
  Col = chol(Rho)
  nGenes = nrow(Rho)
  copular = matrix(rnorm(nGenes*nCells), ncol = nCells)
  copular = t(Col) %*% copular
  copular = pnorm(copular)
  return(copular)
}

#' @export
randcor <- function(ngenes){
  data(puri_data)
  corr = cor(t(puri_data[sample(nrow(puri_data),ngenes),]))
  d <- stats::as.dist((1 - corr)/2)
  h <- stats::hclust(d)
  order <- h$order
  ans <- corr[order, order]
  ans = makespd(ans)
  return(ans)
}

#' @export
rnorm_truc <- function(n, mean, sd, a, b){
  vec1 <- rnorm(n, mean = mean, sd=sd)
  beyond_idx <- which(vec1 < a | vec1 > b)
  if (length(beyond_idx) > 0) { # for each value < rate_2cap_lb
    substi_vec <- sapply(1:length(beyond_idx), function(i){
      while (TRUE){
        temp <- rnorm(1, mean = mean, sd=sd)
        if (temp > a | temp > b) {break}}
      return(temp)} )
    vec1[beyond_idx] <- substi_vec
  }
  return(vec1)
}

# to simulate a pair of different correlation matrix
#' @export
cormatpair<-function(n, p, d, eigenmean = 2, eigengap = 10, eigenprop = 0.3, method = 'gap'){

  # set sd according to signal to noise ratio
  sd = eigenmean/(5*sqrt(log(p)/n))
  ans1 = diag(2*sd,p)
  ans2 = diag(2*sd,p)

  # big gap in first eigen value
  if(method=='gap'){
    eigenvalues1 = sort(abs(rnorm(d, mean=eigenmean, sd=sd)), decreasing = TRUE)
    eigenvalues2 = eigenvalues1
    eigenvalues2[1:max(floor(d*eigenprop),1)] = eigenvalues2[1:max(floor(d*eigenprop),1)] + eigengap

    eigenvectors =  randortho(d, type = "orthonormal")

    ans1[1:d,1:d] = eigenvectors %*% diag(eigenvalues1 + 1e-3) %*% t(eigenvectors)
    ans2[1:d,1:d] = eigenvectors %*% diag(eigenvalues2 + 1e-3) %*% t(eigenvectors)
  }

  # reverse eigen order
  if(method=='order'){
    eigenvalues = sort(abs(rnorm(d, mean=eigenmean, sd=sd)),decreasing = TRUE)
    eigenvectors =  randortho(d, type = "orthonormal")[,1:d]

    ans1[1:d, 1:d] = eigenvectors %*% diag(eigenvalues+ 1e-6) %*% t(eigenvectors)
    ans2[(p-d+1):p, (p-d+1):p] = eigenvectors %*% diag(rev(eigenvalues)+ 1e-6) %*% t(eigenvectors)
  }

  # random generated
  if(method=='rand'){
    load('puri_data.rdata')
    ans1 = cor(t(puri_data[sample(nrow(puri_data),p),]))
    ans1[is.na(ans1)] = 0
    ans1[is.nan(ans1)] = 0
    d = as.dist((1-ans1)/2)
    h = hclust(d)
    ord = h$order
    ans1 = ans1[ord,ord]

    ans2 = cor(t(puri_data[sample(nrow(puri_data),p),]))
    ans2[is.na(ans2)] = 0
    ans2[is.nan(ans2)] = 0
    d = as.dist((1-ans2)/2)
    h = hclust(d)
    ord = h$order
    ans2 = ans2[ord,ord]

  }

  #ans1 = makespd(ans1)
  #ans2 = makespd(ans2)
  return(list(ans1, ans2))
}



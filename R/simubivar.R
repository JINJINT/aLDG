
#' @export
# gaussian mixture
sims.gaussmix<-function(n,k,s){
  if(k==3){
    mean1 = c(0,0)
    mean2 = c(3,3)
    mean3 = c(-4,4)
    sig1 = c(1,1)
    sig2 = c(1,1)
    sig3 = c(1,1)
    corvec = rep(0,k)
    if(s>0)corvec[1:s] = 0.8
    corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
    corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
    corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
    para1list = list(mean1, mean2, mean3)
    para2list = list(sig1, sig2, sig3)
    corrlist = list(corr1, corr2, corr3)
    p = c(0.1,0.1,0.8)
  }
  if(k==4){
    mean1 = c(3,3)
    mean2 = c(-3,3)
    mean3 = c(3,-3)
    mean4 = c(-3,-3)
    sig1 = c(1,1)
    sig2 = c(1,1)
    sig3 = c(1,1)
    sig4 = c(1,1)
    corvec = rep(0,k)
    if(s>0)corvec[1:s] = 0.8
    corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
    corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
    corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
    corr4 = matrix(c(1,corvec[4],corvec[4],1),nrow = 2)
    para1list = list(mean1, mean2, mean3, mean4)
    para2list = list(sig1, sig2, sig3, sig4)
    corrlist = list(corr1, corr2, corr3, corr4)
    p = c(0.25,0.25,0.25,0.25)
  }
  data = sim_mix(distrlist=c('norm','norm','norm','norm'),
                 n = n, d = 2, prob = p,
                 corrlist = corrlist,
                 para1list = para1list,
                 para2list = para2list)
  dat = data[[1]]
  return(list(X =dat[1,],Y = dat[2,]))
}

#' @export
sims.nbmix<-function(n, k, s){
  if(k==3){
    mean1 = c(20,20)
    mean2 = c(1,3)
    mean3 = c(3,1)
    sig1 = c(100,100)
    sig2 = c(1,1)
    sig3 = c(1,1)
    corvec = rep(0,k)
    if(s>0)corvec[1:s] = 0.8
    corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
    corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
    corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
    para1list = list(mean1, mean2, mean3)
    para2list = list(sig1, sig2, sig3)
    corrlist = list(corr1, corr2, corr3)
    p = c(0.1,0.4,0.5)
  }
  if(k==4){
    mean1 = c(10,20)
    mean2 = c(20,10)
    mean3 = c(1,3)
    mean4 = c(3,1)
    sig1 = c(100,100)
    sig2 = c(100,100)
    sig3 = c(1,1)
    sig4 = c(1,1)
    corvec = rep(0,k)
    if(s>0)corvec[1:s] = 0.8
    corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
    corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
    corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
    corr4 = matrix(c(1,corvec[4],corvec[4],1),nrow = 2)
    para1list = list(mean1, mean2, mean3, mean4)
    para2list = list(sig1, sig2, sig3, sig4)
    corrlist = list(corr1, corr2, corr3, corr4)
    p = c(0.25,0.25,0.25,0.25)
  }

  data_nb = sim_mix(distrlist=c('nbinom','nbinom','nbinom','nbinom'),
                    n=n, d = 2, prob = p, corrlist = corrlist,
                    para1list = para1list, para2list = para2list)
  dat_nb = data_nb[[1]]
  return(list(X = dat_nb[1,], Y = dat_nb[2,]))
}

#' @export
sims.gausshardmix<-function(n,k,s){
  if(k==3){
    mean1 = c(0,0)
    mean2 = c(3,3)
    mean3 = c(-4,4)
    sig1 = c(1,1)
    sig2 = c(1,1)
    sig3 = c(1,1)
    corvec = rep(0,k)
    if(s>0)corvec[1:s] = as.numeric(s)*0.2-1
    corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
    corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
    corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
    para1list = list(mean1, mean2, mean3)
    para2list = list(sig1, sig2, sig3)
    corrlist = list(corr1, corr2, corr3)
    p = c(0.1,0.1,0.8)
  }
  if(k==4){
    mean1 = c(3,2)
    mean2 = c(-3,2)
    mean3 = c(2,-3)
    mean4 = c(-2,-3)
    sig1 = c(1,1)
    sig2 = c(1,1)
    sig3 = c(1,1)
    sig4 = c(1,1)
    corvec = rep(0,k)
    if(s>0)corvec[1:s] = as.numeric(s)*0.2-1
    corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
    corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
    corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
    corr4 = matrix(c(1,corvec[4],corvec[4],1),nrow = 2)
    para1list = list(mean1, mean2, mean3, mean4)
    para2list = list(sig1, sig2, sig3, sig4)
    corrlist = list(corr1, corr2, corr3, corr4)
    p = c(0.1,0.1,0.1,0.7)
  }
  data = sim_mix(distrlist=c('norm','norm','norm','norm'),
                 n = n, d = 2, prob = p,
                 corrlist = corrlist,
                 para1list = para1list,
                 para2list = para2list)
  dat = data[[1]]
  return(list(X =dat[1,],Y = dat[2,]))
}

#' @export
sims.nbhardmix<-function(n, k, s){
  if(k==3){
    mean1 = c(20,20)
    mean2 = c(1,3)
    mean3 = c(3,1)
    sig1 = c(100,100)
    sig2 = c(1,1)
    sig3 = c(1,1)
    corvec = rep(0,k)
    if(s>0)corvec[1:s] = 0.8
    corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
    corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
    corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
    para1list = list(mean1, mean2, mean3)
    para2list = list(sig1, sig2, sig3)
    corrlist = list(corr1, corr2, corr3)
    p = c(0.1,0.4,0.5)
  }
  if(k==4){
    mean1 = c(20,30)
    mean2 = c(50,30)
    mean3 = c(1,3)
    mean4 = c(3,1)
    sig1 = c(100,100)
    sig2 = c(100,100)
    sig3 = c(1,1)
    sig4 = c(1,1)
    corvec = rep(0,k)
    if(s>0)corvec[1:s] = 0.8
    corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
    corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
    corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
    corr4 = matrix(c(1,corvec[4],corvec[4],1),nrow = 2)
    para1list = list(mean1, mean2, mean3, mean4)
    para2list = list(sig1, sig2, sig3, sig4)
    corrlist = list(corr1, corr2, corr3, corr4)
    p = c(0.1,0.1,0.4,0.4)
  }
  data_nb = sim_mix(distrlist=c('nbinom','nbinom','nbinom','nbinom'),
                    n=n, d = 2, prob = p, corrlist = corrlist,
                    para1list = para1list, para2list = para2list)
  dat_nb = data_nb[[1]]
  return(list(X = dat_nb[1,], Y = dat_nb[2,]))
}

#' @export
sims.nbharddropmix<-function(n, k, s, dp){

  mean1 = c(20,30)
  mean2 = c(50,30)
  mean3 = c(1,3)
  mean4 = c(3,1)
  sig1 = c(100,100)
  sig2 = c(100,100)
  sig3 = c(1,1)
  sig4 = c(1,1)
  corvec = rep(0,k)
  if(s>0)corvec[1:s] = 0.8
  corr1 = matrix(c(1,corvec[1],corvec[1],1),nrow = 2)
  corr2 = matrix(c(1,corvec[2],corvec[2],1),nrow = 2)
  corr3 = matrix(c(1,corvec[3],corvec[3],1),nrow = 2)
  corr4 = matrix(c(1,corvec[4],corvec[4],1),nrow = 2)
  para1list = list(mean1, mean2, mean3, mean4)
  para2list = list(sig1, sig2, sig3, sig4)
  corrlist = list(corr1, corr2, corr3, corr4)
  p = c(0.15,0.15,0.35,0.35)

  data_nb = sim_mix(distrlist=c('nbinom','nbinom','nbinom','nbinom'),
                    n=n, d = 2, prob = p, corrlist = corrlist,
                    para1list = para1list, para2list = para2list)
  dat_nb = data_nb[[1]]
  label = data_nb[[2]]
  dropx = which(as.integer(rbernoulli(n-4,dp))==1)
  dropy = which(as.integer(rbernoulli(n-4,dp))==1)
  dat_nb[1,dropx] = 0
  dat_nb[2,dropy] = 0
  return(list(X = dat_nb[1,], Y = dat_nb[2,]))
}


#======== simulate different shapes based on unif distribution and normal noise
#' @export
sim_linear <- function(n, d, eps=1, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  ys <- xs%*%w + kappa*eps*nu  # y = xA + nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}

#' @export
sim_exp <- function(n, d, eps=10, ind=FALSE, a=0, b=3) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  ys <- exp(xs%*%w) + kappa*eps*nu  # y = exp(xA) + nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}
#' @export
sim_cubic <- function(n, d, eps=80, ind=FALSE, a=-1, b=1, c.coef=c(-12, 48, 128), s=1/3) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  xw = xs%*%w
  ys <- c.coef[3]*(xw - s)^3 + c.coef[2]*(xw - s)^2 + c.coef[1]*(xw - s) + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=ys))
}
#' @export
sim_step <- function(n, d, eps=1, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- as.numeric(xs%*%w > 0) + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}

#' @export
sim_quad <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- (xs%*%w)^2 + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}
#' @export
sim_2parab <- function(n, d, eps=0.5, a=-1, b=1,prob=0.5) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  u = array(rbinom(n*d,1,prob), dim=c(n,1))
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- ((xs%*%w)^2 + eps*kappa*nu)*(u-0.5)
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}
#' @export
sim_4root <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- (xs%*%w)^{1/4} + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}
#' @export
sim_sine <- function(n, d, period = 1, eps=0.5, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- sin(xs%*%w * period * pi) + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}

#' @export
sim_log <- function(n, d, pirod = 1, eps=0.5, ind=FALSE, a=-1, b=1) {
  xs <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  y <- log(xs%*%w) + eps*kappa*nu
  if (ind) {
    xs <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=xs, Y=y))
}


sim_wshape <- function(n, d, eps=0.5, ind=FALSE, a=-1, b=1) {
  x <- gen.x.unif(n, d, a=a, b=b)
  u <- gen.x.unif(n, d, a=a, b=b)
  w <- gen.coefs(d)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  kappa <- as.numeric(d == 1)
  nu <- rnorm(dim(x)[1], mean=0, sd=1)  # gaussian noise
  y <- 4*(((x%*%w)^2 - 1/2)^2 + u%*%w/500) + kappa*eps*nu
  if (ind) {
    x <- gen.x.unif(n, d, a=a, b=b)
  }
  return(list(X=x, Y=y))
}


sim_spiral <- function(n, d, eps=0.4, a=0, b=5) {
  u <- gen.x.unif(n, 1, a=a, b=b)
  x <- array(cos(pi*u), dim=c(n, d))
  y <- u*sin(pi*u)
  if (d > 1) {
    for (i in 1:(d-1)) {
      x[, i] <- y*x[, i, drop=FALSE]^i
    }
  }
  x[, d] <- u*x[, d, drop=FALSE]
  nu <- rnorm(dim(x)[1], mean=0, sd=1)  # gaussian noise
  y <- y + eps*d*nu
  return(list(X=x, Y=y))
}

sim_circle <- function(n, d, eps=0.05, a=-1, b=1, r=1) {
  u <- gen.x.unif(n, 1, a=a, b=b)
  x <- r*cos(u[,1]*pi)
  y <- sin(u[,1]*pi)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  y <- y + eps*d*nu
  return(list(X=x, Y=y))
}

sim_ellipse <- function(n, d, eps=0.05, a=-1, b=1, r=5) {
  u <- gen.x.unif(n, 1, a=a, b=b)
  x <- r*cos(u[,1]*pi)
  y <- sin(u[,1]*pi)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  y <- y + eps*d*nu
  return(list(X=x, Y=y))
}


sim_diamond <- function(n, d, eps=0.05, a=-1, b=1, period=-pi/4) {
  u <- gen.x.unif(n, 2, a=a, b=b)
  x <- u[,1]*cos(period) + u[,2]*sin(period)
  y <- -u[,1]*sin(period) + u[,2]*cos(period)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  y <- y + eps*d*nu
  return(list(X=x, Y=y))
}

sim_square <- function(n, d, eps=0.05, a=-1, b=1, period=-pi/8) {
  u <- gen.x.unif(n, 2, a=a, b=b)
  x <- u[,1]*cos(period) + u[,2]*sin(period)
  y <- -u[,1]*sin(period) + u[,2]*cos(period)
  nu <- rnorm(n, mean=0, sd=1)  # gaussian noise
  y <- y + eps*d*nu
  return(list(X=x, Y=y))
}

sim_ubern <- function(n, d, eps=0.5, p=0.5) {
  U = rbinom(n=n, size=1, prob=p)
  nu_e1 <- mvrnorm(n=n, mu = array(0, dim=c(d, 1)), Sigma = diag(d))  # gaussian noise
  X = array(rbinom(n=d*n, size=1, prob=p), dim=c(n, d)) + eps*nu_e1
  w <- gen.coefs(d)
  Y <- array(NaN, dim=c(n, 1))
  nu_e2 <- rnorm(n, mean=0, sd=1)
  for (i in 1:n) {
    Y[i] <- (2*U[i] - 1)*w^T %*% X[i,] + eps*nu_e2[i]
  }
  return(list(X=X, Y=Y))
}


sim_multi <- function(n, d) {
  X <- mvrnorm(n=n, mu = array(0, dim=c(d, 1)), Sigma = diag(d))  # gaussian noise
  Y <- mvrnorm(n=n, mu = array(0, dim=c(d, 1)), Sigma = diag(d))
  Y = X*Y

  return(list(X=X, Y=Y))
}

sim_indep <- function(n, d, prob=0.5, sep1=1,sep2=2) {
  X <- mvrnorm(n=n, mu = array(0, dim=c(d, 1)), Sigma = diag(d))  # gaussian noise
  Y <- mvrnorm(n=n, mu = array(0, dim=c(d, 1)), Sigma = diag(d))
  u = rbinom(n*d, 1, prob)
  v = rbinom(n*d, 1, prob)
  X = X/sep1 + sep2 * u -1
  Y = Y/sep1 + sep2 * v -1
  return(list(X=X, Y=Y))
}


#' A helper function for simulating sample labels
#' @param K the number of classes
#' @param class.equal whether the number of samples/class should be equal, with each
#' class having a prior of 1/K, or inequal, in which each class obtains a prior
#' of k/sum(K) for k=1:K. Defaults to \code{TRUE}.
#' @keywords internal
gen.sample.labels <- function(K, class.equal=TRUE) {
  if (isTRUE(class.equal)) {
    priors <- array(1/K, dim=c(K)) # prior is just 1/K for all k
  } else {
    priors <- (1:K)/sum(1:K) # prior is k/sum(K) for all k
  }
  # return class priors
  return(priors)
}


#'
#' A helper function to generate a d-dimensional linear transformation matrix.
#' @param d the number of dimensions.
#' @return A \code{[d]} the coefficient vector.
#' @author Eric Bridgeford
#' @keywords internal
gen.coefs <- function(d) {
  A = as.array(1/1:d, dim=c(d, 1))
  return(A)
}

#' A helper function to generate n samples of a d-dimensional uniform vector.
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param a the lower limit.
#' @param b the upper limit.
#' @param x \code{[n, d]} the simulated data matrix.
#' @author Eric Bridgeford
#' @importFrom stats runif
#' @keywords internal
gen.x.unif <- function(n, d, a=-1, b=1) {
  x <- array(runif(n=(n*d), min=a, max=b), dim=c(n, d))
  return(x)
}

#' GMM Simulate
#'
#' A helper function for simulating from Gaussian Mixture.
#' @param mus \code{[d, K]} the mus for each class.
#' @param Sigmas \code{[d,d,K]} the Sigmas for each class.
#' @param n the number of examples.
#' @param priors \code{K} the priors for each class.
#' @return A list with the following:
#' \item{X}{\code{[n, d]} the simulated data.}
#' \item{Y}{\code{[n]} the labels for each data point.}
#' \item{priors}{\code{[K]} the priors for each class.}
#' @author Eric Bridgeford
#' @importFrom MASS mvrnorm
#' @keywords internal
mgc.sims.sim_gmm <- function(mus, Sigmas, n, priors) {
  K <- dim(mus)[2]
  labs <- sample(1:K, size=n, prob=priors, replace=TRUE)
  ylabs <- as.vector(sort(unique(labs)))
  res <- sapply(ylabs, function(y) mvrnorm(n=sum(labs == y), mus[,y], Sigmas[,,y]), USE.NAMES=TRUE, simplify=FALSE)
  X <- array(0, dim=c(n, dim(Sigmas)[1]))
  for (y in ylabs) {
    X[labs == y,] <- res[[y]]
  }
  return(list(X=X, Y=labs, priors=priors))
}

#' Sample Random Rotation
#'
#' A helper function for estimating a random rotation matrix.
#' @importFrom stats rnorm
#' @param d dimensions to generate a rotation matrix for.
#' @return the rotation matrix
#' @author Eric Bridgeford
#' @keywords internal
mgc.sims.rotation <- function(d) {
  Q <- qr.Q(qr(array(rnorm(d*d), dim=c(d, d))))
  if (det(Q) < -.99) {
    Q[,1] <- -Q[,1]
  }
  return(Q)
}

#' Random Rotation
#'
#' A helper function for applying a random rotation to gaussian parameter set.
#' @param mus means per class.
#' @param Sigmas covariances per class.
#' @param Q rotation to use, if any
#' @author Eric Bridgeford
#' @keywords internal
mgc.sims.random_rotate <- function(mus, Sigmas, Q=NULL) {
  dimm <- dim(mus)
  K <- dimm[2]
  d <- dim(mus)[1]
  if (is.null(Q)) {
    Q <- mgc.sims.rotation(d)
  } else if (!isTRUE(all.equal(dim(Q), c(d, d)))) {
    stop(sprintf("You have specified a rotation matrix with dimensions (%d, %d), but should be (%d, %d).", dim(Q)[1], dim(Q)[2], d, d))
  }

  for (i in 1:K) {
    mus[,i] <- Q %*% mus[,i,drop=FALSE]
    Sigmas[,,i] <- Q %*% Sigmas[,,i] %*% t(Q)
  }
  return(list(mus=mus, S=Sigmas, Q=Q))
}

#' Sample from Unit 2-Ball
#'
#' Sample from the 2-ball in d-dimensions.
#'
#' @importFrom MASS mvrnorm
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param r the radius of the 2-ball. Defaults to \code{1}.
#' @param cov.scale if desired, sample from 2-ball with error sigma. Defaults to \code{NaN},
#' which has no noise.
#' @return the points sampled from the ball, as a \code{[n, d]} array.
#' @examples
#' library(mgc)
#' # sample 100 points from 3-d 2-ball with radius 2
#' X <- mgc.sims.2ball(100, 3, 2)
#' @author Eric Bridgeford
#' @export
mgc.sims.2ball <- function(n, d, r=1, cov.scale=0) {
  Y <- matrix(mvrnorm(n=n, mu=array(0, dim=c(d, 1)), Sigma=diag(d)), nrow=n)
  u <- runif(n)
  r <- r * u^(1/d)
  X <- r * Y/sqrt(apply(Y^2, 1, sum)) + matrix(mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=cov.scale*diag(d)), nrow=n)
  return(X)
}

#' Sample from Unit 2-Sphere
#'
#' Sample from the 2-sphere in d-dimensions.
#'
#' @importFrom MASS mvrnorm ginv
#' @param n the number of samples.
#' @param d the number of dimensions.
#' @param r the radius of the 2-ball. Defaults to \code{1}.
#' @param cov.scale if desired, sample from 2-ball with error sigma. Defaults to \code{0},
#' which has no noise.
#' @return the points sampled from the sphere, as a \code{[n, d]} array.
#' @examples
#' library(mgc)
#' # sample 100 points from 3-d 2-sphere with radius 2
#' X <- mgc.sims.2sphere(100, 3, 2)
#' @author Eric Bridgeford
#' @export
mgc.sims.2sphere <- function(n, d, r, cov.scale=0) {
  u <- matrix(mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=diag(d)), nrow=n)
  unorm <- diag(sqrt(apply(u^2, 1, sum)))
  pts <- r*(ginv(unorm) %*% u)
  pts <- pts + matrix(mvrnorm(n=n, mu=array(0, dim=c(d,1)), Sigma=cov.scale*diag(d)), nrow=n)
  return(pts)
}

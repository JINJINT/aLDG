#' @export
getDiffMatrix<-function(X,Y,methods,thred=-Inf,info=0, norm = FALSE, abs = FALSE,plot=FALSE,extrainfo='',dir='./dat/',ncores=NULL){
  Dhat_list = list()
  if(length(methods)>0){
    CX <- matdep(X, methods = methods,
                          thred = thred, info=paste0(extrainfo,'X'),trial=info,
                          norm = norm, abs = abs, 
                          stat = 'normgap', band = 'fix',
                          wd = 2, qd = 0.05, opt=TRUE, cutoff =1,
                          dir=dir,ncores=ncores)
    CY <- matdep(Y, methods = methods,
                          thred = thred, info=paste0(extrainfo,'Y'),trial=info,
                          norm = norm, abs = abs, 
                          stat = 'normgap', band = 'fix',
                          wd = 2, qd = 0.05, opt=TRUE, cutoff =1,
                          dir=dir,ncores=ncores)
    methods = names(CX)
    for(i in 1:length(methods)){
      CX[[i]][is.infinite(CX[[i]])]=0
      CX[[i]][is.na(CX[[i]])]=0
      CY[[i]][is.infinite(CY[[i]])]=0
      CY[[i]][is.na(CY[[i]])]=0
    }
    Dhat_list = lapply(1:length(methods), function(i)CX[[i]]-CY[[i]])
    names(Dhat_list) = methods
  }
  #if(info==0 & plot==TRUE){
  #    plot_cormat(CX, CY, NULL, ord = NULL, dir='./',
  #                ordmethod = NULL,
  #                extrainfo = extrainfo)
  #    return(list(CX,CY))
  #}
  #
  return(Dhat_list)
}

sLED_one <- function(X, Y, nperm, permidx, thred = -Inf, methods = c('pearson'),
                     norm = FALSE, abs = FALSE,
                     rho=1000, sumabs.seq=0.2,
                     mc.cores=NULL, niter=20, trace=FALSE, extrainfo = ''){
  filename = paste0('./dat/sledstat_', extrainfo,'_',nperm,'.rds')
  if(!file.exists(filename)){
    Z <- cbind(X, Y)
    X = Z[,permidx$i1]
    Y = Z[,permidx$i2]
    
    Dhat_list <-getDiffMatrix(X,Y,methods=methods,
                              thred = thred, norm = norm, info = nperm,
                              abs = abs, extrainfo = extrainfo, ncores=mc.cores)
    Tn_list = list()
    for(method in names(Dhat_list)){
      Tn_list[[method]] <- sLEDTestStat(Dmat=Dhat_list[[method]], rho=rho,
                                        sumabs.seq=sumabs.seq, niter=niter, trace=trace)$stats
    }
    saveRDS(list(Tn_list = Tn_list, permidx = permidx), filename)
  }
}


sLED_all <- function(X, Y, thred = -Inf, methods = c('Pearson'),
                     norm = TRUE, abs = TRUE,
                     rho=1000, sumabs.seq=0.2, start = 0, npermute=100,
                     useMC=FALSE, mc.cores=100, seeds=NULL, verbose=TRUE,
                     niter=20, trace=FALSE, extrainfo = '', save= TRUE,
                     filelist = NULL) {
  
  ## permutation statistics, a list of matrix, each row is sparsity level, each column is of permutation
  n1 <- ncol(X)
  n2 <- ncol(Y)
  Z <- cbind(X, Y) # p*n
  ntest <- length(sumabs.seq)
  
  filenamenull = paste0('./dat/sledstat_', extrainfo,'_0.rds')
  #print(filenamenull)
  if(!file.exists(filenamenull)){
    print(filenamenull)
    Dhat_list <-getDiffMatrix(X,Y,methods=methods, thred = thred,
                              norm = norm, abs = abs, extrainfo = extrainfo)
    Tn_list = list()
    for(method in names(Dhat_list)){
      Tn_list[[method]] <- sLEDTestStat(Dmat=Dhat_list[[method]], rho=rho,
                                        sumabs.seq=sumabs.seq, niter=niter, trace=trace)$stats
    }
    saveRDS(Tn_list, filenamenull)
  }
  Tn = readRDS(filenamenull)
  Tn = Tn$Tn_list
  Tn_permute = list()
  
  if(is.null(filelist)){
    filelist = c()
    for(j in 1:npermute)filelist[j] = paste0('./dat/sledstat_', extrainfo,'_',j,'.rds')
  }
  else{
    npermute = length(filelist)
  }
  for(method in names(Tn)){
    Tn_permute[[method]] = c()
  }
  for(j in 1:npermute){
    filename = filelist[j]
    if(!file.exists(filename)){
      if (!is.null(seeds)){
        set.seed(seeds[i])
      }
      i.permute <- permuteIndex(n1, n2)
      Tn_permute_list <- sLED_one(X, Y, j, i.permute, thred = thred, methods = methods,
                                              norm = norm, abs = abs,
                                              rho = rho, sumabs.seq=sumabs.seq,
                                              mc.cores=mc.cores, niter=niter, trace=trace, extrainfo = extrainfo)
      # test statistic, a list of vector, each vector is of sparsity level
      saveRDS(Tn_permute_list, filename)
    }
    Tn_permute_list = readRDS(filename)
    Tn_permute_list = Tn_permute_list$Tn_list
    for(method in names(Tn)){
      Tn_permute[[method]][j] <- Tn_permute_list[[method]]
    }
    
  }
  pVal_list = sapply(names(Tn), function(method)(mean(Tn_permute[[method]] > Tn[[method]])))
  return(list(pVal = pVal_list, Tn = Tn, Tn_permute = Tn_permute))
}



sLEDTestStat <- function(Dmat, rho=1000, sumabs.seq=0.2,
                         niter=20, trace=FALSE) {
  ndim <- 1 ## only consider the first sparse eigenvector
  p <- nrow(Dmat)
  ntest <- length(sumabs.seq)
  
  results <- list()
  results$sumabs.seq <- sumabs.seq
  results$rho <- rho
  
  results$stats <- rep(NA, ntest)
  results$sign <- rep(NA, ntest)
  results$v <- matrix(NA, nrow=ntest, ncol=p)
  results$leverage <- matrix(NA, nrow=ntest, ncol=p)
  
  ## for each sparsity parameter
  for (i in 1:ntest) {
    sumabs <- sumabs.seq[i]
    pos.out <- symmPMD(Dmat + rho * diag(p),
                       sumabs=sumabs, trace=trace, niter=niter)
    neg.out <- symmPMD(- Dmat + rho * diag(p),
                       sumabs=sumabs, trace=trace, niter=niter)
    
    if (pos.out$d >= neg.out$d) {
      results$sign[i] <- "pos"
      results$stats[i] <- pos.out$d - rho
      results$v[i, ] <- pos.out$v
      results$leverage[i, ] <- (pos.out$v)^2
    } else {
      results$sign[i] <- "neg"
      results$stats[i] <- neg.out$d - rho
      results$v[i, ] <- neg.out$v
      results$leverage[i, ] <- (neg.out$v)^2
    }
  }
  
  return(results)
}

permuteIndex <- function(n1, n2){
  i.sample <- sample(1:(n1+n2), replace=FALSE)
  return(list(i1=i.sample[1:n1], i2=i.sample[(n1+1):(n1+n2)]))
}



symmPMD <- function(x, sumabs=0.3, niter=50, v=NULL, trace=TRUE) {
  if (!isSymmetric(x)) {
    stop("x must be a symmetric matrix.")
  }
  
  x[is.na(x)] <- mean(x[!is.na(x)])
  p <- nrow(x)
  sumabsv <- sqrt(p) * sumabs
  K <- 1
  
  # initial value: first singular vector
  if (is.null(v)) {
    v <- rARPACK::eigs_sym(x, K, "LA")$vectors
  }
  
  # main algorithm (only works for K=1!)
  out <- solvePMD(x, sumabsv=sumabsv, v=v, niter=niter, trace=trace)
  
  return(list(v=out$v, d=out$d, v.init=out$v.init, sumabs=sumabs))
}



solvePMD <- function(x, sumabsv, v, niter=50, trace=TRUE) {
  if (!isSymmetric(x)) {
    stop("x must be a symmetric matrix.")
  }
  
  ## initialize
  v.init <- v
  oldv <- rnorm(ncol(x))
  
  ## iterative updates
  for (iter in 1:niter) {
    if (sum(abs(oldv - v)) <= 1e-07) {
      break
    }
    if (trace) {
      cat(iter, fill=FALSE)
    }
    
    oldv <- v
    
    ## v <- S(X*oldv, lamv) / ||S(X*oldv, lamv)||_2
    ## where S() is soft threshold, lamv >= 0 is such that ||v||_1=sumabsv
    argv <- x %*% v
    lamv <- BinarySearch(argv, sumabsv)
    sv <- soft(argv, lamv)
    v <- matrix(sv / l2n(sv), ncol=1)
  }
  
  ## optimal objective criteria
  d <- as.numeric(t(v) %*% (x %*% v))
  if (trace) {
    cat(fill=TRUE)
  }
  
  return(list(d=d, v=v, v.init=v.init))
}


BinarySearch <- function (argv, sumabsv, maxiter=150) {
  ## no thresholding
  if (l2n(argv) == 0 || sum(abs(argv/l2n(argv))) <= sumabsv) {
    return(0)
  }
  
  ## binary search
  lam1 <- 0
  lam2 <- max(abs(argv)) - 1e-05
  iter <- 1
  while (iter < maxiter) {
    sv <- soft(argv, (lam1 + lam2)/2)
    if (sum(abs(sv/l2n(sv))) < sumabsv) {
      lam2 <- (lam1 + lam2)/2
    }
    else {
      lam1 <- (lam1 + lam2)/2
    }
    if ((lam2 - lam1) < 1e-06)
      return((lam1 + lam2)/2)
    iter <- iter + 1
  }
  return((lam1 + lam2)/2)
}

l2n <- function (vec) {
  a <- sqrt(sum(vec^2))
  if (a == 0)
    a <- 0.05
  return(a)
}


soft <- function (x, d) {
  return(sign(x) * pmax(0, abs(x) - d))
}


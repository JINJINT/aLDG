library(ESCO)
library(VineCopula)
#=============================================================#
#                   simple low rank mixture model             #
#=============================================================#

 
othn<-function(p,d){
  ot = randortho(p, type = "orthonormal")
  ot = ot[,1:d]
  return(ot)
}

 
randcov<-function(p, sigdim, s, sigma){
  mat = matrix(0,nrow=p,ncol=p)
  U = othn(length(sigdim),s)
  mat[sigdim,sigdim] = U%*%t(U)
  diag(mat)[-sigdim] = sigma
  return(mat)
}

 
mixlowrankmv<-function(n, p, dlist, slist, alphalist, mulist, sigma){
  M = length(dlist)
  dat = matrix(rnorm(n*p, 0, sigma), nrow = p, ncol = n)
  cnames = c()
  rnames = rep('None',p)

  nlist = floor(alphalist*n)
  nlist[M] = n-sum(nlist[1:(M-1)])
  nlistdim = c(0,cumsum(nlist))

  for(m in 1:M){
    mu = mulist[[m]]
    U = othn(length(dlist[[m]]),slist[m])
    X = rmvnorm(n = nlist[m], mu, diag(1,slist[m]))
    Y = U%*%t(X)
    dat[dlist[[m]],(nlistdim[m]+1):nlistdim[m+1]] = Y
    cnames = c(cnames,rep(paste0('group ',m),nlist[m]))
    rnames[dlist[[m]]] = paste0('group ',m)
  }
  colnames(dat) = cnames
  rownames(dat) = rnames
  return(dat)
}

 
mixlowranksc<-function(n, p, dlist, slist, alphalist,
                 alist, blist,
                 noise_alist, noise_blist,
                 inside = FALSE,
                 ufrom = 'none', visual=FALSE){
  M = length(dlist)
  rank = length(alist[[1]])
  ranknoise = length(noise_alist[[1]])
  data = matrix(0,p,n)
  gc_z = matrix(0,p,n)
  for(k in 1:ranknoise){
    lambda_z = rgamma(p, shape=noise_alist[k], rate=noise_blist[k])
    L_z = rlnorm(n, meanlog = 1, sdlog = 0.3)
    #L_z = rep(exp(6),n)
    gc_z = gc_z + lambda_z%*%t(L_z)
  }
  for(i in 1:p){
    for(j in 1:n){
      data[i,j] = rpois(1,gc_z[i,j]/ranknoise)
    }
  }
  nlist = floor(alphalist*n)
  if(M>1)nlist[M] = n-sum(nlist[1:(M-1)])
  nlistdim = c(0,cumsum(nlist))

  if(inside)info = paste0(ufrom,'_inside','_rank',rank)
  else info = paste0(ufrom,'_rank',rank)

  cnames = c()
  rnames = rep('None',p)
  for(m in 1:M){
    avec = alist[[m]]
    bvec = blist[[m]]
    if(ufrom!="none"){
      if(ufrom=="othn")U = abs(othn(length(dlist[[m]]),slist[m]))
      if(ufrom=="unif"){
        U = matrix(runif(length(dlist[[m]])*slist[m]), nrow=length(dlist[[m]]))
        U = t(t(U)/sqrt(colSums(U^2)))}
      if(inside){
        X = matrix(0,length(dlist[[m]]),nlist[m])
      }
      else{
        X = matrix(0,slist[m],nlist[m])
      }
      gc = matrix(0,slist[m],nlist[m])

      for(k in 1:rank){
        lambda = rgamma(n = slist[m],shape=avec[k],rate=bvec[k])
        L = rlnorm(nlist[m], meanlog = 1, sdlog = 0.3)
        #L = rep(exp(6),nlist[m])
        gc = gc+ lambda%*%t(L_z[(nlistdim[m]+1):nlistdim[m+1]])
      }
      if(inside)gc = U%*%gc
      for(i in 1:nrow(gc)){
        for(j in 1:nlist[m]){
          X[i,j] = rpois(1,gc[i,j]/rank)
        }
      }
      if(inside) Y = X
      else Y = U%*%X
    }
    else{
      dd = length(dlist[[m]])
      X = matrix(0,dd,nlist[m])
      gc =matrix(0,dd,nlist[m])
      for(k in 1:rank){
        lambda = rgamma(n = dd,shape=avec[k],rate=bvec[k])
        L = rlnorm(nlist[m], meanlog = 1, sdlog = 0.3)
        #L = rep(exp(6),nlist[m])
        gc = gc+ lambda%*%t(L_z[(nlistdim[m]+1):nlistdim[m+1]])
      }
      for(i in 1:dd){
        for(j in 1:nlist[m]){
          X[i,j] = rpois(1,gc[i,j]/rank)
        }
      }
      Y = X
    }

    data[dlist[[m]],(nlistdim[m]+1):nlistdim[m+1]] = Y
    cnames = c(cnames,rep(paste0('group ',m),nlist[m]))
    rnames[dlist[[m]]] = paste0('group ',m)
  }
  colnames(data) = cnames
  rownames(data) = rnames

  if(visual){
    geneinfo = as.factor(rownames(data))
    cellinfo = as.factor(colnames(data))
    rownames(data) = 1:p
    colnames(data) = 1:n
    heatdata(list(data[rev(1:p),]),
             cellinfo = data.frame(cell = 1:ncol(data),
                                   newcelltype =cellinfo),
             geneinfo = data.frame(cell = 1:nrow(data),
                                   newcelltype =geneinfo)[rev(1:p),],
             norm = FALSE, rowv = NA, colv = NA, color = "-YlGnBu:100", dirname = paste0('./plots/data',info))

    normdata = log2(10^3*t(t(data)/colSums(data))+1)
    corr = cor(t(normdata))
    heatgcn(list(corr), ord = 1:p, maxgcn = 1, mingcn = -1,dirname = paste0('./plots/corr',info))
    ccc = 0
    ccn = 0
    dd = length(dlist[[1]])
    for(i in 1:(dd-1)){
      for(j in (i+1):dd){
        if(sum(normdata[i,]>0&normdata[j,]>0)>floor(n*0.6)){
          ccc=ccc+1
          if(ccc<2){
            plot2d(normdata[i,], normdata[j,], index = cnames,
                   type = paste0('./ss_high',info,'_',ccc))}
        }
        if(sum(normdata[i,]>0&normdata[j,]>0)<floor(n*0.4)){
          ccn=ccn+1
          if(ccn<2){
            plot2d(normdata[i,], normdata[j,], index = cnames,
                   type = paste0('./ss_low',info,'_',ccn))}
        }
      }
    }
    ccc = 0
    ccn = 0
    for(i in 1:dd){
      for(j in (dd+1):p){
        if(sum(normdata[i,]>0&normdata[j,]>0)>floor(n*0.6)){
          ccc=ccc+1
          if(ccc<2){
            plot2d(normdata[i,], normdata[j,], index =cnames,
                   type = paste0('./sn_high',info,'_',ccc))}
        }
        if(sum(normdata[i,]>0&normdata[j,]>0)<floor(n*0.4)){
          ccn=ccn+1
          if(ccn<2){
            plot2d(normdata[i,], normdata[j,], index =cnames,
                   type = paste0('./sn_low',info,'_',ccn))}
        }
      }
    }
    ccc = 0
    ccn = 0
    for(i in (dd+1):p-1){
      for(j in (dd+1):p){
        if(sum(normdata[i,]>0&normdata[j,]>0)>floor(n*0.6)){
          ccc=ccc+1
          if(ccc<2){
            plot2d(normdata[i,], normdata[j,], index = cnames,
                   type = paste0('./nn_high',info,'_',ccc))
          }
        }
        if(sum(normdata[i,]>0&normdata[j,]>0)<floor(n*0.4)){
          ccn=ccn+1
          if(ccn<2){
            plot2d(normdata[i,], normdata[j,], index = cnames,
                   type = paste0('./nn_low',info,'_',ccn))}
        }
      }
    }
  }
  return(data)
}


#=============================================================#
#                   simple copula one group model             #
#=============================================================#
 
mixcopula<-function(distrlist=c('norm'), n=100, d = 2,
                 prob = c(0.5, 0.5), corrlist = list(),
                 para1list = list(), para2list = list()){
  if(sum(prob)!=1)cat("warning: the mixture probability does not sum to one!")
  if(length(para1list)!=length(prob))cat('warning: the length of parameter list does not equal to the length of mixtures')
  groupnumber = floor(n*prob)
  groupnumber[length(groupnumber)] = n-sum(groupnumber[1:(length(groupnumber)-1)])
  data = c()
  index = c()
  for(k in 1:length(prob)){
    para1 = para1list[[k]] # k of d
    para2 = para2list[[k]] # k of d
    corr = corrlist[[k]] # k of d*d
    index = c(index, rep(k,groupnumber[k]))
    data = cbind(data, onecopula(Rho=corr, d=d, n=groupnumber[k],
                                distr=distrlist[k],
                                para1=para1, para2=para2))
  }
  return(list(data, index))# d*n
}


# simulate one copula group
 
onecopula<-function(Rho, d,n,distr='norm', para1=NULL, para2=NULL){
  if(!is.positive.definite(Rho)){
    cat('input corr matrix is not positive definite, thus doing correction..')
    Rho = makespd(Rho)
  }
  Col = chol(Rho)
  copula = matrix(rnorm(d*n), ncol = n)
  copula = t(Col) %*% copula
  copula = pnorm(copula)

  single<-function(i){
    if(distr=='unif'){
      if(is.null(para1))para1 = 0
      if(is.null(para2))para2 = 1
      x=qunif(copula[,i], min = para1, max = para2)
    }
    if(distr=='norm'){
      if(is.null(para1))para1 = 0
      if(is.null(para2))para2 = 1
      x=qnorm(copula[,i], mean = para1, sd = para2)
    }
    if(distr=='lnorm'){
      if(is.null(para1))para1 = 0
      if(is.null(para2))para2 = 1
      x=exp(qnorm(copula[,i], mean = para1, sd = para2))
    }
    if(distr=='exp'){
      if(is.null(para1))para1 = 1
      x=qexp(copula[,i], rate = para1)
    }
    if(distr=='pois'){
      if(is.null(para1))para1 = 10
      x=qpois(copula[,i], lambda = para1)
    }
    if(distr=='gamma'){
      if(is.null(para1))para1 = 1
      if(is.null(para2))para2 = 2
      x= qgamma(copula[,i], shape = para1, rate = para2)
    }
    if(distr=='nbinom'){
      if(is.null(para1))para1 = 20
      if(is.null(para2))para2 = 0.3
      x=qnbinom(copula[,i], mu = para1, size = para2)
    }
    if(distr=='cauchy'){
      if(is.null(para1))para1 = 0
      if(is.null(para2))para2 = 1
      x=qcauchy(copula[,i], location = para1, scale = para2)
    }
    return(x)
  }

  data = c()
  for(i in 1:n){
    data = cbind(data,single(i))
  }
  return(data)
}

#===============================================#
#      simple non-parametric kernel model       #
#          (no within dependence)               #
#===============================================#
# alphalist = c(0.1,0.1,0.1)
# mulist = c(1,2,3)
# rlist = c(0.1,0.1,0.1)
# dlist = list(1:10, 5:15, 10:20)
 
mixkernel <-function(name, alphalist,
                                 mulist, rlist,
                                 n, p, dlist, fixed = TRUE){

  n_sig = ceiling(alphalist*n)
  print(n_sig)
  n_noise = n-sum(n_sig)

  kerdens_noise = kerdens(name, 0, 1)

  data_signoise = c()
  labels = c()
  for(i in 1:length(alphalist)){
    kerdens_sig = kerdens(name, mulist[i], rlist[i])
    datasig = matrix(0,p,n_sig[i])

    datasig[dlist[[i]],] = vapply(seq_len(n_sig[i]),
                                  function(x)samp_kerdens(kerdens_sig, length(dlist[[i]]),
                                                          xmin = mulist[i]-rlist[i], xmax = mulist[i]+rlist[i]), numeric(length(dlist[[i]])))
    if(p-length(dlist[[i]])>0){
      datasig[-dlist[[i]],] = vapply(seq_len(n_sig[i]),
                                     function(x)samp_kerdens(kerdens_noise, p-length(dlist[[i]])), numeric(p-length(dlist[[i]])))
    }
    data_signoise = cbind(data_signoise, datasig)
    labels = c(labels,rep(paste0("group",i),n_sig[i]))
  }
  data_purenoise = vapply(seq_len(n_noise),
                          function(x)samp_kerdens(kerdens_noise, p), numeric(p))
  data = cbind(data_signoise,data_purenoise)
  labels = c(labels,rep("noise",n_noise))
  colnames(data) = labels
  return(data)
}

 
kerdens<-function(name, mu, r){
  if(r<=0) warning("r should be bigger than 0!")
  if(name=="box"){
    dens<-function(x)1/2*(abs(x)<=1)
  }
  if(name=="poly2"){
    dens <-function(x)3/4*(1-x^2)*(abs(x)<=1)
  }
  if(name=="poly3"){
    dens <-function(x)70/81*(1-abs(x)^3)^3*(abs(x)<=1)
  }
  dens_final <-function(x)1/r*dens((x-mu)/r)
  return(dens_final)
}

 
samp_kerdens <- function(dens, n, xmin=-1, xmax=1, lower = NULL) {
  if(is.null(lower))lower = xmin

  # assuming dens is symmetric and decreases with |x|
  ymax <- dens((xmax+xmin)/2)
  ymin <- 0
  values <- c()
  nsel <- 0
  while(nsel < n) {
    x <- runif(1e4, xmin, xmax)
    y <- runif(1e4, ymin, ymax)
    sel <- y < dens(x) & x >= lower
    nsel <- nsel + sum(sel)
    values <- c(values, x[sel])
  }
  values <- values[seq_len(n)]
  return(values)
}








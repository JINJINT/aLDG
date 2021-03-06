#' Function to compute the bivariate dependency
#'
#' Estimate simulation parameters for the esco simulation from a real
#' dataset. See the individual estimation functions for more details on how this
#' is done.
#'
#' @param x: a numeric vector
#' @param y: a numeric vector of the same length as x
#' @param sx: a list containing the sorting information, defualt to NULL,
#'            if not NULL, then should have the format
#'            sx = list(rank = vector of value rank, dist = matrix of pairwise distance)
#' @param sy: similar as the definition for sx, but for vector y
#' @param methods: a vector of strings, should be set of the following:
#'            -- 'Pearson': Pearaon's correlation
#'            -- 'Spearman': Spearman's correaltion
#'            -- 'Kendall': Kendall's \eqn{\tau}
#'            -- 'TauStar': Kendall's \eqn{\tau^\star}
#'            -- 'dCor': distance correlation
#'            -- 'HSIC': Hilbert-Schmidt Independence Criterion (normalized to 0,1)
#'            -- 'HoeffD': Hoeffding's D
#'            -- 'HHG': HHG
#'            -- 'MIC': Maximal Information Coefficient
#'            -- 'aLDG': averaged Local Density Gap
#' @param all: logical value, if true then compute for all methods listed
#'
#' @param thred: (aLDG specific parameters) only the data > thred will be used for aLDG
#' @param hx,hy: (aLDG specific parameters) vector consist provided bandwidth for x (or y), default to NULL
#' @param band: (aLDG specific parameters) which method to use for automatically calculating bandwidth (if hx and hy are NULL)
#'              -- 'fix': the window around x is \eqn{B_x = (x-wd \cdot h_n, x + wd \cdot h_n}
#'              -- 'ada': the adaptive window around x is \eqn{P(X \in B_x \mid Y=y) = qd}
#' @param wd: (aLDG specific parameters) the coefficient before the bandwidth h = wd*h_n
#' @param qd: (aLDG specific parameters) the quantile window to compute the adaptive bandwidth \eqn{h_x = sd(x(qd))}
#' @param opt: (aLDG specific parameters) for 'fix' bandwidth calculation,
#'             if TRUE, we use the theoretical optimal rate \eqn{h_n = sd(x) n^{-1/6}};
#'             if FALSE, we use \eqn{h_n = sd}
#' @param stat: (aLDG specific parameters) what statistics to use for \eqn{T_i} in aLDG
#' @param cutoff: (aLDG specific parameters) the coef c in threshold choice for aLDG:
#'             \deqn{aLDG = \frac{1}{n} \sum_{i=1}^n ( T_i > \frac{\Phi^{-1}(1-n^{-c})}{n^{1/3}})}
#' @return a vector
#' @examples
#' # Load example data
#' ans = bidep(runif(10),runif(10),methods='aLDG')
#' @rdname bidep
#' @import TauStar HHG energy dHSIC Hmisc minerva
#' @export
bidep<-function(x, y, sx = NULL, sy = NULL, methods=NULL, all = FALSE,
                      thred=-Inf, hx=NULL, hy=NULL,
                      wd = 1, aldg_thred='error'){

  if(length(x)!=length(y)){
    print('The length of x and y should be the same!')
    return(NULL)
  }

  corr = c()
  n = length(x)

  if(all)methods = c('Pearson','Spearman','Kendall','TauStar','dCor','HSIC','HoeffD','HHG','MIC','MRank','aLDG')

  for(method in methods){
    if(method %in% c('Pearson','Spearman','Kendall')){
      corr[method] = cor(x, y, method = tolower(method))
    }
    if(method=='TauStar'){
      a= TauStar::tStar(x, y)
      # theoretically, taustar lies in range [-1/3,2/3], where 0 represents independence,
      # therefore we do the following to rescale it to [-1,1] with 0 still represents independence
      if(a<=0)corr[method]=a*3
      if(a>0)corr[method]=a*3/2
    }
    if(method=='dCor'){
      corr[method] = dcor(x,y)
    }
    if(method=='HSIC'){
      # we normalize in this way such that HSIC is in [0,1]
      corr[method] = dHSIC::dhsic(x,y)$dHSIC/sqrt(dhsic(x,x)$dHSIC * dhsic(y,y)$dHSIC)
    }
    if(method=='HoeffD'){
      corr[method]= Hmisc::hoeffd(x,y)$D[1,2]
    }
    if(method=='HHG'){
      result = HHG::hhg.univariate.ind.stat(x, y)$statistic
      corr[method] = result[length(result)]
    }
    if(method=='MIC'){
      corr[method]=minerva::mine(x, y)$MIC
    }
    # if(method=='MRank'){
    #   corr[method]=incSubseq(x,y,k=5)/choose(n,5)
    # }
    if(method=='aLDG'){
      ans=aldg(x, y, sx = sx, sy = sy, thred = thred, wd = wd, 
                       trials=max(floor(1000/n),5), chooset = aldg_thred)
      corr[method] = ans$val 
    }
  }
  return(corr)
}




#' Function to compute the multivariate dependency matrix
#' @param data: feature by sample matrix
#' @param methods: a vector of strings, could be subset of c('pearson','kendall','taustar','dcor','hsic','hoeffd','ssd')
#' @param all: logical value, if true then compute for all methods in c('pearson','kendall','taustar','dcor','hsic','hoeffd','ssd')
#' @param thred: only the data > thred will be used
#' @param wd: the coefficient before the bandwidth h = wd*h_n
#' @param qd: the quantile window to compute the adaptive bandwidth h_x = sd(x(qd))
#' @param info: integer denote the trial number
#' @param ncores: number of cores to use
#' @param norm: whether scale the value to [0,1]
#' @param abs: whether take absolute value of the value
#' @param extrainfo: some setting related information
#' @return a list of matrix
#' @examples
#' # Load example data
#' ans = matdep(runif(100,5,20), methods='ssd')
#' @rdname matdep
#' @export
matdep<-function(data, methods=NULL, 
                          all = FALSE,
                          thred = -Inf, 
                          wd = 1, 
                          trial = 0,
                          ncores=NULL, 
                          norm = FALSE, abs = FALSE,
                          info = '', dir='./dat/'){
  
  if(!dir.exists(dir))dir.create(dir)
  
  thisfilename = paste0(dir,'corrmat_',info, trial,'.rds')
  
  if(!file.exists(thisfilename)){
    #print(thisfilename)
    if(all)methods = c('Pearson','Spearman','Kendall','TauStar','dCor','HSIC','HoeffD','MIC','MRank','aLDG')
    p = nrow(data)
    n = ncol(data)
    
    sort_single<-function(gene){
      sortx = sort(data[gene,], index.return = TRUE)
      sx = sortx$x
      rx = sortx$ix
      dx = matrix(0,n,n)
      drx = sapply(1:n, function(i)abs(sx-sx[i]))
      dx[rx,rx] = drx
      return(list(rank=rx,dist=dx))
    }
    
    if (!is.null(ncores)){
      cl <- makeCluster(ncores)
      registerDoSNOW(cl)
      
      pb <- progress_bar$new(
        format = "progress = :letter [:bar] :elapsed | eta: :eta", total = p, width = 60)
      
      progress <- function(n){
        pb$tick(tokens = list(letter = rep("", p)[n]))
      }
      opts <- list(progress = progress)
      
      sortlist <- foreach(i = 1:p, .options.snow = opts)%dopar% {
                            return(sort_single(i))
                          }
      
      
      stopCluster(cl)
    }else{
      sortlist <- lapply(1:p, sort_single) # n*p
    }
    
    index = expand.grid(c(1:p),c(1:p))
    index = index[which(index[,1]<index[,2]),]
    total <- nrow(index)
    
    compute_single<-function(k,verbose = FALSE){
      i = index[k,1]
      j = index[k,2]
      x = data[i,]
      y = data[j,]
      sx = sortlist[[i]]
      sy = sortlist[[j]]
      allvec = bidep(x, y, sx=sx, sy=sy,  wd = wd,  
                            thred=thred, methods = methods)
      return(allvec)
    }
    
    if(!is.null(ncores)){
      cl <- makeCluster(ncores)
      registerDoSNOW(cl)
      
      pb <- progress_bar$new(
        format = "progress = :letter [:bar] :elapsed | eta: :eta", total = total, width = 60)
      
      progress <- function(n){
        pb$tick(tokens = list(letter = rep("", total)[n]))
      }
      
      opts <- list(progress = progress)
      
      allmat <- foreach(i = 1:total, .combine = cbind, .options.snow = opts,
                        .export=c('ldgsingle','kernallist','aldg','bidep'),
                        .packages = c('dHSIC','TauStar','minerva','Hmisc','energy')
      )%dopar% {
        return(compute_single(i))
        }
      stopCluster(cl)
    }
    else{
      allmat <- sapply(1:total,compute_single, verbose=TRUE)
    }
  }
  else{
    print('file already existed...extracting results...')
    allmat <- readRDS(thisfilename)
  }
  
  corrmat = list()
  methods = rownames(allmat)
  for(idx in 1:nrow(allmat)){
    method = methods[idx]
    corrmat[[method]] = matrix(0,p,p)
    corrmat[[method]][upper.tri(corrmat[[method]], diag = FALSE)] = allmat[idx,]
    corrmat[[method]] = t(corrmat[[method]]) + corrmat[[method]]
  }
  methods = names(corrmat)
  for(idx in 1:length(methods)){
    method = methods[idx]
    if(abs){
      corrmat[[method]] = abs(corrmat[[method]])
    }
    if(norm){
      corrmat[[method]] = corrmat[[method]]/max(abs(corrmat[[method]]), na.omit = TRUE, nan.omit = TRUE)
    }
  }
  return(corrmat)
}



#==== helper functions====#
#' Function to compute the aLDG measure
#'
#' @param x: a numeric vector
#' @param y: a numeric vector of the same length as x
#' @param sx: a list containing the sorting information, defualt to NULL,
#'            if not NULL, then should have the format
#'            sx = list(rank = vector of value rank, dist = matrix of pairwise distance)
#' @param sy: similar as the definition for sx, but for vector y
#' @param thred: (aLDG specific parameters) only the data > thred will be used for aLDG
#' @param hx,hy: (aLDG specific parameters) vector consist provided bandwidth for x (or y), default to NULL
#' @param band: (aLDG specific parameters) which method to use for automatically calculating bandwidth (if hx and hy are NULL)
#'              -- 'fix': the window around x is \eqn{B_x = (x-wd \cdot h_n, x + wd \cdot h_n}
#'              -- 'ada': the adaptive window around x is \eqn{P(X \in B_x \mid Y=y) = qd}
#' @param wd: (aLDG specific parameters) the coefficient before the bandwidth h = wd*h_n
#' @param qd: (aLDG specific parameters) the quantile window to compute the adaptive bandwidth \eqn{h_x = sd(x(qd))}
#' @param opt: (aLDG specific parameters) for 'fix' bandwidth calculation,
#'             if TRUE, we use the theoretical optimal rate \eqn{h_n = sd(x) n^{-1/6}};
#'             if FALSE, we use \eqn{h_n = sd}
#' @param stat: (aLDG specific parameters) what statistics to use for \eqn{T_i} in aLDG
#' @param cutoff: (aLDG specific parameters) the coef c in threshold choice for aLDG:
#'             \deqn{aLDG = \frac{1}{n} \sum_{i=1}^n ( T_i > \frac{\Phi^{-1}(1-n^{-c})}{n^{1/3}})}
#' @return a vector
#' @examples
#' # Load example data
#' ans = aldg(runif(10),runif(10))
#' print(ans$aldg)
#' @rdname aldg
#' @export
aldg<-function(x, y, sx = NULL, sy = NULL,
               thred = -Inf, hx=NULL, hy=NULL, 
               wd = 1,  cutoff = 1, chooset = 'error',
               trials=5){
  
  n = length(x)
  
  if(thred > -Inf){
    id = c()
    idx = which(x>thred)
    idy = which(y>thred)
    id = intersect(idx,idy)
  }else{
    id = 1:n
  }
  tinflectlist = c()
  tmiselist = c()
  nt = length(id)
  if(nt>10){
    xt = x[id]
    yt = y[id]
    rid = ((n-nt)+1):n
    result = rep(0,n)
    if(!is.null(sx)&&!is.null(sy)){
      rx = sx[['rank']]
      rx = rx[which(rx %in% rid)]-(n-nt)
      ry = sy[['rank']]
      ry = ry[which(ry %in% rid)]-(n-nt)
      re = ldgsingle(xt, yt, 
                     rx=rx,
                     ry=ry,
                     dx=sx[['dist']][rid,rid],
                     dy=sy[['dist']][rid,rid],
                     hx=hx, hy=hy,
                     wd=wd, kernel = 'box',
                     bandlist = c('fix'),
                     chooset = FALSE,
                     stat = 'normgap', opt=TRUE)
      t_miseold = re[[2]]
      if(trials>0){
        for(i in 1:trials){
          renull = ldgsingle(xt, yt[sample(1:nt,nt)], 
                             hx=hx, hy=hy,
                             wd=wd, kernel = 'box',
                             bandlist = c('fix'),
                             chooset = FALSE,
                             stat = 'normgap', opt=TRUE)
          Tvec = renull[[1]][['fix']]
          tmiselist[i] = max(Tvec)
          sTvec = sort(Tvec)
          sTvec = sTvec[which(sTvec>0)]
          if(length(sTvec)<5)tinflectlist[i]=NA
          else{
            s1gapvec = sTvec[1:(length(sTvec)-2)]-sTvec[2:(length(sTvec)-1)]
            s2gapvec = sTvec[2:(length(sTvec)-1)]-sTvec[3:(length(sTvec))]
            sgapvec = s2gapvec-s1gapvec
            nstar = localMaxima(sgapvec)
            tinflectlist[i] = max(sTvec[nstar])
            if(tinflectlist[i]==-Inf)tinflectlist[i]=NA
          }
        }
      }
    }else{
      re = ldgsingle(xt, yt, 
                     hx=hx, hy=hy,
                     wd=wd, kernel = 'box',
                     bandlist = c('fix'),
                     chooset = FALSE,
                     stat = 'normgap', opt=TRUE)
      t_miseold = re[[2]]
      if(trials>0){
        for(i in 1:trials){
          renull = ldgsingle(xt, yt[sample(1:nt,nt)], 
                             hx=hx, hy=hy,
                             wd=wd, kernel = 'box',
                             bandlist = c('fix'),
                             chooset = FALSE,
                             stat = 'normgap', opt=TRUE)
          Tvec = renull[[1]][['fix']]
          tmiselist[i] = max(Tvec)
          Tvec = Tvec[!is.na(Tvec)]
          sTvec = sort(Tvec)
          sTvec = sTvec[which(sTvec>0)]
          if(length(sTvec)<5)tinflectlist[i]=NA
          else{
            s1gapvec = sTvec[1:(length(sTvec)-2)]-sTvec[2:(length(sTvec)-1)]
            s2gapvec = sTvec[2:(length(sTvec)-1)]-sTvec[3:(length(sTvec))]
            sgapvec = s2gapvec-s1gapvec
            nstar = localMaxima(sgapvec)
            tinflectlist[i] = max(sTvec[nstar])
            if(tinflectlist[i]==-Inf)tinflectlist[i]=NA
          }
        }
      }
    }
    t_norm = qnorm(1-1/(nt)^cutoff)/(wd*sqrt(sd(xt)*sd(yt))*nt^(1/3))
    t_mise = quantile(tmiselist, prob=0.5, na.rm = TRUE)
    t_inflect = quantile(tinflectlist, prob=0.5,na.rm = TRUE)
    
    t = t_norm
    if(chooset=='inflect' & !is.na(t_inflect))t=t_inflect
    if(chooset=='error' & !is.na(t_mise))t=t_mise

    result[id] = re[[1]][['fix']]
    return(list(val = mean(result>t), Tvec = result, 
                t = list(inflect = t_inflect, mise=t_mise, norm=t_norm, miseold=t_miseold)))
  }else{
    return(list(val =0, Tvec = NULL, t = NULL))
  }
}


#' @export
localMaxima <- function(x){
   y <- diff(c(-Inf, x)) >= 0
   y <- cumsum(rle(y)$lengths)
   y <- y[seq.int(1, length(y), 2)]

   if (y[1]==1) {
       y <- y[-c(1)]
   }
   if(length(y)>1){
     if(y[length(y)]==length(x))y <- y[-c(length(y))]
   }
   return(y)
}

#' @export
ldgsingle<-function(x,y, dall = NULL, knn = 100, # the distance matrix of samples using all features
                    dx = NULL, dy = NULL,
                    rx = NULL, ry = NULL,
                    wd = 1, qd = 0.1,
                    kernel='box',
                    hx=NULL,hy=NULL,
                    bandlist=c('fix'),
                    opt=TRUE,
                    chooset = FALSE,
                    stat = c('normgap')){

  k = kernallist(kernel)
  n = length(x)
  if(is.null(rx)){
    sortx = sort(x, index.return = TRUE)
    rx = sortx$ix
  }
  sx = x[rx]
  if(is.null(ry)){
    sorty = sort(y, index.return = TRUE)
    ry = sorty$ix
  }
  sy = y[ry]
  if(is.null(dx)){
    dx = matrix(0,n,n)
    drx = sapply(1:n, function(i)abs(sx-sx[i]))
    dx[rx,rx] = drx
  }
  if(is.null(dy)){
    dy = matrix(0,n,n)
    dry = sapply(1:n, function(i)abs(sy-sy[i]))
    dy[ry,ry] = dry
  }

  gaplist = list()
  t=NULL
  for(band in bandlist){
    if(band=='fix'){
      if(is.null(hx))hx = wd*sd(x)
      if(is.null(hy))hy = wd*sd(y)
    }

    if(band=='ada'){
      kn = ceiling(n*qd/2)
      hx = rep(0,n)
      #print(str(x[ry[max(10-kn,1):min(10+kn,n)]]))
      hx[ry] = wd*sapply(1:n, function(i)sd(x[ry[max(i-kn,1):min(i+kn,n)]]))
      hy = rep(0,n)
      hy[rx] = wd*sapply(1:n, function(i)sd(y[rx[max(i-kn,1):min(i+kn,n)]]))
    }
    #print(str(hx))
    if(opt){
      if(!is.null(dall)){
        hx = hx*n^(-1/6)*knn^(-1/6)
        hy = hy*n^(-1/6)*knn^(-1/6)
      }
      else{
        hx = hx*n^(-1/6)
        hy = hy*n^(-1/6)
      }
    }

    fx = k(dx/hx)/hx
    fy = k(dy/hy)/hy
    if(!is.null(dall)){
      fx = fx*1*(dall<knn)
      fy = fy*1*(dall<knn)
    }
    fxy = fx*fy
    if(is.null(dall)){
      fxyhat =  rowMeans(fxy)
      fxhat = rowMeans(fx)
      fyhat = rowMeans(fy)
    }
    else{
      fxyhat =  rowMeans(fxy)*n/knn
      fxhat = rowMeans(fx)*n/knn
      fyhat = rowMeans(fy)*n/knn
    }
    if(stat=='gap'){
      gap = fxyhat - fxhat*fyhat
    }
    if(stat=='probnormgap'){
      gap = (4*hx*hy*(fxyhat - fxhat*fyhat)*sqrt(length(x)-1))/sqrt(2*hx*fxhat*(1-2*hx*fxhat)*2*hy*fyhat*(1-2*hy*fyhat))
    }
    if(stat=='normgap'){
      gap = (fxyhat - fxhat*fyhat)/sqrt(fxhat*fyhat)
    }
    if(stat=='probmi'){
      gap = hx*hy*fxyhat*log(fxyhat/fxhat*fyhat)
    }
    if(stat=='locratio'){
      gap = fxyhat/fxhat*fyhat
    }
    if(stat=='locgauss'){
      a = localgauss(x,y,b1=sd(x),b2=sd(y),xy.mat=cbind(x,y))
      gap = a$par.est[,5]
    }
    if(stat=='loclinear'){
      Ay = colSums(x*fx)/colSums(fy)
      Ax = colSums(y*fx)/colSums(fx)
      mx = mean(x)
      my = mean(y)
      vx = sd(x)
      vy = sd(y)
      gap = (cor(x,y) + (mx-Ay)*(my-Ax)/vx*vy)/(sqrt(1+(mx-Ay)^2/vx^2)*sqrt(1+(my-Ax)^2/vy^2))
    }
    if(stat=='locdep'){
      gap = (colSums(x*y*fxy) - colSums(x*fxy)*colSums(y*fxy)/colSums(fxy))/(hx*hy/12*colSums(fxy))
    }
    if(stat=='mi'){
      gap = fxyhat*log(fxyhat/fxhat*fyhat)
    }
    if(stat=='joint'){
      gap = fxyhat
    }
    if(stat=='fx'){
      gap = fxhat
    }
    if(stat=='fy'){
      gap = fyhat
    }
    if(stat=='fyfy'){
      gap = fyhat*fyhat
    }
    if(stat=='chi'){
      gap = (fxyhat - fxhat*fyhat)/(fxhat*fyhat)
    }

    gaplist[[band]] = gap

    if(chooset & band=='fix'){
      pmatlist = list(dx, dy)
      hlist = c(hx, hy)
      t = sqrt(mise(pmatlist, hlist, kernel))
    }
  }
  return(list(gaplist, t))
}



kernallist<-function(method){
  if(method=='box'){
    k<-function(d){
      return(0.5*(abs(d)<1))
    }
  }
  if(method=='gauss'){
    k<-function(d){
      return(dnorm(d))
    }
  }
  return(k)
}



# dmatlist: a list distance matrix for x_1, \dots, x_d (d-dimensional random variable)
# kernel: choice of kernel
# hlist: a list bandwidth vectors for x_1, \dots, x_d (d-dimensional random variable)
#' @export
mise<-function(dmatlist, hlist, kernel){
  if(kernel == 'box'){
    sigma = 0.5
  }
  if(kernel == 'gauss'){
    sigma = 1
  }
  d = length(hlist)
  k = kernallist(kernel)

  g <- function(x){
    return(
      integrate(function(z){
        return(k(z)*k(z+x))
      },lower=min(-1-x,-1),upper=max(1-x,1))$value
    )
  }
  n=nrow(dmatlist[[1]])
  pmatlist = lapply(1:d, function(i)dmatlist[[i]]/hlist[i])
  ans1 = lapply(pmatlist, function(mat)apply(mat, 1:2, g))
  ans2 = lapply(pmatlist, function(mat)sapply(diag(mat), g))
  ans1k = lapply(pmatlist, function(mat)apply(mat, 1:2, k))
  ans2k = lapply(pmatlist, function(mat)sapply(diag(mat), k))
  anslist = sapply(1:d, function(i)2*((sum(ans1[[i]])-sum(ans2[[i]]))/(n^2) -(sum(ans1k[[i]]) -sum(ans2k[[i]]))/(n*(n-1)) + sigma/n)/(hlist[i]))

  rans1 = Reduce('*', ans1)
  rans2 = Reduce('*', ans2)
  rans1k = Reduce('*', ans1k)
  rans2k = Reduce('*', ans2k)
  ans = (sum(rans1)-sum(rans2))/(n^2) - (sum(rans1k)-sum(rans2k))/(n*(n-1))
  ans = ans + sigma^d/n
  ansall = 2*ans/(prod(hlist))

  anslist = c(Reduce('*', anslist), ansall)
  return(abs(max(anslist)))
}





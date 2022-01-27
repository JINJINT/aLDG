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
#'            -- 'pearson': Pearaon's correlation
#'            -- 'spearman': Spearman's correaltion
#'            -- 'kendall': Kendall's \eqn{\tau}
#'            -- 'taustar': Kendall's \eqn{\tau^\star}
#'            -- 'dcor': distance correlation
#'            -- 'hsic': Hilbert-Schmidt Independence Criterion (normalized to 0,1)
#'            -- 'hoeffd': Hoeffding's D
#'            -- 'hhg': HHG
#'            -- 'rank': Matching Ranks
#'            -- 'mic': Maximal Information Coefficient
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
#'             if TRUE, we use the theoretical optimal rate \eqn{h_n = sd(x) n^{-1/6};
#'             if FALSE, we use \eqn{h_n = sd}
#' @param stat: (aLDG specific parameters) what statistics to use for \eqn{T_i} in aLDG
#' @param cutoff: (aLDG specific parameters) the coef c in threshold choice for aLDG:
#'             \deqn{aLDG = \frac{1}{n} \sum_{i=1}^n ( T_i > \frac{\Phi^{-1}(1-n^{-c})}{n^{1/3}})}
#' @return a vector
#' @examples
#' # Load example data
#' ans = bidep(runif(10),runif(10),methods='aldg')
#' @rdname bidep
#' @export
bidep<-function(x, y, sx = NULL, sy = NULL, methods=NULL, all = FALSE,
                      thred=-Inf, hx=NULL, hy=NULL, band = 'fix',
                      wd = 0.5, qd = 0.1, opt = TRUE, stat = 'normgap', cutoff = 1){

  if(length(x)!=length(y)){
    print('The length of x and y should be the same!')
    return(NULL)
  }

  corr = c()
  n = length(x)

  if(thred > -Inf){
    id = c()
    idx = which(x>thred)
    idy = which(y>thred)
    id = intersect(idx,idy)
  }else{
    id = 1:n
  }

  if(all)methods = c('pearson','spearman','kendall','taustar','dcor','hsic','hoeffd','hhg','mic','rank','aldg')

  for(method in methods){
    if(method %in% c('pearson','spearman','kendall')){
      corr[method] = cor(x, y, method = method)
    }
    if(method=='taustar'){
      a= tStar(x, y)
      # theoretically, taustar lies in range [-1/3,2/3], where 0 represents independence,
      # therefore we do the following to rescale it to [-1,1] with 0 still represents independence
      if(a<=0)corr['taustar']=a*3
      if(a>0)corr['taustar']=a*3/2
    }
    if(method=='dcor'){
      corr[method] = dcor(x,y)
    }
    if(method=='hsic'){
      # we normalize in this way such that HSIC is in [0,1]
      corr['hsic'] = dhsic(x,y)$dHSIC/sqrt(dhsic(x,x)$dHSIC * dhsic(y,y)$dHSIC)
    }
    if(method=='hoeffd'){
      corr['hoeffd']= hoeffd(x,y)$D[1,2]
    }
    if(method=='hhg'){
      result = hhg.univariate.ind.stat(x, y)$statistic
      corr['hhg'] = result[length(result)]
    }
    if(method=='mic'){
      corr['mic']=mine(x, y)$MIC
    }
    if(method=='rank'){
      corr['rank']=incSubseq(x,y,k=5)/choose(n,5)
    }
    if(method=='aldg'){
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
          re = ldgsingle(xt, yt, NULL,
                         rx=rx,
                         ry=ry,
                         dx=sx[['dist']][rid,rid],
                         dy=sy[['dist']][rid,rid],
                         hx=hx, hy=hy,
                         wd=wd, kernel = 'box',
                         bandlist = c(band),
                         chooset = FALSE,
                         stat = stat, opt=opt)
        }
        else{
          re = ldgsingle(xt, yt, NULL,
                         hx=hx, hy=hy,
                         wd=wd, kernel = 'box',
                         bandlist = c(band),
                         chooset = FALSE,
                         stat = stat,opt=opt)
        }
        if(!is.null(re[[2]]))t = re[[2]]
        else{
          if(stat=='normgap'){
            t = qnorm(1-1/(nt)^cutoff)/(wd*sqrt(sd(xt)*sd(yt))*nt^(1/3))
          }
          if(stat=='probnormgap'){
            t =qnorm(1-1/(nt)^cutoff)
          }
        }
        result[id] = re[[1]][[band]]
        corr['aldg'] = mean(result>t)
      }
      else{
        corr['aldg']=0
      }
    }
  }
  return(corr)
}


#' Function to compute the pairwise dependency matrix for multivariate data
#' @param data: feature (row) by sample (column) matrix
#' @param methods: same as that in function bidep
#' @param all: same as that in function bidep, but omit 'hhg' since it takes too long
#' @param thred,cutoff,wd,qd,opt,band,stat: same as that in function bidep, parameters specific to aLDG method
#' @param trial: integer denote the trial number
#' @param ncores: number of cores to use
#' @param norm: whether scale the value to \eqn{[0,1]}
#' @param abs: whether take absolute value of the value
#' @param info: information related to the settings
#' @param dir: the name of directory to save the results
#' @return a list of matrix
#' @examples
#' Load example data
#' ans = matdep(runif(100,5,20), methods='aldg')
#' @rdname matdep
#' @export
matdep<-function(data, methods=NULL, all = FALSE,
                          thred = -Inf, cutoff = 1,
                          wd = 0.5, qd = 0.1,
                          opt=TRUE, band = 'fix',
                          stat = 'normgap',
                          trial = 0, ncores=NULL,
                          norm = FALSE, abs = FALSE,
                          info = '', dir='./dat/'){

  if(!dir.exists(dir))dir.create(dir)

  thisfilename = paste0(dir,'matdep_',info,'_trial', trial,'.rds')

  if(!file.exists(thisfilename)){
    if(all)methods = c('pearson','spearman','kendall','taustar','dcor','hsic','hoeffd','mic','rank','aldg')
    p = nrow(data)
    n = ncol(data)

    # function to do sorting for a single feature of the data
    sort_single<-function(feature){
      sortx = sort(data[feature,], index.return = TRUE)
      sx = sortx$x
      rx = sortx$ix
      dx = matrix(0,n,n)
      drx = sapply(1:n, function(i)abs(sx-sx[i]))
      dx[rx,rx] = drx
      return(list(rank=rx, dist=dx))
    }

    # sort for all features
    if (!is.null(ncores)){
      cl <- makeCluster(ncores)
      registerDoSNOW(cl)

      pb <- progress_bar$new(
        format = "progress = :letter [:bar] :elapsed | eta: :eta", total = p, width = 60)

      progress <- function(n){
        pb$tick(tokens = list(letter = rep("", p)[n]))
      }
      opts <- list(progress = progress)

      sortlist <- foreach(i = 1:p, .options.snow = opts,
                          .export=c('ldgsingle','kernallist'))%dopar% {
                            return(sort_single(i))
                          }
      stopCluster(cl)
    }else{
      sortlist <- lapply(1:p, sort_single) # list of p lists
    }

    # gather all the possible feature pairs
    index = expand.grid(c(1:p),c(1:p))
    index = index[which(index[,1]<index[,2]),]
    total <- nrow(index)

    # function to compute the bivariate dependence for a single pair of features
    compute_single<-function(k,verbose = FALSE){
      i = index[k,1]
      j = index[k,2]
      x = data[i,]
      y = data[j,]
      sx = sortlist[[i]]
      sy = sortlist[[j]]
      allvec = bidep(x, y, sx=sx, sy=sy,  wd = wd, qd = qd,
                            opt = opt, band=band, stat=stat,
                            thred=thred, methods = methods, cutoff=cutoff)
      return(allvec)
    }

    # compute bivariate dependence for all features pairs
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
                        .export=c('dcor','dhsic','tStar','ldgsingle','hoeffd','compute_corr','kernallist','mine')
      )%dopar% {return(compute_single(i))}
      stopCluster(cl)
    }
    else{
      allmat <- sapply(1:total,compute_single, verbose=TRUE) # a matrix of (# methods) * (# feature pairs)
    }
  }
  else{
    print('file already existed...extracting results...')
    allmat <- readRDS(thisfilename)
  }

  # convert to list of pairwise dependence matrix
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
#' @export
ldgsingle<-function(x,y, id = NULL,
                    dall = NULL, knn = 100, # the distance matrix of samples using all features
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
  t_thred=NULL
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

    if(length(id)>4 & length(id)!=n){
      if(length(hx)>1)hx = hx[id]
      if(length(hy)>1)hy = hy[id]

      if(chooset & band=='fix'){
        pmatlist = list(dx[id,id], dy[id,id])
        hlist = c(hx, hy)
        t_thred =  sqrt(mise(pmatlist, hlist, kernel))
      }

      fx = k(dx[id,id]/hx)/hx
      fy = k(dy[id,id]/hy)/hy
      if(!is.null(dall)){
        fx = fx*1*(dall[id,id]<knn)
        fy = fy*1*(dall[id,id]<knn)
      }
      fxy = fx*fy
      fxyhat =  rowMeans(fxy)
      fxhat = rowMeans(fx)
      fyhat = rowMeans(fy)
      if(stat=='gap'){
        gap = fxyhat - fxhat*fyhat
      }
      if(stat=='probnormgap'){
        gap = (4*hx*hy*(fxyhat - fxhat*fyhat)*sqrt(length(id)-1))/sqrt(2*hx*fxhat*(1-2*hx*fxhat)*2*hy*fyhat*(1-2*hy*fyhat))
      }
      if(stat=='normgap'){
        gap = (fxyhat - fxhat*fyhat)/sqrt(fxhat*fyhat)
      }
      if(stat=='probmi'){
        gap = hx*hy*fxyhat*log(fxyhat/fxhat*fyhat)
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
      gaplist[[paste0(band,'_thred')]] = rep(0,n)
      gaplist[[paste0(band,'_thred')]][id] = gap
    }
  }
  return(list(gaplist, t, t_thred))
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



# dmatlist: a list distance matrix for x and y
# kernel: choice of kernel
# hlist: a list bandwidth vectors for x and y
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

  anslist = c(anslist[[1]]*anslist[[1]], ansall)
  return(abs(max(anslist)))
}





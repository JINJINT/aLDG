# script to run on the server
library(aLDG)

args=(commandArgs(trailingOnly = TRUE))

#' Function to do permutation independence test for simulated data
#'
#' Estimate simulation parameters for the esco simulation from a real
#' dataset. See the individual estimation functions for more details on how this
#' is done.
#'
#' @param n: sample size
#' @param alpha: test level
#' @param eps: noise level
#' @param nperm: number of permutations to compute the null distribution,
#'               if set zero, then only output the value for one single case
#' @param ncores: the number of cores to use to do parallel
#' @param filename: filename related to this experiment
#' @param typelist: a list containing the names of bivariate relationships
#' @param methods: a list containing the names of methods
#' @param all: whether to do the test for all relationships and all methods
#' @param plot: whether plot the data or not, if yes, then return the data
#' @param wd: (aLDG spefici parameters) the coefficient before the bandwidth h = wd*h_n
#' @param qd: (aLDG spefici parameters) the quantile window to compute the adaptive bandwidth h_x = sd(x(qd))
#' @param band:  (aLDG spefici parameters) which method to use for calculating bandwidth
#'              -- 'fix': the theoretical rate optimal one used in the paper: h = wd*h_n
#'              -- 'ada': the adaptive bandwidth calculated as conditional quantile window: h_x = sd(x(qd))
#' @param cutoff:  (aLDG spefici parameters) the coef in threshold choice for aLDG
#' @rdname pairperm
#' @export
pairperm<-function(n=300, alpha=0.05, eps=0.3, nperm =200, ncores=10,
                   filename='', typelist=c('linear'), methods=NULL, all=TRUE, plot=FALSE,
                   wd=0.5,qd=0.1,band='fix',cutoff=1){

  if(all){
    typelist = c('indep','linear',
                 'step','ubern','circle','spiral',
                 'quad','wshape','diamond','multi',
                 'gauss30','gauss31','gauss32','gauss33',
                 'nb30','nb31','nb32','nb33')
    methods = c('pearson', 'spearman', 'kendall', 'taustar', 'dcor', 'hsic','hhg','hoeffd', 'mic','rank','ssd')
  }


  for(type in typelist){
    filenamedata = paste0(filename,type,'_result_n',n,'.rds')
    if(!file.exists(filenamedata)){
      print(filenamedata)
      if(type=='linear')dat = mgc.sims.linear(n,1,eps)
      if(type=='exp')dat = mgc.sims.exp(n,1,eps)
      if(type=='cubic')dat = mgc.sims.cubic(n,1,eps)
      if(type=='gauss')dat = mgc.sims.joint(n,1,eps)
      if(type=='step')dat = mgc.sims.step(n,1,eps)
      if(type=='quad')dat = mgc.sims.quad(n,1,eps)
      if(type=='wshape')dat = mgc.sims.wshape(n,1,eps)
      if(type=='spiral')dat = mgc.sims.spiral(n,1,eps)
      if(type=='ubern')dat = mgc.sims.ubern(n,1,eps)
      if(type=='log')dat = mgc.sims.log(n,1,eps)
      if(type=='4root')dat = mgc.sims.4root(n,1,eps)
      if(type=='sine4pi')dat = mgc.sims.sine(n,1,4,1)
      if(type=='sine16pi')dat = mgc.sims.sine(n,1,16,0.5)
      if(type=='square')dat = mgc.sims.square(n,1,eps)
      if(type=='2parab')dat = mgc.sims.2parab(n,1,eps)
      if(type=='circle')dat = mgc.sims.circle(n,1,eps)
      if(type=='ellipse')dat = mgc.sims.ellipse(n,5,eps)
      if(type=='diamond')dat = mgc.sims.diamond(n,1,eps)
      if(type=='multi')dat = mgc.sims.multi(n,1)
      if(type=='indep')dat = mgc.sims.indep(n,1)
      if(grepl('gauss', type, fixed = TRUE) & !(grepl('hard', type, fixed = TRUE))){
        vec = strsplit(type,split='')[[1]]
        groups = vec[6]
        depgroups = vec[7]
        print(groups)
        print(depgroups)
        dat = sims.gaussmix(n, k=groups, s = depgroups)
      }
      if(grepl('nb', type, fixed = TRUE) & !(grepl('hard', type, fixed = TRUE)) & ! grepl('drop', type, fixed = TRUE)){
        vec = strsplit(type,split='')[[1]]
        groups = vec[3]
        depgroups = vec[4]
        print(groups)
        print(depgroups)
        dat = sims.nbmix(n, k=groups, s = depgroups)
      }
      if(grepl('gauss', type, fixed = TRUE) & grepl('hard', type, fixed = TRUE)){
        vec = strsplit(type,split='')[[1]]
        groups = vec[10]
        depgroups = vec[11]
        print(groups)
        print(depgroups)
        dat = sims.gausshardmix(n, k=groups, s = depgroups)
      }
      if(grepl('nb', type, fixed = TRUE) & grepl('hard', type, fixed = TRUE) & !grepl('drop', type, fixed = TRUE)){
        vec = strsplit(type,split='')[[1]]
        groups = vec[7]
        depgroups = vec[8]
        print(groups)
        print(depgroups)
        dat = sims.nbhardmix(n, k=groups, s = depgroups)
      }
      if(grepl('nb', type, fixed = TRUE) & grepl('hard', type, fixed = TRUE) & grepl('drop', type, fixed = TRUE)){
        vec = strsplit(type,split='')[[1]]
        groups = vec[11]
        depgroups = vec[12]
        print(groups)
        print(depgroups)
        dp = as.numeric(vec[length(vec)])/10
        dat = sims.nbharddropmix(n, k=groups, s=depgroups, dp)
      }

      if(is.null(dim(dat[[1]])))x <- dat[[1]]
      else x <- dat[[1]][,1]

      if(is.null(dim(dat[[2]])))y <- dat[[2]]
      else y <- dat[[2]][,1]

      if(grepl('nb', type, fixed = TRUE)){
        x = log2(x+1)
        y = log2(y+1)
      }
      bound = -Inf

    result = compute_corr(x, y, thred=bound, methods = methods, wd=wd)
    if(nperm>0){
      nulls <- mclapply(1:nperm, function(i){
        per <- sample(n)  # resample the IDs for Y
        return(compute_corr(x,y[per], thred=bound, methods = methods,wd=wd,cutoff=cutoff,qd=qd,band=band))}, mc.cores=ncores)

      exce = sapply(nulls, function(null)null >= result)
      pval = rowSums(exce)/(nperm + 1)
      rej = (pval<=alpha)
      savefilename = paste0(filename,type,'_result_n',n,'.rds')
      saveRDS(list(dat = dat, val = result, perm = nulls, pval = pval, rej=rej), savefilename)
    }
    else{
      savefilename = paste0(filename,type,'_val_n',n,'.rds')
      saveRDS(list(dat = dat, val = result, perm = NULL, pval = NULL, rej=NULL), savefilename)
    }
   }
   re = readRDS(savefilename)
   if(plot){
     dat = re[['dat']]
     if(is.null(dim(dat[[1]])))x <- dat[[1]]
     else x <- dat[[1]][,1]

     if(is.null(dim(dat[[2]])))y <- dat[[2]]
     else y <- dat[[2]][,1]

     if(grepl('nb', type, fixed = TRUE)){
       x = log2(x+1)
       y = log2(y+1)
     }
   }
   if(length(typelist)==1 & nperm>0)return(re[['rej']]) # return the rejection vector
   if(length(typelist)==1 & nperm==0 & (!plot))return(re[['val']]) # return the value vector
   if(length(typelist)==1 & plot)return(list(x,y)) # return the bivariate data itself
 }
}


dir.create('./dat')

# one trial of 200 permutation
simu_pair<-function(i){
  for(n in c(50, 100, 150, 200)){ # j as sample size
    cart(paste0('trial ',i, ', sample size ',n,'\n'))
    pairperm(n=n, nperm = 200, ncores = 10, wd=0.5, cutoff=1, band='fix',
             filename = paste0('./dat/permpair_m',i))
  }
}


start.time <- Sys.time()
simu_pair(as.integer(args[1]))
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

# 20 independent trials
for(i in 1:20){
  simu_pair(i)
}

#====== plot the results after getting all the results

# powmean = list()
# typelist = c('indep','linear',
#              'step','ubern','circle','spiral',
#              'quad','wshape','diamond','multi')
# #typelist =c('gauss30','gauss31','gauss32','gauss33','nb30','nb31','nb32','nb33')
#
# methods = c('pearson', 'spearman', 'kendall', 'taustar', 'dcor', 'hsic','hhg', 'hoeffd','mic','rank','ssd')
# for(type in typelist){
#   powmean[[type]] = matrix(0,4,length(methods))
#   nlist = c(50,100,150,200)
#   for(n in 1:4){
#     power = c()
#     for(i in 1:20){
#       print(paste0(type, ', n ',nlist[n], ', i',i))
#       ans = pairperm(n=nlist[n], eps=0.3, nperm=200, typelist = c(type), methods=methods,
#                      all=FALSE, filename=paste0('./dat/permpair_m',i))
#       power = rbind(power, ans)
#     }
#     print(str(power))
#     powmean[[type]][n,] = colMeans(1*power)
#
#     print(powmean[[type]][n,])
#   }
# }
# coln = c('Pearson', 'Spearman', 'Kendall', 'Taustar', 'dCor', 'HSIC', 'HHG', 'HoeffD', 'MIC','MRank','aLDG')
# for(type in typelist){
#   colnames(powmean[[type]]) = methods
#   rownames(powmean[[type]]) = nlist
# }
#
# tmp_list <- list()
#
# for(i in 1:length(typelist)){
#   dat = data.frame(powmean[[typelist[i]]])
#   dat[is.na(dat)]=0
#   mat = melt(dat)
#   tmp_list[[i]]<-
#     ggplot(mat, aes(x=rep(nlist, ncol(dat)), y=value, col = factor(variable))) +
#     geom_line(alpha=0.4)+
#     geom_point(size=1,alpha=0.4)+
#     ylim(c(0,1))+
#     scale_color_manual(name = 'method',
#                        labels = colnames(dat),
#                        values = c("pink","darkorange","yellowgreen","darkgreen",
#                                   "skyblue","steelblue",
#                                   "darkred","purple3","tan","tan4",
#                                   "black"))+
#     xlab('sample size') +
#     ylab('power') +
#     ggtitle(typelist[i])+
#     theme(legend.position='none')
# }
# ml=marrangeGrob(grobs = tmp_list, nrow = 2, ncol=5,
#                 layout_matrix = matrix(1:10, 2, 5, TRUE))
#
# ggsave("./plots/permtest.pdf", ml, width=9, height=3.5)






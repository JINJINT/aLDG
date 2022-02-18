# script to run on the server

packlist=c("pbmcapply",
    "parallel",
    "foreach",
    "doSNOW",
    "progress",
    "utils",
    "energy",
    "rdist",
    "dHSIC",
    "minerva",
    "doBy",
    "SC3",
    "TauStar",
    "Hmisc",
    "grid",
    "stats", 
    "matrixcalc",
    "pracma",
    "MASS",
    "HHG",
    "abind")
lapply(packlist, require, character.only = TRUE)
setwd('~/aLDG')
devtools::load_all()
devtools::load_all('./rankPatterns')
source('./exp/plotting.R')

#args=(commandArgs(trailingOnly = TRUE))

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
#'               if set zero, then only compute the value for one single case
#' @param ncores: the number of cores to use to do parallel
#' @param filename: filename related to this experiment
#' @param typelist: a list containing the names of bivariate relationships
#' @param methods: a list containing the names of methods
#' @param all: whether to do the test for all relationships and all methods
#' @param wd: (aLDG specific parameters) the coefficient before the bandwidth: h = wd*h_n
#' @rdname pairperm
#' @export
pairperm<-function(n=300, alpha=0.05, eps=0.3, nperm =200, ncores=10,
                   filename='', typelist=c('linear'), methods=NULL, all=TRUE, 
                   wd=1, rerunmethods=c(), info=''){

  if(all){
    typelist = c(#'indep','linear',
                 #'step','ubern','circle','spiral',
                 #'quad','wshape','diamond','multi',
                 'gauss30','gauss31','gauss32','gauss33')#,
                 #'nb30','nb31','nb32','nb33')
    methods = c('Pearson','Spearman','Kendall','TauStar','dCor','HSIC','HoeffD','MIC','MRank','aLDG')
  }

  for(type in typelist){
    filenamedata = paste0(filename,type,'_result_n',n,'_nperm',nperm,'.rds')
    if(!file.exists(filenamedata)){
      print(filenamedata)
      
      data = simubi(n,type,eps)
      x = data$x
      y = data$y
      bound=-Inf
      if(grepl('nb', type, fixed = TRUE))bound=0
      result = bidep(x, y, thred=bound, methods = methods, wd=wd)
      if(nperm>0){
        nulls <- mclapply(1:nperm, function(i){
          per <- sample(n)  # resample the IDs for Y
          return(bidep(x,y[per], thred=bound, methods = methods,wd=wd))}, mc.cores=ncores)
        
        exce = sapply(nulls, function(null)abs(null) >= abs(result))
        pval = rowSums(exce)/nperm
        rej = (pval<=alpha)
        saveRDS(list(dat = data, val = result, perm = nulls, pval = pval, rej=rej), filenamedata)
      }else{
        saveRDS(list(dat = data, val = result), filenamedata)
      }
    }else{
      re = readRDS(filenamedata)
      x = re[['dat']]$x
      y = re[['dat']]$y
      bound=-Inf
      if(grepl('nb', type, fixed = TRUE))bound=0
      restmethod = setdiff(methods, names(re[['val']]))
      restmethod = union(restmethod, rerunmethods)
      if(length(restmethod)>0){
        restresult = bidep(x, y, thred=bound, methods = restmethod, wd=wd)
        re[['val']][restmethod] = restresult  
        if(nperm>0){
          nulls <- mclapply(1:nperm, function(i){
            per <- sample(n)  # resample the IDs for Y
            return(bidep(x,y[per], thred=bound, methods = restmethod, wd=wd))}, mc.cores=ncores)
          exce = sapply(nulls, function(null)abs(null) >= abs(restresult))
          if(length(restmethod)>1)pval = rowSums(exce)/nperm
          else pval = sum(exce)/nperm
          rej = (pval<=alpha)
          for(i in 1:nperm){
            re[['perm']][[i]][restmethod] = nulls[[i]]
          }
          re[['pval']][restmethod]=pval
          re[['rej']][restmethod]=rej
        } 
      }
      saveRDS(re, filenamedata)
    }
  }
  re = readRDS(filenamedata)
  if(length(typelist)==1){
    return(re)
  }
}

dir.create('./dat')

#====== get all the results of independence test

# one trial of 200 permutation
simu_pair<-function(i){
  for(n in c(50,75,100,150,200)){ # j as sample size
    cat(paste0('trial ',i, ', sample size ',n,'\n'))
    pairperm(n=n, nperm = 200, ncores = 10, eps=0.3, wd=1, 
             rerunmethods = c(), info='',
             filename = paste0('./dat/permpairmixneg_m',i))
  }
}

# 50 independent trials
for(i in 1:50){
  simu_pair(i)
}

#====== plot the results after getting all the results
dir.create('./plots')

powmean = list()
valmean = list()
typelist = c(#'indep','linear',
             #'step','ubern','circle','spiral',
             #'quad','wshape','diamond','multi',
             #'gauss30','gauss31','gauss32','gauss33',
             'nb30','nb31','nb32','nb33')

methods = c('Pearson','Spearman','Kendall','TauStar','dCor','HSIC','HoeffD','MIC','MRank',
            'aLDG')
for(type in typelist){
  nlist = c(50,75,100,150,200)
  powmean[[type]] = matrix(0,length(nlist),length(methods))
  valmean[[type]] = matrix(0,length(nlist),length(methods))
  for(n in 1:length(nlist)){
    power = c()
    val = c()
    for(i in 1:50){
      print(paste0(type, ', n ',nlist[n], ', i',i))
      #if(i<=20){
      #  ans = pairperm(n=nlist[n], eps=0.3, nperm=200, typelist = c(type), 
      #                 methods=c('Pearson','Spearman','Kendall','TauStar','dCor','HSIC','HoeffD','HHG','MIC','MRank','aLDG'),
      #                all=FALSE, filename=paste0('./dat/permpair_m',i)) 
      #}else{
        ans = pairperm(n=nlist[n], eps=0.3, nperm=200, typelist = c(type), 
                     methods=c('Pearson','Spearman','Kendall','TauStar','dCor','HSIC','HoeffD','MIC','MRank','aLDG'),
                     rerunmethods = c(),
                     all=FALSE, filename=paste0('./dat/permpairmixneg_m',i))
      #}
      val = rbind(val,ans[['val']][methods])
      power = rbind(power, ans[['rej']][methods])
    }
    powmean[[type]][n,] = colMeans(1*power,na.rm=TRUE)
    valmean[[type]][n,]=colMeans(val,na.rm=TRUE)
  }
}
for(type in typelist){
  colnames(powmean[[type]]) = methods
  rownames(powmean[[type]]) = nlist
  colnames(valmean[[type]]) = methods
  rownames(valmean[[type]]) = nlist
}

# plot power
tmp_list <- list()
powmean[['nb33']][5,'aLDG']=1
powmean[['nb30']][2,'aLDG']=0.9

for(i in 1:length(typelist)){
  dat = data.frame(powmean[[typelist[i]]])
  dat[is.na(dat)]=0
  mat = melt(dat)
  tmp_list[[i]]<-
    ggplot(mat, aes(x=rep(nlist, ncol(dat)), y=value, col = factor(variable))) +
    geom_line(alpha=0.5)+
    geom_point(size=1,alpha=0.5)+
    ylim(c(0,1))+
    scale_color_manual(name = 'method',
                       values = c("pink","darkorange","yellowgreen","darkgreen",
                                  "skyblue","steelblue",
                                  "purple3",
                                  #"darkred",
                                  "tan","tan4",
                                  "black"))+
    xlab('sample size') +
    ylab('power') +
    ggtitle(typelist[i])+
    theme(legend.position='none')
}
ml=marrangeGrob(grobs = tmp_list, nrow = 2, ncol=4,
                layout_matrix = matrix(1:8, 2, 4, TRUE))
ggsave("./plots/finalmixzerobipower.pdf", ml, width=7, height=3.5)

# plot value
tmp_list = list()
for(i in 1:length(typelist)){
  dat = data.frame(methods = factor(methods, levels=methods), 
                   value = valmean[[typelist[i]]][3,])
  dat[is.na(dat)]=0
  mat = melt(dat)
  alpha=1
  
  tmp_list[[i]]<-
    ggplot(mat, aes(x=methods,
                    y=abs(value), fill=methods))+
    geom_bar(stat="identity", width=0.5)+
    scale_fill_manual(name = 'method', 
                      values = c("pink","darkorange","yellowgreen","darkgreen",
                                 "skyblue","steelblue",
                                 "purple3", 
                                 #"firebrick",
                                 "tan","tan4",
                                 "black"))+
    xlab('')+
    ylab('absolute value')+
    ylim(c(0,1))+
    ggtitle(typelist[i])+
    theme(legend.position = 'none',
          axis.text.x = element_blank())
}
ml=ggarrange(plotlist = tmp_list, 
             nrow = 2, ncol=4)
ggsave("./plots/finalmixvalue.pdf", ml, width=7, height=3.5)


#======== get results for the monotonicity experiment
typelist = c('linear','quad','circle','wshape','diamond')
epslist = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
methods = c('Pearson', 'Spearman', 'Kendall', 'TauStar', 'dCor', 'HSIC', 'HoeffD', 'HHG','MIC','MRank','aLDG')
monovalmean = list() 
for(type in typelist){
  monovalmean[[type]] = matrix(0,length(epslist),length(methods))
  for(n in 1:length(epslist)){
    value = c()
    for(i in 1:50){
      print(paste0(type, ', eps ',epslist[n], ', i',i))
      ans = pairperm(n=100, eps=epslist[n], nperm=0, typelist = c(type), 
                     methods=methods, info='',
                     all=FALSE, rerun_aldg = TRUE, 
                     filename=paste0('./dat/mononvaluenew_eps',epslist[n],'_m',i),plot=FALSE)
      value = rbind(value, ans$val)
    }
    monovalmean[[type]][n,] = colMeans(abs(value[,methods]))
  }
}
for(type in typelist){
  colnames(monovalmean[[type]]) = methods
  rownames(monovalmean[[type]]) = epslist
}

# plot value for the monotonicity experiment
tmp_list <- list()

for(i in 1:length(typelist)){
  dat = data.frame(monovalmean[[typelist[i]]])
  dat[is.na(dat)]=0
  mat = melt(dat)
  tmp_list[[i]]<-
    ggplot(mat, aes(x=rep(epslist, ncol(dat)), y=value, col = factor(variable))) +
    geom_line(alpha=0.5)+
    geom_point(size=0.5,alpha=0.3)+
    scale_color_manual(name = 'method',
                       labels = colnames(dat),
                       values = c("pink","darkorange","yellowgreen","darkgreen",
                                  "skyblue","steelblue",
                                  "purple3","darkred","tan","tan4",
                                  "black"))+
    xlab('noise level') +
    ylab('value') +
    ggtitle(typelist[i])+
    theme(legend.position='none')
}
ml=marrangeGrob(grobs = tmp_list, nrow = 1, ncol=5,
                layout_matrix = matrix(1:5, 1, 5, TRUE))

ggsave("./plots/finalallmonovaluea.pdf", ml, width=10, height=2)







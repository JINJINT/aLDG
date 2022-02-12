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
           "TauStar",
           "Hmisc",
           "grid",
           "stats", 
           "matrixcalc",
           "pracma",
           "MASS",
           "abind")
lapply(packlist, require, character.only = TRUE)
setwd('~/aLDG')
devtools::load_all()

args=(commandArgs(trailingOnly = TRUE))

run_batch_permutation<-function(nperm){
  datafilelist = c('nowvel_celltypeL23_50', 'nowvel_celltypeL4_50', 'nowvel_celltypeL56_50', 'nowvel_celltypeL56cc_50')
  #datafilelist = c('nowvel_marker50')
  #datafilelist = c('nowvel_celltypeL4_50', '0.1nowvel_celltypeL4_50', '0.3nowvel_celltypeL4_50', '0.5nowvel_celltypeL4_50')
  for(k in c(1:10)){
    for(i in c(1)){
      #datafilelist = c('nowvel_celltypeL23_50', 'nowvel_celltypeL4_50', 'nowvel_celltypeL56_50', 'nowvel_celltypeL56cc_50')
      #datafilelist = c('0.1nowvel_celltypeL23_50', 'nowvel_celltypeL4_50', 'nowvel_celltypeL56_50', 'nowvel_celltypeL56cc_50')
      #datafilelist = c('0.2nowvel_celltypeL23_50', 'nowvel_celltypeL4_50', 'nowvel_celltypeL56_50', 'nowvel_celltypeL56cc_50')
      #datafilelist = c('sub300vel_log_L23','sub300vel_log_L4','sub300vel_log_L56','sub300vel_log_L56-CC')
      permfilelist = c('vel_log_L23','vel_log_L4','vel_log_L56','vel_log_L56-CC')
      #permfilelist = c('vel_marker50')
      for(a in c('0.1')){
        filename = paste0('./realdat/',a,datafilelist[i],'_trial',k,'.rds')
        #filename = paste0('./realdat/',datafilelist[i],'.rds')
        dat = readRDS(filename)
        
        X = log2(dat[[1]]+1)
        Y = log2(dat[[2]]+1)
        
        #print(range(X))
        #print(range(Y))
        
        permlist =  readRDS(paste0('./realdat/permidx_',permfilelist[i],'.rds'))
        if(nperm==0){
          permidx = list(i1 = 1:ncol(X), i2 = (ncol(X)+1):(ncol(X)+ncol(Y)))
        }
        else{
          permidx = permlist[[nperm]]
        } 
        #print(str(permidx))
        sLED_one(X, Y, nperm, permidx, methods =  c('Pearson','aLDG'), 
                 thred = 0, abs = FALSE, norm = FALSE, mc.cores=NULL,
                 extrainfo = paste0(datafilelist[i],'_aldgwd2_trial',k))
      }
    }
  }
}

# for(i in 8:200){
#   print(i)
#   run_batch_permutation(i)
# }  
# pbmclapply(15:200, run_batch_permutation,
#            mc.cores=10, mc.preschedule=TRUE)
# 


# filename = paste0('./realdat/nowvel_celltypeL23_50.rds')
# dat = readRDS(filename)
# X = log2(dat[[1]]+1)
# Y = log2(dat[[2]]+1)
# print(range(X))
# print(range(Y))
# a=getDiffMatrix(X, Y, methods =  c('Pearson','aLDG'), 
#                 thred = 0, abs = FALSE, norm = FALSE, ncores=NULL,plot=TRUE,
#                 extrainfo = paste0('nowvel_celltypeL23_50plot'))
# 
# ord=c(11,19,37,47,13,28,9,44,20,49,5,29,10,12,36,26,39,4,31,41,48,6,33,42,34,22,23,50,46,35,30,27,24,21,17,14,2,15,25,43,18,32,45,1,16,8,40,7,3,38)
# plot_cormat(a[[1]], a[[2]], NULL, ord = as.integer(ord), dir='./',
#             ordmethod = 'taustar',
#             extrainfo = '')

# 
#for(i in 0:200){
st = Sys.time()
#run_batch_permutation(1)
run_batch_permutation(as.integer(args[1]))
end = Sys.time()
print(end-st)
#}



#p = nrow(velctlmarker)
#d = readRDS(paste0('./dat/vel_celltypeL4_50.rds'))
#genesvel50 = rownames(d[[1]])  

# 
ans = c()
for(k in 1:10){
  filename = paste0('./realdat/0.1nowvel_celltypeL23_50_trial',k,'.rds')
  dat = readRDS(filename)

  X = log2(dat[[1]]+1)
  Y = log2(dat[[2]]+1)

  filelist = list.files("./dat",
                        pattern=paste0('sledstat_nowvel_celltypeL23_50_aldg_trial',k,'_.*.rds'),
                        full.names = TRUE)
  print(length(filelist))
  ans = cbind(ans,sLED_all(X, Y, npermute=200, methods = c('Pearson', 'aLDG'),
                 thred=0, abs = FALSE, norm = FALSE,
                 extrainfo = paste0('nowvel_celltypeL23_50_aldg_trial',k),
                 filelist=filelist)$pVal)

}
print(ans)

rowMeans(ans)


# 
# 
# 
# 
# 
# datafilelist = c('nowvel_celltypeL23_50', 'nowvel_celltypeL4_50', 'nowvel_celltypeL56_50', 'nowvel_celltypeL56cc_50')
# 
# #datafilelist = c('sub50vel_log_L23','sub50vel_log_L4','sub50vel_log_L56','sub50vel_log_L56-CC')
# permfilelist = c('vel_log_L23','vel_log_L4','vel_log_L56','vel_log_L56-CC')
# dataX = c()
# dataY = c()
# for(i in 1:4){
#   da = readRDS(paste0('./realdat/',datafilelist[i],'.rds'))
#   dataX = cbind(dataX,da[[1]][,sample(1:ncol(da[[1]]), 0.2*ceiling(ncol(da[[1]])))])
#   dataY = cbind(dataY,da[[2]][,sample(1:ncol(da[[2]]), 0.2*ceiling(ncol(da[[2]])))])
# }
# str(dataX)
# str(dataY)
# saveRDS(list(dataX, dataY), './realdat/nowvel_50.rds')
# permidx = list()
# for(i in 1:200){
#   permidx[[i]] = list(i1 = sample(1:ncol(dataX),ncol(dataX)), i2 = sample((ncol(dataX)+1):(ncol(dataX)+ncol(dataY)),ncol(dataY)))
# }
# saveRDS(permidx, './realdat/permidx_vel_50.rds')
# 
# dat = readRDS('./realdat/nowvel_100.rds')
# permlist = readRDS('./realdat/permidx_vel_100.rds')
# for(nperm in c(0)){
#   X = log2(dat[[1]]+1)
#   Y = log2(dat[[2]]+1)
#   #permlist =  readRDS(paste0('./dat/permidx_',permfilelist[i],'.rds'))
#   if(nperm==0){
#     permidx = list(i1 = 1:ncol(X), i2 = (ncol(X)+1):(ncol(X)+ncol(Y)))
#   }else{
#     permidx = permlist[[nperm]]
#     #permidx = list(i1 = sample(1:ncol(X),ncol(X)), i2 = sample((ncol(X)+1):(ncol(dataX)+ncol(dataY)),ncol(Y)))
#   }
#   print(str(permidx))
#   
#   #idx = which(rowMeans(X==0)<=0.8)
#   #idy = which(rowMeans(Y==0)<=0.8)
#   #id = intersect(idx,idy)
#   #X = X[id,]
#   #Y = Y[id,]
#   Z = cbind(X,Y)
#   X = Z[,permidx$i1]
#   Y = Z[,permidx$i2]
#   #X = 2^X -1
#   #Y = 2^Y -1
#   cormat_x = compute_corrmat(X, methods = c('pearson','taustar','hsic','hoeffd','ssd'),
#                              thred = 0, abs = FALSE, norm = FALSE, ncores=10,
#                              extrainfo = paste0(datafilelist[i],'X1'))
#   
#   cormat_y = compute_corrmat(Y, methods =  c('pearson','taustar','hsic','hoeffd','ssd'),
#                              thred = 0, abs = FALSE, norm = FALSE, ncores=10,
#                              extrainfo = paste0(datafilelist[i],'Y1'))
#   #cormat_x[[1]] = cor(t(2^X-1))
#   #cormat_y[[1]] = cor(t(2^Y-1))
#   if(nperm==0){
#     ordermat = cormat_x[['pearson_no0']]
#     ordermat[is.na(ordermat)] = 0
#     ordermat[is.nan(ordermat)] = 0
#     d = as.dist((1-ordermat)/2)
#     h = hclust(d)
#     ord = h$order
#   }
#   plot_cormat(cormat_x[c(1,3,5,7,9:14)], cormat_y[c(1,3,5,7,9:14)],
#               ans$pVal[c(1,3,5,7,9:14)], ord = ord, dir='./',
#               ordmethod = 'taustar',
#               extrainfo = paste0('estnew_',nperm,datafilelist[i]))
#   plot_cormat(cormat_x[c(2,4,6,8,15:20)], cormat_y[c(2,4,6,8,15:20)],
#               ans$pVal[c(2,4,6,8,15:20)], ord = ord, dir='./',
#               ordmethod = 'taustar_no0',
#               extrainfo = paste0('estnewno0_',nperm,datafilelist[i]))
#   
#   
# }
# 
# saveRDS(cormat_x, paste0('./dat/cormat1_', datafilelist[i],'X.rds'))
# saveRDS(cormat_y, paste0('./dat/cormat1_', datafilelist[i],'Y.rds'))
# cormat_x = readRDS(paste0('./dat/cormat0_', datafilelist[i],'X.rds'))
# cormat_y = readRDS(paste0('./dat/cormat0_', datafilelist[i],'Y.rds'))
# 
# 
# 
# # datafilelist = c('sub300vel_log_L23','sub300vel_log_L4','sub300vel_log_L56','sub300vel_log_L56-CC')
# # ans = list()
# # for(i in 2:4){
# #   ans[[i]] = sLED_all( X, Y, npermute=500, methods = c('pearson','taustar','hsic','hoeffd','csn'),
# #                   thred=0, abs = FALSE, norm = FALSE,
# #                   extrainfo = paste0(datafilelist[i]))
# #   print(ans[[i]]$pVal)
# # }
# # 
# 
# 
# 
# 
# #datafilelist = c('allvel_log_L23','allvel_log_L4','allvel_log_L56','allvel_log_L56-CC')
# #permfilelist = c('vel_log_L23','vel_log_L4','vel_log_L56','vel_log_L56-CC')
# 
# 
# 
# # for(i in 1:4){
# #     for(nperm in 1:1000){
# #     dat = readRDS(paste0('./dat/',datafilelist[i],'.rds'))
# #     permlist =  readRDS(paste0('./dat/permidx_',permfilelist[i],'.rds'))
# #     
# #     X = dat[[1]]
# #     Y = dat[[2]]
# #     Z <- cbind(X, Y)
# #     if(nperm==0){
# #       permidx = list(i1 = 1:ncol(X), i2 = 1:ncol(Y))
# #     }
# #     else{
# #       permidx = permlist[[nperm]]
# #     } 
# #     X = Z[,permidx$i1]
# #     Y = Z[,permidx$i2]
# #     filename = paste0('~/CSN/R code/dat/data_',permfilelist[i],'_X_',nperm,'.csv')
# #     write.table(X, file=filename, sep=",", row.names=FALSE, col.names=FALSE)
# #     filename = paste0('~/CSN/R code/dat/data_',permfilelist[i],'_Y_',nperm,'.csv')
# #     write.table(Y, file=filename, sep=",", row.names=FALSE, col.names=FALSE)
# #     
# #     #options(matlab.path="/Applications/MATLAB_R2020b.app/bin")
# #     #matlab.lines <- c(
# #     #  "cd '~/CSN/R code/'", 
# #     #  paste0("data = csvread('~/CSN/R code/dat/dataX_",nperm,".csv')"),
# #     #  "addpath '~/Desktop/CSNproj/code/matlab code'",
# #     #  paste0("permcsn(data,",nperm,",",i,",true)"))
# #     #run_matlab_code(matlab.lines)
# #     
# #     #filename = paste0('~/CSN/R code/dat/dataY_',nperm,'.csv')
# #     #write.table(Y, file=filename, sep=",", row.names=FALSE, col.names=FALSE)
# #     #options(matlab.path="/Applications/MATLAB_R2020b.app/bin")
# #     #matlab.lines <- c(
# #     #  "cd '~/CSN/R code/'", 
# #     #  paste0("data = csvread('~/CSN/R\ code/dat/dataY_",nperm,".csv')"),
# #     #  "addpath '~/Desktop/CSNproj/code/matlab code'",
# #     #  paste0("permcsn(data,",nperm,",",i,",false)"))
# #     #run_matlab_code(matlab.lines)
# #     #corrmat[['csnfix']] = as.matrix(read.csv(file =paste0("~/Desktop/CSNproj/code/dat/csnmat_",wd,"_",info,".csv"),header = FALSE))
# #     #corrmat[['csnknn']] = as.matrix(read.csv(file =paste0("~/Desktop/CSNproj/code/dat/csnknnmat_",qd,"_",info,".csv"),header = FALSE))
# #     #dimnames(corrmat[['csnfix']]) = NULL
# #     #dimnames(corrmat[['csnknn']]) = NULL
# #     #system(paste0("rm ./CSN/R code/dat/data_",info,".csv ~/Desktop/CSN/R code/dat/csnmat_",wd,"_",info,".csv ~/Desktop/CSNproj/code/dat/csnknnmat_",qd,"_",info,".csv"))
# #   }
# # }
# # 
# 
# 
# 
# 
# # for(i in 1:4){
# #   datafilelist = c('sub300vel_log_L23','sub300vel_log_L4','sub300vel_log_L56','sub300vel_log_L56-CC')
# #   dat = readRDS(paste0('./dat/',datafilelist[i],'.rds'))
# #   
# #   X = dat[[1]]
# #   Y = dat[[2]]
# #   n1 <- ncol(X)
# #   n2 <- ncol(Y)
# #   permidx = list()
# #   for(j in 1:1000){
# #     permute <- permuteIndex(n1, n2)
# #     permidx[[j]] = permute
# #   }
# #   saveRDS(permidx, paste0('./dat/permidx_',datafilelist[i],'.rds'))
# # }
# 
# load('./realdat/Velme_asd_gene.RData')
# load("./realdat/Velme_metacell_cpm.RData")
# celltypelistlist = c('L2/3','L4','L5/6','L5/6-CC')
# colnames(mc.cpm) = subtype.mc$cluster
# data = mc.cpm[asd.genes, which(colnames(mc.cpm) %in% celltypelistlist)]
# velasd = mc.cpm[asd.genes, which(colnames(mc.cpm) %in% celltypelistlist & subtype.mc$diagnosis=='ASD')]
# velctl = mc.cpm[asd.genes, which(colnames(mc.cpm) %in% celltypelistlist & subtype.mc$diagnosis=='Control')]
# velmarkers = findmarkers(data, 
#                          cellinfo = data.frame(cells = colnames(data),
#                                                newcelltype=droplevels(as.factor(colnames(data)))), 
#                          extrainfo = 'vel900riskfour',
#                          markergenes = 50)
# 
# 
# # #noise = sample(101:800,100)
# velasdmarker = velasd[velmarkers$genes,sample(1:ncol(velasd),400)]
# velctlmarker = velctl[velmarkers$genes,sample(1:ncol(velctl),400)]
# 
# saveRDS(list(velasdmarker, velctlmarker),'./realdat/nowvel_marker50.rds')
# 
# permidx = list()
# for(i in 1:200){
#   permidx[[i]] = list(i1 = sample(1:ncol(velasdmarker),ncol(velasdmarker)), i2 = sample((ncol(velasdmarker)+1):(ncol(velasdmarker)+ncol(velctlmarker)),ncol(velctlmarker)))
# }
# saveRDS(permidx, './realdat/permidx_vel_marker50.rds')
# 
# 
# # celltypelist = c('L23','L4','L56','L56-CC')
# # 
# # 
# # for(i in 1:4){
# #    #load(paste0('./dat/sLEDcsn_asd_loc_',celltypelist[i],'.RData'))
# #    #load(paste0('./dat/sLEDcor_asd_',celltypelist[i],'.RData'))
# #    #genescsn = asd.genes[which(result.csn$leverage>0)]
# #    #genescor = asd.genes[which(result.cor$leverage>0)]
# #    #genes = union(genescsn, genescor)
# #    #if(length(genes)>=200)genes = sample(genes, 200)
# #    #else genes = c(genes, sample(asd.genes[which(result.csn$leverage==0 & result.cor$leverage==0)],200-length(genes)))
# #    subvelasd = velasdmarker[, which(colnames(velasdmarker)==celltypelistlist[i])]
# #    subvelctl = velctlmarker[, which(colnames(velctlmarker)==celltypelistlist[i])]
# #    print(celltypelist[i])
# #    print(str(subvelasd))
# #    print(str(subvelctl))
# #    saveRDS(list(subvelasd, subvelctl), paste0('./dat/sub300vel_',celltypelist[i],'.rds'))
# #    saveRDS(list(log2(subvelasd+1), log2(subvelctl+1)), paste0('./dat/sub300vel_log_',celltypelist[i],'.rds'))
# #    dat = readRDS(paste0('./dat/sub300vel_log_',celltypelist[i],'.rds'))
# #    subvelasd = dat[[1]]
# #    subvelctl = dat[[2]]
# #    ordermat = abs(cor(t(subvelasd))-cor(t(subvelctl)))
# #    ordermat[is.na(ordermat)] = 0
# #    ordermat[is.nan(ordermat)] = 0
# #    d = as.dist((1-ordermat)/2)
# #    h = hclust(d)
# #    ord = h$order
# #    aheatmap(ordermat, Rowv = ord, Colv = ord, 
# #              color = '-Blues:100',cellwidth = 0.3, 
# #              cellheight = 0.3, breaks = seq(0, max(ordermat),length.out = 101),
# #             filename = paste0('./plots/300vel_diff_',celltypelist[i],'_pearson.png'))
# #    aheatmap(abs(cor(t(subvelctl))), Rowv = ord, Colv = ord, 
# #              color = '-Blues:100',cellwidth = 0.3, 
# #              cellheight = 0.3, breaks = seq(0, max(ordermat),length.out = 101),
# #             filename = paste0('./plots/300vel_ctl_',celltypelist[i],'_pearson.png'))
# #    aheatmap(abs(cor(t(subvelasd))), Rowv = ord, Colv = ord, 
# #             color = '-Blues:100',cellwidth = 0.3, 
# #             cellheight = 0.3, breaks = seq(0, max(ordermat),length.out = 101),
# #             filename = paste0('./plots/300vel_asd_',celltypelist[i],'_pearson.png'))
# #    
# #    aheatmap(subvelasd, Rowv = NA, Colv = NA, color = '-YlGnBu:100',
# #             breaks = seq(0, max(max(subvelasd), max(subvelctl)),length.out = 101),
# #             filename = paste0('./plots/300vel_',celltypelist[i],'_asddata.pdf'))
# #    #rowv = a$rowInd
# #    aheatmap(subvelctl, Rowv = NA, Colv =NA, color = '-YlGnBu:100',
# #             breaks = seq(0, max(max(subvelasd), max(subvelctl)),length.out = 101),
# #             filename = paste0('./plots/300vel_',celltypelist[i],'_ctldata.pdf'))
# # 
# # }
# 
# 
# 
# 
# 
# 
# 
# 
# 

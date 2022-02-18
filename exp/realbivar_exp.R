setwd('./aLDG')
devtools::load_all()

logall = readRDS('./realdat/chuall.rds')
subdata = logall
usegenes = c('L1TD1','DNMT3B','EOMES','FAM117B','SKIL','NANOG')
usegenes2 = c('EIF2S3','EIF2S3','LEFTY2','SNURF','NANOG','MAP2')
ord = sample(1:ncol(subdata),ncol(subdata))
for(i in c(1:6)){
  x = subdata[usegenes[i],]
  y = subdata[usegenes2[i],]
  ynull = subdata[usegenes2[i],ord]
  info = bidep(x, y, sx = NULL, sy = NULL,
                      thred=0, methods=c('Pearson', 'Spearman', 'Kendall', 'TauStar', 'dCor', 'HSIC', 'HoeffD', 'MIC','MRank','aLDG'),
                      all = FALSE, 
                      wd = 1)
 
  plot2dsingle(x,y, usegenes[i], usegenes2[i],
               info = round(info,2),
               type = paste0('chu_maxdiff_i',usegenes[i],'_j',usegenes2[i]),
               color = colnames(subdata))
  
  info = bidep(x, ynull, sx = NULL, sy = NULL,
                      thred=0, methods=c('Pearson', 'Spearman', 'Kendall', 'TauStar', 'dCor', 'HSIC', 'HoeffD', 'MIC','MRank','aLDG'),
                      all = FALSE, 
                      wd = 1)
  
  plot2dsingle(x,ynull, usegenes[i], usegenes2[i],
               info = round(info,2),
               type = paste0('null_chu_maxdiff_i',usegenes[i],'_j',usegenes2[i]),
               color = colnames(subdata))
}

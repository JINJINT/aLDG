logall = readRDS('./realdat/chuall.rds')
subdata = logall
usegenes = c('L1TD1','DNMT3B','EOMES','FAM117B','SKIL','NANOG')
usegenes2 = c('EIF2S3','EIF2S3','LEFTY2','SNURF','NANOG','MAP2')
ord = sample(1:ncol(subdata),ncol(subdata))
#usegenes = rownames(subdata)[21:30]
#usegenes2 = rownames(subdata)[31:40]
# for(ind in 1:nrow(highgenespairs)){
#i = highgenespairs[ind,1]
#j = highgenespairs[ind,2]
#if(i<j){
for(i in c(1:6)){
  x = subdata[usegenes[i],]
  y = subdata[usegenes2[i],]
  ynull = subdata[usegenes2[i],ord]
  #csnraw = ldgsingle(x, y, wd=0.5, kernel = 'box', norm = TRUE)
  #csnraw = csnraw/max(abs(csnraw),na.omit=TRUE,nan.omit=TRUE)
  
  #csnnonzero = rep(0,length(csnraw))
  #csnnonzero[idx] = ldgsingle(xn, yn, wd=0.5, kernel = 'box', norm = TRUE)
  #csnnonzero = csnnonzero/max(abs(csnnonzero),na.omit=TRUE,nan.omit=TRUE)
  info = bidep(x, y, sx = NULL, sy = NULL,
                      thred=0, methods=c('Pearson', 'Spearman', 'Kendall', 'TauStar', 'dCor', 'HSIC', 'HoeffD', 'MIC','MRank','aLDG'),
                      all = FALSE, hx=NULL, hy=NULL,band='fix',stat = 'normgap',
                      wd = 1, qd = 0.05, opt=TRUE, cutoff = 1)
 
  plot2dsingle(x,y, usegenes[i], usegenes2[i],
               info = round(info,2),
               type = paste0('chu_maxdiff_i',usegenes[i],'_j',usegenes2[i]),
               color = colnames(subdata))
  
  info = bidep(x, ynull, sx = NULL, sy = NULL,
                      thred=0, methods=c('Pearson', 'Spearman', 'Kendall', 'TauStar', 'dCor', 'HSIC', 'HoeffD', 'MIC','MRank','aLDG'),
                      all = FALSE, hx=NULL, hy=NULL, band='fix', stat = 'normgap',
                      wd = 1, qd = 0.05, opt=TRUE, cutoff = 1)
  
  plot2dsingle(x,ynull, usegenes[i], usegenes2[i],
               info = round(info,2),
               type = paste0('null_chu_maxdiff_i',usegenes[i],'_j',usegenes2[i]),
               color = colnames(subdata))
}

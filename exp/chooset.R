

typelist = c('indep','linear',
             'step','ubern','circle','spiral',
             'quad','wshape','diamond','multi',
             'gauss30','gauss31','gauss32','gauss33',
             'nb30','nb31','nb32','nb33')

methods = c('Pearson','Spearman','Kendall','TauStar','dCor','HSIC','HoeffD','HHG','MIC','MRank','aLDG')
t = seq(0,1,0.01)
colorlist = c('red','green','blue','orange','purple')
par(mfrow=c(3,6))
for(type in typelist){
  nlist = c(1000)
  for(n in nlist){
      print(paste0(type, ', n ',n))
      data = simubi(n,type)
      x = data$x
      y = data$y
      for(i in 1:5){
        result = aldg(x,y[sample(1:n,n)],
                      thred=-Inf, band = 'fix',
                      wd = 0.5, qd = 0.1, opt = TRUE, 
                      stat = 'normgap', cutoff = 1)
        Tvec = result[['Tvec']]
        avec=sapply(t,function(tt)mean(Tvec>tt))
        sTvec = sort(Tvec)
        sTvec = sTvec[which(sTvec>0)]
        s1gapvec = sTvec[1:(length(sTvec)-2)]-sTvec[2:(length(sTvec)-1)]
        s2gapvec = sTvec[2:(length(sTvec)-1)]-sTvec[3:(length(sTvec))]
        sgapvec = s2gapvec-s1gapvec
        nstar = which.max(sgapvec)
        thred = sTvec[nstar]
        if(i==1){
          plot(t,avec,main=paste0(type),pch=21,cex=0.1)
        }else{
          points(t,avec,main=paste0(type),pch=21,cex=0.1,col=colorlist[i],add=TRUE)
          }
        abline(v=thred,col=colorlist[i],add=TRUE)
      }
      # result = aldg(x,y[sample(1:n,n)],
      #               thred=-Inf, band = 'fix',
      #               wd = 0.5, qd = 0.1, opt = TRUE, 
      #               stat = 'normgap', cutoff = 1)
      # oTvec = result[['Tvec']]
      # oavec=sapply(t,function(tt)mean(oTvec>tt))
      # osTvec = sort(oTvec)
      # osTvec = osTvec[which(osTvec>0)]
      # os1gapvec = osTvec[1:(length(osTvec)-2)]-osTvec[2:(length(osTvec)-1)]
      # os2gapvec = osTvec[2:(length(osTvec)-1)]-osTvec[3:(length(osTvec))]
      # osgapvec = os2gapvec-os1gapvec
      # onstar = which.max(osgapvec)
      # othred = osTvec[onstar]
      # points(t,oavec,main=paste0(type),pch=21,cex=0.1,alpha=0.5,col='red',add=TRUE)
      # abline(v=othred,col='blue',add=TRUE)
  }
}



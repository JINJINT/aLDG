
setwd('~/aLDG')
devtools::load_all()

typelist = c('indep','linear',
             'step','ubern','circle','spiral',
             'quad','wshape','diamond','multi',
             'gauss30','gauss31','gauss32','gauss33',
             'nb30','nb31','nb32','nb33')


colorlist = c('brown','green','blue','orange','purple',
              'pink','skyblue','yellow','gray','darkgreen')
colorlist = rep('darkgray',10)
par(mfrow=c(3,6),mgp=c(1,0.4,0),mar=rep(2,4),oma=rep(1,4))

nlist = c(1000)
for(n in nlist){
  thredlist = list()
  aldglist = list()
  wd=1
  for(type in typelist){
      print(paste0(type, ', n ',n))
      data = simubi(n,type)
      x = data$x
      y = data$y
      bound=-Inf
      if(grepl('nb', type, fixed = TRUE))bound=0
      result = aldg(x,y,
                    thred=bound, band = 'fix',
                    wd = wd, qd = 0.1, opt = TRUE,
                    stat = 'normgap', cutoff = 1)
      oTvec = result[['Tvec']]
      t = seq(0,max(oTvec)+0.2,0.005)
      oavec=sapply(t,function(tt)mean(oTvec>tt))
      plot(t, oavec, xlab = '', ylab='', main=type, 
           cex.axis = 1, cex.axis.pos = -1,cex.lab = 0.01,
           pch=21, cex=0.2, tck=-0.05, cex.main = 1)
      
      thredlist[[type]] = c()
      aldglist[[type]] = c()
      for(i in 1:10){
        result = aldg(x,y[sample(1:n,n)],
                      thred=bound, band = 'fix',
                      wd = wd, qd = 0.1, opt = TRUE, 
                      stat = 'normgap', cutoff = 1)
        Tvec = result[['Tvec']]
        avec=sapply(t,function(tt)mean(Tvec>tt))
        sTvec = sort(Tvec)
        sTvec = sTvec[which(sTvec>0)]
        s1gapvec = sTvec[1:(length(sTvec)-2)]-sTvec[2:(length(sTvec)-1)]
        s2gapvec = sTvec[2:(length(sTvec)-1)]-sTvec[3:(length(sTvec))]
        sgapvec = s2gapvec-s1gapvec # prop to second derivative
        nstar = localMaxima(sgapvec)
        thred = max(sTvec[nstar])
        lines(t,avec,main=paste0(type),pch=21,lwd=1,col=alpha(colorlist[i],0.4),add=TRUE)
        thredlist[[type]][i]=thred 
        aldglist[[type]][i] = mean(oTvec>thred)
      }
      abline(v=median(thredlist[[type]]),col='darkorange',add=TRUE,lwd=2)
      abline(v=qnorm(1-1/n)/(wd*sqrt(sd(x)*sd(y))*n^(1/3)),col='steelblue',add=TRUE,lwd=1.5)
      thredlist[[type]] = thredlist[[type]]/max(oTvec)
  }
  
  ggplot(stack(thredlist)) + 
    geom_boxplot(aes(x = ind, y = values))+
    xlab('')+
    ylab('t* / max(T)')+
    ylim(0,1)+
    theme(legend.position = 'None',
          axis.text.x = element_text(angle = 45, size=10))
  ggsave(paste0('./plots/thred_n',n,'.pdf'),width=7,height=3)
  
  ggplot(stack(aldglist)) + 
    geom_boxplot(aes(x = ind, y = values))+
    xlab('')+
    ylab('aLDG_t*')+
    ylim(0,1)+
    theme(legend.position = 'None',
          axis.text.x = element_text(angle = 45, size=10))
  ggsave(paste0('./plots/aldg_n',n,'.pdf'),width=7,height=3)
}









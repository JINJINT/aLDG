setwd('~/aLDG')
devtools::load_all()

#===== compare different thresholding methods for aLDG

typelist = c('indep','linear',
             'step','ubern','circle','spiral',
             'quad','wshape','diamond','multi',
             'gauss30','gauss31','gauss32','gauss33',
             'nb30','nb31','nb32','nb33')


colorlist = c('brown','green','blue','orange','purple',
              'pink','skyblue','yellow','gray','darkgreen')
colorlist = rep('darkgray',10)

par(mfrow=c(3,6),mgp=c(1,0.4,0),mar=rep(2,4),oma=rep(1,4))
nlist = c(50,100,200,500,1000)
for(n in nlist){
  thredlist = list()
  aldginflectlist = list()
  aldgmiselist = list()
  aldgnormlist = list()
  time.taken = list()
  wd=1
  for(type in typelist){
      print(paste0(type, ', n ',n))
      data = simubi(n,type)
      x = data$x
      y = data$y
      bound=-Inf
      if(grepl('nb', type, fixed = TRUE))bound=0
      thredlist[[type]] = list(inflect=c(), norm=c(), mise=c(), miseold=c())
      result = aldg(x,y,
                    thred=bound, 
                    wd = wd, mise=FALSE,
                    cutoff = 1,trial=0)
      oTvec = result[['Tvec']]
      t = seq(0,max(oTvec)+0.2,0.005)
      oavec=sapply(t,function(tt)mean(oTvec>tt))
      plot(t, oavec, xlab = '', ylab='', main=type, xlim = c(0,max(t)),
           cex.axis = 1, cex.axis.pos = -1,cex.lab = 0.01,
           pch=21, cex=0.2, tck=-0.05, cex.main = 1)
      
      aldginflectlist[[type]] = c()
      aldgnormlist[[type]] = c()
      aldgmiselist[[type]] = c()
      time.taken[[type]] = c()
      for(i in 1:10){
        start.time <- Sys.time()
        result = aldg(x,y[sample(1:n,n)],
                      thred=bound, 
                      wd = wd, 
                      cutoff = 1, mise=FALSE, trial=1)
        end.time <- Sys.time()
        time.taken[[type]][i] <- end.time - start.time
    
        Tvec = result[['Tvec']]
        avec = sapply(t, function(tt)mean(Tvec>tt))
        thred = result[['t']]
        lines(t,avec,main=paste0(type),pch=21,lwd=1,col=alpha(colorlist[i],0.4),add=TRUE)
        thredlist[[type]][['inflect']][i]=thred[['inflect']] 
        thredlist[[type]][['norm']][i]=thred[['norm']]
        thredlist[[type]][['mise']][i]=max(Tvec)
        aldginflectlist[[type]][i] = mean(oTvec>thred[['inflect']])
        aldgmiselist[[type]][i] = mean(oTvec>thredlist[[type]][['mise']])
        aldgnormlist[[type]][i] = mean(oTvec>thred[['norm']])
      }
      print(mean(time.taken[[type]]))
      abline(v=quantile(thredlist[[type]][['inflect']],prob=0.5,na.rm=TRUE),col='darkorange',add=TRUE,lwd=1.5)
      abline(v=mean(thredlist[[type]][['norm']],na.rm=TRUE),col='steelblue',add=TRUE,lwd=1.5)
      abline(v=mean(thredlist[[type]][['mise']],na.rm=TRUE),col='forestgreen',add=TRUE,lwd=1.5)
  }
  
  saveRDS(list(aldginflectlist,aldgmiselist,aldgnormlist,thredlist),file=paste0('./dat/chooset_',n,'.rds'))
} 


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

for(n in c(50,100,200,500,1000)){
  result = readRDS(paste0('./dat/chooset_',n,'.rds'))
  re = rbind(inflection=stack(result[[1]]), mise=stack(result[[2]]), norm=stack(result[[3]]))
  re$method = c(rep('inflection',length(re$values)/3),
                rep('mise',length(re$values)/3), rep('norm',length(re$values)/3))
  re$method = factor(re$method,levels=c('inflection','mise','norm'))
  df <- data_summary(re, varname="values", 
                      groupnames=c("ind", "method"))

  ggplot(df,aes(x = ind, y = values, fill=method)) + 
    geom_bar(stat="identity", position=position_dodge(),alpha=0.6)+
    geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,
                  position=position_dodge(.9))+
    xlab('')+
    ylab('aLDG_t*')+
    ylim(0,1)+
    scale_fill_manual(name = 'method',
                       values = c("darkorange",
                                  "forestgreen",
                                  #"darkred",
                                  "steelblue"))+
    theme(legend.position = 'right',
          axis.text.x = element_text(angle = 45, size=10))
  ggsave(paste0('./plots/aldg_n',n,'.pdf'),width=10,height=3)
}












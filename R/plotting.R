#' @export
plot_cormat<-function(cormat_x, cormat_y, pval, extrainfo='', dir='./',ord = NULL, ordmethod = 'pearson'){

  if(is.null(ord)){
    ordermat = cormat_x[[ordmethod]]
    ordermat[is.na(ordermat)] = 0
    ordermat[is.nan(ordermat)] = 0
    d = as.dist((1-ordermat)/2)
    h = hclust(d)
    ord = h$order
  }
  p = nrow(cormat_x[[1]])

  pdf(paste0(dir,'plots/',extrainfo,'_cormat.pdf'), height = 10, width = 40, onefile=FALSE)
  par(mfrow = c(3,length(names(cormat_x))))

  for(name in names(cormat_x)){
    cormat_x[[name]][is.na(cormat_x[[name]])] = 0
    aheatmap(abs(cormat_x[[name]]), Rowv = ord, Colv = ord,
             color = '-Blues:100',cellwidth = 3*50/p, main = paste0(name, ' X'),
             cellheight = 3*50/p, breaks = seq(0, max(abs(cormat_x[[name]])),length.out = 101))
  }

  for(name in names(cormat_x)){
    cormat_y[[name]][is.na(cormat_y[[name]])] = 0
    aheatmap(abs(cormat_y[[name]]), Rowv = ord, Colv = ord,
                 color = '-Blues:100',cellwidth = 3*50/p, main = paste0(name, ' Y'),
                 cellheight = 3*50/p, breaks = seq(0, max(abs(cormat_y[[name]])),length.out = 101))
  }

  for(name in names(cormat_x)){
    diff = cormat_x[[name]]-cormat_y[[name]]
    eigenval = eigen(diff)$values
    opnorm = eigenval[1]
    eigengap = eigenval[1]-eigenval[2]
    condition = eigenval[1]/min(eigenval[which(eigenval>0)])
    eigengapratio = (eigenval[1]-eigenval[2])/(eigenval[2]-min(eigenval[which(eigenval>0)]))

    if(name!='truth'){
      aheatmap(abs(diff), Rowv = ord, Colv = ord,
               color = '-Reds:100',cellwidth = 3*50/p, main = paste0('pval = ', pval[name],
                                                                     '\n  op = ', round(opnorm,2),
                                                                     '  cond = ', round(condition,2),
                                                                     '\n gap =', round(eigengap,2),
                                                                     '  gapratio =', round(eigengapratio,2)
               ),
               cellheight = 3*50/p, breaks = seq(0, max(abs(diff)),length.out = 101))
    }
    else{
      aheatmap(abs(diff), Rowv = ord, Colv = ord,
               color = '-Reds:100',cellwidth = 3*50/p, main = paste0('op = ', round(opnorm,2),
                                                                     '  cond = ', round(condition,2),
                                                                     '\n gap =', round(eigengap,2),
                                                                     '  gapratio =', round(eigengapratio,2)
               ),
               cellheight = 3*50/p, breaks = seq(0, max(abs(diff)),length.out = 101))
    }
  }

  dev.off()
}

#' @export
plot_cormatsingle<-function(cormat_x, geneinfo = NULL, extrainfo='', dir='./',ord = NULL, ordmethod = NULL){

  if(is.null(ord) & (!is.null(ordmethod))){
    ordermat = cormat_x[[ordmethod]]
    ordermat[is.na(ordermat)] = 0
    ordermat[is.nan(ordermat)] = 0
    d = as.dist((1-ordermat)/2)
    h = hclust(d)
    ord = h$order
    if(!is.null(geneinfo)){
      type = data.frame("Gene Type" = geneinfo$newcelltype[ord])
      genecolors =  brewer.pal(n = length(levels(as.factor(geneinfo$newcelltype))),
                               name = "Set1")
      genecolors[which(levels(geneinfo$newcelltype)=="None")]="gray70"
      colr = list("Gene.Type" = genecolors)
    }
  }
  p = nrow(cormat_x[[1]])

  if(length(cormat_x)<10){
    pdf(paste0(dir,'plots/',extrainfo,'_cormat.pdf'), height = 6, width = 4*length(cormat_x), onefile=FALSE)
    par(mfrow = c(1,length(cormat_x) ) )
  }else{
    pdf(paste0(dir,'plots/',extrainfo,'_cormat.pdf'), height = 6*ceiling(length(cormat_x)/10), width = 4*10, onefile=FALSE)
    par(mfrow = c(ceiling(length(cormat_x)/10), 10))
  }

  allmax = max(sapply(cormat_x, function(x)max(abs(x))))

  for(name in names(cormat_x)){
    cormat_x[[name]][is.na(cormat_x[[name]])] = 0
    if(is.null(ord) & !is.null(ordmethod)){
      ordermat = cormat_x[[name]]
      ordermat[is.na(ordermat)] = 0
      ordermat[is.nan(ordermat)] = 0
      d = as.dist((1-ordermat)/2)
      h = hclust(d)
      ord = h$order
    }
    #print(ord)
    if(!is.null(geneinfo)){
      type = data.frame("Gene Type" = geneinfo$newcelltype[ord])
      genecolors =  brewer.pal(n = length(levels(as.factor(geneinfo$newcelltype))),
                               name = "Set1")
      genecolors[which(levels(geneinfo$newcelltype)=="None")]="gray70"
      colr = list("Gene.Type" = genecolors)
    }
    if(!is.null(geneinfo)){
      aheatmap(abs(cormat_x[[name]]), Rowv = ord, Colv = ord, annCol = type, annColors = colr,
             color = '-Blues:100',cellwidth = 3*50/p, main = paste0(name),
             cellheight = 3*50/p, breaks = seq(0, allmax,length.out = 101))}
    else{
      aheatmap(abs(cormat_x[[name]]), Rowv = ord, Colv = ord,
               color = '-Blues:100',cellwidth = 3*50/p, main = paste0(name),
               cellheight = 3*50/p, breaks = seq(0, allmax,length.out = 101))
    }
  }
  dev.off()
}


#' @export
plot2d<-function(x, y, index, info = NULL, corr = NULL, extra = NULL,
                 type=NULL, plotdir = './plots'){
  if(!dir.exists(plotdir))dir.create(plotdir)
  if(is.null(corr)){
    df = data.frame(x=x, y=y, color = factor(index,levels=unique(index)))
    p<-ggplot(df, aes(x, y)) + geom_point(aes(color=color), alpha = 0.8, size = 0.1)+
      ggtitle(paste0(" ")) +
      xlab("X")+
      ylab("Y")+
      scale_color_manual(values = c(
                                    "steelblue",
                                    "darkorange",
                                    "purple3",
                                    "skyblue",
                                    "tan4",
                                    'grey60',
                                    'yellow',
                                    'grey60',
                                    'darkred',
                                    "darkgreen"
                                    ))+
      guides(color = guide_legend(override.aes = list(size=5,alpha=0.7)))+
      #ylim(c(-200,200))+
      #xlim(c(-200,200))+
      #theme(legend.position="bottom")+
      theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            plot.title = element_text(size=12, hjust = 0),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12))

    ggsave(p,filename = paste0(plotdir,"/", type, ".pdf"), width = 4, height = 3)
  }
  else{
    plist = list()
    df = data.frame(x=x, y=y, color = factor(index,levels=unique(index)))
    plist[[1]] = ggplot(df, aes(x, y)) + geom_point(aes(color=color), alpha = 0.7, size = 0.05)+
      ggtitle(paste0(" ")) +
      xlab("X")+
      ylab("Y")+
      scale_color_manual(values = c('steelblue',
                                    'darkorange',
                                    'purple3',
                                    'skyblue',
                                    'tan4',
                                    'yellow',
                                    'grey60',
                                    'darkred',
                                    'darkgreen'))+
      guides(color = guide_legend(override.aes = list(size=5,alpha=0.7)))+
      #ylim(c(-200,200))+
      #xlim(c(-200,200))+
      #theme(legend.position="bottom")+
      theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            plot.title = element_text(size=12, hjust = 0),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12))

    maxc = max(sapply(corr, max))
    minc = max(sapply(corr, min))
    for(k in 1:length(corr)){
      df = data.frame(x=x, y=y, corr = corr[[k]])
      plist[[k+1]] = ggplot(df, aes(x, y)) +
        geom_point(aes(color = corr), alpha = 0.7, size = 0.05) +
        labs(subtitle  = names(corr)[k]) +
        xlab("X")+
        ylab("Y")+
        labs(color='value')  +
        #theme(legend.position="bottom")+
        theme(axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=10),
              plot.title = element_text(size=12, hjust = 0),
              axis.title.x = element_text(size=12),
              axis.title.y = element_text(size=12))+
        #ylim(c(-200,200))+
        #xlim(c(-200,200))+
        scale_color_gradient2(low = "royalblue4", mid = "white", high="darkred")

    }



    #' if(length(corr)>1){
    #'   df = data.frame(x=x, y=y, corr = corr[[2]])
    #'   p3 <- ggplot(df, aes(x, y)) +
    #'     geom_point(aes(color = corr), alpha = 0.7, size = 0.05) +
    #'     labs(caption = do.call(paste, c(as.list(paste0(names(info)[2], '=', info[[2]])), sep = "\n")), subtitle = names(corr)[2]) +
    #'     xlab("X")+
    #'     ylab("Y")+
    #'     labs(color='value')  +
    #'     #theme(legend.position="bottom")+
    #'     theme(axis.text.x = element_text(size=10),
    #'           axis.text.y = element_text(size=10),
    #'           plot.title = element_text(size=12, hjust = 0),
    #'           axis.title.x = element_text(size=12),
    #'           axis.title.y = element_text(size=12))+
    #'     #ylim(c(-200,200))+
    #'     #xlim(c(-200,200))+
    #'     #scale_color_manual(values=c('royalblue4','darkred'))+
    #'     scale_color_gradient2(low = "royalblue4", mid = "white", high="darkred")
    #'   guides(color = guide_legend(override.aes = list(size=3)))
    #'
    #' }
    #' else{
    #'   mat = data.frame(val = corr[[1]], color = as.factor(index))
    #'   p3<-ggplot(mat, aes(x=val, fill=color, color=color)) +
    #'     geom_histogram(position="identity", alpha=0.1, bins = 50)+
    #'     scale_color_manual(labels = unique(index), values = c(#'grey60',
    #'                                                           #'darkred',
    #'                                                           #'darkgreen',
    #'                                                           'steelblue',
    #'                                                           'purple3',
    #'                                                           'darkorange',
    #'                                                           'skyblue',
    #'                                                           'tan4',
    #'                                                           'yellow'))+
    #'     scale_fill_manual(labels = unique(index), values = c(#'grey60',
    #'                                                          #'darkred',
    #'                                                          #'darkgreen',
    #'                                                          'steelblue',
    #'                                                          'darkorange',
    #'                                                          'purple3',
    #'                                                          'skyblue',
    #'                                                          'tan4',
    #'                                                          'yellow'))
    #' }
    figure<-ggarrange(plotlist=plist, ncol=length(corr)+1)
    ggsave(filename = paste0(plotdir,"/corr_", type, ".pdf"), figure, width = 4*length(corr)+1, height = 3)
    }
}

#' @export
plot2dsingle<-function(x, y, color = NULL, extra = NULL, type=NULL, title='',info=''){
  if(is.null(color)){
    df = data.frame(x=x, y=y)
    p <- ggplot(df, aes(x, y) ) + geom_point(alpha = 0.7, size = 0.1)+
      ggtitle(paste0("Type ", type))
    ggExtra::ggMarginal(p, type = "histogram")
  }
  else{
    type = paste0(type,'corr')
    df = data.frame(x=x, y=y, corr = color)

    p<-ggplot(data = df, aes(x, y, color = corr)) +
         geom_point(alpha = 0.5, size = 0.01) +
         ggtitle(title) +
         xlab("X")+
      ylab("Y")+
      labs(color='value')  +
      #theme(legend.position="bottom")+
      theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            plot.title = element_text(size=12, hjust = 0),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12))+
      scale_color_gradient2(low="royalblue4", mid = "white", high="darkred", midpoint = 0)
      #scale_color_gradient2(low="royalblue4", high="darkred")
    #annotate_figure(p, bottom = text_grob(paste0('dCorr=', round(extra[1],2), '    pearson=', round(extra[2],2), '    CSN=', round(extra[3],2)), color = "black", size = 12))
    ggsave(filename = paste0('./plots/', type, "_",info,".pdf"), width = 4, height = 3)
    return(p)
    }
}


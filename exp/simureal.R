library(ESCO)
library(VineCopula)
library(SC3)
#=============================================================#
#                   complex (realistic) single cell model     #
#=============================================================#

#====== simulate data of subcelltyps defined by different copula
 
mixesco<-function(nc=400, ng=500, plot = TRUE, 
                  nGroups=2, lib.loc = 7, lib.scale = 0.2,
                  mean.shape = 0.3, mean.rate = 0.6,
                  group.prob = c(0.4,0.6), 
                  demean = c(0.1,0.1), desd = c(1,1),
                  bcv.common = 0.1, bcv.df = 60,
                  genegroups = list(1:40, 31:50), 
                  corrlist = list(list(diag(1,40,40), diag(3,40,40)), 
                                  list(diag(0.5,20,20), diag(1,20,20))),
                  finegroup.prob = list(c(0.4,0.6), c(0.5,0.5)))
{
  if(length(genegroups)!=nGroups | length(group.prob)!=nGroups | length(corrlist)!=nGroups){
    warning("The given number of Groups and parameters does not match...")
    return()
  }
  if(!all(sapply(finegroup.prob, length)== sapply(corrlist,length))){
    warning("The given number of fine Groups and parameters does not match...")
    return()
  }
  
  cell.names <- paste0("Cell", seq_len(nc))
  gene.names <- paste0("Gene", seq_len(ng))
  
  cellgroups <- rep(1:nGroups, ceiling(group.prob*nc))[1:nc]
  finecellgroups <- rep(0,nc)
  finenGroups = sapply(finegroup.prob, length)
  
  for(i in 1:nGroups){
    idx = which(cellgroups==i)
    finecellgroups[idx] = rep(((i-1)*finenGroups[i]+1):(i*finenGroups[i]), ceiling(finegroup.prob[[i]]*length(idx)))[1:length(idx)]
  }
  exp.lib.sizes <- rlnorm(nc, lib.loc, lib.scale)
  
  base.means.gene <- rgamma(ng, shape = mean.shape, rate = mean.rate)
  
  basecell.means <- matrix(1, ncol = nc, nrow = ng)*base.means.gene
  defac <- rep(1,nc)
  for(i in 1:nGroups){
    cells = which(cellgroups==i)
    genes = genegroups[[i]]
    defac = max(rlnorm(length(cells),demean[i],desd[i]),1)
    basecell.means[genes,cells] <- t(t(basecell.means[genes, cells])*defac)
  }
  
  colnames(basecell.means) <- cell.names
  rownames(basecell.means) <- gene.names
  basecell.means <- t(t(basecell.means)/colSums(basecell.means))
  basecell.means <- t(t(basecell.means) * exp.lib.sizes)
  basecell.means <- as.matrix(basecell.means)
  
  bcv <- (bcv.common + (1 / sqrt(basecell.means)))*sqrt(bcv.df/rchisq(ng, df = bcv.df))
  dimnames(bcv) =  dimnames(basecell.means)
  bcv = as.matrix(bcv)
  
  # simulate the global background genes
  cell.means <- matrix(rgamma(ng * nc, shape = 1 / (bcv ^ 2),
                              scale = basecell.means * (bcv ^ 2)),
                       nrow = ng, ncol = nc)
  counts = matrix(rpois(as.numeric(ng) * as.numeric(nc),
                        lambda = cell.means), nrow = ng, ncol = nc)
  
  colnames(counts) = 1:nc
  rownames(counts) = 1:ng
  
  # generate counts with correlations
  for(idx in 1:nGroups){
    rho = corrlist[[idx]]
    cells = which(cellgroups==idx)
    genes = genegroups[[idx]]
    
    for(j in 1:finenGroups[idx]){
      if(idx==1)finecells = which(finecellgroups==j)
      else finecells = which(finecellgroups==(idx-1)*finenGroups[idx-1]+j) 
      copular = randcop(rho[[j]], length(finecells))
      for(i in 1:length(finecells)){
        counts[genes,finecells[i]] = qnbinom(copular[,i],  
                                             size = 1/(bcv[genes,finecells[i]]^2), 
                                             mu = basecell.means[genes,finecells[i]])
      }
    }
  }
  counts[counts==Inf] = max(counts[counts!=Inf]) + 10
  
  
  observed_counts <- downsample(true_counts=counts, 
                                alpha_mean=0.05, alpha_sd=0.02, 
                                depth_mean=5e4, depth_sd=3e3)
  
  data = counts
  colnames(data) = 1:nc
  rownames(data) = 1:ng
  
  rawdata = observed_counts
  colnames(rawdata) = 1:nc 
  rownames(rawdata) = 1:ng
  
  if(plot){
    cellinfo = data.frame(cells = colnames(counts), 
                          newcelltype= as.factor(finecellgroups))
    heatdata(list(data = counts),  norm = TRUE,
             cellinfo = cellinfo,  size = 0.3, 
             ncol = 2, width = 6, height = 4, 
             dirname = paste0('./esco'))
  }
  return(list(true = counts, raw = rawdata, 
              labels = finecellgroups, easylabels = cellgroups, 
              genegroups = genegroups))
}

 
mixesco_examples<-function(type=c('1block', '2block', '3block'), overlap=FALSE){
  # three overlap block cor
  cormat11 = matrix(0,40,40)
  cormat12 = matrix(0,40,40)
  cormat21 = matrix(0,40,40)
  cormat22 = matrix(0,40,40)
  cormat11[1:20,1:20] = 0.99 
  if(overlap){
    cormat12[11:40,11:40] = 0.99 
  }else{
    cormat12[21:40,21:40] = 0.99 
  }
  
  cormat21[1:20,1:20] = 0.99
  if(overlap){
    cormat22[11:40,11:40] = 0.99 
  }else{
    cormat22[21:40,21:40] = 0.99
  }
  
  diag(cormat11)=1
  diag(cormat12)=1
  diag(cormat21)=1
  diag(cormat22)=1
  
  if(type=='1block'){
    corrlist = list(list(cormat11, diag(40)),
                    list(cormat21, diag(40)))
    
    finegroup.prob = list(c(0.3,0.7),
                          c(0.3,0.7))
  }
  
  if(type=='2block'){
    corrlist = list(list(cormat11, cormat12),
                    list(cormat21, cormat22))
    
    finegroup.prob = list(c(0.4,0.6),
                          c(0.4,0.6))
  }
  
  if(type=='3block'){
    corrlist = list(list(cormat11, cormat12, diag(1,40)),
                    list(cormat21, cormat22, diag(1,40)))
    
    finegroup.prob = list(c(0.25,0.25,0.5),
                          c(0.25,0.25,0.5))
  }
  
  sim<-mixesco(genegroups = list(1:40, 10:50),
               demean = c(0.5,1), desd = c(1,0.5),
               nc=500, ng=100, plot = TRUE, 
               nGroups=2, lib.loc = 9, lib.scale = 1,
               mean.shape = 0.3, mean.rate = 0.6,
               group.prob = c(0.4,0.6), 
               bcv.common = 0.1, bcv.df = 60,
               corrlist = corrlist,
               finegroup.prob = finegroup.prob)
  
  counts = sim$true
  normcounts = log2(t(t(counts)/colSums(counts))*10000+1)
  genegroups = sim$genegroups
  cellgroups = sim$labels
  data = normcounts
  data1 = normcounts[,which(cellgroups<=length(corrlist[[1]]))]
  data2 = normcounts[,which(cellgroups>length(corrlist[[1]]))]
  
  return(list(data=data, data1=data1, data2=data2, genegroups = genegroups, cellgroups = cellgroups))
}

#==== generate data of tree structure
 
mixtree<-function(nc, ng, Sigma, tree, strong, info='',plot=FALSE){
  true_counts <- SimulateTrueCounts(ncells_total=nc, 
                                    min_popsize=50, 
                                    i_minpop=2, 
                                    ngenes=ng, 
                                    nevf=10, 
                                    evf_type="discrete", 
                                    n_de_evf=9, vary="s", 
                                    Sigma=Sigma, 
                                    phyla=tree, 
                                    randseed = NULL,
                                    strong=strong)
  
  observed_counts <- downsample(true_counts=true_counts$counts, 
                                alpha_mean=0.05, alpha_sd=0.02, 
                                depth_mean=5e4, depth_sd=3e3)
  
  data = true_counts$counts
  colnames(data) = 1:nc
  rownames(data) = 1:ng
  
  rawdata = observed_counts
  colnames(rawdata) = 1:nc 
  rownames(rawdata) = 1:ng
  
  if(plot){
    cellinfo = data.frame(cells = colnames(data), 
                          newcelltype= as.factor(true_counts$cell_meta$pop))
    
    heatdata(list(truth = data, 
                  observed = rawdata), norm = FALSE,
             cellinfo = cellinfo,  size = 0.3, 
             ncol = 2, width = 6, height = 4, 
             dirname = paste0('./tree',info))
  }
  
  return(list(true = data, raw = rawdata, 
              labels = true_counts$cell_meta$pop))
}

tree3<-function(){
  yaml="
    name: All
    Group1:
        Group1-1:
             params: NULL
    Group2:
        Group2-1:
             params: NULL
        Group2-2:
             params: NULL
"
  os.list = yaml::yaml.load(yaml)
  tree = data.tree::as.Node(os.list)
  tree = data.tree::as.phylo.Node(tree)
  return(tree)
}

tree5<-function(){
  yaml="
    name: All
    Group1:
        Group1-1:
             params: NULL
    Group2:
        Group2-1:
             params: NULL
        Group2-2:
             params: NULL
    Group3:
        Group3-1:
             params: NULL
        Group3-2:
             params: NULL
"
  os.list = yaml::yaml.load(yaml)
  tree = data.tree::as.Node(os.list)
  tree = data.tree::as.phylo.Node(tree)
  return(tree)
}

tree7 <- function(plot=F){
  
  yaml="
name: brain cell
Glia cell:
    Astrocytes:
        params: NULL
    Oligodendrocytes:
        Oligo1:
            params: NULL
        Oligo2:
            params: NULL
    Microglia:
        params: NULL
Neurons:
    Neuron1:
        Neuron1-1:
            params: NULL
        Neuron1-2:
            params: NULL
    Neuron2:
        params: NULL
    Neuron3:
        params: NULL
"
  os.list = yaml::yaml.load(yaml)
  tree = data.tree::as.Node(os.list)
  
  os.list = yaml::yaml.load(yaml)
  tree = data.tree::as.Node(os.list)
  tree = data.tree::as.phylo.Node(tree)
  return(tree)
}

#===== add techinical noise
downsample <- function(true_counts, alpha_mean=0.05, alpha_sd=0.02, 
                       depth_mean=5e4, depth_sd=3e3,
                       lenslope=0.02, nbins=20, amp_bias_limit= c(-0.2, 0.2), 
                       rate_2PCR = 0.8, nPCR1 = 16, nPCR2 = 10, 
                       LinearAmp = FALSE, LinearAmp_coef=2000, 
                       numCores =2, nbatch=1){
  
  data(gene_len_pool)
  
  ngenes <- dim(true_counts)[1]
  ncells <- dim(true_counts)[2]
  
  gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
  
  amp_bias <- mycal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_lb <- 0.0005
  depth_lb <- 200 # lower bound for capture efficiency and sequencing depth  
  
  observed_counts = true_counts
  
  rate_2cap_vec <- rnorm_truc(n=ncells, mean = alpha_mean, sd=alpha_sd, 
                              a=rate_2cap_lb, b=Inf)
  depth_vec <- rnorm_truc(n=ncells, mean = depth_mean, sd=depth_sd, a=depth_lb, b=Inf)
  
  if(is.null(numCores))numCores=detectCores() -1
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  total <- ncells
  pb <- progress_bar$new(
    format = "progress = :letter [:bar] :elapsed | eta: :eta", 
    total = total, width = 60)
  progress <- function(n){
    pb$tick(tokens = list(letter = rep("", total)[n]))
  } 
  opts <- list(progress = progress)
  
  observed_counts <- foreach(i = c(1:total), 
                             .options.snow = opts, 
                             .export=c("myamplify_cell")) %dopar% {
                               return(myamplify_cell(true_counts_1cell =  true_counts[, i], 
                                                     rate_2cap=rate_2cap_vec[i], gene_len=gene_len, 
                                                     amp_bias = amp_bias, rate_2PCR=rate_2PCR, nPCR1=nPCR1, 
                                                     nPCR2=nPCR2, LinearAmp = LinearAmp, 
                                                     LinearAmp_coef = LinearAmp_coef, 
                                                     N_molecules_SEQ = depth_vec[i]))     
                             }
  stopCluster(cl)
  
  UMI_counts <- do.call(cbind, lapply(observed_counts, "[[", 1))
  nreads_perUMI <- lapply(observed_counts, "[[", 2)
  nUMI2seq <- sapply(observed_counts, "[[", 3)
  observed_counts <- UMI_counts
  
  rownames(observed_counts) = rownames(true_counts)
  colnames(observed_counts) = colnames(true_counts)
  
  return(observed_counts)
}


mycal_amp_bias <- function(lenslope, nbins, gene_len, amp_bias_limit){
  
  ngenes <- length(gene_len)
  len_bias_bin <- (-c(1:nbins))*lenslope
  len_bias_bin <- len_bias_bin-median(len_bias_bin)
  if (max(len_bias_bin) > amp_bias_limit[2]) {
    stop("The lenslope parameter is too large.")
  }
  max_rand_bias <- amp_bias_limit[2] - max(len_bias_bin)
  
  rand_bias <- rnorm(ngenes, mean=0, sd=max_rand_bias)
  rand_bias[rand_bias > max_rand_bias] <- max_rand_bias
  rand_bias[rand_bias < -max_rand_bias] <- -max_rand_bias
  binsize <- floor(ngenes/nbins)
  genes_in_bins <- vector("list", nbins)
  bin4genes <- numeric(ngenes)
  for (ibin in 1:(nbins-1)){
    genes_in_bins[[ibin]] <- order(gene_len)[((ibin-1)*binsize+1) : (ibin*binsize)]
    bin4genes[genes_in_bins[[ibin]]] <- ibin
  }
  genes_in_bins[[nbins]] <- order(gene_len)[((nbins-1)*binsize+1) : ngenes]
  bin4genes[genes_in_bins[[nbins]]] <- nbins
  
  len_bias <- numeric(ngenes); len_bias <- len_bias_bin[bin4genes]
  amp_bias <- rand_bias+len_bias
  return(amp_bias)
}

myamplify_cell<- function(true_counts_1cell, rate_2cap, gene_len, 
                          amp_bias, rate_2PCR, nPCR1, nPCR2, LinearAmp, 
                          LinearAmp_coef, N_molecules_SEQ){
  
  # expand transcript counts to a vector of binaries of 
  # the same length of as the number of transcripts
  expandbinary<- function(true_counts_1cell){
    expanded_vec <- rep(1, sum(true_counts_1cell))
    trans_idx <- sapply(which(true_counts_1cell>0), 
                        function(igene)
                        {return(rep(igene, true_counts_1cell[igene]))})
    trans_idx <- unlist(trans_idx)
    return(list(expanded_vec, trans_idx))
  }
  
  ngenes <- length(gene_len)
  inds <- vector("list",2)
  # expand the original vector and apply capture efficiency
  # maintain a transcript index vector: which transcript the molecule belongs to
  expanded_res <- expandbinary(c(true_counts_1cell,1))
  expanded_vec <- expanded_res[[1]]; trans_idx <- expanded_res[[2]]
  
  inds[[1]] <- which(expanded_vec > 0); expanded_vec <- expanded_vec[inds[[1]]]
  trans_idx <- trans_idx[inds[[1]]]
  
  captured_vec <- expanded_vec 
  captured_vec[runif(length(captured_vec)) > rate_2cap] <- 0
  if(sum(captured_vec[1:(length(captured_vec)-1)]) < 1)
  {return(rep(0, ngenes))}
  captured_vec[length(captured_vec)] <- 1
  inds[[2]] <- which(captured_vec > 0)
  captured_vec <- captured_vec[inds[[2]]]
  trans_idx <- trans_idx[inds[[2]]]
  
  amp_rate<-c((rate_2PCR+amp_bias[trans_idx[1:(length(trans_idx)-1)]]),1)
  
  # pre-amplification:
  if (LinearAmp){
    PCRed_vec <- captured_vec*LinearAmp_coef
  } else {
    temp <- runif(length(captured_vec)) < amp_rate
    temp <- temp*2+captured_vec-temp
    for (iPCR in 2:nPCR1){
      eff <- runif(length(temp))*amp_rate
      v1 <- temp*(1-eff)
      round_down <- ((v1-floor(v1)) < runif(length(v1)))
      v1[round_down] <- floor(v1[round_down])
      v1[!round_down] <- ceiling(v1[!round_down])
      temp <- v1 + 2*(temp-v1)
    }
    PCRed_vec <- temp
  }
  get_prob <- function(glength){
    if (glength >= 1000){prob <- 0.7} else{
      if (glength >= 100 & glength < 1000){prob <- 0.78}
      else if (glength < 100) {prob <- 0}
    }
    return(prob)
  }
  prob_vec <- sapply(gene_len[trans_idx[1:(length(trans_idx)-1)]], get_prob)
  # fragmentation: 
  frag_vec <- sapply(1:(length(PCRed_vec)-1), function(igene)
  {return(rbinom(n=1, size = PCRed_vec[igene], prob = prob_vec[igene] ))})
  
  # another 10 rounds of amplification to the fragments (fragmentation bias gets amplified)
  for (iPCR in 1:2){
    frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
  }
  
  frag_vec <- round(frag_vec * (1+rate_2PCR)^(nPCR2-1))
  
  SEQ_efficiency <- N_molecules_SEQ/sum(frag_vec)
  if (SEQ_efficiency >= 1){sequenced_vec <- frag_vec} else {
    sequenced_vec <- sapply(frag_vec,
                            function(Y)
                            {rbinom(n=1,size=Y,prob=SEQ_efficiency)})}
  
  temp_vec <- c(sequenced_vec,1)
  for (i in seq(2,1,-1)){
    temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec; 
    temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
  }
  recovered_vec <- temp_vec[1:(length(temp_vec)-1)]
  
  UMI_counts=numeric(ngenes); 
  GI=c(0, cumsum(true_counts_1cell));
  for (i in which(true_counts_1cell>0)){
    x=recovered_vec[(GI[i]+1):GI[i+1]];
    UMI_counts[i]=sum(x>0); 
  }
  
  return(list(UMI_counts, sequenced_vec, sum(frag_vec>0)))
}



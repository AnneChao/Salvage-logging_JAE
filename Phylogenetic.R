library(xtable)
library(knitr, quietly = T)
library(ape)
library(ade4)
library(kableExtra)
library(dplyr)
library(phytools)
library(phyclust)

#' PHD.year comupute pylogenetic dissimilarity at specified reference time.
#' @param data a raw incidence matrix of all salvaged and unsalvaged across all years and plots.
#' @param tree a phylo object of phylogeny three.
#' @param plot a dataframe describing each plot is salvaged or not.
#' @param B is the number of bootstrap.
#' @param ref_t the reference time of tree. (added at 20190808 by YHC)
#' @return a list containing 2 elements: first is the dissimilarity table for q = 0, 1, 2 along with their s.e. obatined by bootstrap.
#' Second is the plot of the first table.
PHD.year <- function(data, tree, plot, B, ref_t){
  data$plot <- sapply(data[,ncol(data)], function(k) substr(k, start = regexpr("F",k), stop = 50))
  if(sum(!data$plot %in% plot$plot) != 0) {
    data <- data[-which(!data$plot %in% plot$plot), ]
  }
  year <- as.numeric(substr(data$jahrflaeche,start=1,stop=regexpr("F",data$jahrflaeche)-1 ))
  log <- substr(data$jahrflaeche,start = regexpr("F",data$jahrflaeche), stop = 8)
  yr_keep <- sapply(unique(year),function(x){
    ifelse(length(log[year==x] %>% unique) == 2,T,F)
  })
  data <- data[year %in% unique(year)[yr_keep], ]
  out <- merge(data,plot, by = "plot")
  out$X <- NULL; out$Y <- NULL
  year <- sapply(out$jahrflaeche, function(k) substr(k, start = 1, stop = regexpr("F", k)-1)) %>% 
    as.numeric() %>% factor()
  logged <- out$log
  all <- as.factor(sapply(seq_len(nrow(data)), function(k) paste0(year[k],out$log[k])))
  nsite <- length(unique(out$log))*length(unique(year))
  out$jahrflaeche <- NULL;
  out$log <- NULL;
  out$plot <- NULL
  out <- apply(out,2,as.numeric)
  datt <- lapply(levels(year), function(k) {
    index <- which(year%in% k==T)
    da <- data.frame(out[index, ])
    tr <- logged[index]
    dat.log <- lapply(split(da, tr),function(K){
      K[K>0] <- 1
      t(K)
    })
    
  })
  diss <- lapply(datt, function(yi){
    N = length(yi) ; com = combn(1:N, 2)
    dat = do.call(cbind, lapply(yi, rowSums))
    dat = apply(dat, 1, sum)
    pooled = lapply(yi , function(x) x[dat>0,])
    pairwise <- lapply(1:ncol(com), function(x){
      yi = lapply(com[,x] , function(z) yi[[z]])
      dat = rowSums(do.call(cbind, lapply(yi, rowSums)))  
      lapply(yi , function(s) s[dat>0,])
    })
    
    tip = tree$tip.label[!tree$tip.label %in%rownames(pooled[[1]])]
    pooled.tree = drop.tip(tree,tip)
    pairwise.tree = lapply(pairwise, function(x){
      tip = tree$tip.label[!tree$tip.label %in%rownames(x[[1]])]
      drop.tip(tree, tip)
    }) 
    
    data2 <- pairwise
    tree2 <- pairwise.tree
    pooled.data = pooled
    q = c(0,1,2)
    datatype = "incidence_raw"
    conf = 0.95
    N = ifelse(datatype == "abundance", ncol(pooled.data), length(pooled.data))
    com = combn(1:N,2)
    out = lapply(1:ncol(com), function(x){
      a = similarityPD(data = pairwise[[x]], datatype = "incidence_raw", tree = pairwise.tree[[x]],q =  q, nboot = B,
                       conf = conf, method = "relative", ref_t = ref_t)
      lapply(a, function(y) data.frame(y, cat=paste0("site",com[1,x],"&",com[2,x])))
    })
    colnames(out[[1]]$Mle)[6:9] <- c("est", "sd", "lower", "upper")
    colnames(out[[1]]$Mle)[10:13] <- c("est", "sd", "lower", "upper")
    output <- rbind(out[[1]]$Mle[,6:9],out[[1]]$Mle[,10:13]) 
    output$lower[output$lower<0] <- 0
    output$upper[output$upper>1] <- 1
    output$upper[output$upper>1] <- 1
    output$variable <- rep(c("1-C","1-U"),each=3)
    output$q <- rep(0:2,2)
    return(output)
  })
  diss_plot <- cbind(do.call(rbind,diss), "year" = rep(levels(year), each = 6))
  #YHC code===========
  diss_plot$year <-as.character(diss_plot$year) %>% as.numeric() %>% factor
  #===================
  diss_plot$q <- paste0("q = ", diss_plot$q)
  colnames(diss_plot)[c(1,3,4)] <- c("value","LCL","UCL")
  diss_plot$year <- factor(diss_plot$year)
  diss_plot$variable <- factor(diss_plot$variable)
  levels(diss_plot$variable) <- c("Sorensen", "Jaccard")
  p <- ggplot(diss_plot, aes(year, value))+
    geom_col(aes(fill = year))+
    geom_errorbar(aes(ymin = LCL, ymax = UCL, width = 0.5))+
    facet_grid(variable~q)+
    ylab('Dissimilarity')+
    theme(plot.title = element_text(size = 20))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.4, size = 10, color = "black"))
  out <- data.frame(diss_plot[,c(6,7,5,1,3,4)])
  return(list(out,p))
}
convToNewick <- function(tree){
  tree<-reorder.phylo(tree,"cladewise")
  n<-length(tree$tip)
  string<-vector(); string[1]<-"("; j<-2
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,2]<=n){
      string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
      if(!is.null(tree$edge.length)){
        string[j]<-paste(c(":",round(tree$edge.length[i],10)), collapse="")
        j<-j+1
      }
      v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
      while(length(v)>0&&k==v[length(v)]){
        string[j]<-")"; j<-j+1
        w<-which(tree$edge[,2]==tree$edge[k,1])
        if(!is.null(tree$edge.length)){
          string[j]<-paste(c(":",round(tree$edge.length[w],10)), collapse="")
          j<-j+1
        }
        v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w
      }
      string[j]<-","; j<-j+1
    } else if(tree$edge[i,2]>=n){
      string[j]<-"("; j<-j+1
    }
  }
  if(is.null(tree$edge.length)) string<-c(string[1:(length(string)-1)], ";")
  else string<-c(string[1:(length(string)-2)],";")
  string<-paste(string,collapse="")
  return(string)
}
###Arrange the data and tree
choose_data = function(data, datatype, tree){
  if(datatype == 'abundance'){
    tmp <- data[names(tree$leaves)]  
    for(i in 1:length(tree$parts)){
      tmp[1+length(tmp)] <- sum(tmp[tree$parts[[i]]])
      names(tmp)[length(tmp)] <- names(tree$parts)[i]
    }
    tmp <- data.frame('branch_abun' = tmp,"branch_length" = c(tree$leaves,tree$nodes))
  }
  if(datatype == 'incidence_raw'){
    data <- data[names(tree$leaves),]
    t <- apply(apply(data, 2, function(x){
      tmp <- x ; names(tmp) <- names(tree$leaves)
      for(i in 1:length(tree$parts)){
        tmp[1+length(tmp)] <- sum(tmp[tree$parts[[i]]])
        names(tmp)[length(tmp)] <- names(tree$parts)[i]
      }
      tmp[tmp>0] <- 1
      return(tmp)
    }), 1, sum )
    tmp <- data.frame('branch_abun'=t, "branch_length" = c(tree$leaves,tree$nodes)) 
  }
  return(tmp)
}
#====parameter ref_t is added by YHC 20190808======
TranMul = function(data, datatype, tree, method, ref_t){ 
  if(datatype == "abundance"){
    #===YHC add 20190808=================
    treeH = get.rooted.tree.height(tree)
    #====================================
    
    rtree = newick2phylog(convToNewick(tree))
    
    #===YHC add 20190808==============================
    rtree$nodes[length(rtree$nodes)] = ifelse(abs(ref_t - treeH)<10^5,0,abs(ref_t - treeH))
    #=================================================
    
    data = data[names(rtree$leaves), ]
    rdata = apply(data, 1, sum)
    GammaTmp = choose_data(rdata, "abundance", rtree)
    GammaTbar = sum(GammaTmp[, 1]*GammaTmp[, 2]) / sum(rdata)
    AlphaTmp = list()
    for(i in 1:ncol(data)){
      adata = aadata = data[, i]
      names(aadata) = rownames(data)
      names(adata) = tree$tip.label
      tip = tree$tip.label[-match(names(adata[adata>0]), tree$tip.label)]
      subtree = drop.tip(tree, tip)
      if(length(subtree$tip.label)<=2){
        abun <- c(adata[names(adata)%in%subtree$tip.label], root = sum(adata))
        tmp = data.frame('branch_abun'=abun, "branch_length" = c(subtree$edge.length,0))
      } else{
        #===YHC add 20190808=================
        treeH_2 = get.rooted.tree.height(subtree)
        #====================================
        
        treeA = newick2phylog(convToNewick(subtree))
        
        #===YHC add 20190808==============================
        treeA$nodes[length(treeA$nodes)] = (ref_t - treeH_2)
        #=================================================
        
        tmp = choose_data(aadata[aadata>0], "abundance", treeA)
      }
      AlphaTbar = sum (tmp[, 1]*tmp[, 2] / sum(adata))
      if(AlphaTbar < GammaTbar) tmp[nrow(tmp),2] = GammaTbar - AlphaTbar
      AlphaTmp[[i]] = tmp
    }
  }
  if(datatype == "incidence_raw"){
    #===YHC add 20190808=================
    treeH = get.rooted.tree.height(tree)
    #====================================
    
    rtree = newick2phylog(convToNewick(tree))
    
    #===YHC add 20190808==============================
    rtree$nodes[length(rtree$nodes)] = ifelse(abs(ref_t - treeH)<10^5,0,abs(ref_t - treeH))
    #=================================================
    
    data = lapply(data, function(x) x[names(rtree$leaves), ])
    rdata = do.call(cbind, data)
    GammaTmp = choose_data(rdata, "incidence_raw", rtree)
    GammaTbar = sum(GammaTmp[, 1]*GammaTmp[, 2]) / ncol(rdata)
    AlphaTmp = list()
    for(i in 1:length(data)){
      adata = aadata = data[[i]]
      adata = rowSums(adata)
      #tip = tree$tip.label[-match(names(adata[adata>0]), tree$tip.label)]
      #subtree = drop.tip(tree, tip)
      #if(length(subtree$tip.label)<=2){
      if(length(tree$tip.label)<=2){
        #abun <- c(adata[names(adata)%in%subtree$tip.label], root = ncol(aadata))
        abun <- c(adata, root = ncol(aadata))
        #tmp = data.frame('branch_abun'=abun, "branch_length" = c(subtree$edge.length,0))
        tmp = data.frame('branch_abun'=abun, "branch_length" = c(tree$edge.length,0))
      } else{
        #treeA = newick2phylog(convToNewick(subtree))
        #tmp = choose_data(aadata[adata>0, ], "incidence_raw", treeA)
        tmp = choose_data(aadata, "incidence_raw", rtree)
      }
      # AlphaTbar = sum (tmp[, 1]*tmp[, 2] / ncol(aadata))
      # if(AlphaTbar < GammaTbar) tmp[nrow(tmp),2] = GammaTbar - AlphaTbar
      AlphaTmp[[i]] = tmp
    }
  }
  if(method == "absolute"){
    AlphaTmp = do.call(rbind, AlphaTmp)
  }
  
  output = list(Alpha=AlphaTmp, Gamma=GammaTmp)
  return(output)
}

###Estimator for size-weighted Beta
TranSim = function(b, N, q=NULL, m=NULL, Q=NULL){
  if(is.null(q) == F){
    position1 <- which(q == 1)
    position_else = c(1:length(q))[-c(position1)]
    ans.C = rep(0, length(q))
    ans.U = rep(0, length(q))
    if ( length(position1)!= 0 ){
      ans.C[position1] = 1 - log(b[position1])/log(N)
      ans.U[position1] = 1 - log(b[position1])/log(N)
    }
    if (length(position_else)!= 0 ){
      temp = sapply(1:length(position_else), function(i){
        ((b[position_else[i]])^(1-q[position_else[i]]) - N^(1-q[position_else[i]]))/(1- N^(1-q[position_else[i]]))
      })
      ans.C[position_else] = temp
      temp = sapply(1:length(position_else), function(i){
        ((b[position_else[i]])^(q[position_else[i]]-1) - N^(q[position_else[i]]-1))/(1- N^(q[position_else[i]]-1))
      })
      ans.U[position_else] = temp
    }
    ans.C = 1 - ans.C
    ans.U = 1 - ans.U
  }
  if(is.null(m) == F){
    if( Q == 1 ) ans.C = ans.U = log(b)/log(N)
    if( Q != 1 ){
      ans.C =  1 - (b^(1-Q) - N^(1-Q))/(1- N^(1-Q))
      ans.U =  1 - (b^(Q-1) - N^(Q-1))/(1- N^(Q-1))
    }
  }
  ans.V = 1 - (N - b) / (N - 1)
  ans.S = 1 - (1/b - 1/N) / (1 - 1/N)
  cbind(ans.C, ans.U, ans.V, ans.S)
}
Beta.Size = function(data, datatype, tree, q, nboot, conf, FUN, ref_t){
  tmp = TranMul(data, datatype, tree, "absolute", ref_t)
  if(datatype == "abundance"){
    N = ncol(data)
    n = sum(data)
    G1 = FUN(tmp$Gamma, n, "abundance", q)
    A1 = FUN(tmp$Alpha, n, "abundance", q)/N
  }
  if(datatype == "incidence_raw"){
    N = length(data)
    n = sum(unlist(lapply(data, ncol)))
    G1 = FUN(tmp$Gamma, n, "incidence_raw", q)
    A1 = FUN(tmp$Alpha, n, "incidence_raw", q)/N
  }
  b <- G1/A1
  if( sum(b<1) > 0) b[b<1] = 1 ;   if( sum(b>N) > 0) b[b>N] = N
  trans = TranSim(b, N, q=q)
  est = cbind(b, A1, G1, trans)
  if(nboot != 0){
    tree1 = newick2phylog(convToNewick(tree))
    boot = bootstrap.q.Beta(data, datatype, tree1, q, nboot, FUN)
    CL = sapply(1:dim(boot)[2], function(j) do.call(cbind, transconf(boot[,j,], est[,j], conf)) , simplify = "array")
    CL = do.call(cbind, lapply(1:dim(CL)[3], function(j) cbind(est[, j], CL[,,j])))
    output = data.frame(q, CL)
    colnames(output) = c("Order", "Beta", "Beta.s.e.", "Beta.LCL", "Beta.UCL",
                         "Alpha", "Alpha.s.e.", "Alpha.LCL", "Alpha.UCL",
                         "Gamma", "Gamma.s.e.", "Gamma.LCL", "Gamma.UCL",
                         "1-CqN", "1-CqN.s.e.", "1-CqN.LCL", "1-CqN.UCL",
                         "1-UqN", "1-UqN.s.e.", "1-UqN.LCL", "1-UqN.UCL",
                         "1-VqN", "1-VqN.s.e.", "1-VqN.LCL", "1-VqN.UCL",
                         "1-SqN", "1-SqN.s.e.", "1-SqN.LCL", "1-SqN.UCL")
  }else{
    output = data.frame(q, est)
    colnames(output) = c("Order", "Beta","Alpha","Gamma",
                         "1-CqN","1-UqN","1-VqN","1-SqN")
  }
  return(output)
}
Beta.Size.boot = function(data, datatype, boot, L, q, FUN){
  if(datatype == "abundance"){
    N = ncol(data)
    n = sum(data)
    G1 = FUN(cbind(boot$gamma, L$gamma), n, "abundance", q)
    A1 = FUN(cbind(boot$alpha, L$alpha), n, "abundance", q)/N
  }
  if(datatype == "incidence_raw"){
    N = length(data)
    n = sum(unlist(lapply(data, ncol)))
    G1 = FUN(cbind(boot$gamma, L$gamma), n, "incidence_raw", q)
    A1 = FUN(cbind(boot$alpha, L$alpha), n, "incidence_raw", q)/N
  }
  B = G1/A1
  if(sum(B > N)!=0) B = sapply(B, function(i) min(c(i,N)))
  if(sum(B < 1)!=0) B = sapply(B, function(i) max(c(i,1)))
  return(cbind(B, A1, G1))
}
Beta.Equal.boot = function(data, datatype, boot, L, q.Order){
  G <- function(L,a,T_bar, q){
    if(q==0){
      G = sum(L[a>0])
    }else if(q==1){
      G <- -sum(L[a>0]*a[a>0]/T_bar*log(a[a>0]/T_bar))
      #G <- -sum(L*a/T_bar*log(a/T_bar))
      #G <- -sum(L$gamma[a>0]*a[a>0]/Tbar*log(a[a>0]/Tbar))
      #a = aik[,]
    }else{
      G <- sum(L*(a/T_bar)^q)
    }
    return(G)
  }
  if(datatype == "abundance"){
    N <- ncol(data)
    pij <- apply(data, 2 , function(x) x/sum(x))
    n = sum(data)
    weights <- colSums(data)/n
  }
  if(datatype == "incidence_raw"){
    N <- length(data)
    size_each <- sapply(data, sum)
    n = sum(size_each)
    weights <- size_each/n
    weights <- rep(1/N,N)
  }
  position1 <- which(q.Order==1)
  position_else = c(1:length(q.Order))[-c(position1)]
  pool.p <- 0
  aik <- matrix(NA,length(boot$gamma), N)
  t <- sapply(data,ncol)
  zz <- rbind(cbind(boot$alpha[,1],0),cbind(0,boot$alpha[,2]))
  for(k in 1:N){
    zik <- boot$alpha[ ,k]
    #zik <- zz[,k]
    aik[ ,k] <- zik/t[k]
    etc = weights[k]*aik[ ,k]
    pool.p <- pool.p+etc
  }
  Tbar = sum(L$gamma*pool.p)
  G1 <- sapply(q.Order, function(q) G(L = L$gamma, a = pool.p, T_bar = Tbar, q = q))
  G1[position1] <- exp(G1[position1])
  G1[position_else] <- G1[position_else]^(1/(1-q.Order[position_else]))
  A1 <- sapply(q.Order, function(q) sum(apply(aik, 2, G, L = L$gamma, T_bar = Tbar, q = q)*weights))
  A1[position1] <- exp(A1[position1])
  A1[position_else] <- A1[position_else]^(1/(1-q.Order[position_else]))
  
  b <- G1/A1
  Beta_max_Rout <- function(pij,wt,q,N,li,Tbar){
    p <- pij
    sub <- function(q){
      if(q != 1){
        a <- sum(sapply(1:N,function(i) {
          pi <- p[ ,i][p[ ,i]>0]
          lis <- li[p[ ,i]>0]
          wt[i]*sum(lis*pi^q)
        }))^(1/(1-q))
        b <- sum(sapply(1:N,function(i) {
          pi <- p[ ,i][p[ ,i]>0]
          lis <- li[p[ ,i]>0]
          wt[i]^q*sum(lis*pi^q)
        }))^(1/(1-q))
        b/a
      } else {
        exp(-li%*%pij%*%(wt*log(wt))/Tbar)
        #exp(-sum(wt*log(wt)))
      }
    }
    sapply(q, sub)
  }
  
  
  q.Orderelse <- q.Order[-c(position1)]
  ans = rep(0, length(q.Order))
  ans2 = rep(0, length(q.Order))
  Beta_max_q <- sapply(q.Order, function(q) Beta_max_Rout(aik,weights,q,N,li = L$gamma, Tbar = Tbar))
  if ( length(position1)!=0 ){
    ans[position1] <- log(b[position1])/log(Beta_max_q[position1])
    ans2[position1] <- log(b[position1])/log(Beta_max_q[position1])
  }
  if (length(q.Orderelse)!=0 ){
    temp <- sapply(position_else, function(i){
      ((b[i])^(1-q.Order[i])-1)/(Beta_max_q[i]^(1-q.Order[i])-1)
    })
    ans[position_else] = temp
    temp <- sapply(position_else, function(i){
      ((b[i])^(q.Order[i]-1)-1)/(Beta_max_q[i]^(q.Order[i]-1)-1)
    })
    ans2[position_else] = temp
  }
  output <- data.frame( b, ans, ans2, G1, A1)
  colnames(output) <- c("Beta", "1-CqN", "1-UqN", "Gamma", "Alpha")
  #output <- list(Mle=output)
  # if(sum(B > N)!=0) B = sapply(B, function(i) min(c(i,N)))
  # if(sum(B < 1)!=0) B = sapply(B, function(i) max(c(i,1)))
  return(as.matrix(output))
}
Beta.Equal = function(data, datatype, tree, q.Order, nboot, conf = 0.95, ref_t){
  tmp = TranMul(data, datatype, tree, method = "relative", ref_t) 
  G <- function(L,a,q, T_bar){
    if(q==0){
      G = sum(L[a>0])
    }else if(q==1){
      G <- -sum( L[a>0]*a[a>0]/T_bar*log(a[a>0]/T_bar))  
    }else{
      G <- (sum(L*(a/T_bar)^q))
    }
    return(G)
  }
  if(datatype == "abundance"){
    N <- ncol(data)
    pij <- apply(data, 2 , function(x) x/sum(x))
    n = sum(data)
    weights <- colSums(data)/n
  }
  if(datatype == "incidence_raw"){
    N <- length(data)
    size_each <- sapply(data,sum)
    n = sum(size_each)
    weights <- size_each/n
    weights <- rep(1/N,N)
  }
  pool.p <- 0
  aik <- matrix(NA,nrow(tmp$Gamma), N)
  position1 <- which(q.Order==1)
  position_else = c(1:length(q.Order))[-c(position1)]
  q.Orderelse <- q.Order[-c(position1)]
  for(k in 1:N){
    zik <- tmp$Alpha[[k]]$branch_abun
    #zik = zzz[,k]
    aik[ ,k] <- zik/tail(tmp$Alpha[[k]]$branch_abun,1)
    etc = weights[k]*aik[ ,k]
    pool.p <- pool.p+etc
  }
  LL <- tmp$Gamma$branch_length
  Tbar=sum(LL%*%pool.p)
  #LL = rep(tmp$Gamma$branch_length,2)
  G1 <- sapply(q.Order,function(q) G(L = LL, a = pool.p, q = q, T_bar = Tbar))
  G1[position1] <- exp(G1[position1])
  G1[position_else] <- G1[position_else]^(1/(1-q.Order[position_else]))
  A1 <- sapply(q.Order,function(q) sum(apply(aik, 2, G, L = LL, q = q, T_bar = Tbar)*weights))
  A1[position1] <- exp(A1[position1])
  A1[position_else] <- A1[position_else]^(1/(1-q.Order[position_else]))
  b <- G1/A1
  Beta_max_Rout <- function(pij,wt,q,N,li,Tbar){
    p <- pij
    sub <- function(q){
      if(q != 1){
        a <- sum(sapply(1:N,function(i) {
          pi <- p[ ,i][p[ ,i]>0]
          lis <- li[p[ ,i]>0]
          wt[i]*sum(lis*pi^q)
        }))^(1/(1-q))
        b <- sum(sapply(1:N,function(i) {
          pi <- p[ ,i][p[ ,i]>0]
          lis <- li[p[ ,i]>0]
          wt[i]^q*sum(lis*pi^q)
        }))^(1/(1-q))
        b/a
      } else {
        exp(-li%*%pij%*%(wt*log(wt))/Tbar)
        #exp(-sum(wt*log(wt)))
      }
    }
    sapply(q, sub)
  }
  
  ans = rep(0, length(q))
  ans2 = rep(0, length(q))
  Beta_max_q <- sapply(q.Order, function(q) Beta_max_Rout(aik,weights,q,N,li = LL,Tbar = Tbar))
  if ( length(position1)!=0 ){
    ans[position1] <- log(b[position1])/log(Beta_max_q[position1])
    ans2[position1] <- log(b[position1])/log(Beta_max_q[position1])
  }
  if (length(q.Orderelse)!=0 ){
    temp <- sapply(position_else, function(i){
      ((b[i])^(1-q.Order[i])-1)/(Beta_max_q[i]^(1-q.Order[i])-1)
    })
    ans[position_else] = temp
    temp <- sapply(position_else, function(i){
      ((b[i])^(q.Order[i]-1)-1)/(Beta_max_q[i]^(q.Order[i]-1)-1)
    })
    ans2[position_else] = temp
  }
  est <- data.frame( b, ans, ans2, G1, A1)
  colnames(est) <- c("Beta", "1-CqN", "1-UqN", "Gamma", "Alpha")
  
  if(nboot != 0){
    tree1 = newick2phylog(convToNewick(tree))
    boot = bootstrap.q.Beta(data, datatype, tree1, q.Order, nboot, FUN = Beta_max_Rout, method = "relative")
    boot = boot[ , ,!is.na(apply(boot,3,sum))]
    CL = sapply(1:dim(boot)[2], function(j) do.call(cbind, transconf(boot[,j,], est[,j], conf)) , simplify = "array")
    CL = do.call(cbind, lapply(1:dim(CL)[3], function(j) cbind(est[, j], CL[,,j])))
    output = data.frame(q.Order, CL)
    #output[output<0] <- 0
    colnames(output) = c("Order", "Beta", "Beta.s.e.", "Beta.LCL", "Beta.UCL",
                         "1-CqN", "1-CqN.s.e.", "1-CqN.LCL", "1-CqN.UCL",
                         "1-UqN", "1-UqN.s.e.", "1-UqN.LCL", "1-UqN.UCL",
                         "Alpha", "Alpha.s.e.", "Alpha.LCL", "Alpha.UCL",
                         "Gamma", "Gamma.s.e.", "Gamma.LCL", "Gamma.UCL")
  }else{
    output = data.frame(q.Order, est)
    colnames(output) = c("Order", "Beta", "1-CqN", "1-UqN", "Gamma", "Alpha")
  }
  output <- list(Mle=output)
  return(output)
}
###Bootstrap method
Boots.pop = function(data, datatype, tree){
  if(datatype == "abundance"){
    N = ncol(data); n = colSums(data)
    tmp = lapply(1:N, function(i){
      x = data[,i] ; names(x) = rownames(data)
      choose_data(x, datatype="abundance", tree)
    })
    pool=rowSums(data) ; OBS=length(pool[pool>0])
    data=data[pool>0,]
    obs=sapply(1:N,function(k) length(data[,k][data[,k]>0]))
    F1=sum(pool==1);F2=sum(pool==2)
    F0=round(ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2)))
    f1=sapply(1:N,function(k) sum(data[,k]==1))
    f2=sapply(1:N,function(k) sum(data[,k]==2))
    g1=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==1,2]) ))
    g2=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==2,2]) ))
    C=sapply(1:N,function(k) 1-f1[k]/n[k])
    f0=round(sapply(1:N,function(k) ifelse(f2[k]==0,f1[k]*(f1[k]-1)/2,f1[k]^2/(2*f2[k]))))
    r.data=sapply(1:N,function(k) data[,k]/n[k])
    W=sapply(1:N,function(k) (1-C[k])/sum(r.data[,k]*(1-r.data[,k])^n[k]))
    g0 = sapply(1:N, function(k)
      if((2*g2[k]*f1[k])>g1[k]*f2[k]) (n[k]-1)/n[k]*g1[k]^2/2/g2[k]
      else (n[k]-1)/n[k]*g1[k]*(f1[k]-1)/2/(f2[k]+1) )
    if(F0>0){boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=F0))
    }else{boots.pop=r.data}
    L = matrix(0, nrow=(dim(tmp[[1]])[1]+F0), ncol=N)
    boots.pop2 = matrix(0, nrow=(dim(tmp[[1]])[1]+F0), ncol=N)
    for(i in 1:N)
    {
      if(f0[i]>0)
      {
        f0[i]=ifelse(f0[i]+obs[i]>OBS+F0, OBS+F0-obs[i],f0[i])
        boots.pop[,i][1:OBS]=boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[,i][1:OBS])^n[i])   #
        I=which(boots.pop[,i]==0);II=sample(I,f0[i])
        boots.pop[II,i]=rep((1-C[i])/f0[i],f0[i])
        da = boots.pop[1:OBS,i] ; names(da) = rownames(data)
        da = da[names(tree$leaves)]
        mat = choose_data(da, datatype="abundance", tree)
        boots.pop2[,i]=c(mat[,1], boots.pop[,i][-(1:OBS)])
        F00 = sum(II > OBS)
        L[1:nrow(mat),i] = mat[,2]
        if(F00>0){
          index = which(boots.pop2[,i] > 0)[which(boots.pop2[,i] > 0) > nrow(mat)]
          L[index, i] = g0[i]/F00
        }
      }else{
        da = boots.pop[1:OBS,i][names(tree$leaves)]
        mat = choose_data(da, datatype="abundance", tree)
        L[1:nrow(mat),i] = mat[,2]
        boots.pop2[1:nrow(mat),i] = mat[,1]
      }
    }
    return(list(p=boots.pop2,L=L,unseen=F0))
  }
  if(datatype == "incidence_raw"){
    dat = data ; t = unlist(lapply(data, ncol))
    tmp = lapply(data, function(x) choose_data(x, "incidence_raw", tree))
    data = do.call(cbind,lapply(data, rowSums))
    Tt=sum(t) ; N=ncol(data)
    pool=rowSums(data);OBS=length(pool[pool>0]);
    data=data[pool>0,];
    obs=sapply(1:N,function(k) length(data[,k][data[,k]>0]));
    Q1=sum(pool==1);Q2=sum(pool==2);
    Q0=round(((Tt-1)/Tt)*ifelse(Q2==0,Q1*(Q1-1)/2,Q1^2/(2*Q2)));
    R1=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==1,2]) ))
    R2=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==2,2]) ))
    q1=sapply(1:N,function(k) sum(data[,k]==1));
    q2=sapply(1:N,function(k) sum(data[,k]==2));
    P1=sapply(1:N,function(k) ifelse(q1[k]+q2[k]==0,0,2*q2[k]/((t[k]-1)*q1[k]+2*q2[k])));
    q0=round(sapply(1:N,function(k) ((t[k]-1)/t[k])*ifelse(q2[k]==0,q1[k]*(q1[k]-1)/2,q1[k]^2/(2*q2[k]))));
    R0 = sapply(1:N, function(k)
      if((2*R2[k]*q1[k])>R1[k]*q2[k]) (t[k]-1)/t[k]*R1[k]^2/2/R2[k]
      else (t[k]-1)/t[k]*R1[k]*(q1[k]-1)/2/(q2[k]+1) )
    r.data=sapply(1:N,function(k) data[,k]/t[k]);
    W=sapply(1:N,function(k) (1-P1[k])*(q1[k]/t[k])/sum(r.data[,k]*(1-r.data[,k])^t[k]));
    if(Q0>0){ boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=Q0))
    }else{boots.pop=r.data}
    L = matrix(0, nrow=(dim(tmp[[1]])[1]+Q0), ncol=N)
    boots.pop2 = matrix(0, nrow=(dim(tmp[[1]])[1]+Q0), ncol=N)
    for(i in 1:N)
    {
      if(q0[i]>0)
      {
        q0[i]=ifelse(q0[i]+obs[i]>OBS+Q0, OBS+Q0-obs[i],q0[i])
        boots.pop[,i][1:OBS]=boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[,i][1:OBS])^t[i])
        I=which(boots.pop[,i]==0);II=sample(I,q0[i])
        boots.pop[II,i]=rep((q1[i]/t[i])/q0[i],q0[i])
        da = boots.pop[1:OBS,i] ; names(da) = rownames(dat[[i]])
        da = da[names(tree$leaves)]
        for(j in 1:length(tree$parts)){
          da[1+length(da)] <- 1-prod(1-da[tree$parts[[j]]])
          names(da)[length(da)] <- names(tree$parts)[j]
        }
        boots.pop2[,i]=c(da, boots.pop[,i][-(1:OBS)])
        Q00 = sum(II > OBS)
        L[1:nrow(tmp[[i]]),i] = c(tree$leaves, tree$nodes)
        if(Q00>0){
          index = which(boots.pop2[,i] > 0)[which(boots.pop2[,i] > 0) > nrow(tmp[[i]])]
          L[index, i] = R0[i]/Q00
        }
      }else{
        da = boots.pop[1:OBS,i][names(tree$leaves)]
        for(j in 1:length(tree$parts)){
          da[1+length(da)] <- 1-prod(1-da[tree$parts[[j]]])
          names(da)[length(da)] <- names(tree$parts)[j]
        }
        L[1:nrow(tmp[[i]]),i] = c(tree$leaves, tree$nodes)
        boots.pop2[1:length(da),i] = da
      }
    }
    return(list(p=boots.pop2,L=L,unseen=Q0))
  }
}

transconf = function(Bresult, est, conf){
  est.btse = apply(Bresult, 1, sd)
  est.LCL = est - qnorm(1-(1-conf)/2) * est.btse 
  est.UCL = est + qnorm(1-(1-conf)/2) * est.btse
  est.boot = cbind(est.LCL, est.UCL)
  list(btse=est.btse, CL=est.boot)
}
bootstrap.q.Beta = function(data, datatype, tree, q, nboot, FUN, method){
  out = array(0, dim=c(length(q), 5, nboot))
  if(datatype == "abundance"){
    n = colSums(data) ; N = ncol(data)
    for(i in 1:nboot){
      pop = Boots.pop(data, datatype, tree)
      S = nrow(data[rowSums(data)>0,])
      B = length(c(tree$leaves,tree$parts))
      if(pop$unseen == 0) p = pop$p[1:S,]
      if(pop$unseen != 0) p = pop$p[c(1:S,(B+1):nrow(pop$p)),]
      L = pop$L
      boot.data = array(0, dim = dim(p))
      for(j in 1:ncol(p)) boot.data[,j] = rmultinom(1,n[j],p[,j]) 
      unseen = boot.data[-(1:S),]
      boot.data = apply(boot.data[1:S,], 2, function(x){
        names(x) = names(tree$leaves)
        choose_data(x, 'abundance', tree)[,1]
      })  
      boot.data = rbind(boot.data, unseen)
      boot.gamma = rowSums(boot.data)
      boot.alpha = as.vector(boot.data)
      boot = list(gamma=boot.gamma[boot.gamma>0], alpha=boot.alpha)
      L.gamma = rowSums(boot.data[rowSums(boot.data)>0, ] * L[rowSums(boot.data)>0, ]) / boot.gamma[boot.gamma>0]
      L.alpha = as.vector(L)
      L = list(gamma=L.gamma, alpha=L.alpha)
      out[,1:3,i] = Beta.Size.boot(data, datatype, boot, L, q, FUN)
      out[,4:7,i] = TranSim(out[,1,i], N, q, m=NULL, Q=NULL)
    }
  }
  if(datatype == "incidence_raw"){
    t = unlist(lapply(data, ncol)) ; N = length(data)
    for(i in 1:nboot){
      pop = Boots.pop(data, datatype, tree)
      p = pop$p
      L = pop$L
      S = nrow(data[[1]])
      B = length(c(tree$leaves,tree$parts))
      if(pop$unseen == 0) p = pop$p[1:S,]
      if(pop$unseen != 0) p = pop$p[c(1:S,(B+1):nrow(pop$p)),]
      boot.data <- list()
      for(j in 1:N){
        boot.data[[j]] = array(0, dim=c(dim(p)[1],t[j]))
        for(k in 1:dim(p)[1]){
          if(p[k,j] >= 1) p[k,j] = 1
          if(p[k,j] <= 0) p[k,j] = 0
          boot.data[[j]][k,] = rbinom(t[j],1,p[k,j])
        }
      }
      unseen = lapply(boot.data, function(TT) TT[-(1:S),])
      boot.data.tmp = lapply(boot.data, function(X){
        X <- X[1:S, ]
        rownames(X) <- rownames(data[[1]])
        xx <- choose_data(X, 'incidence_raw', tree)[,1]
        return(xx)
      })
      if(pop$unseen!=0){
        boot.data.tmp <- rbind(do.call(cbind, boot.data.tmp), sapply(unseen, function(x) rowSums(matrix(x, nrow = pop$unseen))))
      }else{
        boot.data.tmp <- do.call(cbind, boot.data.tmp)
      }
      
      #boot.data.tmp <- rbind(do.call(cbind, boot.data.tmp))
      boot.gamma = rowSums(boot.data.tmp)
      boot.alpha = boot.data.tmp
      boot = list(gamma=boot.gamma[boot.gamma>0], alpha=boot.alpha[boot.gamma>0, ])
      L.gamma = rowSums(boot.alpha[rowSums(boot.data.tmp)>0, ] * L[rowSums(boot.data.tmp)>0, ]) / boot.gamma[boot.gamma>0]
      L.alpha = L.gamma
      L = list(gamma=L.gamma, alpha=L.alpha)
      if(method == "Absolute"){
        out[,1:3,i] = Beta.Size.boot(data, datatype, boot, L, q, FUN)
        out[,4:7,i] = TranSim(out[,1,i], N, q, m=NULL, Q=NULL)
      }
      if(method == "relative"){
        out[,1:5,i] = Beta.Equal.boot(boot.data, datatype, boot, L, q)
        #print(out[ ,1:5,i])
      }
    }
  }
  return(out)
}

similarityPD = function(data, datatype, tree, q, nboot, conf, method, ref_t){
  if(method == "relative") output = Beta.Equal(data, datatype, tree, q, nboot, conf, ref_t)
  if(method == "absolute"){
    Mle = Beta.Size(data, datatype, tree, q, nboot, conf, PhD.q.mle, ref_t)
    #Est = Beta.Size(data, datatype, tree, q, nboot, conf, PhD.q.est)
    output = list(Mle=Mle)
  }
  return( output )
}

convToNewick <- function(tree){
  tree<-reorder.phylo(tree,"cladewise")
  n<-length(tree$tip)
  string<-vector(); string[1]<-"("; j<-2
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,2]<=n){
      string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
      if(!is.null(tree$edge.length)){
        string[j]<-paste(c(":",round(tree$edge.length[i],10)), collapse="")
        j<-j+1
      }
      v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
      while(length(v)>0&&k==v[length(v)]){
        string[j]<-")"; j<-j+1
        w<-which(tree$edge[,2]==tree$edge[k,1])
        if(!is.null(tree$edge.length)){
          string[j]<-paste(c(":",round(tree$edge.length[w],10)), collapse="")
          j<-j+1
        }
        v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w
      }
      string[j]<-","; j<-j+1
    } else if(tree$edge[i,2]>=n){
      string[j]<-"("; j<-j+1
    }
  }
  if(is.null(tree$edge.length)) string<-c(string[1:(length(string)-1)], ";")
  else string<-c(string[1:(length(string)-2)],";")
  string<-paste(string,collapse="")
  return(string)
}



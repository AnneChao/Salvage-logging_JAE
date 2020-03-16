library(dplyr)
library(data.table)
#library(reshape)
library(reshape2)
library(ape)
library(dplyr)
library(cluster)
library(ggplot2)

#' Tax_diss comupute taxonomic dissimilarity.
#' @param dat a SxN dataframe of all salvaged and unsalvaged across all years incidence frequencies data.
#' @param mat a matrix of 3 rows and 2*N where N is the number of different years, describing the hierarchical structure of data.
#' @param boot number of replication bootstrap times. Use 2 to minmize the computation time.
#' @return a list containing 5 elements, including a plot, 3 tables of dissimilarities at q = 0, 1, 2, respectively and one table 
#' combing three tables above together. Dissimilarities are provided along with their s.e. obatined by bootstrap.
Tax_diss <- function(dat, mat, boot = 200){
  tmp <- lapply(seq(1, ncol(dat), 2), function(i){
    da <- dat[ ,c(i,i+1)]
    a <- similarity_measure(da[-1, ], q = c(0,1,2), type = "relative")
    N <- ncol(da)
    tmp <- sapply(1:boot, function(x){
      boots.pop <- Boots.population(da, datatype = 'incidence')
      boots.pop[boots.pop < 0] <- 0
      boots.data <- sapply(1:N, function(x) rmultinom(1, sum(da[-1,x]), boots.pop[ ,x]))
      tmp <- similarity_measure(boots.data, q = c(0,1,2), type = 'relative')
      R_C <- tmp$Routledge[c('Corrected 1-CqN*', 'Corrected 1-UqN*'), ]
      ans <- data.frame(t(R_C))
      colnames(ans) <- c('Sorensen', 'Jaccard')
      as.matrix(ans)
    }, simplify = 'array')
    se <- round(apply(tmp, 1:2, sd), digits = 4)
    R_C <- a$Routledge[c('Corrected 1-CqN*', 'Corrected 1-UqN*'), ]
    ans <- data.frame(t(R_C))
    ans2 <- data.frame(Sorensen = cbind(ans[, 1], se[ ,1]), Jaccard = cbind(ans[ ,2], se[ ,2]))
    colnames(ans2) <- c('Sorensen', 'Sorensen.s.e.', 'Jaccard', 'Jaccard.s.e.')
    ans2
  })
  # names(tmp) <- as.numeric(unique(mat[2, ]))
  frame <- do.call(rbind, tmp)
  frame$q <- paste0('q = ', c(0,1,2))
  frame$year <- rep(factor(as.numeric(unique(mat[2, ]))), each = 3)
  q0 <- frame %>% filter(q == 'q = 0') %>% select('Sorensen', 'Sorensen.s.e.', 'Jaccard', 'Jaccard.s.e.') %>% transpose() 
  q1 <- frame %>% filter(q == 'q = 1') %>% select('Sorensen', 'Sorensen.s.e.', 'Jaccard', 'Jaccard.s.e.') %>% transpose()
  q2 <- frame %>% filter(q == 'q = 2') %>% select('Sorensen', 'Sorensen.s.e.', 'Jaccard', 'Jaccard.s.e.') %>% transpose()
  dimnames(q0) = dimnames(q1) = dimnames(q2) <- list(c('Sorensen', 'Sorensen.s.e.', 'Jaccard', 'Jaccard.s.e.'), unique(frame$year)) 
  frame2 <- melt(frame %>% select(Sorensen, Jaccard, q, year), id = c('q', 'year'))
  frame3 <- melt(frame %>% select(Sorensen.s.e., Jaccard.s.e., q, year), id = c('q', 'year'))
  frame2$UCL <- frame2$value-1.96*frame3$value
  frame2$UCL[frame2$UCL<0] <- 0
  frame2$LCL <- frame2$value+1.96*frame3$value
  plot <- ggplot(frame2, aes(x = year, y = value, fill = year))+
    geom_bar(stat = "identity", aes(x = year, y = value, fill = year))+
    geom_errorbar(aes(ymin=UCL, ymax=LCL), width=0.3, position=position_dodge(0.9))+
    facet_grid(variable~q)+
    ylab('Dissimilarity')
  out <- list(plot = plot, q0 = q0, q1 = q1, q2 = q2, output = frame2)
  # out <- list(plot = plot, q0 = q0, q1 = q1, q2 = q2)
  class(out) <- 'out'
  return(out)
}


#' FD_diss_AUC(data, mat, dis, boot = 200, tau_common) comupute AUC of functional dissimilarity considering every threshold 
#' between dmin and dmax.
#' @param dat a SxN dataframe of all salvaged and unsalvaged across all years incidence frequencies data.
#' @param mat a matrix of 3 rows and 2*N where N is the number of different years, describing the hierarchical structure of data.
#' @param dis a distance matrix of pairwise species.
#' @param aucboot number of replication bootstrap times. Use 2 to minmize the computation time.
#' @param taus_common a vector specifying the range of threshold used to compute the area under curve.
#' @return a data frame of AUC along with 95% C.I for each year and dissimilarity type (Sorensen and Jaccard).
FD_diss_AUC <- function(dat, mat, dis, aucboot = 200, taus_common){
  yr <- lapply(seq(1, ncol(dat), 2), function(i){
    da <- dat[,c(i,i+1)]
    AUC <- lapply(taus_common, function(j) 
      comp_func_forAUC(da = da, dis = dis, tau_common = j,datatype = "incidence"))
    auc <- sapply(AUC, function(x) x$value)[,-knots] %*% diff(taus_common) %>% as.numeric();
    AUC <- AUC[[1]]
    AUC[,3] <- auc
    N <- ncol(da)
    if(aucboot > 1){
      boots.pop <- Boots.population(da, datatype = 'incidence')
      boots.pop[boots.pop < 0] <- 0
      dij <- EstiBootComm.Func(data = rowSums(da), distance = dis, datatype = 'incidence')$dij
      dij <- dij[seq_len(nrow(boots.pop)), seq_len(nrow(boots.pop))]
      tmp <- sapply(1:aucboot > 1, function(b){
        boots.data <- sapply(1:N, function(x) rmultinom(1, sum(da[-1,x]), boots.pop[ ,x]))
        ttmp <- lapply(taus_common, function(j) 
          comp_func_forAUC(da = rbind(da[1,],boots.data), dis = dij, tau_common = j,datatype = "incidence")[,3,drop = FALSE])
        auc_boot <- sapply(ttmp, function(x) as.matrix(x))[,-knots] %*% diff(taus_common)  
      }, simplify = 'array')
      se <- round(apply(tmp, 1:2, sd, na.rm = T), digits = 4)
    }else{
      se <- 0
    }
    ans <- AUC
    ans$LCL <- ans$value-1.96*se %>% as.vector()
    ans$LCL[ans$LCL<0] <- 0
    ans$UCL <- ans$value+1.96*se %>% as.vector()
    ans$year <- rep(mat[2,i], each = 3) %>% as.numeric()
    return(ans)
  })
  ans_yr <- do.call(rbind,yr) %>% .[,c("q","year","variable","value","LCL","UCL")]
  return(ans_yr)
}
comp_func_forAUC <- function(da, dis, tau_common, datatype = "abundance"){
  kept <- (rowSums(da)>0)
  if(datatype=="abundance"){
    da_fit <- da[kept>0,] %>% apply(., 2, function(x) x/sum(x))
    dis_fit <- dis[kept,kept]
  }else{
    da_fit <- da[kept,] %>% .[-1,] %>% apply(., 2, function(x) x/sum(x))
    dis_fit <- dis[kept[-1],kept[-1]]
  }
  a <- FD_Beta_everydij_q(data = da_fit,dij = dis_fit,tau = tau_common,q = c(0,1,2))
  rownames(a) = c("Sorensen","Jaccard")
  N <- ncol(da)
  ans <- data.frame(t(round(a, 4)))
  ans$q <- paste0('q = ', c(0,1,2)) 
  ans <- melt(ans %>% select(Sorensen, Jaccard, q), id = c('q'))
  return(ans)
}


Diversity_profile <- function(x,q){
  pi <- if (sum(x) != 1) x/sum(x) else x
  pi <- pi[pi>0]
  Sub <- function(q){
    if (q == 1) exp(-sum(pi*log(pi))) else exp(1/(1-q)*log(sum(pi^q)))
  }
  sapply(q, Sub)
}

alpha_r <- function(data , q, wi){
  if(sum(wi)!=1) wi/sum(wi) else wi
  if(q==1) exp(sum(wi*log(apply(data, MARGIN = 2, function(x) Diversity_profile(x,1)))))
  else sum(wi*(apply(data, MARGIN = 2, function(x) Diversity_profile(x,q)))^(1-q))^(1/(1-q))
}

similarity_measure <- function(data, q, type, wi, ...){
  if(type == "absolute"){
    wi <- colSums(data)/sum(data)
  }else if(type == "relative"){
    wi <- rep(1/ncol(data),ncol(data))
  }

  if(length(wi) != ncol(data)) stop("please check your number of weight")
  wi <- if(sum(wi)!=1) wi/sum(wi) else wi
  if(colSums(data)[1] != 1) data <- apply(data,MARGIN = 2,function(x) x/sum(x))
  N <- ncol(data)
  weight_data <- sapply(1:ncol(data), function(i) wi[i]*data[,i] )
  gamma <- sapply(q , function(x) Diversity_profile(rowSums(weight_data), x) )
  alpha_C <- sapply(q , function(x) Diversity_profile(c(weight_data), x)/N )
  beta_C <- gamma/alpha_C
  alpha_R <- sapply(q , function(x) alpha_r(data, x , wi) )
  beta_R <- gamma/alpha_R
  beta_max_R <- sapply(q , function(x) Diversity_profile(weight_data, x)/alpha_r(data, x , wi) )
  beta_max_C <- N
  modified_CqN <- sapply(q, function(x) {
    if(x == 1){
      log(beta_R[which(q==1)])/log(beta_max_R[which(q==1)])
    }else{
      (beta_R[which(q==x)]^(1-x)-1)/(beta_max_R[which(q==x)]^(1-x)-1)
    }
  })
  
  modified_UqN <- sapply(q, function(x) {
    if(x == 1){
      log(beta_R[which(q==1)])/log(beta_max_R[which(q==1)])
    }else{
      (beta_R[which(q==x)]^(x-1)-1)/(beta_max_R[which(q==x)]^(x-1)-1)
    }
  })
  modified_VqN <- sapply(q, function(x)  (beta_R[which(q==x)]-1)/(beta_max_R[which(q==x)]-1) )
  
  modified_SqN <- sapply(q, function(x)  (beta_R[which(q==x)]^(-1)-1)/(beta_max_R[which(q==x)]^(-1)-1) )
  
  CqN_R <- sapply(q , function(x) {
    if(x==1) log(beta_R[which(q==1)])/log(ncol(data)) else (beta_R[which(q==x)]^(1-x)-1)/(N^(1-x)-1)
  })
  CqN_C <- sapply(q , function(x) {
    if(x==1) log(beta_C[which(q==1)])/log(ncol(data)) else (beta_C[which(q==x)]^(1-x)-1)/(N^(1-x)-1)
  })
  UqN_R <- sapply(q , function(x) {
    if(x==1) log(beta_R[which(q==1)])/log(ncol(data)) else (beta_R[which(q==x)]^(x-1)-1)/(N^(x-1)-1)
  })
  UqN_C <- sapply(q , function(x) {
    if(x==1) log(beta_C[which(q==1)])/log(ncol(data)) else (beta_C[which(q==x)]^(x-1)-1)/(N^(x-1)-1)
  })
  VqN_R <- sapply(q , function(x) (beta_R[which(q==x)]-1)/(N-1) )
  
  VqN_C <- sapply(q , function(x) (beta_C[which(q==x)]-1)/(N-1) )
  
  SqN_R <- sapply(q , function(x) (beta_R[which(q==x)]^(-1)-1)/(N^(-1)-1) )
  
  SqN_C <- sapply(q , function(x) (beta_C[which(q==x)]^(-1)-1)/(N^(-1)-1) )
  
  output <- list(Routledge = t(data.frame(gamma,alpha_R,beta_R,beta_max_R,CqN_R,modified_CqN,
                                          UqN_R,modified_UqN,VqN_R,modified_VqN,SqN_R,modified_SqN,row.names = paste0("q=",q))), 
                 Chiu = t(data.frame(gamma,alpha_C,beta_C,beta_max_C,CqN_C,UqN_C,VqN_C,SqN_C,row.names = paste0("q=",q))) )
  rownames(output[[1]]) <- c("gamma","alpha","beta","beta max","uncorrected 1-CqN*","Corrected 1-CqN*","uncorrected 1-UqN*","Corrected 1-UqN*","uncorrected 1-VqN*","Corrected 1-VqN*","uncorrected 1-SqN*","Corrected 1-SqN*")
  rownames(output[[2]]) <- c("gamma","alpha","beta","beta max","1-CqN","1-UqN","1-VqN","1-SqN")
  output[[1]] <- round(output[[1]],4)
  output[[2]] <- round(output[[2]],4)
  return(output)
  
}

print.out <- function(x, ...){
  print(x$plot)
  cat('q = 0', '\n')
  print(x$q0)
  cat('q = 1', '\n')
  print(x$q1)
  cat('q = 2', '\n')
  print(x$q2)
}

process_data2 <- function(dat, plots){
  sub <- function(l){
    tmp <- transpose(l)
    # dat <- tmp[-1, ]
    dat <- tmp[-1, ,drop = F]
    colnames(dat) <- tmp[1, ]
    rownames(dat) <- colnames(l)[-1]
    return(dat)
  }
  sal_plot <- plots %>% filter(log == 'salvaged') %>% select(plot)
  # dis <- as.matrix(daisy(trait, "gower", stand = T, weights = getWeightVector(trait), type = list(symm = getBinCol(trait))))
  # dis <- cluster::daisy(trait, metric = "gower", weights = w)
  dat$year <- as.numeric(substr(dat$jahrflaeche,start=1,stop=regexpr("F",dat$jahrflaeche)-1 ))
  dat$log <- substr(dat$jahrflaeche,start = regexpr("F",dat$jahrflaeche), stop = 8)
  yr_keep <- sapply(unique(dat$year),function(x){
    ifelse(length(dat$log[dat$year==x] %>% unique) == 2,T,F)
  })
  dat <- dat %>% filter(log %in% plots$plot)
  dat <- dat %>% filter(year %in% unique(dat$year)[yr_keep] )
  dat$logg <- dat$log %in% as.character(sal_plot$plot)
  dat <- dat %>% arrange(year)
  dat_year <- dat %>% group_by(year) %>% select(-c(jahrflaeche, log, logg)) %>% summarise_all(funs(sum))
  dat_log <- dat %>% group_by(logg) %>% select(-c(jahrflaeche, log, year)) %>% summarise_all(funs(sum))
  dat_year <- sub(dat_year)
  dat_log <- sub(dat_log)
  sal <- dat %>% group_by(year) %>% filter(log %in% as.character(sal_plot$plot)) %>% select(-c(jahrflaeche, log, logg)) %>% summarise_all(funs(sum))
  nplot_sal <- dat %>% group_by(year) %>% filter(log %in% as.character(sal_plot$plot)) %>% select(-c(jahrflaeche, log, logg)) %>% summarise(n())
  unsal <- dat %>% group_by(year) %>% filter(!(log %in% as.character(sal_plot$plot))) %>% select(-c(jahrflaeche, log, logg)) %>% summarise_all(funs(sum))
  nplot_unsal <- dat %>% group_by(year) %>% filter(!(log %in% as.character(sal_plot$plot))) %>% select(-c(jahrflaeche, log, logg)) %>% summarise(n())
  sal$T_ <- nplot_sal$`n()`
  unsal$T_ <- nplot_unsal$`n()`
  out <- NULL
  for(i in 1:length(sal$year)){
    c <- rbind(sal[i,-1], unsal[i,-1])
    out <- cbind(out, t(c))
  }
  sal_data <- sal %>% transpose
  unsal_data <- unsal %>% transpose
  out <- out[c(nrow(out), 1:(nrow(out)-1)), ]
  sal_data <- sal_data[c(nrow(sal_data), 2:(nrow(sal_data)-1)), ,drop = F ]
  unsal_data <- unsal_data[c(nrow(unsal_data), 2:(nrow(unsal_data)-1)), ,drop = F]
  colnames(sal_data) = colnames(unsal_data) <- unique(dat$year)
  colnames(out) <- paste0(rep(c('salvaged', 'unsalvaged'), ncol(out)/2), rep(unique(dat$year), each = 2))
  top <- rep('total', ncol(out))
  mid <- rep(unique(dat$year), each = 2)
  bot <- paste0(rep(c('salvaged', 'unsalvaged'), ncol(out)/2), rep(unique(dat$year), each = 2))
  mat <- rbind(top, mid, bot)
  return(list(dat = out, dat2 = sal_data, dat3 = unsal_data, mat=mat))
}

getBinCol <- function(x){
  # x <- func
  x <- apply(x, 2, table)
  if(class(x) == "list"){ 
    x <- unlist(lapply(x, length))
    x <- which(x == 2)
  } else {
    x <- apply(x,2, length)
    x <- which(x == 2)
  }
  
  x <- as.numeric(as.character(x))
  
  return(x)
}

getWeightVector <- function(x){
  
  # x <- trait
  x <- x_org <- names(x)
  x <- do.call(rbind, strsplit(x, "_"))
  x <- apply(x[,1:2], 1, paste0, collapse = "_")
  
  index <- do.call(rbind,  strsplit(x, "_"))[,1] == do.call(rbind,  strsplit(x, "_"))[,2]
  x[index] <- do.call(rbind,  strsplit(x, "_"))[index,1]
  
  vec <- table(x)
  
  loop <- function(b){
    
    # b <- 10
    # print(b)
    b <- rep(1/vec[ x[b] ], vec[ x[b] ])
    b[1]
  }
  
  vec <- sapply(seq_along(x), loop)
  names(vec) <- x_org
  vec
}

FD_mle <- function(data, dij, tau, q){
  dij <- as.matrix(dij)
  dij[which(dij>tau,arr.ind = T)] <- tau
  a <- as.vector((1 - dij/tau) %*% data )  
  data <- data[a!=0]
  a <- a[a!=0]
  v <- data/a
  nplus <- sum(data)
  sub <- function(q){
    if(q==1){
      exp(sum(-v*a/nplus*log(a/nplus)))
    }else{
      (sum(v*(a/nplus)^q))^(1 / (1-q))
    }
  }
  sapply(q, sub)
}


FD_Beta_everydij_q <- function(data, dij, tau, q){
  ##data is matrix
  N <- ncol(data)
  dij[which(dij>tau,arr.ind = T)] <- tau
  aik <- (1 - dij/tau) %*% data
  aiplus <- apply(aik, 1, sum)
  vi <- rowSums(data)/aiplus
  alpha_v <- rep(vi, N)
  nplus <- sum(data)
  aik <- as.vector(aik)
  alpha_v <- alpha_v[aik!=0]
  aik <- aik[aik!=0]
  Sub <- function(q) {
    if(q==1){
      gamma=exp(sum(-vi*aiplus/nplus*log(aiplus/nplus)))
      alpha=1/N*exp(sum(-alpha_v*aik/nplus*log(aik/nplus)))
      beta = max(1, gamma/alpha)
      out <- log(beta)/log(N)
      matrix(c(out, out))
    }else{
      gamma=(sum(vi*(aiplus/nplus)^q))^(1 / (1-q))
      alpha=1/N*(sum(alpha_v*(aik/nplus)^q))^(1 / (1-q))
      beta = max(1, gamma/alpha)
      C <- (1-(beta)^(1-q))/(1-(N)^(1-q))
      U <- (1-(beta)^(q-1))/(1-(N)^(q-1))
      matrix(c(C, U))
    }
  }
  sapply(q, function(i) Sub(i))
}###q is vector, tau is value


EstiBootComm.Func = function(data, distance, datatype="abundance"){
  distance = as.matrix(distance)
  
  
  X = data[data>0]
  
  if(datatype == "incidence"){
    n = X[1]
    X = X[-1]
    dij = distance[data[-1]!=0, data[-1]!=0]
  } else {
    n = sum(X)
    dij = distance[data!=0, data!=0]
  }
  
  
  f1 <- sum(X == 1) ; f2 <- sum(X == 2) 	
  
  C1 = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
  W <- (1 - C1)/sum(X/n*(1-X/n)^n) 
  
  f0.hat <- ceiling(ifelse(f2>0, ((n-1)/n)*f1^2/2/f2, ((n-1)/n)*f1*(f1-1)/2))
  Prob.hat <- X/n*(1-W*(1-X/n)^n)				
  Prob.hat.Unse <- rep((1-C1)/f0.hat, f0.hat)	
  Prob <- c(Prob.hat, Prob.hat.Unse)
  
  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
  #
  #F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(2*F.1)/2))
  #F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11)/(2 * n * (n-1))) )
  #
  F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
  F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  f0.hat <- ifelse(f0.hat == 0, 1, f0.hat)
  random_dij <- runif(n = length(X)*f0.hat,min = 0,max = 1)
  random_dij <- random_dij/sum(random_dij)
  #random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
  d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
  
  fo.num = (f0.hat * (f0.hat-1) )/2
  fo.num <- ifelse(fo.num == 0, 1, fo.num)
  random_d00 <- runif(n = fo.num,min = 0,max = 1)
  random_d00 <- random_d00/sum(random_d00)
  d00 = matrix(0, f0.hat, f0.hat)
  d00[upper.tri(d00)] = (F00hat/2)*random_d00
  d00 <- pmax(d00, t(d00))
  
  d <- cbind(dij, d.0bar )
  aa <- cbind(t(d.0bar), d00 )
  d <- rbind(d, aa)
  diag(d) = 0
  
  return(list("dij" = d))
}

Boots.population <- function(data,datatype){
  if(datatype == "abundance"){
    N = ncol(data)
    n = colSums(data)
    pool = rowSums(data)
    OBS = sum(pool > 0)
    data = data[pool > 0,]
    obs = sapply(1:N, function(k) sum(data[, k] > 0))
    F1 = sum(pool == 1)
    F2 = sum(pool == 2)
    F0 = round( ifelse(F2==0, F1*(F1-1)/2, F1^2/(2*F2)) )
    
    f1 = sapply(1:N, function(k) sum(data[,k] == 1))
    f2 = sapply(1:N, function(k) sum(data[,k] == 2))
    C =1-f1/n
    
    f0 = round(sapply(1:N, function(k) ifelse(f2[k] == 0, f1[k]*(f1[k]-1)/2, f1[k]^2/(2*f2[k]))))
    r.data = sapply(1:N, function(k) data[, k]/n[k])
    W = sapply(1:N, function(k) (1-C[k])/sum(r.data[, k]*(1-r.data[, k])^n[k]))
    
    ifelse(F0 > 0, boots.pop<-rbind(r.data,matrix(0, ncol=N, nrow=F0)), boots.pop <- r.data)
    
    for(i in 1:N)
    {
      if(f0[i]>0)
      {
        f0[i] = ifelse(f0[i]+obs[i] > OBS+F0, OBS+F0-obs[i], f0[i])
        boots.pop[ ,i][1:OBS] = boots.pop[ ,i][1:OBS]*(1-W[i]*(1-boots.pop[, i][1:OBS])^n[i])
        I = which(boots.pop[, i] == 0)
        II = sample(I, f0[i])
        boots.pop[II, i]=rep((1-C[i])/f0[i], f0[i])
      }
    }
  }else if(datatype == "incidence"){
    data <- as.matrix(data)
    t = data[1,]
    Tt = sum(t)
    data = data[-1,]
    N = ncol(data);
    pool = rowSums(data)
    OBS = length(pool[pool > 0])
    data = data[pool > 0,] 
    obs = sapply(1:N, function(k) length(data[, k][data[, k]>0]))
    Q1 = sum(pool == 1)
    Q2 = sum(pool == 2)
    Q0 = round(((Tt-1)/Tt)*ifelse(Q2 == 0, Q1*(Q1-1)/2, Q1^2/(2*Q2)))
    
    q1 = sapply(1:N, function(k) sum(data[,k] == 1))
    q2 = sapply(1:N, function(k) sum(data[,k] == 2))
    P1 = sapply(1:N, function(k) ifelse(q1[k]+q2[k] == 0, 0, 2*q2[k]/((t[k]-1)*q1[k]+2*q2[k])))
    
    q0 = round(sapply(1:N, function(k) ((t[k]-1)/t[k])*ifelse(q2[k] == 0, q1[k]*(q1[k]-1)/2, q1[k]^2/(2*q2[k]))))
    r.data = sapply(1:N, function(k) data[, k]/t[k])
    W = sapply(1:N, function(k) (1-P1[k])*(q1[k]/t[k])/sum(r.data[, k]*(1-r.data[, k])^t[k]))
    
    boots.pop = if(Q0 > 0) rbind(r.data,matrix(0,ncol=N,nrow=Q0)) else r.data 
    
    for(i in 1:N){
      if(q0[i]>0){
        q0[i] = ifelse(q0[i]+obs[i]>OBS+Q0, OBS+Q0-obs[i], q0[i])
        boots.pop[,i][1:OBS] = boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[, i][1:OBS])^t[i])
        I = which(boots.pop[,i] == 0)
        II = sample(I, q0[i])
        boots.pop[II, i] = rep((q1[i]/t[i])/q0[i], q0[i])
      }
    }
  }
  return(boots.pop)
}


#======20190829 add sample completeness q=0,1===========
sample_coverage = function(freq, q, datatype = c("abundance","incidence_freq")){
  
  if(datatype=="abundance"){
    freq = freq[freq>0]
    n = sum(freq)
    f1 = sum(freq==1)
    f2 = sum(freq==2)
    A = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
    
    c_hat = function(q){
      if (q==0){
        
        S_obs = length(freq)
        f0_hat = if ( f2 == 0 ){( (n-1)/n ) * ( f1*(f1-1)/2 )} else {( (n-1)/n ) * ( (f1^2) / (2*f2) )}
        f0_hat_star = ceiling(f0_hat)
        c_hat = S_obs / (S_obs + f0_hat_star)
        return(c_hat)
        
      } else if (q==1){  
        
        c_hat = 1 - (f1/n)*(1-A)
        return(c_hat)
        
      } else if (q==2){
        
        x = freq[freq>=2]
        c_hat = 1 - (f1/n)*( (A*(1-A))/sum( x*(x-1) / (n*(n-1)) ) )
        return(c_hat)
        
      } else {
        
        r <- 0:(n-1)
        sort.data = sort(unique(freq))
        tab = table(freq)
        term = sapply(sort.data,function(z){
          k=0:(n-z)
          sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
        })
        lambda_hat =  sum(tab*term) + ifelse(f1==0|A==1,0,f1/n*(1-A)^(1-n)*(A^(q-1)-sum(choose(q-1,r)*(A-1)^r)))
        c_hat = 1 - ((f1/n)*(A^(q-1))*(1-A)/lambda_hat)
        return(c_hat)
        
      }
    }
  } else {
    
    t = freq[1]
    freq = freq[-1]; freq = freq[freq>0]
    u = sum(freq)
    Q1 = sum(freq==1)
    Q2 = sum(freq==2)
    B = ifelse(Q2>0,2*Q2/((t-1)*Q1+2*Q2),ifelse(Q1>0,2/((t-1)*(Q1-1)+2),1))
    
    c_hat = function(q){
      if (q==0){
        
        S_obs = length(freq)
        Chao2 = S_obs + ceiling(if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )})
        c_hat = S_obs / Chao2
        return(c_hat)
        
      } else if (q==1){  
        
        c_hat = 1 - (Q1/u)*(1-B)
        return(c_hat)
        
      } else if (q==2){
        
        x = freq[freq>=2]
        c_hat = 1 - (t-1)*Q1*( (B*(1-B))/sum( x*(x-1) ) )
        return(c_hat)
        
      } else {
        
        r <- 0:(t-1)
        sort.data = sort(unique(freq))
        tab = table(freq)
        term = sapply(sort.data,function(z){
          k=0:(t-z)
          sum(choose(k-q,k)*exp(lchoose(t-k-1,z-1)-lchoose(t,z)))
        })
        phi_hat = sum(tab*term) + ifelse(Q1==0|B==1,0,Q1/t*(1-B)^(1-t)*(B^(q-1)-sum(choose(q-1,r)*(B-1)^r)))
        c_hat = 1 - ((Q1/t)*(B^(q-1))*(1-B)/phi_hat)
        return(c_hat)
      }
    }
    
  }
  
  sapply(q,c_hat)
  
}




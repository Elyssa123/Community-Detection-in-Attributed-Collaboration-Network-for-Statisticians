#' ECV
#' We thank the code provided by
#' @article{li2020network, title={Network cross-validation by edge sampling}, author={Li, Tianxi and Levina, Elizaveta and Zhu, Ji}, journal={Biometrika}, volume={107}, number={2}, pages={257--276}, year={2020}, publisher={Oxford University Press}}


library(AUC)
library(softImpute)
library(Matrix)
library(irlba)

ECV.undirected.Rank.weighted <- function(A,max.K,B=3,holdout.p=0.1,fast=fast){
  
  # input: 
  #    A                -  adjacency matrix of a graph
  #    max.K            -  the largest number of communities
  #    B                - the times of sampling
  #    holdout.p        -  the proportion of hold out nodes 
  # output: SSE matrix
  
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  
  holdout.index.list <- list()
  holdout.n <- floor(holdout.p*edge.n)
  
  for(j in 1:B){
    holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)  # sample edges
  }
  result <- lapply(holdout.index.list,missing.undirected.Rank.weighted.fast.all,A=A,max.K=max.K,fast=fast,p.sample=1-holdout.p)
  sse.mat <- roc.auc.mat <- matrix(0,nrow=B,ncol=max.K)
  
  for(b in 1:B){
    sse.mat[b,] <- result[[b]]$sse
  }
  
  sse.seq <- colMeans(sse.mat)
  sse.sd <- apply(sse.mat,2,sd)/sqrt(B)
  return(list(sse=sse.seq,sse.sd=sse.sd))
}


missing.undirected.Rank.weighted.fast.all <- function(A,max.K,holdout.index,fast=fast,p.sample=1){
  
  # input: 
  #    A                -  adjacency matrix of a graph
  #    max.K            -  the largest number of communities
  #    B                -  the times of sampling
  #    holdout.index    -  the indices of hold out nodes 
  #    p.sample         -  the proportion of sampling nodes
  # output: imputed adjacency metrix, the positions of hold out nodes, SSE
  
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  A.new <- matrix(0,n,n)
  A.new[upper.tri(A.new)] <- A[edge.index]
  A.new[edge.index[holdout.index]] <- NA
  A.new <- A.new + t(A.new)
  diag(A.new) <- diag(A)
  degrees <- colSums(A.new,na.rm=TRUE) 
  no.edge <- 0
  no.edge <- sum(degrees==0)
  
  Omega <- which(is.na(A.new))
  imputed.A <- list()
  sse <- roc.auc <- rep(0,max.K)
  SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)
  for(k in 1:max.K){
    tmp.est <- SVD.result[[k]]
    A.approx <- tmp.est$A
    response <- A[Omega]
    predictors <- A.approx[Omega]
    sse[k] <- mean((response-predictors)^2)
    imputed.A[[k]] <- A.approx
  }
  return(list(imputed.A=imputed.A,Omega=Omega,sse=sse))
}


iter.SVD.core.fast.all <- function(A,Kmax,sparse=TRUE,tau=0,fast=FALSE,p.sample=1){
  
  # input: 
  #    A                -  adjacency matrix of a graph
  #    Kmax             -  the largest number of communities
  #    p.sample         -  the proportion of sampling nodes
  # output: SVD result
  
  if(sparse) A <- Matrix(A,sparse=TRUE)
  avg.p <- mean(as.numeric(A),na.rm=TRUE)
  cap <- 1
  A[which(is.na(A))] <- 0
  A <- A/p.sample
  svd.new <- irlba(A,nu=Kmax,nv=Kmax)
  result <- list()
  for(K in 1:Kmax){
    print(K)
    if(K==1){
      A.new <- svd.new$d[1]*matrix(svd.new$u[,1],ncol=1)%*%t(matrix(svd.new$v[,1],ncol=1))
    }else{
      A.new <- A.new + svd.new$d[K]*matrix(svd.new$u[,K],ncol=1)%*%t(matrix(svd.new$v[,K],ncol=1))
    }
    A.new.thr <- A.new
    A.new.thr[A.new < 0+tau] <- 0+tau
    A.new.thr[A.new >cap] <- cap
    
    tmp.SVD <- list(u=svd.new$u[,1:K],v=svd.new$v[,1:K],d=svd.new$d[1:K])
    result[[K]] <- list(iter=NA,SVD=tmp.SVD,A=A.new,err.seq=NA,A.thr=A.new.thr)
  }
  return(result)
}

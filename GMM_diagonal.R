library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(RcppMLPACK)


f_gmm<- function(samp,g,threshold=1e-6){
  m <- ncol(samp)
  n <- nrow(samp) 
  prob <- matrix(rep(0,g*n),nrow=n)
  weight <- matrix(rep(0,g*n),nrow=n)
  #####inicialize proportion,mean vector and covariance matrix
  kc <- kmeans(samp,g,nstart=10,iter.max = 50)  #k-means
  pie <- kc$size/n      # proportion
  miu <- t(as.matrix(kc$centers)) #mean vector
  #covariance matrix for each cluster
  sig <- list()
  vari <- rep(0,m)
  for(j in 1:g){
    for(w in 1:m) vari[w] <- var(samp[kc$cluster==j,w])
    sig[[j]] <- diag(vari)
  }
  
  # initialised all the parameters (log)
  for(j in 1:g){
    for(i in 1:n) prob[i,j]=dmvnorm(samp[i,],miu[,j],sig[[j]])
    weight[,j] <- pie[j]*prob[,j]
  }
  row_sum <- rowSums(weight)
  prob <- weight/row_sum
  oldlog <- sum(log(rowSums(weight))) 
  BIC <- -2*oldlog+((m+m+1)*g-1)*log(n)  ###diagonal matrix for covariance
  newlog <- oldlog
  s <- 0
  
  out <- iteration_gmm(step=200,g,n,m,samp,prob,miu,sig,pie,threshold,oldlog,BIC,newlog,s)
  
  ##Label 
  label <- rep(0,n)
  for(i in 1:n) label[i] <- which.max(out$prob[i,])
  return(list(step=out$step,pie=out$pie,mu=out$miu,sigma=out$sig,BIC=out$BIC,label=label,log=out$log))
}

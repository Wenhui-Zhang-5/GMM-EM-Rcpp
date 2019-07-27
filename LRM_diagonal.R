library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(RcppMLPACK)


f_lrm <- function(samp,g,covariates,threshold=1e-6){
  m <- ncol(samp)
  n <- nrow(samp)
  X <- covariates
  k <- ncol(as.matrix(X[[1]]))  ### the number of covariates 
  for(w in 1:m) X[[w]] <- cbind(rep(1,n),X[[w]]) ## m design matrixes(n*k+1) in a list
  beta <- list()
  for(j in 1:g) beta[[j]] <- matrix(rep(0,(k+1)*m),ncol=m)
  
  prob <- matrix(rep(0,g*n),nrow=n)
  weight <- matrix(rep(0,g*n),nrow=n)
  t <- matrix(rep(0,m*n),nrow=n)
  l <- list()
  
  #####inicialize proportion,mean vector and covariance matrix
  kc <- kmeans(samp,g,nstart=20,iter.max = 50)  #k-means
  pie <- kc$size/n      # proportion
  miu <- t(as.matrix(kc$centers)) #mean vector m*g
  ####initialize beta if there are more than one covariate:using LM function
  for (j in 1:g){
    for(w in 1:m){
      fit <- lm(samp[kc$cluster==j,w]~X[[w]][kc$cluster==j,-1])
      beta[[j]][,w] <- fit$coefficients
    }
  }
  
  #covariance matrix for each cluster
  sig <- list()
  vari <- rep(0,m)
  for(j in 1:g){
    for(w in 1:m) vari[w] <- var(samp[kc$cluster==j,w])
    sig[[j]] <- diag(vari)
  }
  
  ##Generate list for updating sigma
  store <- list()
  for(j in 1:g) store[[j]] <- matrix(rep(0,n*m),nrow=n)
  
  
  # initialised all the parameters (log)
  #E step
  for(j in 1:g){
    for(i in 1:n) {
      for(w in 1:m) t[i,w] <-X[[w]][i,]%*%beta[[j]][,w] 
      prob[i,j] <- dmvnorm(samp[i,],t[i,],sig[[j]])
    }
    store[[j]] <- t
    weight[,j] <- pie[j]*prob[,j]
  }
  
  row_sum <- rowSums(weight)
  prob <- weight/row_sum
  oldlog <- sum(log(rowSums(weight))) 
  BIC <- -2*oldlog+(((k+1)*m + m +1)*g-1)*log(n)
  s <- 0
  
  out <- iteration_LRM2(step=200,g,n,m,k,samp,prob,sig,X,beta,store,pie,threshold=1e-6,oldlog,BIC,s)
  
  ##Label 
  label <- rep(0,n)
  for(i in 1:n) label[i] <- which.max(out$prob[i,])
  return(list(step=out$step,pie=out$pie,beta=out$beta,sigma=out$sig,BIC=out$BIC,label=label))
}



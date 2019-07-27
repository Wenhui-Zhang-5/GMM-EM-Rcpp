##Simulation of m measurement g cluster k covariates
# m,g,k,N are numbers
# pie should be a vector with length of g and sum up to 1
library(mvtnorm)
data_lrm <- function(m,g,k,N,pie,error_mean=0.5,error_cov=0.05,range_mu=10,range_int=5,range_beta=5){
  N <- N #sample size
  pie <- pie # proportion
  g <- g
  m <- m
  k <- k
  n <- rep(0,g) # vector of group size 
  for(j in 1:g) n[j] <- N*pie[j] # calculate the group size
  mu <- list() # mean of covariates
  ## generate the mean vector of covariates at random for each group
  for(j in 1:g) {
    mu[[j]] <- matrix(rep(0,m*k),nrow=m)
    for(co in 1:k) mu[[j]][,co] <-sample(-range_mu:range_mu,size=m) ##sample() -- random integer
  }
  
 
  ##generate the covariance matrix of covariates at random for each group
  ##for simplicity, the covariance matrix for different covariates within a group is set to be the same 
  sigma <- list() #covariance of covariates
  for(j in 1:g){
    s <- runif(1,0,0.5)
    sigma[[j]] <- diag(s,m)
  }
  
 
  ##generate data for covariates(k) for different groups
  data_x <- list()
  for(j in 1:g){
    data_c <- list()
    for(co in 1:k) data_c[[co]] <- rmvnorm(n[j],mu[[j]][,co],sigma[[j]])
    data_x[[j]] <- data_c
  }
  
  ###combine all the data into formal form
  covariates <- list()
  for(w in 1:m){
    temp <- list()
    temp2 <- list()
    for(j in 1:g){
      for(co in 1:k) temp[[co]]<- data_x[[j]][[co]][,w]
      temp2[[j]] <- do.call("cbind",temp)
    }
    covariates[[w]] <- do.call("rbind",temp2)
  }
  
  ####generate error&intercept term
  #we control 0.05 for error term
  #for simplicity, we set the intercept term within group to be the same across m measurements
  e <- list()
  int <- list()
  intt <- list()
  for(j in 1:g) {
    random <- rep(0,m)
    e[[j]] <- rmvnorm(n[j],rep(error_mean,m),sigma=diag(error_cov,m))
    random<- sample(-range_int:range_int,size=m)
    intt[[j]] <- random
    int[[j]] <- matrix(rep(0,n[j]*m),ncol=m)
    for(w in 1:m)int[[j]][,w] <- runif(n[j],min=random[w]-0.02,max=random[w]+0.02)
    
  }
   
  
  ##beta
  beta <- list()
  betaa <- list()
  tbeta <- matrix(rep(0,k*m),nrow=k)
  for(j in 1:g){
    for(w in 1:m) tbeta[,w]<- sample(-range_beta:range_beta,size=k)
    beta[[j]] <- tbeta
    betaa[[j]] <- rbind(t(intt[[j]]),tbeta)
  }
  
  
  ###generate measurements
  samp <- matrix(rep(0,N*m),nrow=N,ncol=m)
  te <- list()
  X <- list()
  for(j in 1:g){
    te2 <- matrix(rep(0,n[j]*k),ncol=k)
    for(w in 1:m){
      for(co in 1:k) te2[,co] <- data_x[[j]][[co]][,w]
      te[[w]] <- te2
    }
    X[[j]] <- te
  }
  
  
###do.call cbind

  y <- list()
  tem <- list()
  for(j in 1:g){
    for(w in 1:m){
      tem[[w]] <- X[[j]][[w]]%*%beta[[j]][,w]+int[[j]][,w]
    }
    y[[j]] <- do.call("cbind",tem)+e[[j]]
  }
  Y <- do.call("rbind",y)
  
  return(list(measurement=Y,covariates=covariates,beta=betaa))

} 


set.seed(870)
d <- data_lrm(m=2,g=3,k=1,N=1000,pie=c(0.3,0.4,0.3))
system.time(result2 <- f2(d$measurement,g=3,d$covariates))
system.time(result3 <- f_lrm(d$measurement,g=3,d$covariates))

result2$step
result3$step

result2$pie
result3$pie

d$beta
result2$beta
result3$beta
result3$sigma



set.seed(854)
d <- data_lrm(m=2,g=3,k=1,N=1000,pie=c(0.3,0.4,0.3))
system.time(result2 <- f2(d$measurement,g=3,d$covariates))
system.time(result3 <- f_lrm(d$measurement,g=3,d$covariates))

result2$step
result3$step

result2$pie
result3$pie

d$beta
result2$beta
result3$beta


set.seed(873)
d <- data_lrm(m=6,g=3,k=2,N=1000,pie=c(0.3,0.4,0.3))
system.time(result2 <- f2(d$measurement,g=3,d$covariates))
system.time(result3 <- f_lrm(d$measurement,g=3,d$covariates))

result2$step
result3$step

result2$pie
result3$pie

d$beta
result2$beta
result3$beta






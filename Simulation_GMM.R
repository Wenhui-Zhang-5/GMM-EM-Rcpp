##Simulation of m measurement g cluster k covariates
# m,g,k,N are numbers
# pie should be a vector with length of g and sum up to 1
library(mvtnorm)

data_gmm <- function(m,g,N,pie,error=0.3,range=20){
  N <- N #sample size
  pie <- pie # proportion
  g <- g
  m <- m
  n <- rep(0,g) # vector of group size 
  for(j in 1:g) n[j] <- N*pie[j] # calculate the group size
  mu <- matrix(rep(0,g*m),nrow=m) # mean of measurements
  
  ## MODIFY::generate the mean vector 
 
  for(j in 1:g) mu[,j] <- sample(-range:range,size=m)
  
  
  ##generate the covariance matrix each group
  sigma <- list()
  for(j in 1:g){
    s <- runif(1,1,2)    
    sigma[[j]] <- diag(s,m)
  }
  
  ####generate error
  #we control 0.3 for error term
 
  e <- list()
  for(j in 1:g) {
    e[[j]] <- rmvnorm(n[j],rep(0.1,m),sigma=diag(error,m))
  }
  
 
  ###generate measurements
  samp <- matrix(rep(0,N*m),nrow=N,ncol=m)
  temp <- list()
  for(j in 1:g){
    temp[[j]] <- rmvnorm(n[j],mu[,j],sigma[[j]])+e[[j]]
  }
  samp <- do.call("rbind",temp)
  
  return(list(mu=mu,sigma=sigma,measurement=samp))
  
} 



set.seed(459)
d <- data_gmm(m=2,g=3,N=5000,pie=c(0.3,0.3,0.4),error=0.6)
system.time(result2 <- f2(d$measurement,g=3))
system.time(result3 <- f_gmm(d$measurement,g=3))

result2$step
result3$step

result2$pie
result3$pie

d$mu
result2$mu
result3$mu
result3$sigma
####As the number of measurements gets bigger,it is easier to cluster. However, the error term affects more in that case.

set.seed(476)
d<- data_gmm(m=10,g=3,N=5000,pie=c(0.3,0.3,0.4),error=0.01)
system.time(result2 <- f2(d$measurement,g=3))
system.time(result3 <- f_gmm(d$measurement,g=3))
result2$step
result3$step

result2$pie
result3$pie

d$mu
result2$mu
result3$mu


####When doing variable selection

set.seed(407)
d2<- data_gmm(m=10,g=3,N=5000,pie=c(0.3,0.3,0.4),error=0.01)
d2$mu



library(mclust)
mod1 <- Mclust(d2$measurement)
summary(mod1)




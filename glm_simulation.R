##Simulation of m measurement g cluster k covariates
# m,g,k,N are numbers
# pie should be a vector with length of g and sum up to 1

data_glm <- function(m,g,k,N,pie,range_int=5,range_beta=5){
  N <- N #sample size
  pie <- pie # proportion
  g <- g
  m <- m
  k <- k
  n <- rep(0,g) # vector of group size 
  for(j in 1:g) n[j] <- N*pie[j] # calculate the group size

  ### all the covariates have standard normal distribution. norm
  
  
  ##generate data for covariates(k) for different groups at different measurements
  data_xg <- list()  ### g elements
  
   for(j in 1:g){
     data_gc <- list()  ###m elements  each element is a matrix n[j]*k
     for(w in 1:m) {
       data_gc[[w]] <- matrix(rep(0,n[j]*k),ncol=k)
       for(co in 1:k) data_gc[[w]][,co] <- rnorm(n[j])
     }
     data_xg[[j]] <- data_gc
   }  
  
  
  
  ###combine all the data into formal form
  covariates <- list()
   for(w in 1:m){
     temp <- list()
     for(j in 1:g) temp[[j]] <- data_xg[[j]][[w]]
     covariates[[w]] <- do.call("rbind",temp)
   }
  
  
  ####generate intercept term
  int <- list()
  intt <- list()
  random <- rep(0,m)
  for(j in 1:g) {
    int[[j]] <- matrix(rep(0,n[j]*m),ncol=m)
    for(w in 1:m){
      random[w]<- sample(-range_int:range_int,size=1)
      int[[j]][,w] <- runif(n[j],min=random[w]-0.01,max=random[w]+0.01)
    }
    intt[[j]] <- random
  }
  
  

  ##beta
  beta <- list()
  tbeta <- matrix(rep(0,k*m),nrow=k)
  for(j in 1:g){
    for(w in 1:m) tbeta[,w]<- sample(-range_beta:range_beta,size=k)
    beta[[j]] <- tbeta
  }
  
  
  ###generate Z (mean)
  samp <- matrix(rep(0,N*m),nrow=N,ncol=m)
  te <- list()
  X <- list()
  for(j in 1:g){
    te <- list()
    for(w in 1:m) te[[w]] <- data_xg[[j]][[w]]
    X[[j]] <- te
  }
  
  ###do.call cbind
  
  y <- list()
  tem <- list()
  for(j in 1:g){
    for(w in 1:m){
      tem[[w]] <- X[[j]][[w]]%*%beta[[j]][,w]+int[[j]][,w]
      tem[[w]] <- exp(tem[[w]])/(1+exp(tem[[w]]))
      tem[[w]] <- rbinom(n[j],1,tem[[w]])
    }
    y[[j]] <- do.call("cbind",tem)
  }
  Y <- do.call("rbind",y)
  
  return(list(measurement=Y,covariates=covariates,beta=beta,intercept=intt))
  
} 

set.seed(231)
d <- data_glm (m=2,g=2,k=1,N=1000,pie=c(0.4,0.6),range_int=5,range_beta=5)
d$beta
d$intercept
set.seed(237)
d1 <- data_glm (m=2,g=3,k=2,N=1000,pie=c(0.3,0.3,0.4),range_int=5,range_beta=5)

set.seed(239)
d <- data_glm (m=8,g=2,k=1,N=1000,pie=c(0.4,0.6),range_int=5,range_beta=5)
d$intercept
d$beta









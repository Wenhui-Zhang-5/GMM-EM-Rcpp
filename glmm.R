library(mvtnorm)
 ####Simulation: 2 measurement with 1 covariate --- 2groups
n <- 1000
g <- 2
m <- 2
k <- 1
pie1 <- 0.4
pie2 <- 0.6
n1 <- n*pie1
n2 <- n*pie2

###generate covariate and measurement for measurement 1 seperately for each group
x11 <- rnorm(n1)
z11 <- 2*x11+2
x21 <- rnorm(n2)
z21 <- -3*x21-1
x1<- c(x11,x21)


y11 <- exp(z11)/(1+exp(z11))
y11 <- rbinom(n1,1,y11)
y21 <- exp(z21)/(1+exp(z21))
y21 <- rbinom(n2,1,y21)
y1 <- c(y11,y21)


###generate covariate and measurement for measurement 1 seperately for each group
x12 <- rnorm(n1)
z12 <- 4*x12+1
x22 <- rnorm(n2)
z22 <- -2*x22-2
x2<- c(x12,x22)


y12 <- exp(z12)/(1+exp(z12))
y12 <- rbinom(n1,1,y12)
y22 <- exp(z22)/(1+exp(z22))
y22 <- rbinom(n2,1,y22)
y2 <- c(y12,y22)



samp <- cbind(y1,y2)
X <- list(x1,x2)


f_glm <- function(samp,g,covariates,step=1000){
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
  tt <- matrix(rep(0,m*n),nrow=n)
  l <- list()
  
  #####inicialize proportion,mean vector and covariance matrix
  kc <- kmeans(samp,g)  #k-means
  pie <- kc$size/n      # proportion
  miu <- t(as.matrix(kc$centers)) #mean vector m*g

  ###initialize beta
  for (j in 1:g){
    for(w in 1:m){
      
      beta[[j]][,w] <- runif(k+1,min=-5,max=5)
    }
  }
  
  
  ###initial E step
  for(j in 1:g){
    for(i in 1:n) {
      for(w in 1:m) {
        t[i,w] <-exp(X[[w]][i,]%*%beta[[j]][,w])/(1+exp(X[[w]][i,]%*%beta[[j]][,w]))
        tt[i,w] <- dbinom(samp[i,w],1,prob=t[i,w])}
      prob[i,j] <- prod(tt[i,])
    }
    weight[,j] <- pie[j]*prob[,j]
  }
  
  row_sum <- rowSums(weight)
  prob <- weight/row_sum
  old.logdata <- sum(log(rowSums(weight))) 
  

  
  for(s in 1:step){
    
    #M step
    
    for(j in 1:g){
      sum1 <- sum(prob[,j])
      pie[j] <- sum1/n
      W <-prob[,j]
      
      for(w in 1:m) {
        fit <- glm(samp[,w]~X[[w]][,-1],family=binomial,weights = W)
        beta[[j]][,w] <- fit$coefficients
      }
    }
    
    ##E step
    for(j in 1:g){
      for(i in 1:n) {
        for(w in 1:m) {
          t[i,w] <-exp(X[[w]][i,]%*%beta[[j]][,w])/(1+exp(X[[w]][i,]%*%beta[[j]][,w]))
          tt[i,w] <- dbinom(samp[i,w],1,prob=t[i,w])}
        prob[i,j] <- prod(tt[i,])
      }
      weight[,j] <- pie[j]*prob[,j]
    }
    
    row_sum <- rowSums(weight)
    prob <- weight/row_sum
    logdata <- sum(log(rowSums(weight))) 
    BIC <- -logdata+((k+1)*m*g+g-1)*log(n)
    
    #set threshold
    threshold <- 1e-6
    if(abs(logdata-old.logdata)<threshold) 
      break
    cat('step',s,'pie',pie,'logdata',logdata,'BIC',BIC,'\n')
    old.logdata <- logdata
  }
  
  ##Label
  label <- rep(0,n)
  for(i in 1:n) label[i] <- which.max(prob[i,])
  return(list(pie=pie,beta=beta,BIC=BIC,label=label,s=s))
  
}


r <- f_glm(samp=d$measurement,g=2,covariates=d$covariates)
r$s
system.time(r <- f_glm(samp=d$measurement,g=2,covariates=d$covariates))
r$pie
table(r$label)
la <- c(rep(1,400),rep(2,600))
table(la,r$label)

d$beta
d$intercept

r$pie
r$beta

r1 <- f_glm(samp=d1$measurement,g=3,covariates=d1$covariates)

d$beta
d$intercept

r$beta


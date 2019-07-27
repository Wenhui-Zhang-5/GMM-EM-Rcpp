#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double Mahalanobis(arma::vec x, arma::vec center, arma::mat cov){
  arma::vec h=(x-center).t()*arma::inv(cov)*(x-center);
  double r=sum(h);
  return r;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double dmvnorm ( arma::vec x,  arma::vec mean,  arma::mat sigma){
  
  int mm=x.n_rows;
  double distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double log2pi = 1.8378770664093454835606594728112352797227949472755668;
  double  logretval = -( (mm * log2pi + logdet + distval)/2  ) ;
  return(exp(logretval));
  
}



//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List iteration_gmm(int step,int g,int n,int m,arma::mat samp, arma::mat prob, arma::mat miu, List sig,NumericVector pie,double threshold,double oldlog, double BIC,double newlog,int s){
  for( s=0;s<step;s++){
    oldlog=newlog;
    // M step
    for(int j=0; j<g;j++){
      double sum1=sum(prob.col(j));
      pie(j)=sum1/n;
      arma::mat t(n,m);
      NumericMatrix tem(m,m);
      arma::mat temp= as<arma::mat>(tem);
      //Mstep for miu
      for(int i=0; i<n; i++){
        t.row(i)=prob(i,j)*samp.row(i);
      }
      miu.col(j)=(sum(t,0)/sum1).t(); 
      
      //Mstep for covariance
      for(int i=0;i<n;i++){
        arma::vec ve=samp.row(i).t()-miu.col(j);
        arma::vec vee=ve%ve;
        temp=temp+prob(i,j)*diagmat(vee);
      }
      sig[j]=temp/sum1;
    }
    
    
    //Estep 1
    
    
    arma::mat store(n,g);
    arma::mat weight(n,g);
    
    for(int j=0; j<g ;j++){
      arma::mat cova=sig[j];
      for(int i=0 ; i<n; i++){
        store(i,j)= dmvnorm(samp.row(i).t(),miu.col(j),cova);
      }
      weight.col(j)=pie(j)*store.col(j);
    }
    
    
    
    
    //Estep 2:normalization & calculate loglikelihood
    
    arma::vec row_sum=sum(weight,1);
    for(int j=0;j<g;j++) prob.col(j)=weight.col(j)/row_sum;
    newlog=sum(log(row_sum));
    BIC= -2*oldlog+((m+m+1)*g-1)*log(n);
    
    
    // Convergence or not
    if(std::abs(newlog-oldlog)<threshold) break;
    
    
  }
  List outcome=List::create(Named("pie")=pie,_["miu"]=miu,_["sig"]=sig,_["prob"]=prob,_["log"]=newlog,_["BIC"]=BIC,_["step"]=s+1);
  return(outcome);
}


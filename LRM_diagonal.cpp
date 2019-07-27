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

List iteration_LRM2(int step,int g,int n,int m,int k,arma::mat samp, arma::mat prob, List sig,List X, List beta,List store,NumericVector pie,long double threshold,long double oldlog,long double BIC,int s){
  for(s=0;s<step;s++){
    
    // M step
    List beta_new=beta;
    for(int j=0; j<g;j++){
      double sum1=sum(prob.col(j));
      pie(j)=sum1/n;
      arma::mat W=diagmat(prob.col(j));
      arma::mat bet((k+1),m);
      //Mstep for beta
      for(int w=0;w<m; w++){
        arma::mat x=X[w];
        arma::vec be=inv(x.t()*W*x)*x.t()*W*samp.col(w);
        bet.col(w)=be;
      }
      beta_new[j]=bet;
    }
    
    // M step for covariance matrix
    
    for(int j=0;j<g;j++){
      double sum1=sum(prob.col(j));
      NumericMatrix tem(m,m);
      arma::mat temp= as<arma::mat>(tem);
      arma::mat te=store[j];
      for(int i=0;i<n;i++){
        arma::vec ve=samp.row(i).t()-te.row(i).t();
        arma::vec vee=ve%ve;
        temp=temp+prob(i,j)*diagmat(vee);
      }
      sig[j]=temp/sum1;
    }
    
    
    
    //Estep 1
    arma::mat weight(n,g);
    for(int j=0; j<g ;j++){
      arma::mat t(n,m);
      arma::mat be=beta[j];
      arma::mat cova=sig[j];
      for(int i=0 ; i<n; i++){
        for(int w=0;w<m;w++) {
          arma::mat cv=X[w];
          NumericVector re=wrap(cv.row(i)*be.col(w));
          t(i,w)=as<long double>(re);
        }
        prob(i,j)=dmvnorm(samp.row(i).t(),t.row(i).t(),cova);
      }
      store[j]=t;
      weight.col(j)=pie(j)*prob.col(j);
    }
    
    //Estep 2:normalization & calculate loglikelihood
    
    arma::vec row_sum=sum(weight,1);
    for(int j=0;j<g;j++) prob.col(j)=weight.col(j)/row_sum;
    long double newlog=sum(log(row_sum));
    BIC = -2*newlog+(((k+1)*m + m +1)*g-1)*log(n);
    
    
    // Convergence or not
    if(std::abs(newlog-oldlog)<threshold) break;
     oldlog=newlog;
    
  }
  
  List outcome=List::create(Named("pie")=pie,_["beta"]=beta,_["sig"]=sig,_["prob"]=prob,_["BIC"]=BIC,_["step"]=s+1);
  return(outcome);
  
}

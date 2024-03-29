// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
#include <RcppParallel.h>
#include <iostream>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;
using namespace Eigen;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double weighted_median(const Eigen::VectorXd & beta, const Eigen::VectorXd & weight){
  VectorXi ind = VectorXi::LinSpaced(beta.size(),0,beta.size()-1);
  double weight_half_sum = 0.5*weight.sum();
  double weight_cum = 0;
  auto rule=[beta](int i, int j)->bool{return beta(i)<beta(j);};
  std::sort(ind.data(),ind.data()+ind.size(),rule);
  Eigen::VectorXd sorted_beta(beta), sorted_weight(weight);
  int i=0;
  while(weight_cum < weight_half_sum){
    weight_cum += weight[ind(i)];
    i++;
  }
  
  return beta(ind(i-1));
}

Eigen::VectorXd QCD(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd initial,
                    Eigen::VectorXd lambda, const double thresh, const int maxin){
  
  int p = X.cols();
  int n = X.rows();
  Eigen::VectorXd pre_value_final(n+1), beta(p), xj(n), beta_j(p-1);
  Eigen::MatrixXd X_j(n,p-1);
  Eigen::MatrixXd Weight_mat(n+1,p);
  Weight_mat << X.cwiseAbs()/n, lambda.transpose();
  beta = initial;
  Eigen::VectorXd pre_value = y - X*beta;
  int countj, count = 0;
  double temp1, temp2, err = 1;
  while((err > thresh) & (count < maxin)){
    for(countj = 0; countj < p; countj++){
      xj = X.col(countj);
      pre_value_final << pre_value.array()/xj.array() + beta[countj], 0.0;
      temp1 = weighted_median(pre_value_final,Weight_mat.col(countj));
      temp2 = fabs(temp1)>thresh? temp1: 0.0;
      pre_value -= (temp2-beta[countj])*xj;
      beta[countj] = temp2;
    }
    err = (beta - initial).norm();
    initial = beta;
    count++ ;
  }
  return beta;
}


struct QICD_init : public Worker
{
  // design matrix and outcome
  const Eigen::MatrixXd& X;
  const Eigen::VectorXd& y;
  const Eigen::MatrixXd& lambda_list;
  Eigen::VectorXd& initial;
  const double thresh;
  const int maxin;
  const int maxout; 
  
  // output
  Eigen::MatrixXd& beta_final;
  
  // initialize with input and output
  QICD_init(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::MatrixXd& lambda_list,
            Eigen::VectorXd& initial, const double thresh, const int maxin, 
            const int maxout, Eigen::MatrixXd& beta_final)
    : X(X), y(y), lambda_list(lambda_list),initial(initial), thresh(thresh),
      maxin(maxin), maxout(maxout), beta_final(beta_final) {}
  
  void operator()(std::size_t begin, std::size_t end) {  
    Eigen::VectorXd beta_temp; 
    int count;
    double err;  
    for(unsigned int i = begin; i < end; ++i)
    {
      count = 0;
      err = 1;
      while((err > thresh) & (count < maxout)){
        beta_temp = QCD(X, y, initial, lambda_list.col(i), thresh, maxin);
        err = (beta_temp - initial).norm();
        initial = beta_temp;
        count++ ;
      }
      beta_final.row(i) = beta_temp;
    }
  }
};


// [[Rcpp::export]]
Eigen::MatrixXd QICD(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::MatrixXd lambda_list, 
                     Eigen::VectorXd initial, const double thresh, const int maxin, const int maxout){
  
  int p = X.cols();
  int nlam = lambda_list.cols();  
  Eigen::MatrixXd beta_final(nlam,p); 
  
  // pass input and output
  QICD_init qicd_init(X,y,lambda_list,initial,thresh,maxin,maxout,beta_final);
  
  // parallelFor to do it
  parallelFor(0, nlam, qicd_init); 
  return beta_final;
}

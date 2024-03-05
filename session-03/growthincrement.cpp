#include <TMB.hpp>
#include <iostream>

template <class Type> 
Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Minimal example */
  DATA_VECTOR(len1);
  DATA_VECTOR(len2);
  DATA_VECTOR(dt);
  int n = len1.size();
  
  PARAMETER_VECTOR(u_age);
  PARAMETER(log_sigma_inf);
  PARAMETER(log_sigma);
  PARAMETER(log_mu_linf);
  PARAMETER(log_k);
  PARAMETER(log_shape);
  PARAMETER(log_scale);
  
  Type nll=Type(0.0);
  
  //random effects
  /*
  nll -= sum(dnorm(u_age,Type(0),Type(1),true));       // Assign N(0,1) distribution u 
  vector<Type> v = pnorm(u_age,Type(0),Type(1));  // Uniformly distributed variables (on [0,1])
  //std::cout << u_age << std::endl;
  //std::cout << v << std::endl;
  vector<Type> w =   qgamma(v,exp(log_shape),exp(log_scale));  
  //std::cout << w << std::endl;
  */
  
  nll -= sum(dnorm(u_age,log_shape,exp(log_scale),true));       // Assign N(0,1) distribution u 
  
  // predictions
  vector<Type> sigma1(n);
  vector<Type> sigma2(n);
  vector<Type> f1(n);
  vector<Type> f2(n);
  vector<Type> mu1(n);
  vector<Type> mu2(n);
  vector<Type> rho(n);
  vector<Type> q(n);
  vector<Type> h(n);
  
  for (int i=0;i<n;i++) {
    //f1(i) = 1.0-exp(-exp(log_k)*w(i));
    //f2(i) = 1.0-exp(-exp(log_k)*(w(i)+dt(i)));
    f1(i) = 1.0-exp(-exp(log_k)*exp(u_age(i)));
    f2(i) = 1.0-exp(-exp(log_k)*(exp(u_age(i))+dt(i)));
    
    mu1(i) = exp(log_mu_linf)*f1(i);
    mu2(i) = exp(log_mu_linf)*f2(i);
    sigma1(i) = sqrt(square(exp(log_sigma_inf))*square(f1(i)) 
                       + square(exp(log_sigma)));
    sigma2(i) = sqrt(square(exp(log_sigma_inf))*square(f2(i)) 
                       + square(exp(log_sigma)));
    /*rho(i) = square(exp(log_sigma_inf))*f1(i)*f2(i)/(sigma1(i)*sigma2(i));
    q(i) = square((len1(i)-mu1(i))/sigma1(i));
    q(i) -= Type(2.0)*rho(i)*(len1(i)-mu1(i))*(len2(i)-mu2(i))/(sigma1(i)*sigma2(i));
    q(i) += square((len2(i)-mu2(i))/sigma2(i));
    h(i) = -log(sigma1(i))-log(sigma2(i))-log(sqrt(Type(1.0)-square(rho(i))));
    h(i) -= Type(0.5)*q(i)/(Type(1.0)-square(rho(i)));
    */
    
    matrix<Type> Sigma(2, 2);
    Sigma.row(0) << square(sigma1(i)), square(exp(log_sigma_inf))*f1(i)*f2(i);
    Sigma.row(1) << square(exp(log_sigma_inf))*f1(i)*f2(i), square(sigma2(i));
    
    using namespace density;
    MVNORM_t<Type> neg_log_dmvnorm(Sigma);
    
    vector<Type> residual(2);
    residual(0) = len1(i) - mu1(i);
    residual(1) = len2(i) - mu2(i);
    
    nll += neg_log_dmvnorm(residual);
  }
  
  Type sigma_inf = exp(log_sigma_inf);
  Type sigma = exp(log_sigma);
  Type k = exp(log_k);
  Type mu_age = exp(log_shape);
  Type sd_age = exp(log_scale);
  REPORT(exp(u_age));
  ADREPORT(sigma_inf);
  ADREPORT(sigma);  
  ADREPORT(k);
  ADREPORT(mu_age);
  ADREPORT(sd_age);
    
  return nll;

  /* Quick Reference
     ===============

     ** Macros to read data and declare parameters:

     _Template_Syntax_              _C++_type_                     _R_type_
     DATA_VECTOR(name)              vector<Type>                   vector
     DATA_MATRIX(name)              matrix<Type>                   matrix
     DATA_SCALAR(name)              Type                           numeric(1)
     DATA_INTEGER(name)             int                            integer(1)
     DATA_FACTOR(name)              vector<int>                    factor
     DATA_SPARSE_MATRIX(name)       Eigen::SparseMatrix<Type>      dgTMatrix
     DATA_ARRAY(name)               array<Type>                    array
     PARAMETER_MATRIX(name)         matrix<Type>                   matrix
     PARAMETER_VECTOR(name)         vector<Type>                   vector
     PARAMETER_ARRAY(name)          array<Type>                    array
     PARAMETER(name)                Type                           numeric(1)

     ** Macro to report intermediate expressions back to R:

     REPORT(x)
     ADREPORT(x)

     ** Basic constructors:

     vector<Type> v(n1);
     matrix<Type> m(n1,n2);
     array<Type> a(n1,n2,n3)

     ** Basic operations:

     v+v,v-v,v*v,v/v                Pointwise binary operations
     m*v                            Matrix-vector multiply
     a.col(i)                       R equivalent of a[,,i]
     a.col(i).col(j)                R equivalent of a[,j,i]
     a(i,j,k)                       R equivalent of a[i,j,k]
     exp(v)                         Pointwise math
     m(i,j)                         R equivalent of m[i,j]
     v.sum()                        R equivalent of sum(v)
     m.transpose()                  R equivalent of t(m)

     ** Distributions:

     Type dnbinom2(const Type &x, const Type &mu, const Type &var, int give_log=0)
     Type dpois(const Type &x, const Type &lambda, int give_log=0)
     Type dlgamma(Type y, Type shape, Type scale, int give_log=0)
     Type dnorm(Type x, Type mean, Type sd, int give_log=0)

     ** Parallel accumulator declaration (only methods "+=" and "-="):
     
     parallel_accumulator<Type> res(this);

  */

}

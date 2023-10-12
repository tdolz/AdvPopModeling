#include <TMB.hpp>

//Gavin's way

//Some notes:
  //Data: L1, L2, dt, for each fish(i), so L1_i, L2_i, dt_i
  //Model: L1hat_i, L2hat_i, A1_i
  //Parameters: L_inf, mu_A, sigma_Linf, sigma_obs, mu_A1, sigma_A1, K
  //sigma obs is the observation error. 
//relationships
  //length is a function of age, so: A1_i has a dist described by mu_A1, sigma_A1
  //the predicted L1hat_i is described by the growth model of Linf, K, A1_i
  //L2hat_i is described by Linf, K, A1_i, and dt_i for time at liberty.  
//the likelihood function describes how similar the predictions are to the data
  //the likelihood function has an inherent distribution. 
  //in this case the distribution is MVnormal.
  //we used MVnorm because there is covariance in the length. 
  //the params of a mvnorm are a vector of means (mu_i)
  //and a var-covar matrix which is going to be a 2x2 matrix because two timepoints.
  //the var-covar matrix will contain some correlation coefficient rho. 
  //varcovar[i] will have sigma_Linf, sigma_obs, A1_i, K
  //varcovar will look like VarL1, covar(L1,L2) then second row covar(l1,l2), varL2
//the random effects update the likelihood function 
  //in this case, A1_i with params MuA1 and SigmaA1 
  //mu_i will be the contribution to the likelihood function
//everthing subscripted by i is unique for each fish. 
  //TMB thinks about random effects as params. 
//The param section
  //inputs passed from r: L1_i, L2_i, dt_i
  

template <class Type>Type square(Type x){return x*x;}
template <class Type>Type objective_function<Type>::operator()()
{
  /* Minimal example */
  //DATA
  DATA_VECTOR(len1); // lenMark, L1, length at first capture
  DATA_IVECTOR(len2); // lenRecap, length at 2nd capture -not sure why this is an ivector
  DATA_VECTOR(dt); // delta time, t2-t1
  int n = len1.size(); //index
  
  //PARAMS
  PARAMETER_VECTOR(u_age); 
  PARAMETER(log_sigma_inf);
  PARAMETER(log_sigma);
  PARAMETER(log_mu_linf);
  PARAMETER(log_k); 
  PARAMETER(log_mu_a1); 
  PARAMETER(log_sigma_a1);
  
  // define the likelihood function
  Type nll = 0.0;
  
  //contribution to the A_i random effects
  nll-=sum(dnorm(age_1,log_mu_age1,exp(log_sigma_a1),true)); //true=log
  
  //create storage vectors
  //length = n = length of original dataframe
  vector<Type> sigma1(n);
  vector<Type> sigma2(n);
  vector<Type> f1(n);
  vector<Type> f2(n);
  vector<Type> mu1(n);
  vector<Type> mu2(n);
  vector<Type> rho(n);
  

  for (int i=0;i<n;i++) {
    //equations 1.3, 1.4 in Laslett et al. 2002
    f1=1-exp(-exp(log_k)*exp(u_age(i))); //growth before first encounter
    f1=1-exp(-exp(log_k)*(exp(u_age(i))+dt(i))); //growth between encounters
    
    //equations 3.2-3.6 in Laslett et al. 2002
    mu1(i)=exp(log_mu_linf)*f1(i); //adding Linf to complete the vonB for first growth increment
    mu1(i)=exp(log_mu_linf)*f2(i); //adding Linf to complete the vonB for second growth increment
    sigma1(i)=sqrt(square(exp(log_sigma_inf))*square(f1(i))
                     +square(exp(log_sigma)));
    sigma2(i)=sqrt(square(exp(log_sigma_inf))*square(f2(i))
                     +square(exp(log_sigma)));
    
    
    //populate the variance-covariance matrix
    matrix<Type> Sigma(2,2);
    Sigma.row(0) << square(sigma1(i)), square(exp(log_sigma_inf))*f1(i)*f2(i); //variance covariance
    Sigma.row(1) << square(exp(log_sigma_inf))*f1(i)*f2(i), square(sigma2(i)); //covariance variance
    
    using namespace density;
    MVNORM_t<Type> neg_log_dmvnorm(Sigma)
    
    vector<Type> residual(2)
      residual(0)=len1(i)-mu1(i);
      residual(1)=len2(i)-mu1(i);
  }
  nll += sum(residuals);
  AIC = -2*(-1*nll)+2*6;
  ADREPORT(Age);
  REPORT(AIC);
  return nll;
}  




//Derived quantities

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


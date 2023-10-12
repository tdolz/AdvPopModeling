#include <TMB.hpp>

//Gary's way
template <class Type>Type square(Type x){return x*x;}
template <class Type>Type objective_function<Type>::operator()()
{
  /* Minimal example */
  //DATA
  DATA_VECTOR(l1); // lenMark, L1, length at first capture
  DATA_IVECTOR(l2); // lenRecap, length at 2nd capture -not sure why this is an ivector
  DATA_VECTOR(dt); // delta time, t2-t1
  int n = l1.size(); //index
  
  //Recall that the vonB growth equation is L(a)=Linf*(1-exp(-k(a-t0)))
  //We have to split this up into two growth increments:
  //growth t0-t1 (aka "A") - random effect
  //growth t1-t2 
  //k is a param we will estimate
  //Linf is a param we will estimate
  //there is an overall Linf with a mean and variance (Linf_mu & sigma2_Linf)
  //Laslett et al. 2002 specify t at first encounter (A) 
  //is a random effect (A_i) for each fish
  //with a mean (mu_A) and variance (sigma2_A)

  //QUESTIONS
  //I don't know what sigma2_resid is for
  //why is K just one parameter? (fixed effect?)
  //why is the function f1 and f2 squared in the variance on A calculation? 
  //...this must have something to do with the variance covariance matrix. 
  
  //PARAMS
  PARAMETER(K); //from the vonB growth model
  PARAMETER(Linf_mu); //Linf est. for all fish
  PARAMETER(sigma2_Linf); // variance on Linf
  PARAMETER(sigma2_resid); //variance (on what??) for all fish
  
  PARAMETER_VECTOR(A_i); //random Age at first encounter for each fish, A = t1-t0
  PARAMETER(sigma2_A); //variance on A for all fish
  PARAMETER(mu_A); //mean A for all fish

  // specify object to store the likelihood function
  Type nll = 0.0;
  
  //contribution to the A_i random effects
  nll-=sum(dnorm(A_i,0,sqrt(sigma2_A),true)); //true=log
  
  //DECLARE TEMP VARIABLES
  //length = n = length of original dataframe
  vector<Type>mu1(n); //mean length at t1
  vector<Type>mu2(n); //mean length at t2
  vector<Type>sigma2_l1(n); //variance on l1
  vector<Type>sigma2_l2(n); //variance on l2
  vector<Type>cov(n); //covariance
  vector<Type>rho(n); //not even sure where this is used
  vector<Type>residuals(n); //
  vector<Type>pairs(2); //length 2 because done for each fish
  Type f1=0.0; //not sure why this is
  Type f2=0.0; //not sure why you have to initiate this function. 
  vector<Type>Age(n); //storing the ages of each fish
  //empty var-covar matrix to populate with the loop
  //mean vector for each fish has a length 2
  //thus var-covar matrix is 2x2
  //the random effect is in the variance on A
  matrix<Type> big_sigma(2,2);
  Type AIC=0.0; //not sure why declaring this but could be derived quantity
  //"namespace" a collection of functions, in this case related to the multivariate gaussian.
  using namespace density; 


  //this bit is not in Gary's code
  // Var-covar matrix is for the neg_log_dmvnorm likelihood
 //MVNORM_t<Type> neg_log_dmvnorm(sigma2_A);
    
  //evaluate the likelihood for the observed lengths given their expected values
  //by calculating the multivariate normal likelihood and 
  //adding the contribution to the objective function
  for (int i=0;i<n;i++) {
    //each age of fish at first encounter is some deviation from the mean Age
    Age(i)=mu_A+A_i(i); 
    
    //equations 1.3, 1.4 in Laslett et al. 2002
    f1=1-exp(-K*Age(i)); //growth before first encounter
    f2=1-exp(-K*Age(i)+dt(i)); //growth between encounters
    
    //equations 3.2-3.6 in Laslett et al. 2002
    mu1(i)=Linf_mu*f1; //adding Linf to complete the vonB for first growth increment
    mu2(i)=Linf_mu*f2; //adding Linf to complete the vonB for second growth increment
    sigma2_l1(i)=sigma2_resid+sigma2_Linf*square(f1); //
    sigma2_l2(i)=sigma2_resid+sigma2_Linf*square(f2);//
    cov(i)=sigma2_Linf*f1*f2; // equation 
    
    //populate the variance-covariance matrix
    big_sigma.row(0) << sigma2_l1(i), cov(i); //variance covariance
    big_sigma.row(1) << cov(i), sigma2_l2(i); //covariance variance
    
    //pairs - what is this doing?? something about the mean change in length?
    //recalculated each fish.
    pairs(0)=l1(i)-mu1(i); //length of each fish at t1 MINUS growth t0-t1
    pairs(1)=l2(i)-mu2(i); //length of each fish at t2 MINUS growth t1-t2
      //multiply the variance covariance matrix*growth increments
      //the MVNORM function is a multivariate normal distribution with a 
      //user-supplied covar matrix.
    residuals=MVNORM(big_sigma)(pairs); 
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



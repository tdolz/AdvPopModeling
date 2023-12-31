#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Minimal example  */  //declaring what type of variables things are
  DATA_VECTOR(count);
  DATA_IVECTOR(month); //adding the month
  PARAMETER_VECTOR(beta);  //because were doing months, now, beta is a vector of parameters. 
  // PARAMETER(beta);
  int nobs = count.size();
  int nmonth = month.size();
  
  vector<Type> lamda(nobs);
  
  Type nll = 0;
  
  for (int iobs=0; iobs < nobs; iobs++) {
    for(int imonth=0; imonth < nmonth; imonth++){
    lamda(iobs) = exp(beta(imonth));   //fill vector of predictions for means
    }}
    
  nll -= sum(dpois(count,lamda,true));          //calculate the objective function

  vector<Type> lamda_coef = exp(beta);
  //Type lamda_coef = exp(beta);
  ADREPORT(lamda_coef);
    
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

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data */
  DATA_VECTOR(catches);
  DATA_VECTOR(index);
  DATA_IVECTOR(index_years);
  int nyrs = catches.size();
  int nobs = index.size();
  DATA_VECTOR(selex);
  DATA_VECTOR(weight);
  DATA_VECTOR(maturity);
  int nages = maturity.size();
  DATA_INTEGER(k);

  /* PARAMETERS */
  //PARAMETER(dummy);
  PARAMETER(logR0);
  PARAMETER(logfracR1);
  PARAMETER(logith);
  PARAMETER(logM);
  PARAMETER(logSigmaC);
  PARAMETER(logSigmaI);
  PARAMETER_VECTOR(logF);

  Type Rzero = exp(logR0);
  Type R1 = Rzero*exp(logfracR1);
  Type h = 0.2 + 0.8*exp(logith)/(1.0+exp(logith));
  Type M = exp(logM);
  Type SigmaC = exp(logSigmaC);
  Type SigmaI = exp(logSigmaI);

  Type neglogL = 0.;
  
  // F & Z calculations
  matrix<Type> F(nyrs,nages);
  matrix<Type> Z(nyrs,nages);
  F.setZero();
  Z.setZero();
  for (int iyr=0; iyr<nyrs; iyr++)
    F.row(iyr) = exp(logF(iyr))*selex;
  Z = M + F.array();
  
    
  /* POPULATION DYNAMICS */
  vector<Type> Bio(nyrs+1);
  vector<Type> SSB(nyrs+1);
  vector<Type> R(nyrs+1);  
  matrix<Type> N(nyrs+1,nages);
  
  //initialize
  N.setZero();
  SSB.setZero();
  R.setZero();
  Bio.setZero();
  
  //   Initial numbers at age & 
  //   unfished Spawning biomass calculation
  N(0,0) = R1;
  for (int age=1;age<nages;age++) {
    N(0,age) = N(0,age-1)*exp(-M);
  }
  N(0,nages-1) = N(0,nages-1)/(1-exp(-M));
  
  Type SSB0 = 0.0;
  vector<Type> nb0(nages);
  nb0(0) = Rzero;
  for (int age=1;age<nages;age++) 
    nb0(age) = nb0(age-1)*exp(-M);
  nb0(nages-1) = nb0(nages-1)/(1-exp(-M));
  for (int age=1;age<nages;age++) 
    SSB0 += 0.5*nb0(age)*weight(age)*maturity(age);
   
  // 1st year stock quantities
  for (int age=0; age<nages; age++) {
    Bio(0) += selex(age)*N(0,age)*weight(age)*exp(-0.5*Z(0,age));
    if (age != 0) SSB(0) += 0.5*N(0,age)*maturity(age)*weight(age);  
  }
  R(0) = N(0,0);
  
  //   Begin Loop over time
  for (int iyr=1; iyr<=nyrs; iyr++) {
    
  //     Numbers at age update 
  for (int age=1; age<nages; age++)
    N(iyr,age) = N(iyr-1,age-1)*exp(-Z(iyr-1,age-1));
  N(iyr,nages-1) += N(iyr-1,nages-1)*exp(-Z(iyr-1,nages-1));
  
  //Biomass predictions
  Bio(iyr) = 0.;  //bio is survey available biomass
  SSB(iyr) = 0.;
  for (int age=0; age<nages; age++) {
   if (iyr != nyrs) Bio(iyr) += selex(age)*N(iyr,age)*weight(age)*exp(-0.5*Z(iyr,age));
   if (age != 0) SSB(iyr) += 0.5*N(iyr,age)*maturity(age)*weight(age);  
  }

  // Recruitment for this year given annual stock size
  R(iyr) = 4.*h*Rzero*SSB(iyr-k)/(SSB0*(1.-h)+SSB(iyr-k)*(5*h-1.));
  N(iyr,0) = R(iyr);
  
  
  //       End Loop over time
  }
  
  //       Contribution to the Likelihood function
  
  vector<Type> index_pred(nobs);
  Type q = exp(log(index/Bio(index_years)).sum()/nobs);
  index_pred = q*Bio(index_years); 

  // //estimate sigmaI analytically
  // Type SigmaI = 0.;
  // for (int iobs=0; iobs<nobs; iobs++) {
  //   SigmaI += pow(log(index(iobs))-log(index_pred(iobs)),2);
  // }
  // SigmaI = sqrt(SigmaI/nobs);
  
  vector<Type> catch_pred(nyrs);
  for (int iyr=0;iyr<nyrs;iyr++) {
    catch_pred(iyr) = 0.;
    for (int age=0; age<nages; age++)
     catch_pred(iyr) += weight(age)*F(iyr,age)*N(iyr,age)*(1.0-exp(-1.*Z(iyr,age)))/Z(iyr,age);
  }
  
  // Likelihood
  //CATCHES
  neglogL -= sum(dnorm(log(catches),log(catch_pred),SigmaC, true));
  //SURVEY
  neglogL -= sum(dnorm(log(index),log(index_pred),SigmaI, true));
  
  //neglogL = dummy*dummy;
  
  REPORT(F);
  REPORT(N);
  REPORT(Z);
  REPORT(catch_pred);
  REPORT(index_pred);
    
  ADREPORT(h);
  ADREPORT(Rzero);
  ADREPORT(SigmaI);
  ADREPORT(q);
  ADREPORT(SSB0);
  ADREPORT(SSB);
  ADREPORT(R);
  ADREPORT(Bio);
  
  return neglogL;

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

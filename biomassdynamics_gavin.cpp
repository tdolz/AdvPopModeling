#include <TMB.hpp>

//function to implement the schaefer production model for a single time step
template <class Type> Type schaefer(Type b1, Type r, Type K, Type m, Type u_i){
    Type u_use = exp(u_i)/(1.0+exp(u_i)); //logit transforming exploitation rate Ui
    Type b2 = b1 + b1*r*(1.0-pow(b1/K,m-1.0)) - u_use*b1; //literally just the shaefer model equation.
    return b2;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* DATA */
  DATA_VECTOR(catches);
  DATA_VECTOR(index);
  DATA_IVECTOR(index_years); //integer vector
  int nyrs = catches.size(); //number of unique observations of catch
  int nobs = index.size(); // number of unique observations of the survey biomass. 
  
  /* PARAMETERS */
  PARAMETER(logK);
  PARAMETER(logr);
  PARAMETER(m);
  PARAMETER(logSigmaC);
  PARAMETER(logSigmaI);
  //PARAMETER(logq);
  PARAMETER_VECTOR(upar); //exploitation rate in each year. 

  Type neglogL = 0.;
  
  /* POPULATION DYNAMICS */
  vector<Type> biomass(nyrs+1); //a vector where we will store predicted true biomass, adding an extra year because we start at 0. 
  biomass(0) = exp(logK); //the first year of biomass is carrying capacity. Looks like we are storing this in unlogged space. 
  for (int iyr=1;iyr<=nyrs;iyr++) {
    //schaefer biomass upodate goes here!
    //Type u_use = exp(upar(iyr-1))/(1.0+exp(upar(iyr-1)));
    biomass(iyr) = schaefer(biomass(iyr-1),exp(logr),exp(logK),m,upar(iyr-1)); //giving it all the logged params
    //biomass(iyr) = biomass(iyr-1)*(1.0+exp(logr)*(1.0-pow(biomass(iyr-1)/exp(logK),m-1.0)) - u_use);

    //std::cout << iyr << " " << biomass(iyr) << std::endl;
  }
    
  //model predictions
  //predicted surveys

  
  vector<Type> index_pred(nobs); //store predicted survey biomass (different from true biomass)
  //estimate q analytically
  Type logq = 0.;
  vector<Type> predbio(nobs); //store predicted survey biomass (different from true biomass)
  for (int iobs=0; iobs<nobs; iobs++) {
    predbio(iobs) = 0.5*(biomass(index_years(iobs))+biomass(index_years(iobs) + 1));
    logq += log(index(iobs)/predbio(iobs));   
  }
  logq = logq/nobs;
  index_pred = exp(logq)*predbio; //*exp(logq);
  
  //predicted catches
  vector<Type> catch_pred(nyrs);
  for (int iyr=0;iyr<nyrs;iyr++)
    catch_pred(iyr) = biomass(iyr)*exp(upar(iyr))/(1.0+exp(upar(iyr)));
  //std::cout << catch_pred << std::endl;
  
  // Likelihood
   //CATCHES
   neglogL -= sum(dnorm(log(catches),log(catch_pred),exp(logSigmaC), true));
   //SURVEY
   neglogL -= sum(dnorm(log(index),log(index_pred),exp(logSigmaI), true));
  
  ADREPORT(logq);
  ADREPORT(biomass);
  ADREPORT(logSigmaI); //not estimated analytically
  ADREPORT(m); //estimated in the model
  ADREPORT(catch_pred); //predicted catch
  ADREPORT(index_pred); //predicted index
  ADREPORT(logSigmaC);
  //ADREPORT(upar);
  
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

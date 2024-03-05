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
  
  
  //PROJECTIONS
  //SIMULATE{
    //storage vectors
    vector<Type> storageF1(30);
    vector<Type> storageF2(30);
    vector<Type> storageF3(30);
    vector<Type> proj_catchF1(30);
    vector<Type> proj_catchF2(30);
    vector<Type> proj_catchF3(30);
    //exploitation
    Type F1=0.1;
    Type F2=0.2;
    Type F3=0.3;
    
    SIMULATE{
    //bridge year values equal to terminal years of dataset
    storageF1(0) = schaefer(biomass(nyrs), exp(logr), exp(logK), m, F1);
    storageF2(0) = schaefer(biomass(nyrs), exp(logr), exp(logK), m, F2);
    storageF3(0) = schaefer(biomass(nyrs), exp(logr), exp(logK), m, F3);
    proj_catchF1(0)= biomass(nyrs)*exp(F1)/(1.0+exp(F1)); //exploitation rate is fixed
    proj_catchF2(0)= biomass(nyrs)*exp(F2)/(1.0+exp(F2));
    proj_catchF3(0)= biomass(nyrs)*exp(F3)/(1.0+exp(F3));
    //biomass projection
    for (int i = 1; i<30; i++){
      storageF1(i) = schaefer(storageF1(i-1), exp(logr), exp(logK), m, F1);
      storageF2(i) = schaefer(storageF2(i-1), exp(logr), exp(logK), m, F2);
      storageF3(i) = schaefer(storageF3(i-1), exp(logr), exp(logK), m, F3); 
    //catch projection
    proj_catchF1(i)= storageF1(i)*exp(F1)/(1.0+exp(F1)); //exploitation rate is fixed
    proj_catchF2(i)= storageF2(i)*exp(F2)/(1.0+exp(F2));
    proj_catchF3(i)= storageF3(i)*exp(F3)/(1.0+exp(F3));
    }
  }
  //Report projections
  REPORT(storageF1);
  REPORT(storageF2);
  REPORT(storageF3);
  REPORT(proj_catchF1);
  REPORT(proj_catchF2);
  REPORT(proj_catchF3);
  
  //likelihood
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

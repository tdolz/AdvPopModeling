#include <TMB.hpp>

template<class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* DATA */
  DATA_SCALAR(M);// 
  DATA_VECTOR(catches);
  DATA_VECTOR(index);
  DATA_VECTOR(selectivity); //selectivity at age
  DATA_VECTOR(WAA); //weights at age
  DATA_VECTOR(maturity); //maturity
  DATA_IVECTOR(index_years); //integer vector
  int nyrs = catches.size(); //number of unique observations of catch
  int nobs = index.size(); // number of unique observations of the survey biomass. 
  int nages =10; //number of ages - 
  
  ///////////////* PARAMETERS *//////////////////
  //process model - pop dy params
  //PARAMETER(R1); //recruitment in the first year of the series - not doing this at first. 
  PARAMETER(log_R0); //unfished recruitment 
  PARAMETER(logit_h); //steepness parameter for BH recruitment
  //PARAMETER_VECTOR(eps); //recruitment process noise
  //PARAMETER(logit_m);//Natural Mortality - WE WILL MAP THIS = 0.4!
  PARAMETER_VECTOR(Fpar); //fishing mortality 
  //observation model params
  PARAMETER(logSigmaC);
  //PARAMETER(logSigmaI); //estimating analytically.
  //PARAMETER(logq); //estimating analytically
  

  Type neglogL = 0.;
  
  ///TRANSFORM PARAMETERS//
  //Type M = 3*exp(logit_m)/(1.0+exp(logit_m))+1.0001; // 0<m<4
  //Type F = 2/(1.0+exp(logit_F))+1.0001; //0<F<2 - not sure if this needs to be here. 
  Type h = exp(logit_h)/(1+exp(logit_h));
  Type F = exp(Fpar);
 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //* POPULATION DYNAMICS - THE PROCESS MODEL *//
 
 
 //Calculate S0
 Type S0 =0.0;
 Type SPR =0.0;
 Year0(0) = R0;//initial recruitment
 vector <Type> Year0(nages);
 for (int jages=1; jages<nages; jages++){ //ages 1-9 aka 2-10
   //ages 1-8 aka 2-9
   if(jages<=(nages-2)) Year0(jages)=Year0(jages-1)*exp(-M);
   //age 9 a, somka 10
   if(jages==(nages-1)) Year0(jages)=(Year0(jages-1)*exp(-M))/(1-exp(-M));
   SPR += Year0(jages)*WAA(jages)*maturity(jages);
 }
 S0 =R0*SPR;
 
 /////////////////////NUMBERS AT AGE MATRIX//////////////////////////////////
 matrix <Type> NAA(nyrs,nages); //plus 1 because we start at zero
 matrix <Type> F(nyrs+1,nages); // fishing mortality at age each year. 
 vector <Type> SSB(nyrs);
 vector<Type> catch_pred(nyrs);
 vector<Type> biomass(nyrs);

for(iyrs=0;iyrs<nyrs;iyrs++){ //rows starting at year 0 - catch series -1 (aka end of catch series);
  for (int jages=0; jages<nages; jages++){  //columns  = ages 0-9 aka 1-10
    
    ///////FIRST YEAR of the NAA///////
    if (iyrs == 0){
      //age 0 aka 1, time step 0, aka 1
      if(jages==0) NAA[iyrs][jages]=R0;
      //ages 1-8 aka 2-9
      if(jages<=(nages-2)) NAA[iyrs][jages]=NAA[iyrs][jages-1]*exp(-M);
      //age 9 aka 10 - plus group
      if(jages==(nages-1)) NAA[iyrs][nages]=NAA[iyrs][jages-1]*exp(-M))/(1-exp(-M));
      //calculate SSB and catch in the same loop. 
      SSB(iyrs)+=NAA[iyrs][jages]*WAA(jages)*maturity(jages); //maturity here because only spawners.
      biomass(iyrs)+=NAA[iyrs][jages]*WAA(jages)/1000; //metric tons. maturity not used here because all biomass
      //catch = the ratio of fishing mortality to all mortality * function of all mortality applied to biomass
      catch_pred(iyrs)+=(selectivity(jages)*F(iyrs))/(selectivity(jages)*F(iyrs))+M)*(1-exp(-selectivity(jages)*F(iyrs))+M)))*biomass(iyrs);
    }
  
  ///////OTHER YEARS OF THE NAA ////////
      //RECRUITMENT//
      if(jages==0){ 
        //BEVERTON-HOLT
        NAA[iyrs][jages] = (4*h*R0*SSB(iyrs-1))(S0*(1-h)+SSB(iyrs-1)*(5*h-1));//*exp(eps(iyrs));
      } 
      //POPULATION PROCESS AGES 2-9 - the diagonal
      if (jages>=1&&jages<nages-2){
      NAA[iyrs][jages]  = NAA[iyrs-1][jages-1] *exp(-(F(iyrs-1)*selectivity(jages-1)+M));
    }
      //POPULATION PROCESS 10 PLUS GROUP 
      if (jages==nages-1){
     NAA[iyrs][jages]= NAA[iyrs-1][jages-1] *exp(-(F(iyrs-1)*selectivity(jages-1)+M))+ NAA[iyrs-1][jages]*exp(-(F(iyrs-1)*selectivity(jages)+M)); 
    }
    //Calculate catch and biomass here//
    SSB(iyrs)+=NAA[iyrs][jages]*WAA(jages)*maturity(jages); //maturity here because only spawners.
    biomass(iyrs)+=NAA[iyrs][jages]*WAA(jages)/1000; //metric tons. maturity not used here because all biomass
    //catch = the ratio of fishing mortality to all mortality * function of all mortality applied to biomass
    catch_pred(iyrs)+=(selectivity(jages)*F(iyrs))/(selectivity(jages)*F(iyrs))+M)*(1-exp(-selectivity(jages)*F(iyrs))+M)))*biomass(iyrs);
  }

  
  /////PROJECT FORWARD BEYOND CATCH SERIES///////
  Type biomNext =0.0;
  //I really don't understand why this is necessary?
  for (int jages=0; jages<nages; jages++){
    //recruits
    if(jages==0) biomNext(iyrs)=(4*h*R0*SSB(iyrs-1))(S0*(1-h)+SSB(iyrs-1)*(5*h-1))*WAA(jages)/1000;
    //ages 1-8 aka 2-9
      if(jages>0 && jages<=(nages-2))  biomNext(iyrs)+=NAA[iyrs-1][jages-1]*exp(-M-selectivity(jages-1)*F(iyrs-1))*WAA(jages)/1000;
      //age 9 aka 10
      if(jages==(nages-1)) biomNext(iyrs)+=(NAA[iyrs-1][jages-1]*exp(-M-selectivity(jages-1)*F(iyrs-1))*WAA(jages)+NAA[iyrs-1][jages]*exp(-M-selectivity(jages)*M(iyrs-1))*WAA(jages))/1000;
  }
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////  
 *// OBSERVATION MODEL //*
   
//SURVEY
  vector<Type> index_pred(nobs); //store predicted survey biomass (different from true biomass)
  //Implementing the analytical estimator for q
  Type logq = 0.;
  vector<Type> predbio(nobs); //store predicted survey biomass (different from true biomass)
  for (int iobs=0; iobs<nobs; iobs++) {
    predbio(iobs) = 
    logq += log(index(iobs)/predbio(iobs));   
  }
  logq = logq/nobs;
  index_pred = exp(logq)*predbio; //*exp(logq);
  //Implementing the analytical estimator for logSigmaI
  Type logSigI = 0;
  Type resid; //
  for (int iobs=0; iobs<nobs; iobs++) {
    resid += square(log(index(iobs))-log(exp(logq)*biomass(index_years(iobs))));
  }
  logSigI = sqrt(resid/nobs);
  
  
  ///////////////////////////// LIKELIHOOD/////////////////////////////////////////////////////
 
   //CATCHES
   neglogL -= sum(dnorm(log(catches),log(catch_pred),exp(logSigmaC), true));
   //SURVEY
   neglogL -= sum(dnorm(log(index),log(index_pred),logSigI, true));
 
 //////////////////////////////REPORT///////////////////////////////////////////////////////
  
  ADREPORT(logq);//analytically estimated
  ADREPORT(logSigI); //analytically estimated
  ADREPORT(biomass);//predicted biomass
  ADREPORT(M); //estimated in the model
  ADREPORT(catch_pred); //predicted catch
  ADREPORT(index_pred); //predicted index
  
  return neglogL;

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

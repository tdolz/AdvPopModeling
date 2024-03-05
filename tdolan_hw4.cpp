#include <TMB.hpp>
#include <iostream>
#include <Eigen/Dense>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data */
  DATA_MATRIX(age_comp); //age proportions as a matrix
  DATA_VECTOR(catches); //catch starting in 1963
  DATA_VECTOR(index); //survey index starting in 1963
  int nyrs = catches.size(); //catches should be the same length as index since 1963 start. 
  DATA_VECTOR(weight); //weights at age for 10 age classes
  DATA_VECTOR(maturity); //proportion mature for 10 age classes
  int nages = maturity.size(); //10
  int nage_prop = 5; //5 age classes in the age_comp data
  
  /* PARAMETERS */
  PARAMETER(logR0); //unfished recruitment
  PARAMETER(logfracR1); //proportion relating r1 to r0
  PARAMETER(logSigmaI); //survey observation model variance
  PARAMETER_VECTOR(logF); //fishing mortality
  //new params
  PARAMETER(logSigmaR); //std dev annual recr variation. (process error variance)
  PARAMETER_VECTOR(EtaT); //mean recruitment
  //mapped params
  PARAMETER(logith);
  PARAMETER(logSigmaC); //fishery observation model variance
  PARAMETER(logM);
  //selectivity params
  PARAMETER(logA50s);
  PARAMETER(logA95s);
  PARAMETER(logA50f);
  PARAMETER(logA95f);

  Type Rzero = exp(logR0); //we are estimating it in log space but converting back to normal
  Type R1 = Rzero*exp(logfracR1); //first year recruiment is some fraction of virgin recruitment
  Type h = 0.2 + 0.8*exp(logith)/(1.0+exp(logith));
  Type M = exp(logM);
  Type SigmaC = exp(logSigmaC);
  Type SigmaI = exp(logSigmaI);
  Type SigmaR =exp(logSigmaR);
  Type A50s =exp(logA50s);
  Type A95s =exp(logA95s);
  Type A50f =exp(logA50f);
  Type A95f =exp(logA95f);
  
  using namespace density;
  
  Type neglogL = 0.;
  
  ////////////////////////////// Selectivity //////////////////////////////////////////
  vector<Type> selex_s(nages); //survey selectivity
  vector<Type> selex_f(nages); //fishery selectivity 
  selex_s.setZero();
  selex_f.setZero();
  for (int age=0;age<nages;age++) { 
    selex_s(age)= 1/(1+exp(-log(19)*(age-logA50s)/logA95s));
    selex_f(age)= 1/(1+exp(-log(19)*(age-logA50f)/logA95f));
  }
  
  /////////////////////// F & Z calculations/////////////////////////////////////////
  matrix<Type> F(nyrs,nages); //fishing mortality
  matrix<Type> Z(nyrs,nages); //total mortality
  F.setZero();
  Z.setZero();
  //Loop to iterate over age as well. 
  for (int iyr=0; iyr<nyrs; iyr++){
    for (int age=0;age<nages;age++) {
    F(iyr,age) = exp(logF(iyr))*selex_f(age); //populate the F matrix: age influenced by selectivity, 
                                        //year from the estimated param vector. 
    Z(iyr,age) = F(iyr,age)+M;
    }}
  
  ////////////////////////////////////* POPULATION DYNAMICS *///////////////////////////
  vector<Type> Bio(nyrs+1); //vector for storing predicted total biomass
  vector<Type> SSB(nyrs+1); //vector for storing biomass scaled by proportion mature
  vector<Type> R(nyrs+1);  //the first age class for all years. 
  matrix<Type> N(nyrs+1,nages); //numbers at age matrix
  N.setZero();
  SSB.setZero();
  R.setZero();
  Bio.setZero();
  
  //////////////////   Initial numbers at age Matrix //////////////////////////////////////
  //   **unfished** Spawning biomass calculation
  //   fill the first year such that each age class is deprecated by natural mortality only. 
  N(0,0) = R1;
  for (int age=1;age<nages-1;age++) { //from age 1 to age 9
    N(0,age) = N(0,age-1)*exp(-M); //year 0, age i = previous age deprecated by M. 
  }
  N(0,nages-1) = N(0,nages-2)/(1-exp(-M)); 
 
 
 ////////////////// SSB0 //////////////////////////////////////////////////////////////////// 
  //unfished spawning stock biomass vector --> we are trying to calculate the SSB0 value
  // This is repeating the process in the loop above, but this time only to calculate SSB0
  Type SSB0 = 0.0;
  vector<Type> nb0(nages);
  nb0(0) = Rzero; //Rzero is an estimated parameter
  for (int age=1;age<nages-1;age++) //for ages 1-10
    nb0(age) = nb0(age-1)*exp(-M); //year 0, age i = previous age deprecated by M.
  
  nb0(nages-1) = nb0(nages-2)/(1-exp(-M)); //changed to nages-2
  for (int age=1;age<nages;age++) 
    SSB0 += 0.5*nb0(age)*weight(age)*maturity(age); //SSB0 is the sum of all mature age classes at year0. 
  //0.5 because it's female SSB
 
   
  /////////////////// 1st year stock quantities////////////////////////////////////////////
  //for the time loop we need separate values for the first year. 
  for (int age=0; age<nages; age++) {
    Bio(0) += N(0,age)*weight(age)*exp(-0.5*Z(0,age)); //predicted biomass
    //selectivity at age moderates fishing mortality, which is contained within z. 
    //number of fish in year 0 from the NAA matrix*WAA to get biomass.
    //deprecated by total mortality (Z) But why by -0.5? because it's mid-year. 
    SSB(0) += 0.5*N(0,age)*maturity(age)*weight(age);  
    // SSB is the NAA*maturity*weight to get the biomass of mature females (0.5 because females?)
  }
  // first year recruitment
  R(0) = N(0,0); //R0 = N(0,0) = log(nb0(0))
  
  
  //////////////// BEGIN LOOP OVER TIME TO FILL THE REST OF YEARS /////////////////////////////////////
  for (int iyr=1; iyr<=nyrs; iyr++) {
    
  /////    Numbers at age update 
  for (int age=1; age<nages-1; age++) //for ages 1-9?
    // the N this year is the N from the previous age class last year, deprecated by total mortality.
    N(iyr,age) = N(iyr-1,age-1)*exp(-Z(iyr-1,age-1)); 
  //then for the plus group
  N(iyr,nages-1) += N(iyr-1,nages-2)*exp(-Z(iyr-1,nages-2));
  
  //Biomass predictions
  Bio(iyr) = 0.;
  SSB(iyr) = 0.;
  for (int age=0; age<nages; age++) {
    // for all years other than the last year
    // biomass is the sum of numbers at all age classes * weight
    // deprecated by total mortality, which is changing over time because fishing mortality is 
    // which is multiplied by -0.5, WHY??
    // which is moderated by selectivity - why is the selectivity not with the mortality??
   if (iyr != nyrs) Bio(iyr) += N(iyr,age)*weight(age)*exp(-0.5*Z(iyr,age));
   if (age != 0) SSB(iyr) += 0.5*N(iyr,age)*maturity(age)*weight(age);  
  }
  
  ////// RECRUITMENT OVER TIME /////////////
  //moved this section from above biomass predictions
  // Recruitment for this year given annual stock size from the recruitment equation
  // h is fixed at 0.6
  //recruitment is the BH + process error variance that is annual deviations & annual variance. 
 //ssb needs to be lagged to produce recruitment in the following year. 
  R(iyr) = 4.*h*Rzero*SSB(iyr-1)/(SSB0*(1.-h)+SSB(iyr-1)*(5*h-1.))*exp(EtaT(iyr-1)-0.5*logSigmaR*logSigmaR);
  N(iyr,0) = R(iyr);
   
  //////////// End Loop over time ////////////////
  }
  
  
  /////////////////////// PRED CATCH AND SURVEY ////////////////////////////////////
  
  ///// PRED SURVEY BIO AND CATCH ////
  //survey index pred
  vector<Type> index_pred(nyrs);
    Type q = exp(log(index/Bio(nyrs)).sum()/nyrs); //biomass (kg)
    //index_pred = q*Bio(index_years); 
    for (int iyr=0;iyr<nyrs;iyr++) {
      index_pred(iyr) = 0.;
      for (int age=0; age<nages; age++)
        //not sure if i have correctly implemented q
        index_pred(iyr) += q*N(iyr,age)*weight(iyr,age)*selex_s(age)*exp(-0.5*Z(iyr,age));
    }
    
  //catch index pred
  vector<Type> catch_pred(nyrs);
  for (int iyr=0;iyr<nyrs;iyr++) {
    catch_pred(iyr) = 0.;
    for (int age=0; age<nages; age++)
    catch_pred(iyr) += N(iyr,age)*weight(iyr,age)*(selex_f(age)*F(iyr,age)/Z(iyr,age))*(1.0-exp(-1.*Z(iyr,age))); //catch (mt)
  }
 
  //fishery age proportion prediction - following format sent to Maya from Gavin
  matrix<Type> CAA_pred(nyrs,nage_prop); //survey at age predictions
  CAA_pred.setZero();
  for (int iyr=0;iyr<nyrs;iyr++) {
    for (int age=0; age<nages; age++) {
      int jage = age;
      if (age>=nage_prop) jage = nage_prop-1;  //ages 6-10 correspond to age 5 of age data
      CAA_pred(iyr,jage) += ((N(iyr,age)*selex_f(age)*F(iyr,age))/Z(iyr,age))*(1.0-exp(-1.*Z(iyr,age))); 
    }
    CAA_pred.row(iyr) = (CAA_pred.row(iyr)/CAA_pred.row(iyr).sum());
  }  

  //Survey age proportion prediction - following format sent to Maya from Gavin
  matrix<Type> IAA_pred(nyrs,nage_prop); //survey at age predictions
  IAA_pred.setZero();
  for (int iyr=0;iyr<nyrs;iyr++) {
    for (int age=0; age<nages; age++) {
      int jage = age;
      if (age>=nage_prop) jage = nage_prop-1;  //ages 6-10 correspond to age 5 of age data
      IAA_pred(iyr,jage) += selex_s(age)*N(iyr,age)*exp(-0.5*Z(iyr,age));
    }
    IAA_pred.row(iyr) = (IAA_pred.row(iyr)/IAA_pred.row(iyr).sum()); 
  }  
  
  //////// renormalize the predicted proportions which must sum to 1.
  //add a small constant then divide by the rowsum to renormalize to 1.  
  for (int iyr=0;iyr<nyrs;iyr++) {
    for (int age=0; age<nage_prop; age++) {
      IAA_pred(iyr,age) = IAA_pred(iyr,age)+0.0001;
      CAA_pred(iyr,age) = CAA_pred(iyr,age)+0.0001;
    }
    IAA_pred.row(iyr) = (IAA_pred.row(iyr)/IAA_pred.row(iyr).sum()); 
    CAA_pred.row(iyr) = (CAA_pred.row(iyr)/CAA_pred.row(iyr).sum()); 
  }  
  

  /////////////////////// Contribution to the Likelihood function ////////////////////////////////////

  
  //catches and survey aggregate time series
  neglogL -= sum(dnorm(log(catches),log(catch_pred),SigmaC, true));
  neglogL -= sum(dnorm(log(index),log(index_pred),SigmaI, true));
  
  //AGE proportions for catches and survey
  for (int iyr=0;iyr<nyrs;iyr++) {
    neglogL -= dmultinom(age_comp.row(iyr),CAA_pred.row(iyr),true);
    neglogL -= dmultinom(age_comp.row(iyr),IAA_pred.row(iyr), true);
    } 
  
  // RECRUITMENT DEVIATIONS
  for (int iyr=0;iyr<nyrs;iyr++) {
  //neglogL -= sum(dnorm(R(iyr),EtaT(iyr),exp(logSigmaR), true));
  neglogL -= sum(dnorm(EtaT(iyr),0,exp(logSigmaR), true));
  }

  
  //////////////////////////////// REPORT SECTION ////////////////////////////////////////  
  REPORT(F); //time series of fishing mortality. 
  REPORT(N);
  REPORT(Z);
  REPORT(catch_pred);
  REPORT(index_pred);
    
  //ADREPORT(h);
  ADREPORT(logSigmaR);//process error variance - std. dev of recruitment variation. 
  ADREPORT(Rzero); //unfished recruitment
  ADREPORT(logfracR1); //initial year of recruitment relative to unfished recruitment. 
  ADREPORT(SigmaI);
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

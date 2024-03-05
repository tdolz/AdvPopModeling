#include <TMB.hpp>

//function to implement the PT production model for a single time step
template <class Type> Type schaefer(Type b1, Type r, Type K, Type m, Type u_i){
    Type u_use = exp(u_i)/(1.0+exp(u_i)); //logit transforming exploitation rate Ui
    Type b2 = b1 + b1*r*(1.0-pow(b1/K,m-1.0)) - u_use*b1; //literally just the shaefer model equation.
    return b2;
}

template<class Type> Type square(Type x){return x*x;}

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
  PARAMETER(logit_m);
  PARAMETER(logSigmaC);
  //PARAMETER(logSigmaI); //estimating analytically.
  //PARAMETER(logq);
  PARAMETER_VECTOR(upar); //exploitation rate in each year. 

  Type neglogL = 0.;
  
  //transform some params
  Type m = 3*exp(logit_m)/(1.0+exp(logit_m))+1.0001; // 0<m<4
  
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
  
  //Implementing the analytical estimator for q
  Type logq = 0.;
  vector<Type> predbio(nobs); //store predicted survey biomass (different from true biomass)
  for (int iobs=0; iobs<nobs; iobs++) {
    predbio(iobs) = 0.5*(biomass(index_years(iobs))+biomass(index_years(iobs) + 1));
    logq += log(index(iobs)/predbio(iobs));   
  }
  logq = logq/nobs;
  index_pred = exp(logq)*predbio; //*exp(logq);
  
  //Implementing the analytical estimator for logSigmaI
  Type logSigI = 0;
  Type resid =0; //
  for (int iobs=0; iobs<nobs; iobs++) {
    resid += square(log(index(iobs))-log(exp(logq)*predbio(iobs)));
  }
  logSigI = sqrt(resid/nobs);
  
  //predicted catches
  vector<Type> catch_pred(nyrs);
  for (int iyr=0;iyr<nyrs;iyr++)
    catch_pred(iyr) = biomass(iyr)*exp(upar(iyr))/(1.0+exp(upar(iyr)));
  //std::cout << catch_pred << std::endl;
  
  // Likelihood
   //CATCHES
   neglogL -= sum(dnorm(log(catches),log(catch_pred),exp(logSigmaC), true));
   //SURVEY
   neglogL -= sum(dnorm(log(index),log(index_pred),logSigI, true));

  
  ADREPORT(logq);//analytically estimated
  ADREPORT(logSigI); //analytically estimated
  ADREPORT(biomass);//predicted biomass
  ADREPORT(m); //estimated in the model
  ADREPORT(catch_pred); //predicted catch
  ADREPORT(index_pred); //predicted index

  //When I put the simulate section here - i get new predictions each time,
  //but predicted catches are orders of magnitude too high. 
  // Gavin says it needs to be inside the objective function. 
  SIMULATE{
    catches = exp(rnorm(log(catch_pred), exp(logSigmaC)));//simulate catch time series
    REPORT(catches);
    index=exp(rnorm(log(index_pred),logSigI));//simulate survey time series
    REPORT(index);
   } 

  return neglogL;
 
 //When I put the simulate section here predicted catches are reasonable
 //but each run gives me the exact same output for predicted catch and index.
 //SIMULATE{
  // catches = exp(rnorm(log(catch_pred), exp(logSigmaC)));//simulate catch time series
  // REPORT(catches);
  // index=exp(rnorm(log(index_pred),logSigI));//simulate survey time series
  // REPORT(index);
 //} 
  
}

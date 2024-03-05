#include <TMB.hpp>
#include <iostream>
template<class Type> Type square(Type x){return x*x;}
template<class Type>Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  int n = y.size(); // add a minus 1?
  PARAMETER(logit_alpha); //constrained to range from 0 to 2. 
  PARAMETER(log_beta);
  PARAMETER(log_sigma_obs);
  PARAMETER(logN0);
  PARAMETER_VECTOR(u); //unobserved state, actual number of juveniles. 
  
  using namespace density;
  
  //back-transform parameters
  Type beta = exp(log_beta);
  Type alpha =2.0*exp(logit_alpha)/(1.0+exp(logit_alpha));
  Type N0 = exp(logN0);

  //results that we want to populate
  Type nll = 0.0;
  Type lambda = 0.0;
  Type lmean = 0.0;
  
  //Process model
  nll-=dpois(u[0],alpha*N0*exp(-beta*N0),true);//define likelihood for first value in the state variable.
  for(int i=1;i<n;i++){ //year 2 onward
   lambda=alpha*u[i-1]*exp(-beta*u[i-1]);
   nll-=dpois(u[i],lambda,true);
  }
  //Observation model
  for(int i=0;i<n;i++){
    lmean=log(u[i])-square(log_sigma_obs)/2;
    nll-=dnorm(log(y[i]),lmean,log_sigma_obs,true);
  }
  ADREPORT(beta);
  ADREPORT(alpha);
  ADREPORT(N0);
  REPORT(beta);
  REPORT(alpha);
  REPORT(N0);
  REPORT(log_sigma_obs);
  return nll;
}

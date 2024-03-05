#include <TMB.hpp>
#include <vector>
#include <cmath>
#include <iostream>

template<class Type> Type square(Type x){return x*x;}
template<class Type>
Type objective_function<Type>::operator() ()
{
  //DATA//
  DATA_VECTOR(pups_data); //pups
  DATA_VECTOR(nonpups_data); //nonpups
  DATA_VECTOR(year); //time step
  int n = year.size();
  
  //PARAMETERS//
  PARAMETER_VECTOR(pups); //unobserved state variable, actual number of pups
  PARAMETER_VECTOR(nonpups); //unobserved state variable, actual number of pups
  PARAMETER(f); //fecundity
  PARAMETER(phi_np); //nonpup survival
  PARAMETER(sigma2_p); //process error on obs model pups
  PARAMETER(sigma2_np); //process error on obs model nonpups
  PARAMETER(tau2); //variance of the process error
   
//DECLARE TEMP VARIABLES
//everything has to be either data or a parameter or a temp variable.
 using namespace density;

vector<Type> nll =0;
std::vector<double> ypups(79, log(1000)); //sample mean pupcount in log space. 
std::vector<double> ynp(79, log(1300)); //sample mean adultcount in log space. 
Type psi_pups = 0.0; 
Type psi_np = 0.0;
Type eps_p = 1000;
Type eps_np = 2000; 
Type q = 0.3;
Type phi_p = 0.6;

//backtransform_params
Type ln_psi_pups = log(psi_pups);
Type ln_psi_np = log(psi_np);
Type ln_eps_p = log(eps_p);
Type ln_eps_np = log(eps_np);

// PROCESS MODEL  
// there is a process model output for every year, not for every observation. 
for (int i=0; i<n;i++) { //looping over years
  pups(i+1)=nonpups(i)*f + ln_psi_pups;
  nonpups(i+1)=pups(i)*phi_p + nonpups(i)*phi_np + ln_psi_np;
  nll-=dnorm(psi_pups,0,tau2,true);
  nll-=dnorm(psi_np,0,tau2,true);
  nll-=sum(dnorm(log(pups(i+1)),log(f*nonpups(i)),psi_pups,true));
  nll-=sum(dnorm(log(nonpups(i+1)),log(phi_p*pups(i))+log(phi_np*nonpups(i)),psi_np,true));
}

// OBSERVATION MODEL 
for (int i=0; i<n;i++){
  ypups(i) = pups(i)+ln_eps_p; //observations of pups are the state variable pups plus error. 
  nll-=sum(dnorm(log(pups_data(i)),log(ypups),eps_p,true)); //dnorm(conditioned, expectation,variance)
  nll-=dnorm(eps_p,0,sigma2_p,true); //observation error
  ynp(i) = q*nonpups(i)+ln_eps_np; //observations of nonpups are the state variable nonpups plus error. 
  nll-=sum(dnorm(log(nonpups_data(i)),log(ynp(i)*q),eps_np,true)); //dnorm(conditioned, expectation,variance)
  nll-=dnorm(eps_np,0,sigma2_p,true); //observation error
}


//OUTPUTS  
// put ADREPORT and REPORT here.
ADREPORT(phi_np);
ADREPORT(sigma2_p); 
ADREPORT(sigma2_np); 
ADREPORT(tau2);
  return nll;
}
  



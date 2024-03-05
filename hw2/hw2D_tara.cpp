#include <TMB.hpp>
#include <vector>
#include <cmath>
#include <iostream>

//PROCESS ERROR ONLY MODEL//

template<class Type> Type square(Type x){return x*x;}
template<class Type>
Type objective_function<Type>::operator() ()
{
  //DATA//
  DATA_VECTOR(pups_data); //pups
  DATA_VECTOR(nonpups_data); //nonpups
  DATA_VECTOR(year); //time step
  DATA_VECTOR(pup_years); //years in which pup counts were made
  DATA_VECTOR(adult_years); //years in which adult counts were made
  int year = year.size(); //list years including replicate counts (??)
  int pup_obs = pup_years; // year associated with every pup count, even with replicates. 
  int adult_obs = adult_years; // year associated with every adult count, even with replicates.
  
  //PARAMETERS//
  PARAMETER_VECTOR(pups); //unobserved state variable, actual number of pups
  PARAMETER_VECTOR(nonpups); //unobserved state variable, actual number of pups
  PARAMETER(f); //fecundity
  PARAMETER(phi_p); //pup survival - could use a constant from the pilot experiment at 0.6
  PARAMETER(phi_np); //nonpup survival
  PARAMETER(sigma2_p); //process error on obs model pups
  PARAMETER(sigma2_np); //process error on obs model nonpups
  PARAMETER(q); //observability of nonpups - could use a constant from the pilot experiment at 0.3
  PARAMETER(tau2); //variance of the process error
   
//DECLARE TEMP VARIABLES
//everything has to be either data or a parameter or a temp variable.
using namespace density

Type nll =0.0;
std::vector<double> ypups(79, log(1000)); //sample mean pupcount in log space. 
std::vector<double> ynp(79, log(1300)); //sample mean adultcount in log space. 
Type psi_pups = 0.0;
Type psi_np = 0.0;

//backtransform_params
Type ln_psi_pups = log(psi_pups);
Type ln_psi_np = log(psi_np);

// PROCESS MODEL  
// there is a process model output for every year, not for every observation. 
for (int i=0; i<uniqueYear;i++) { //looping over unique years.
  pups(i+1)= nonpups(i)*f + ln_psi_pups;
  nonpups(i+1)= pups(i)*phi_p + nonpups(i)*phi_np + ln_psi_np;
  nll-=dnorm(psi_pups,0,tau2,true);
  nll-=dnorm(psi_np,0,tau2,true);
  nll-=sum(dnorm(log(pups(i+1)),log(f*nonpups(i)),psi_pups,true));
  nll-=sum(dnorm(log(nonpups(i+1)),log(phi_p*pups(i)+phi_np*nonpups(i)),psi_np,true));
}

// OBSERVATION MODEL - I decided to split this into two for loops, so that we could loop over correct n observations for each.
// observation model for pups
for(int i=0;i<pup_obs;i++){ //looping over observations of pups, even though there are multiple observations of year.
  ypups(i) = pups(i) //observations of pups are the state variable pups plus error. 
  nll-=sum(dnorm(log(pups_data(i)),log(ypups(i)),0,true)); //dnorm(conditioned, expectation,variance)
}

// observation model for nonpups
for(int i=0;i<adult_obs;i++){ //looping over number of observations of nonpups, even though there are multiple observations of year.
  ynp(i) = q*nonpups(i); //observations of nonpups are the state variable nonpups plus error. 
  nll-=sum(dnorm(log(nonpups_data(i)),ynp(i)*q,0,true)); //dnorm(conditioned, expectation,variance)
}


//likelihood penalties
  nll-=dnorm(phi_p,0.6,0.05,true); // penalty on pup survival(phi_p) mean=0.6, sd=0.05
  nll-=dnorm(log(q),log(0.3),log(0.1),true);//penalty on *fraction* of nonpups hauled out. mean=0.3, sd=0.1

//OUTPUTS  
// put ADREPORT and REPORT here.
  ADREPORT(pups);
  ADREPORT(nonpups);
  ADREPORT(f);
  ADREPORT(phi_p);
  ADREPORT(phi_np);
  ADREPORT(sigma2_p); 
  ADREPORT(sigma2_np); 
  ADREPORT(q); 
  ADREPORT(tau2);
  return nll;
}
  


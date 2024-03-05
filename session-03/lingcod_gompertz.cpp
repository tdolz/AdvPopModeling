// State-space Gompertz model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_VECTOR(obs);
  
  // parameters:
  PARAMETER(logr); // population growth rate parameter
  PARAMETER(logit_beta); // density dependence parameter
  PARAMETER(log_sigma_proc); // log(process SD)
  PARAMETER(log_sigma_obs); // log(observation SD)
  PARAMETER_VECTOR(u); // unobserved state vector
  
  // procedures: (transformed parameters)
  Type sigma_proc = exp(log_sigma_proc);
  Type sigma_obs = exp(log_sigma_obs);
  Type beta = -2.0 + 2.0*exp(logit_beta)/(1.0+exp(logit_beta));
  Type r = exp(logr);
  
  // reports on transformed parameters:
  ADREPORT(sigma_proc)
  ADREPORT(sigma_obs)
  ADREPORT(beta)
  ADREPORT(r)
    
    int n = obs.size(); // get time series length
    int n_u = n-1;
    
  Type nll = 0.0; // initialize negative log likelihood
  
  // process model:
  Type m = -r/beta;   //initial state
  nll -= dnorm(u[0], m, sigma_proc, true);
  for(int i = 1; i < n_u; i++){
    Type m = r + (1+beta) * u[i-1]; // Gompertz
    nll -= dnorm(u[i], m, sigma_proc, true);
  }
  
  // observation model:
  vector<Type> b(n);   //storage vector
  b[0] = -r/beta;
  for(int i = 1; i < n; i++){
   b[i] = u[i-1];
  }
  nll -= sum(dnorm(obs, b, sigma_obs, true));
  
  ADREPORT(b)
  
  return nll;
}


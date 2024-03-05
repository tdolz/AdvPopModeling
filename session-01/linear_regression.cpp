#include <TMB.hpp> //header

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() () //objective function
{
  //data section
  DATA_VECTOR(x);    //data vector, passed in list object from R
  DATA_VECTOR(y);    //data vector, passed in list object from R
  int n = y.size();    //define an integer n that is the length of the y vector
  
  //parameter section
  PARAMETER(b0);    //model parameter to be estimated, real number
  PARAMETER(b1);    //model parameter to be estimated, real number
  PARAMETER(logSigma);    //model parameter to be estimated, real number
  //sigma is on the log scale to force it to be positive. 
  
  vector<Type> ypred(n);   //define object to store our model predictions
  // note this is a function of model parameters, 
  // therefore needs to be differntiable
  // we also specify the size of the vector
  
  Type neglogL = 0.0;    //Initialize the objective function to 0.
  // specify the model & objective function
  
  ypred = b0 + b1*x;    // model predictions
  
  //negative log likelihood. aka the objective function
  //dnorm(observations, mean, standard deviation, true=force to log likelihood)
  //recall, we are summing the log likelihood over all predictions. so over all preds,
  //the probability density (continuous dist) for each value y, where it's coming
  //from a normal distribution with a mean of yhat pred and a standard deviation of sigma
  //
  neglogL = -sum(dnorm(y, ypred, exp(logSigma), true));  
   
   //derived quantities
   SIMULATE {
     y = rnorm(ypred,exp(logSigma));
     REPORT(y); //get results back into R
   }
  
  //compute standard errors for derived quantities. Then use sdreport(model) on the R side. 
  Type sigma2 = square(exp(logSigma));
  ADREPORT(sigma2);  
  
  //adding priors
  //vector<Type> mu= log(0.2);
  //vector<Type> obs= Type(logSigma);
  //neglogL -= dnorm(obs,mu,0.1);
  //neglogL -= square((logSigma-log(0.2))/0.1);
  
  return neglogL;  //return the objective function
}

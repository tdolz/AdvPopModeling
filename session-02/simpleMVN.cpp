#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator()()
{
  //import data matrix X from R and extract its dimensions. 
  DATA_MATRIX(X);
  int nobs = X.cols();
  int n = X.rows();
  
  //declare the parameters- total of 7. 
  PARAMETER_VECTOR(mu); //vector of 3 means
  PARAMETER_VECTOR(logSigma); //vector of 3 standard deviations
  PARAMETER(rhopar); //correlation parameter

  // we need the standard deviations to be positive 
  // and the correlations to be between -1 and 1.
  // estimate standard deviation on the log scale. 
  vector<Type> sd = exp(logSigma);
  //shifted logistic transformation for the correlations. 
  Type rho = 2.0 / (1.0 + exp(-rhopar)) - 1.0;
  
  //create a matrix that is a mirror of the one in R and fill it in using <<
  matrix<Type> Sigma(n,n);
  Sigma.row(0) << sd[0]*sd[0], sd[0]*sd[1]*rho, sd[0]*sd[2]*rho;
  Sigma.row(1) << sd[1]*sd[0]*rho, sd[1]*sd[1], sd[1]*sd[2]*rho;
  Sigma.row(2) << sd[2]*sd[0]*rho, sd[2]*sd[1]*rho, sd[2]*sd[2];
  
  vector<Type> residual(n);
  
  Type nll = 0.0;
  
  //evaluate multivariate normal density
  //"using namespace" is like "library" in R. 
  //MVNORM function in package/namespace density
  //mean is always 0. We pass it the covariance matrix Sigma we constructed. 
  using namespace density;
  MVNORM_t<Type> neg_log_dmvnorm(Sigma);
  
  //loop over the observations to calculate the negative log likelihood. 
  for (int iobs=0;iobs<nobs;iobs++) {
    //for each row x, we define a residual vector that has a mean of 0
    //by subtracting mu from x values.
    //we must make it compatible with the row of mu by typecasting that row as a vector. 
    residual = vector<Type>(X.row(iobs))-mu;
    //evaluate multivariate normal density of the residuals. 
    //note we are returning the negative log density. 
    nll += neg_log_dmvnorm(residual);
    // another way to do this would be:
    // nll += MVNORM(Sigma)(residual); 
    // negative log density of the covariance matrix evaluated at the residuals.
    // you would not need lines 38 or 48 in that case. 
  }  
  
  REPORT(Sigma);
  REPORT(sd);
  REPORT(rho);
  
  return nll;
}

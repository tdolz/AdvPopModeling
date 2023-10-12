#include <TMB.hpp>
template<class Type> Type square(Type x){return x*x;}
template<class Type>Type objective_function<Type>::operator()()
{

  DATA_VECTOR(l1);
  DATA_IVECTOR(l2);
  DATA_VECTOR(dt);
  int n = l1.size();
  
  PARAMETER(K);
  PARAMETER(Linf_mu);
  PARAMETER(sigma2_Linf);
  PARAMETER(sigma2_resid);
  PARAMETER_VECTOR(A_i);
  PARAMETER(sigma2_A);
  PARAMETER(mu_A);
  
  Type nll = 0.0;
  
  //contribution to the A_i random effects
  nll -= sum(dnorm(A_i,0,sqrt(sigma2_A),true)); //true = log

  vector<Type>mu1(n);
  vector<Type>mu2(n);
  vector<Type>sigma2_l1(n);
  vector<Type>sigma2_l2(n);
  vector<Type>cov(n);
  vector<Type>rho(n);
  vector<Type>residuals(n);
  vector<Type>pairs(2);
  Type f1=0.0;
  Type f2=0.0;
  vector<Type>Age(n);
  matrix<Type>big_sigma(2,2);
  Type AIC=0.0;
  using namespace density;
  
  for (int i=0;i<n;i++){
   Age(i)=mu_A+A_i(i);
   f1=1-exp(-K*(Age(i)));
   f2=1-exp(-K*(Age(i)+dt(i)));
   mu1(i) = Linf_mu*f1;
   mu2(i) = Linf_mu*f2;
   sigma2_l1(i)=sigma2_resid+sigma2_Linf*square(f1);
   sigma2_l2(i)=sigma2_resid+sigma2_Linf*square(f2);
   cov(i)=sigma2_Linf*f1*f2;
   big_sigma.row(0) << sigma2_l1(i), cov(i);
   big_sigma.row(1) << cov(i), sigma2_l2(i);
   pairs(0)=l1(i)-mu1(i);
   pairs(1)=l2(i)-mu2(i);
   residuals(i)=MVNORM(big_sigma)(pairs);
 }
   nll+=sum(residuals);
  AIC = -2*(-1*nll)+2*6;
  ADREPORT(Age);
  REPORT(AIC);
    return nll;

}

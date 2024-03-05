#include <TMB.hpp>
#include <iostream>

template<class Type> Type square(Type x){return x*x;}
template <class Type> Type recruits(Type r0, Type h1, Type S,Type s0){
  Type rr = (4*h1*r0*S)/(s0*(1-h1)+S*(5*h1-1)); 
  return rr;
}
template<class Type>Type objective_function<Type>::operator()()
{
  DATA_VECTOR(catches);
  DATA_VECTOR(index);
  DATA_IVECTOR(index_years);
  int nyrs = catches.size();
  int niyrs = index.size();
  DATA_VECTOR(WAA);
  DATA_VECTOR(select);
  DATA_VECTOR(maturity);
  DATA_SCALAR(NM);
  DATA_INTEGER(nages);
  
  PARAMETER(log_R0);
  PARAMETER(log_sigma_catch);
  PARAMETER(logit_h);
  PARAMETER_VECTOR(u);//F parameters

  Type R0 = exp(log_R0);
  Type h = exp(logit_h)/(1+exp(logit_h));
  
  //Storage Arrays
  vector<Type> FM(nyrs);
  matrix<Type> N(nyrs,nages);
  vector<Type> SSB(nyrs);
  vector<Type> B(nyrs+1);
  vector<Type> pred_catch(nyrs);
  vector<Type> pred_index(niyrs);
  vector<Type> ntemp(nages);
  
  //Counters
  int y;
  int iy;
  int a;
  
  //Define real variables
  Type q=0.0;
  Type s1=0.0;
  Type log_sigma_index=0.0;
  Type nll = 0.0;
  Type S0=0.0;
  Type SPR=0.0;
  
  using namespace std;
  using namespace density;
  
  // Full Fishing Mortality
  for(y=0;y<nyrs;y++) FM(y)=exp(u[y]);
  //-------------Calculate SPR F=0 for S0 -----------------------------------------------
  ntemp(0)=1;
  SPR=ntemp(0)*WAA(0)*maturity(0);
  for(a=1;a<nages;a++){
    if(a<=(nages-2)) ntemp(a)=ntemp(a-1)*exp(-NM);
    if(a==(nages-1)) ntemp(a)=(ntemp(a-1)*exp(-NM))/(1-exp(-NM));
    SPR+=ntemp(a)*WAA(a)*maturity(a);
  } //checked calculation in EXCEL
  //Calculate S0
  S0=R0*SPR;
  
  //---------------------------- Abundance at Age Matrix ----------------------------------

  for(y=0;y<nyrs;y++){
    B(y)=0.0;
    pred_catch(y)=0.0;
    SSB(y)=0.0;
    for(a=0;a<nages;a++){
       if(y==0){ //first row of N matrix
         if(a==0) N(y,a)=R0;
         if(a>0)  N(y,a)=N(y,a-1)*exp(-NM);
          SSB(y)+=N(y,a)*WAA(a)*maturity(a);
          B(y)+=N(y,a)*WAA(a)/1000;
          pred_catch(y)+=(select(a)*FM(y))/(select(a)*FM(y)+NM)*(1-exp(-select(a)*FM(y)-NM))*N(y,a)*WAA(a)/1000;
       }
      if(y>0){ //remaining rows of N matrix
         if(a==0) N(y,a)=recruits(R0,h,SSB(y-1),S0);
         if(a>0 && a<=(nages-2))  N(y,a)=N(y-1,a-1)*exp(-NM-select(a-1)*FM(y-1));
         if(a==(nages-1))  N(y,a)=N(y-1,a-1)*exp(-NM-select(a-1)*FM(y-1))+N(y-1,a)*exp(-NM-select(a)*FM(y-1));
         SSB(y)+=N(y,a)*WAA(a)*maturity(a);
         B(y)+=N(y,a)*WAA(a)/1000;
         pred_catch(y)+=(select(a)*FM(y))/(select(a)*FM(y)+NM)*(1-exp(-select(a)*FM(y)-NM))*N(y,a)*WAA(a)/1000;
       }
    }
  }
  
  //Biomass last year+1
  B(nyrs)=0.0;
  for(a=0;a<nages;a++){
    if(a==0) B(nyrs)=recruits(R0,h,SSB(nyrs-1),S0)*WAA(a)/1000;
    if(a>0 && a<=(nages-2))  B(nyrs)+=N(nyrs-1,a-1)*exp(-NM-select(a-1)*FM(nyrs-1))*WAA(a)/1000;
    if(a==(nages-1)) B(nyrs)+=(N(nyrs-1,a-1)*exp(-NM-select(a-1)*FM(nyrs-1))*WAA(a)+N(nyrs-1,a)*exp(-NM-select(a)*FM(nyrs-1))*WAA(a))/1000;
   }
  
  //--------------Observation Model--------------------------------------
  for(iy=0;iy<niyrs;iy++){
    q+=log(index(iy))-log(B(index_years(iy)));
  }
  q=exp(q/niyrs); //catchability for index
  
  for(iy=0;iy<niyrs;iy++){
    s1+=square(log(index(iy))-log(q*B(index_years(iy))));
  }
  log_sigma_index=sqrt(s1/niyrs); //catchability for index
  
  for(iy=0;iy<niyrs;iy++){
    pred_index(iy)=q*(B(index_years(iy))+B(index_years(iy)+1))/2;
    nll-=dnorm(log(index(iy)),log(pred_index(iy)),log_sigma_index,true);
  } 
  //Catch 
  for(y=0;y<nyrs;y++){
    nll-=dnorm(log(catches(y)),log(pred_catch(y)),exp(log_sigma_catch),true);
  }
 
  //----------------------Output-----------------------------------------

  ADREPORT(FM);
  REPORT(B);
  ADREPORT(R0);
  //ADREPORT(frac);
  ADREPORT(h);
  REPORT(log_sigma_index);
  REPORT(pred_index);
  REPORT(pred_catch);
  REPORT(FM);
  REPORT(SSB);
  REPORT(N);
  REPORT(index);
  REPORT(catches);
  return nll;
 
}

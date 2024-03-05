#include <TMB.hpp>

//function to implement the schaefer production model for a single time step 
//It should be the same as the biomass dynamics exercise
template <class Type> Type schaefer(Type b1, Type r, Type K, Type m, Type u_i){
  Type u_use = exp(u_i)/(1.0+exp(u_i)); //logit transforming exploitation rate Ui
  Type b2 = b1 + b1*r*(1.0-pow(b1/K,m-1.0)) - u_use*b1; //literally just the shaefer model equation.
  return b2;
}

/* DATA */
DATA_VECTOR(catches);
DATA_VECTOR(index);
DATA_IVECTOR(index_years); //integer vector
int nyrs = catches.size(); //number of unique observations of catch
int nobs = index.size(); // number of unique observations of the survey biomass. 

/* PARAMETERS */
PARAMETER(logK); //carrying capacity in log space
PARAMETER(logr); // intrinsic growth rate in log space
PARAMETER(m); // m is the shape param that determine location of BMSY must be in logit space.
PARAMETER(logSigmaC); // observation error on catch
PARAMETER(logSigmaI); // observation error on survey
PARAMETER_VECTOR(upar); //we are estimating exploitation rate for each year.
PARAMETER(logq); //we weren't given an estimator for q in this model, so i assume we estimate it in the model.

//temporary variables//
vector<Type> biomass(nyrs+1); //a vector where we will store predicted true biomass in the fishing series, adding an extra year because we start at 0. 
vector<Type> catches(iyr)
vector<Type> index_pred(nobs); //store predicted survey biomass (different from true biomass)
vector<Type> predbio(nobs); //store predicted survey biomass (different from true biomass)

//backtransform parameters//
Type m_use = exp(m)/(1.0+exp(m)); //m is in logit space, so you have to un-logit it?

//initialize the nll 
Type neglogL = 0.; 

//* EFFECTS OF FISHING ON BIOMASS *//
//the shaefer model describes how fishing is affecting the population starting at virgin biomass. 
biomass(0) = exp(logK); //the first year of biomass in the fishing series is carrying capacity. 
catch_pred(0)= biomass(0)*exp(upar(0))/(1.0+exp(upar(0)));
for (int iyr=1;iyr<=nyrs;iyr++) { //iterating over the observations of catch
  biomass(iyr) = schaefer(biomass(iyr-1),exp(logr),exp(logK),m_use,upar(iyr-1)); //giving it all the logged params
  catch_pred(iyr) = biomass(iyr)*exp(upar(iyr))/(1.0+exp(upar(iyr))); //in the schaefer model the u parametr is already logit transformed. here it isn't. 
//contribution to the nll
}
//I am putting this outside the for loop even though the arguments "catches" and "catch_pred" are vectors
// whereas the argument "logSigmaC" is not.
neglogL -= sum(dnorm(log(catches),log(catch_pred),exp(logSigmaC), true));


//*SURVEY MODEL*//
for (int iobs=0; iobs<nobs; iobs++) {
  //the index in year t is the average of the biomass in year t
  predbio(iobs) = 0.5*(biomass(index_years(iobs))+biomass(index_years(iobs) + 1));//why is this plus 1?? not  -1??
  logq += log(index(iobs)/predbio(iobs));   
}
logq = logq/nobs;
index_pred = exp(logq)*predbio; //*exp(logq);
neglogL -= sum(dnorm(log(index),log(index_pred),exp(logSigmaI), true));


//std::cout << catch_pred << std::endl;


ADREPORT(logq);
ADREPORT(biomass);

return neglogL;

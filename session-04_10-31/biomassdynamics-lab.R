## surplus production model
## 2023-10-31

## Load TMB TRUE
library(TMB)
library(tidyverse)
library(readxl)


## Data and parameters
catch <- read_xlsx("session-04_10-31/flounder_index_data.xlsx","catch")
index <- read_xlsx("session-04_10-31/flounder_index_data.xlsx","survey")

index_years <- index$Year - min(catch$Year)

data <- list(catches = catch$Catch,
               index = index$Index,
               index_years = index_years)
data

parameters <- list(logK=10, #in log space because carrying capacity must be positive
                   logr=log(0.4), 
                   m=2, #shape param
                   logSigmaC=log(0.05), #observation error for catch, must be positive
                   logSigmaI=2, #observation error for survey
                   #logq=-2,
                   upar=rep(-1.6, length(data$catches))) #logit space, must be between 0-1 and 
                    #it's the time series of exploitation rates. one for every year. 
  
parameters


## Make C++ file
#TMB::template("session-04_10-31/biomassdynamics.cpp")

## Compile and load the model
dll_name <- "biomassdynamics"
if(is.loaded(dll_name)) {
  dyn.unload(dynlib(dll_name))
}
compile("biomassdynamics.cpp")
dyn.load(dynlib(dll_name))


## Make a function object
obj <- MakeADFun(data, parameters, DLL= dll_name,
                 map = list( ______ ),
                 control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr, control= list(eval.max = 10000,
                                                     iter.max = 10000,
                                                     rel.tol = 1e-15))


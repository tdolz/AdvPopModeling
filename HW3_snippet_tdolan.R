#HW 3 R snippet
# T Dolan 11/19/23

#Load required packages
library(TMB)
library(kableExtra)
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", force=TRUE) # helper function to get AIC
library(TMBhelper)
library(tidyverse) 
library(readxl) 
library(parallel)

#catch data from 1935-2016
catch <- read_xlsx("/Users/tdolan/Documents/R-Github/AdvPopModeling/session-04_10-31/flounder_index_data.xlsx","catch") 
#survey index from 1963-2016
index <- read_xlsx("/Users/tdolan/Documents/R-Github/AdvPopModeling/session-04_10-31/flounder_index_data.xlsx","survey")

#set year 0 = 1935 the beginning of the catch series,such that the survey starts 28 years after the first catch year. 
index_years <- index$Year - min(catch$Year)

#Data for input into all models. 
data <- list(catches = catch$Catch,
             index = index$Index,
             index_years = index_years)
data

parameters <- list(logK=9, 
                   logr=log(0.56), 
                   logit_m=-0.9, 
                   logSigmaC=log(0.05),   
                   upar = rep(-1.0,length(data$catches)))
parameters

## Compile and load the model
dll_name <- "Simulate_HW3"
if(is.loaded(dll_name)) {
  dyn.unload(dynlib(dll_name))
}
compile("Simulate_HW3.cpp")
dyn.load(dynlib(dll_name))

######### Run the simulation multiple times

simdat <-list()
for (i in 1:5){ #only doing 5 for now
  
  ## Make a function object - whether I index the function object with [[i]] doesn't seem to make a difference.
  obj_sim[[i]] <- MakeADFun(data, parameters, DLL= dll_name,
                            map = list(logSigmaC=factor(NA)),
                            control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))
  
  sim_data <-obj_sim[[i]]$simulate(complete=TRUE)
  simdat[[i]] <-sim_data
}

#check if catches are too high (in the 1000s when they are supposed to be in the 10s) 
#and check if simulation runs are identical
simdat[[1]]
simdat[[2]]


############### Try simulating in R ##########################################

## Make a function object - whether I index the function object with [[i]] doesn't seem to make a difference.
obj_sim <- MakeADFun(data, parameters, DLL= dll_name,
                          map = list(logSigmaC=factor(NA)),
                          control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

simdat2 <-list()
for (i in 1:5){ #only doing 5 for now
  
  ## Make a function object - whether I index the function object with [[i]] doesn't seem to make a difference.
  list_sims <- simulate(obj_sim, nsim=5, seed=123) #simulate not working with this list
  simdat[[i]] <-list_sims
}

#check if catches are too high (in the 1000s when they are supposed to be in the 10s) 
#and check if simulation runs are identical
simdat2[[1]]
simdat2[[2]]
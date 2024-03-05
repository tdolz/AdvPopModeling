## surplus production model
## 11/6/23

## generalized Pella-Tomlinson model:
## B_(y+1)= B_y + B_y*r*(1-(B_y/K)^(m-1))-U_y*B_y
#
#B_y = biomass in year y (from our survey series)
#r = intrinsic growth rate
#K = carrying capacity
#U_y = exploitation rate in year y. (from our catch series)
#m = shape parameter that determines the location of Bmsy. When m=2, we get the Schaefer model. 

#The schaefer model: 
#assume the initial biomass is at K
#assume survey index obs are proportional to biomass and normally distributed. 
# Ihat_y = q*(B_y+B_(y+1)/2) #survey index is a moving average of past and present biomass, scaled by catchability q. 
# logs of the catch are normally distributed
# we have an analytical MLE for q. 

#Estimable parameters: 
# K,r,m, sigmaI, sigmaC. 
# std. dev. of log catches =0.05 and m=2. 

## Load TMB TRUE
library(TMB)
library(tidyverse)
library(readxl)

setwd("/Users/tdolan/Documents/R-Github/AdvPopModeling")
## Data and parameters

#catch data from 1935-2016
catch <- read_xlsx("session-04_10-31/flounder_index_data.xlsx","catch") 
#survey index from 1963-2016
index <- read_xlsx("session-04_10-31/flounder_index_data.xlsx","survey")
#set year 0 = 1935 the beginning of the catch series, 
#such that the survey starts 28 years after the first catch year. . 
index_years <- index$Year - min(catch$Year)

data <- list(catches = catch$Catch,
               index = index$Index,
               index_years = index_years)
data

#a function from Gary to find the initial value of M in logit space if it's 2 normally
ss<-function(m){
  ee<-(3*exp(m)/(1+exp(m)))+1 #applying our same modified logit transformation from the cpp side
  return((ee-2)^2)#assume m=2, return the squared difference
}
dd <-optimize(ss,c(-10,10)) #find the likelihood minimum over a range -10,10
dd$minimum

#providing initial values for params based on param values from the Schaefer model. 
parameters <- list(logK=9, 
                      #list of previous tries
                            #10
                            #6.432107,
                            #9
                   logr=log(0.56), 
                      #list of previous tries
                            #0.4, but can't do that because it's negative.
                            #exp(-0.578463)= 0.5607596? 
                            #0.56
                   logit_m=-0.9, 
                      #list of previous tries
                            #-0.6931638, #logit transformed value of 2.
                            #
                   logSigmaC=log(0.05),  
                      #list of previous tries
                            #exp(-2.99) = 0.05005004,
                            #previously set to 0.05, param value was -2.99, but cant be negative
                   #logSigmaI=log(0.23), #now estimating analytically. 
                      #list of previous tries
                            #previously set to 0.3, param value was -1.217, but can't be negative.
                            #what if we try exp(-1.217081) = 0.2960932, so set it back to that. 
                   #logq = -2, 
                      #list of previous tries
                   upar = rep(-1.0,length(data$catches))) #previously -1.6 
                      #list of previous tries
                      #mean of the upar estimates from the schaefer model = -0.9
                      #Gary said to try -1.
                      #try back transforming -0.9 from logit

#check to make sure there are no NAN in params!
parameters


## Make C++ file
#TMB::template("HW3_biomassdynamics_gavin.cpp")

## Compile and load the model
dll_name <- "HW3_biomassdynamics_gavin"
if(is.loaded(dll_name)) {
  dyn.unload(dynlib(dll_name))
}
compile("HW3_biomassdynamics_gavin.cpp")
dyn.load(dynlib(dll_name))


## Make a function object
obj <- MakeADFun(data, parameters, DLL= dll_name,
                 #map = list(m=factor(NA), #turning this off so that we are estimating m
                            map = list(logSigmaC=factor(NA)),
                 control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr, control= list(eval.max = 10000,
                                                       iter.max = 10000))

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

pl <- obj$env$parList(opt$par) 

q <- summary(sdr, type = "report") %>% 
  as.data.frame() %>% 
  janitor::clean_names() %>% 
  rownames_to_column(var="type") %>% 
  filter(str_detect(type, "logq"))

summary(sdr, type = "report") %>% 
  as.data.frame() %>% 
  janitor::clean_names() %>% 
  rownames_to_column(var="type") %>% 
  filter(str_detect(type, "biomass")) %>% 
  separate(type, c("type","yr")) %>% 
  mutate(yr = as.numeric(ifelse(is.na(yr),0,yr))) %>% 
  ggplot() +
  aes(x = yr, y = estimate) +
  geom_line() +
  geom_point(data = index, aes(x = Year-1935, y = Index/exp(q[1,2])),
             col = "darkgreen")

  





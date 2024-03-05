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

parameters <- list(logK=10, 
                   logr=log(0.4), 
                   m=2, 
                   logSigmaC=log(0.05), 
                   logSigmaI=log(0.2),
                   #logq = -2, 
                   upar = rep(-1.6,length(data$catches)))
parameters


## Make C++ file
#TMB::template("biomassdynamics.cpp")

## Compile and load the model
dll_name <- "biomassdynamics_gavin"
if(is.loaded(dll_name)) {
  dyn.unload(dynlib(dll_name))
}
compile("biomassdynamics_gavin.cpp")
dyn.load(dynlib(dll_name))


## Make a function object
obj <- MakeADFun(data, parameters, DLL= dll_name,
                 map = list(m=factor(NA),logSigmaC=factor(NA)),
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

  





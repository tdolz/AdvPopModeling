#Hw 4 snippet
#Tara Dolan 12/3/23 UPDATED 12/5/23

library(TMB)
library(tidyverse)
library(readxl)

##################################### load and format data ######################################################
#catch data from 1935-2016
catch <- read_xlsx("session-04_10-31/flounder_index_data.xlsx","catch") 
catch1963 = filter(catch, Year > 1962) # for this first exercise, we're setting year 0 to 1963
#survey index from 1963-2016
index <- read_xlsx("session-04_10-31/flounder_index_data.xlsx","survey")
#set year 0 = 1963 the beginning of the survey series, 
index_years <- index$Year - min(catch1963$Year)

## age comps ##
age <- read_xlsx("hw4/floundah_agecomps.xlsx") 
age_comp = sapply(age[,-1], '*', 100)%>%as.matrix()

##biodata
maturity <- c(0.1,0.5,0.9,1,1,1,1,1,1,1)
WAA <-c(0.040,0.160,0.331,0.511,0.675,0.812,0.922,1.006,1.070,1.117)

####################################### set up data and params ##################################################################
#data
data <- list(age_comp=age_comp,
             catches = catch1963$Catch,
             index = index$Index,
             weight=WAA, 
             maturity=maturity)
data

#parameters
parameters <- list(logR0= log(1000), 
                   logfracR1= log(1),  
                   logSigmaI= log(0.2), 
                   logF= rep(-1.6,length(data$catches)),
                   #Mapped params
                   logith=0.6,
                   logM= log(0.4), 
                   logSigmaC= log(0.05),
                   #New params
                   logEtaT=rep(0,length(data$catches)), #based on aspm R.
                   logSigmaR=log(0.6), #from Gavin's email.  
                   logA50s=log(7.389), #about 2 
                   A95s=log(2.71828), # about 1
                   A50f=log(7.389), #about 2 
                   A95f=log(2.71828) # about 1
)  

parameters

#################################### compile model ####################################################################
## Compile and load the model
dll_name <- "tdolan_hw4" # # changed 
if(is.loaded(dll_name)) {
  dyn.unload(dynlib(dll_name))
}
compile("tdolan_hw4.cpp")
dyn.load(dynlib(dll_name))


## Make a function object
obj2 <- MakeADFun(data, parameters, DLL= dll_name,
                  map = list(
                    logM=factor(NA),
                    logith=factor(NA),
                    logSigmaC=factor(NA)),
                  control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, control= list(eval.max = 10000,
                                                         iter.max = 10000))

###################################### inspect outputs #########################################################

## Get parameter uncertainties and convergence diagnostics
sdr2 <- sdreport(obj2)
sdr2

p2 <- obj2$env$parList(opt2$par) 
p2
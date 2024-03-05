## Load TMB
library(TMB)
library(tidyverse)
library(FSAdata)

## Make C++ file
#TMB::template("growthincrement.cpp")

## 
#data(BluegillIL)
#str(BluegillIL)
head(BluegillIL)
#plot((lenRecap-lenMark)~deltaTime,data=BluegillIL)
data <- tibble(len1 = BluegillIL$lenMark,
               len2 = BluegillIL$lenRecap,
               dt = BluegillIL$deltaTime)
parameters <- list(age_1 = rep(0,nrow(data)),
                   log_sigma_inf = 0, 
                   log_sigma = 0,
                   log_mu_linf = 0,
                   log_k = 0,
                   log_mu_a1 = 0,
                   log_sigma_a1 = 0)
                   # log_shape,
                   # log_scale);


## Compile and load the model
compile("growthincrement.cpp")
dyn.load(dynlib("growthincrement"))


## Make a function object
obj <- MakeADFun(data = data, parameters = parameters,
                 random = "age_1", DLL = "growthincrement")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace=1))

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
summary(sdr)
#summary(sdr,select="random")



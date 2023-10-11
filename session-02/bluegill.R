#growth estimation from lengths at release and recapture, bluegill
#10/10/23

library("FSA")
library("FSAdata")

head(FSAdata::BluegillIL)


##### to the TMB Robin! ########

## Load TMB TRUE
library(TMB)

## Make C++ file
TMB::template("session-02/Bluegill.cpp")

## Compile and load the model
compile("session-02/Bluegill.cpp")
dyn.load(dynlib("session-02/Bluegill"))

## Data and parameters
data <- list(density = Streams$Density,
             stream = Streams$Stream)
parameters <- list(mu=0, logSigma=0, logSigmaS=0,StreamDevs=rep(0,length(unique(Streams$Stream))))

## Make a function object
obj <- MakeADFun(data, parameters, random = "StreamDevs", DLL="streams")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
#derived quantities too
summary(sdr)

#random effects or fixed effects alone. 
summary(sdr,select="random")
summary(sdr,select="fixed")


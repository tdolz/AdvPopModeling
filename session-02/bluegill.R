#growth estimation from lengths at release and recapture, bluegill
#10/10/23
library("FSA")
library("FSAdata")

head(FSAdata::BluegillIL)
BG = FSAdata::BluegillIL


##### to the TMB Robin! ########

## Load TMB TRUE
library(TMB)

## Make C++ file - run only once to generate the file
#TMB::template("session-02/Bluegill.cpp")

## Compile and load the model
compile("session-02/Bluegill.cpp")
dyn.load(dynlib("session-02/Bluegill"))

## Data and parameters - gary's way
#data <- list(l1 = BG$lenMark, l2=BG$lenRecap,dt=BG$deltaTime)
#parameters <- list(K=0, Linf_mu=0, sigma2_Linf=0,sigma2_resid=0,sigma2_A=0,mu_A=0,
 #                  A_i=rep(0,length(unique(BG$tag))))

## Data and parameters - gavin's way
data <- tibble(len1 = BG$lenMark, len2=BG$lenRecap,dt=BG$deltaTime)
parameters <- list(log_sigma_inf=0,
                   log_sigma=0,
                   log_mu_linf=0,
                   log_k=0,
                   mu_a1=0,
                   sigma_a1=0,
                   age1=rep(0,nrow(BG)))#last line is random effects. 
#estimate everything in log space to constrain to positive numbers. 


## Make a function object
obj <- MakeADFun(data, parameters, random = "A_i", DLL="Bluegill")

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


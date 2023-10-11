### first TMB example
library(TMB)
library(tidyverse) 

#TMB:::setupRStudio()


#linear regression
# Generate some data
set.seed(8675309)
data <- tibble(x = 1:20,
               y = 0.5 + 2*x + rnorm(20,0,2))
# view the data set
data
ggplot(data) +
  geom_point(aes(x=x,y=y)) +
  geom_smooth(aes(x=x,y=y),method="lm") +
  theme_minimal()

#  yhat(i) = b0 + b1*x(i)
lm1 <- lm(y~x, data = data)
lm1

# model parameters
parameters <- list(b0=0, b1=0, logSigma=0)
print(parameters)

#compile the TMB model
#require(TMB)
compile("session-01/linear_regression.cpp")

# load the compiled model object
dyn.load(dynlib("session-01/linear_regression"))

################################################################################
# make the TMB model
model <- MakeADFun(data = data, 
                   parameters = parameters,
                   #map = list(b0=factor(NA)),
                   DLL="linear_regression")
#control=list(eval.max=10000,iter.max=1000,rel.tol=1e-15),silent=T)

# look at the built model
print(attributes(model))

# fit the model using nlminb
fit <- nlminb(model$par, model$fn, model$gr) #estimator nlminb - TMB has created the gradient
fit # at the minimum the gradient should be ~0
#$par is the params, 
#$objective function value
#$converge yest or no
#These outputs are explained in the help file for nlminb

# get variance estimates for model parameters
rep <- sdreport(model) #asymptotic variances for parameters
rep #intercept slope and estimated variance in the log scale. 

# Sumamrize ALL
print(summary(rep,p.value=T)) #but we also looked for the standard deviation in real space. 
#sigma has to be a positive number so he logged it so it could be constrained as like a hack so 
#he didn't have to put on bounds. then he just unlogged. 

#summary(rep,select="random")
summary(rep,select="fixed") #the summary object can get us the derived variable
#plot the predictions
ggplot(data) +
  aes(x=x,y=y) +
  geom_point()+
  geom_smooth(method=lm)+
  geom_abline(intercept = fit$par["b0"],
              slope = fit$par["b1"]) +
  theme_minimal()



#Use of the SIMULATE function
# uncomment out the SIMULATE section of the code then re-compile
compile("session-01/linear_regression.cpp")
dyn.load(dynlib("session-01/linear_regression"))
# make the TMB model
model <- MakeADFun(data = data, 
                   parameters = parameters,
                   DLL="linear_regression")
fit <- nlminb(model$par, model$fn, model$gr)

#generate 10 simulated data sets based on the estimated parameters
sim_data <- sapply(1:10,function(x)model$simulate(complete=TRUE))
sim_data[,1]$y
sim_data[,2]$y




### TMB with STAN
library(shinystan)
library(tmbstan)
# Run a single chain in serial with defaults
fit <- tmbstan(model, chains=1)

## Run in parallel with a init function
cores <- parallel::detectCores()-1
options(mc.cores = cores)
 init.fn <- function()
   list(u=rnorm(114), beta=rnorm(2), logsdu=runif(1,0,10), logsd0=runif(1,0,1))
fit <- tmbstan(model, chains=cores, open_progress=FALSE, init="par")

## To explore the fit use shinystan
launch_shinystan(fit)




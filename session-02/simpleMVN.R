#simpleMVN
## Load TMB TRUE
library(TMB)

p <- 3 # three groups
n <- 10000 #10000 individuals per group
mu <- c(10, 20, 30) #mean of each group
sd <- c(3, 2, 1) #sd of each group
rho <- 0.5 #covariance
#VCV matrix
# recall, variance is sd^2, and covariance is rho = sdxy/sdx*sdy, solve for sdxy if you have rho. 
Sigma <- matrix(c(sd[1]*sd[1], sd[1]*sd[2]*rho, sd[1]*sd[3]*rho,
                  sd[2]*sd[1]*rho,sd[2]*sd[2], sd[2]*sd[3]*rho,
                  sd[3]*sd[1]*rho,sd[3]*sd[2]*rho, sd[3]*sd[3]),
                nrow=p, ncol = p, byrow=TRUE)
set.seed(42)
X <- MASS::mvrnorm(n=n, mu=mu, Sigma) #multivariate normal dist
#pairs(X, lower.panel = NULL)

data <- list(X=X) #X is a mvnorm dist, with n draws where n =1000
parameters <- list(mu=rep(0,3),logSigma=rep(0,3),rhopar=1) #initialize empty param vector


## Make C++ file
#TMB::template("session-02/simpleMVN")

## Compile and load the model
compile("session-02/simpleMVN.cpp")
dyn.load(dynlib("session-02/simpleMVN"))

## Data and parameters

## Make a function object
obj <- MakeADFun(data, parameters, DLL="simpleMVN")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

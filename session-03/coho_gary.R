## univariate coho salmon example, section 3.2.1

set.seed(2345)
num.times <- 40 #40 years
Ricker.a <- 1.5 #from fig 3.2, implicitly includes survival between t-1 and t, and survival from egg to juv.
# a < 2.69 to avoid chaotic behavior, and a , 2 to avoid limit cycles. we must contstrain. 
Ricker.b <- 0.0003 #from fig 3.2, a measure of density dependence. 
equilibrium <-  log(Ricker.a)/Ricker.b
cat("equilibrium=", round(equilibrium),"\n")
initial.n <- round(0.1*equilibrium)
juveniles <- estimates <- numeric(num.times)

# 
juveniles[1] <- rpois(1,lambda = Ricker.a*initial.n*exp(-Ricker.b*initial.n))
CV.obs <- .30
obs.sd <- sqrt(log(CV.obs^2+1))
estimates[1] <- round(rlnorm(1,meanlog=log(juveniles[1]),sdlog=obs.sd))
for(i in 2:num.times) {
  juveniles[i] <- rpois(1,lambda=Ricker.a*juveniles[i-1]*exp(-Ricker.b*juveniles[i-1]))
  estimates[i] <-  round(rlnorm(1,meanlog=log(juveniles[i])-obs.sd^2/2,sdlog=obs.sd))
}
my.ylim <- range(c(juveniles,estimates))
plot(1:num.times,juveniles,xlab="Year",ylab="",ylim=my.ylim,type="l")
lines(1:num.times,estimates,lty=2,col=2)
legend(2,0.95*max(my.ylim),legend=c("State","Observation"),lty=1:2,col=1:2)

## Compile and load the model
## 
## Load TMB TRUE
library(TMB)
compile("session-03/coho_gary.cpp")
dyn.load(dynlib("session-03/coho_gary"))

data <- tibble(y = estimates)
parameters <- list(logit_alpha=0,
                   log_beta=0,
                   log_sigma_obs=0,
                   logN0=0,
                   u=rep(0,40))

## Make a function object
obj <- MakeADFun(data, parameters, DLL="coho_gary")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr


## Load TMB
library(TMB)
library(tidyverse)
library(FSAdata)

## Make C++ file
#TMB::template("coho.cpp")


## univariate coho salmon example, section 3.2.1

set.seed(2345)
num.times <- 40
Ricker.a <- 1.5
Ricker.b <- 0.0003
equilibrium <-  log(Ricker.a)/Ricker.b
cat("equilibrium=", round(equilibrium),"\n")

initial.n <- round(0.1*equilibrium)
juveniles <- estimates <- numeric(num.times)
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


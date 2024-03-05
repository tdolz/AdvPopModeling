### first TMB example
library(TMB)
library(tidyverse) 


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

data <- list(x = data$x,
             y = data$y)

#compile the TMB model
#require(TMB)
compile("session_11_2/linear_regression.cpp")

# load the compiled model object
dyn.load(dynlib("session_11_2/linear_regression"))

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
fit <- nlminb(model$par, model$fn, model$gr)
fit
# get variance estimates for model parameters
rep <- sdreport(model)
rep

# Sumamrize ALL
print(summary(rep,p.value=T))

#summary(rep,select="random")
summary(rep,select="fixed")
#plot the predictions
ggplot(data) +
  geom_point(aes(x=x,y=y)) +
  geom_abline(intercept = fit$par["b0"],
              slope = fit$par["b1"]) +
  theme_minimal()



#Use of the SIMULATE function
#kind of like a parametric bootstrap. 
# uncomment out the SIMULATE section of the code then re-compile
# can be used for anything stochastic. 
# could be process or observation errors. 
compile("session_11_2/linear_regression.cpp")
dyn.load(dynlib("session_11_2/linear_regression"))
# make the TMB model
model <- MakeADFun(data = data, 
                   parameters = parameters,
                   DLL="linear_regression")
fit <- nlminb(model$par, model$fn, model$gr)

#om values
om <- sdreport(model) %>% 
  summary() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "parameter") %>% 
  janitor::clean_names() %>%
  rename(om_val = estimate) %>% 
  select(-std_error)

#generate 10 simulated data sets based on the estimated parameters
sim_data <- map(1:10,function(x)model$simulate(complete=TRUE))
str(sim_data, max.level = 1)
sim_data[[1]]  #list of simulated datasets
sim_data[[1]]$y
sim_data[[2]]$y

#make a function to fit the model to a dataset
#refit the model to the simulated data. this is a function here to make the fitting and parameter extraction in one go. 
refit <- function(data, parameters, dll_name) {
  model_sim <- MakeADFun(data = data, 
                     parameters = parameters,
                     DLL=dll_name)
  fit <- nlminb(model_sim$par, model_sim$fn, model_sim$gr)
  sdr <- sdreport(model_sim) %>% 
    summary() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "parameter")
  return(sdr)
}

#re-fit the model to each simulated dataset
#now we are using the refit function to fit all the simulated datasets. 
sim_estimates <- map_dfr(sim_data, refit, 
                         parameters = parameters, 
                         dll_name = "linear_regression", .id = "isim") %>% 
  janitor::clean_names()

#plot the estimates
sim_estimates %>% 
  ggplot() +
  aes(x = parameter, y = estimate) +
  geom_boxplot()
  
#plot the relative error of the estimates
sim_estimates %>% 
  left_join(om) %>%
  mutate(ree = (estimate-om_val)/om_val) %>% 
  ggplot() +
  aes(x = parameter, y = ree) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2) +
  NULL




# ### TMB with STAN
# library(shinystan)
# library(tmbstan)
# # Run a single chain in serial with defaults
# fit <- tmbstan(model, chains=1)
# 
# ## Run in parallel with a init function
# cores <- parallel::detectCores()-1
# options(mc.cores = cores)
# # init.fn <- function()
# #   list(u=rnorm(114), beta=rnorm(2), logsdu=runif(1,0,10), logsd0=runif(1,0,1))
# fit <- tmbstan(model, chains=cores, open_progress=FALSE, init=fit$par)
# 
# ## To explore the fit use shinystan
# launch_shinystan(fit)




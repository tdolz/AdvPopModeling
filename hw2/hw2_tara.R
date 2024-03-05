## Tara Dolan
## hw2 adv pop modelling
## 10/23/23

library(tidyverse)
############## EXPLORE DATA #################
sslData <- readRDS("hw2/ssl-data.rds")
pups <-as.data.frame(sslData$forrester_pups)%>%rename("icount_p"="icount")
nonpups <-as.data.frame(sslData$forrester_nonpups)%>%rename("icount_np"="icount")
ssl <-full_join(pups, nonpups, by="Year")%>%arrange(Year)
ssl_piv <-select(ssl, -icount_p, -icount_np)%>%
  pivot_longer(c("pupcount","adultcount"), names_to = "stage",values_to = "count")
ssl_sum <-ssl_piv%>%group_by(Year,stage)%>%
  summarize(av_count=mean(count),sd_count=sd(count), cv_count=sd(count)/mean(count),
            n_years=n_distinct(Year),.groups="keep")
sslsum_wide <-select(ssl_sum, av_count, Year, stage)%>%
  pivot_wider(names_from = stage, values_from = av_count)%>%
  mutate(prod=pupcount/adultcount)

#In some years there are more than one count of either adults or pups. 
#However, there are two few cases with multiple counts to include count as a factor, so we will summarize. 
#Gary solved this with averaging them. 
  ggplot()+
  geom_point(aes(Year,count, col=stage),data=ssl_piv)+
    #geom_line(aes(Year,count, col=stage),data=ssl_piv)+
    scale_color_manual(values = c("red", "blue"))+
  #geom_line(aes(Year, av_count,col=stage),data=ssl_sum)+
  labs(title="Count of Stellar Sea Lions")+
  ylab("count of individuals")+
  theme_minimal()
  
#average count
ssl_sum%>%
  ggplot()+
  geom_point(aes(Year,av_count, col=stage))+
  scale_color_manual(values = c("red", "blue"))+
  labs(title="Count of Stellar Sea Lions")+
  ylab("count of individuals")+
  theme_minimal()

#pups per adult
sslsum_wide%>%
  ggplot()+
  geom_point(aes(Year,prod))+
  labs(title="Stellar sea lion productivity")+
  ylab("pups per adult")+
  theme_minimal()

ssl%>%
  ggplot()+
  geom_point(aes(adultcount,pupcount))+
  labs(title="Stellar sea lion productivity")+
  ylab("pup count")+xlab("adult count")+
  theme_minimal()

ssl%>%
  ggplot()+
  geom_point(aes(adultcount,pupcount))+
  labs(title="Stellar sea lion productivity")+
  ylab("pup count")+xlab("adult count")+
  theme_minimal()

### Prepare the data for TMB ###
Year = data.frame(Year=seq(min(ssl$Year),max(ssl$Year),1))
ssltmb = right_join(ssl,Year)%>%arrange(Year)


################ TMB INTERFACE ################################
## Load TMB TRUE
library(TMB)

## Make C++ file - run only once to generate the file
#TMB::template("hw2/hw2_tara.cpp")

#compile the script
compile("hw2/hw2_tara.cpp")
dyn.load(dynlib("hw2/hw2_tara"))

#feed in data
data <- tibble(pups_data = ssltmb$pupcount,
               nonpups_data = ssltmb$adultcount,
               year =ssltmb$Year)
#feed in params
params <-list(pups=rep(0,length(data$year)),
              nonpups=rep(0,length(data$year)),
              f=0,
              phi_p=0,
              phi_np=0,
              logsigma2_p=0,
              logsigma2_np=0,
              q=0,
              tau2=0)

#make a function object and declare random variables. 
obj<- MakeADFun(data, 
                parameters, 
                DLL="hw2_tara",
                random="tau2")#not really sure what the random effects would be here.

#fit the model(?)
fit <- nlminb(model$par, model$fn, model$gr)

#Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
summary(sdr)

# Question 5: once we have our parameter values, we can run a stochastic projection #

# maybe with the simulate function? but how to extend out to 2023
    #generate 10 simulated data sets based on the estimated parameters
    sim_data <- sapply(1:10,function(x)obj$simulate(complete=TRUE))
    
# maybe using the predict() function in R, though I am not sure if that's deterministic?

    

    

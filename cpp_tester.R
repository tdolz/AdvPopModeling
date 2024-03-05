### Run C++ test file tester.cpp

library(TMB)
library(tidyverse)
library(readxl)

setwd("/Users/tdolan/Documents/R-Github/AdvPopModeling")

##################################### input data ############################################################
#catch data from 1935-2016
catch <- read_xlsx("session-04_10-31/flounder_index_data.xlsx","catch") 
catch1963 = filter(catch, Year > 1962)
index <- read_xlsx("session-04_10-31/flounder_index_data.xlsx","survey")
## age comps ##
age <- read_xlsx("hw4/floundah_agecomps.xlsx") 
##biodata 
maturity <- c(0.1,0.5,0.9,1,1,1,1,1,1,1)
WAA <-c(0.040,0.160,0.331,0.511,0.675,0.812,0.922,1.006,1.070,1.117)


####################################### format data and params ################################################
### inputs
data <- list(age_comp=age,
             catches = catch1963$Catch,
             index = index$Index,
             weight=WAA, 
             maturity=maturity)
data

########################################### compile and load the model ################################################
## Compile and load the model
dll_name <- "tester" # # changed 
if(is.loaded(dll_name)) {
  dyn.unload(dynlib(dll_name))
}
compile("tester.cpp")
dyn.load(dynlib(dll_name))


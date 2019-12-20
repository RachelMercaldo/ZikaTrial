## Code to initiate hpc cluster and run analyzefxn() over all  possible parameter sets
################################################################
## Code accompanies:
## 
## Mercaldo, RA, Bellan, SE. Evaluation of Alternative Endpoints for Zika Virus Vaccine Efficacy Trials. 2019.
##
## Rachel Mercaldo, 2019
## License at bottom.


#!/usr/bin/Rscript
library(foreach)
library(doParallel)
library(iterators)
library(gsDesign)
library(tidyverse)
library(data.table)
library(mefa)
library(survival)


# The code to run  the full simulation is provided below in lines 60-84. Un-comment to run. The following code, lines
#    29:55, runs the simulation for a small portion of the parameter space. This will still take a few minutes to run.


####### Sample run #######

### Calculate number of cores:
no_cores <- 2

### Initiate cluster:

cl <- makeCluster(no_cores)
registerDoParallel(cl)

clusterEvalQ(cl, c(library(tidyverse), library(data.table), library(mefa), library(gsDesign), library(survival)))

### Source function files and load rate data:
source('simulationFXNs.R')
source('analyzeFXN.R')

load('paho.Rdata')

### create parameter value set
params<-makeParms()
params <- params[params$regSize == 500,][1:9,] #keep a small portion of total parameter possibilities
params$iter <- 50 #reduce number of iterations to 50, from 250, to speed up computation.

results <- foreach(parms = iter(params, by='row')) %dopar% analyzeTrial(parms)

stopCluster(cl)

trialOut <- do.call('rbind', results)
trialOut <- cbind(params,trialOut)

view(trialOut)


######### Full simulation #########

# ### Calculate number of cores:
# 
# no_cores <- 2
# 
# ### Initiate cluster:
# 
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# 
# clusterEvalQ(cl, c(library(tidyverse), library(data.table), library(mefa), library(gsDesign), library(survival)))
# 
# source('simulationFXNs.R')
# source('analyzeFXN.R')
# 
# load('paho.Rdata')
# 
# params<-makeParms()
# 
# results <- foreach(parms = iter(params, by='row')) %dopar% analyzeTrial(parms)
# 
# stopCluster(cl)
# 
# trialOut <- do.call('rbind', results)
# trialOut <- cbind(params,trialOut)
# write.csv(trialOut, file='trialOut.csv')




## LICENSE
##
## This code is made available under a Creative Commons Attribution 4.0
## International License. You are free to reuse this code provided that you
## give appropriate credit, provide a link to the license, and indicate if
## changes were made.
## You may do so in any reasonable manner, but not in any way that suggests
## the licensor endorses you or your use. Giving appropriate credit includes
## citation of the above publication *and* providing a link to this repository:
##
## https://github.com/RachelMercaldo/ZikaTrial
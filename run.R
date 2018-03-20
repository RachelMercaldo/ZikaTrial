#!/usr/bin/env Rscript
library(foreach)
library(doParallel)
library(iterators)
library(gsDesign)
library(tidyverse)
library(data.table)
library(mefa)
library(survival)


### Calculate number of cores:

no_cores <- #

### Initiate cluster:

cl <- makeCluster(no_cores)
registerDoParallel(cl)

clusterEvalQ(cl, c(library(tidyverse), library(data.table), library(mefa), library(gsDesign), library(survival)))

source('simulationFXNs.R')
source('analyzeFXN.R')
load('paho.Rdata')

params<-makeParms()

results <- foreach(parms = iter(params, by='row')) %dopar% analyzeTrial(parms) 

stopCluster(cl)

trialOut <- do.call('rbind', results)
trialOut <- cbind(params,trialOut)
#write.csv(trialOut, file='trialOut.csv')


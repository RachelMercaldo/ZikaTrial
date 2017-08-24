library(tidyverse)
library(magrittr)
library(data.table)
library(mefa)
library(gmodels)
library(survival)

set.seed(628496)

PAHOdata<-read.csv("PAHOdata.csv") #Andersen lab data from 1/3/2016 thru 5/21/2017, 73 wks, https://github.com/andersen-lab/Zika-cases-PAHO
PAHOdata$date<-as.Date(PAHOdata$date,"%m/%d/%Y")
popSize<-c(472000,396000,3681000,8190000,128624000,6184000,209553000,48650000,276000,31518000)
PAHOdata[,2:11]<-(PAHOdata[,2:11]*4.5)/popSize

region <- colnames(PAHOdata[,-1])


makePop <- function(trial,regSize){
  trial %<>% group_by(region) %<>% nest(.key = 'pop')
  pop<-list()
  for(i in trial$region){
    ID = paste(substr(i, 1, 2),1:regSize,sep='')
    pop[i]<-data.frame(ID)
  }
  names(pop)<-paste(trial$region)
  trial$pop<-pop
  trial
}

#trial <- data.frame(region)
#trial<-makePop(trial,300)

getIndRR <- function(trial,regSize){
  for(reg in 1:length(trial$pop)){
    trial$pop[[reg]]<-as.data.frame(trial$pop[[reg]])
    trial$pop[[reg]]$indRR<-rlnorm(regSize,0,1)
  }
  trial
}

#trial<-getIndRR(trial,300)

randomize <- function(trial, regSize){
  for(reg in 1:length(trial$pop)){
    trial$pop[[reg]]$immuneStatus <- factor(rbinom(regSize, 1, prob = 0.500), labels = c('vaccine','control'))
    colnames(trial$pop[[reg]])<-c("ID","indRR","studyArm")
  }
  trial <- suppressWarnings(as.data.frame(unnest(trial))) #suppress bind_rows warning (coerce factor to character)
  trial
}

#trial<-randomize(trial,300)

mergeData<-function(paho, trial,regSize){   #merge PAHO rates and trial df
  paho<-rep(PAHOdata,regSize)
  paho<-gather(paho,'region','regRate',2:11)
  paho$region<-as.factor(paho$region)
  trial<-as.data.frame(rep(trial,each=nrow(PAHOdata)))
  trial<-cbind(trial,paho)
}

#trial<-mergeData(paho,trial,300)

totalRate<-function(trial,immuneDate = '2016-02-28', vaccEff){
  trial$totalRate<-ifelse(trial$date < as.Date(immuneDate),trial$indRR*trial$regRate,
                          trial$indRR*trial$regRate*ifelse(trial$studyArm=='vaccine',1-vaccEff,1))
  trial
}

simInf<-function(trial){
  trial<-trial[trial$totalRate != 0,]
  trial$infectTime<-rexp(nrow(trial),rate=trial$totalRate)
  trial
}


getSurvTime<-function(trial,immuneDate = '2016-02-28',endDate=max(trial$date)){
  trial<-trial[,c(5,1,2,3,4,7,8,9)] #reordering columns for pretty-ness and removing duplicate 'region'
   #remove rows with totalRate = 0
  
  preImmune<-trial[trial$date < as.Date(immuneDate),]
  preImmune<-preImmune[preImmune$infectTime<=1,]
  preImmune<-preImmune[!duplicated(preImmune$ID),]
  preImmune$survt<-(preImmune$date - as.Date(immuneDate)) + 7*preImmune$infectTime
  
  trial<-trial[!(trial$date < as.Date(immuneDate)) & !(trial$ID %in% preImmune$ID),] 
  
  infected<-trial[trial$infectTime<=1,]
  infected<-infected[!duplicated(infected$ID),] 
  infected$survt<-(infected$date - as.Date(immuneDate)) + 7*infected$infectTime
  
  notInfected<-trial[trial$infectTime>1 & !(trial$ID %in% infected$ID),]
  notInfected<-notInfected[!duplicated(notInfected$ID),] 
  notInfected$survt<-(as.Date(endDate)-as.Date(immuneDate))
 
  trial<-rbind(infected,notInfected)
  trial<-trial[order(trial$date,trial$region),]
  trial$survt<-as.numeric(trial$survt)
  
  preAndTrial<-list(pre=preImmune,trial=trial)
  preAndTrial
}


getSurvStatus<-function(trial){
  trial$status <- ifelse(trial$infectTime <= 1, 1, 0)
  trial
}


trial <- data.frame(region)
trial<-makePop(trial,300)
trial <- getIndRR(trial,300)
trial<-randomize(trial,300)
trial<-mergeData(paho,trial,300)
trial<-totalRate(trial,vaccEff=0.8)
trial<-simInf(trial)  

both<-getSurvTime(trial)
preImmune<-both$pre
trial<-both$trial
trial<-getSurvStatus(trial)

survdat<-trial[,c(3,5,9,10)]


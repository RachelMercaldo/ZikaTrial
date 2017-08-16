library(tidyverse)
library(magrittr)
library(GGally)
library(purrr)
library(data.table)
library(mefa)
library(gmodels)

PAHOdata<-read.csv("PAHOdata.csv") #Andersen lab data from 1/3/2016 thru 5/21/2017, 73 wks, https://github.com/andersen-lab/Zika-cases-PAHO
PAHOdata$date<-as.Date(PAHOdata$date,"%m/%d/%Y")
popSize<-c(472000,396000,3681000,8190000,128624000,6184000,209553000,48650000,276000,31518000)
PAHOdata[,2:11]<-(PAHOdata[,2:11]*4.5)/popSize

region <- colnames(PAHOdata[,-1])



makePop <- function(trial,regSize){
  pop<-list()
  for(i in trial$region){
    ID = paste(substr(i, 1, 2),1:regSize,sep='')
    pop[i]<-data.frame(ID)
  }
  names(pop)<-paste(trial$region)
  pop
}



getIndRR <- function(trial,regSize){
  for(reg in 1:length(trial$pop)){
    trial$pop[[reg]]<-as.data.frame(trial$pop[[reg]])
    trial$pop[[reg]]$indRR<-rlnorm(regSize,0,1)
  }
  trial
}


getImmunity <- function(trial, regSize){
  for(reg in 1:length(trial$pop)){
    trial$pop[[reg]]$immuneStatus <- factor(rbinom(regSize, 1, prob = 0.500), labels = c('vaccine','control'))
    colnames(trial$pop[[reg]])<-c("ID","indRR","studyArm")
  }
  trial
}


mergeData<-function(paho, trial,regSize){
  paho<-rep(PAHOdata,regSize)
  paho<-gather(paho,'region','rate',2:11)
  paho$region<-as.factor(paho$region)
  trial<-as.data.frame(rep(trial,each=nrow(PAHOdata)))
  trial<-cbind(trial,paho)
}


rexpRate<-function(trial,vaccEff){
  trial$rexpRates<-ifelse(trial$date<'2016-01-31',trial$indRR*trial$rate,
                          trial$indRR*trial$rate*ifelse(trial$studyArm=='vaccine',1-vaccEff,1))
  trial
}

getInfTimes<-function(trial){
  trial$infTime<-rexp(nrow(trial),rate=trial$rexpRates)
  trial
}


trial <- data.frame(region)
trial %<>% group_by(region) %<>% nest(.key = 'pop')
trial$pop<-makePop(trial,300)
trial <- getIndRR(trial,300)
trial<-getImmunity(trial,300)
trial <- as.data.frame(unnest(trial))
trial<-mergeData(paho,trial,300)
trial<-rexpRate(trial,0.8)
trial<-getInfTimes(trial)

infected<-na.omit(trial[trial$infTime<=1,])
infected<-infected[!duplicated(infected[,2]),]
CrossTable(infected$studyArm)
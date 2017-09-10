library(RCurl)
library(tidyverse)
library(magrittr)
library(data.table)
library(mefa)
library(gmodels)
library(survival)

#get the Andersen lab data from 1/3/2016 thru 5/21/2017, https://github.com/andersen-lab/Zika-cases-PAHO
fileURLs<-c("https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/caribbean.csv",
            "https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/central_america.csv",
            "https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/south_america.csv")
PAHOdata <- as.data.frame(unlist(lapply(fileURLs,fread), recursive=FALSE)) 

#selecting top-10 places to be to get Zika
PAHOdata <- subset(PAHOdata, select = c('French.Guiana','Brazil','Colombia',
                                        'Venezuela','Martinique','Puerto.Rico','Guadeloupe','Honduras',
                                        'Mexico','Nicaragua', 'susp.con.ZIKV.cases','V2')) 

#Andersen data is awesome but is terribly formatted:
PAHOdata <- setnames(PAHOdata[-1,], 11:12, c('year','date'))
PAHOdata[1:52,11]<-'2016'; PAHOdata[53:nrow(PAHOdata),11]<-2017  
PAHOdata$date <- as.Date(paste(PAHOdata$date,'-',PAHOdata$year, sep=""),"%d-%b-%Y")

#These are the population sizes for the top-10 countries, taken from PAHO 
##http://www.paho.org/hq/index.php?option=com_content&view=article&id=12390&Itemid=42090&lang=en 
popSize<-c(276000,209553000,48650000,31518000,396000,3681000,472000,8190000,128624000,6184000)

for(i in 1:10){
  PAHOdata[,i]<-PAHOdata[,i]*4.5/popSize[i]
}


PAHOdata<-PAHOdata[,c(12,1:10)] #getting rid of the extra "year" column

region <- colnames(PAHOdata[,2:11])

set.seed(628496)

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


getIndRR <- function(trial,regSize){
  for(reg in 1:length(trial$pop)){
    trial$pop[[reg]]<-as.data.frame(trial$pop[[reg]])
    trial$pop[[reg]]$indRR<-rlnorm(regSize,0,1)
  }
  trial
}


randomize <- function(trial, regSize){
  for(reg in 1:length(trial$pop)){
    trial$pop[[reg]]$immuneStatus <- factor(rbinom(regSize, 1, prob = 0.500), labels = c('vaccine','control'))
    colnames(trial$pop[[reg]])<-c("ID","indRR","studyArm")
  }
  trial <- suppressWarnings(as.data.frame(unnest(trial))) #suppress bind_rows warning (coerce factor to character)
  trial
}



mergeData<-function(paho, trial,regSize){   #merge PAHO rates and trial df
  paho<-rep(PAHOdata,regSize)
  paho<-gather(paho,'region','regRate',2:11)
  paho$region<-as.factor(paho$region)
  trial<-as.data.frame(rep(trial,each=nrow(PAHOdata)))
  trial<-cbind(trial,paho)
}


totalRate<-function(trial,startDate = min(trial$date), vaccEff){
  immuneDate = as.Date(startDate) + 4
  trial$totalRate<-ifelse(trial$date < as.Date(immuneDate),trial$indRR*trial$regRate,
                          trial$indRR*trial$regRate*ifelse(trial$studyArm=='vaccine',1-vaccEff,1))
  trial
}


simInf<-function(trial){
  trial<-trial[trial$totalRate != 0,] #removing weeks with rates = 0 to avoid NaN with rexp() in next step
  trial$infectTime<-rexp(nrow(trial),rate=trial$totalRate)
  trial
}


getSurvTime<-function(trial,startDate = min(trial$date), endDate=max(trial$date)){
  immuneDate = as.Date(startDate) + 28
  trial<-trial[,c(5,1,2,3,4,7,8,9)] #reordering columns for pretty-ness and removing duplicate 'region'
  
  preImmune<-trial[trial$date < as.Date(immuneDate),]
  preImmune<-preImmune[preImmune$infectTime<=1,]
  preImmune<-preImmune[!duplicated(preImmune$ID),]
  
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
  
  trial
}


getSurvStatus<-function(trial){
  trial$status <- ifelse(trial$infectTime <= 1, 1, 0)
  trial
}



simPower<-function(trial = trial, regSize=15,startDate = min(trial$date), vaccEff=0.80, iter=99,binMod=FALSE,coxMod=FALSE){
  immuneDate = as.Date(startDate) + 28
  pvec<-rep(NA,iter)
  for(i in 1:iter){
    trial<-data.frame(region)
    trial<-makePop(trial,regSize)
    trial<-getIndRR(trial,regSize)
    trial<-randomize(trial,regSize)
    trial<-mergeData(paho,trial,regSize)
    trial<-totalRate(trial,startDate, vaccEff)
    trial<-simInf(trial)  
    trial<-getSurvTime(trial,startDate)
    trial<-getSurvStatus(trial)
    
    if(binMod==TRUE){
      mod<-glm(status~studyArm,family=binomial(link='logit'),data=trial)
      pvec[i]<-coef(summary(mod))[2,4]
    }
    
    if(coxMod==TRUE){
      coxphMod<-coxph(Surv(survt,status)~studyArm, data=trial)
      pvec[i]<-summary(coxphMod)$logtest["pvalue"]
    }
  }
  return(mean(pvec<.05))
}


#Testing, testing:

simByPopSize<-function(regSize, startDate,binMod=FALSE,coxMod=FALSE){
  bySizeDF<-data.frame(RegionSize=rep(NA,length(regSize),Power=NA))
  for(r in 1:length(regSize)){
    bySizeDF$RegionSize[r] <- regSize[r]
    if(binMod==TRUE){
      bySizeDF$PowerBin[r] <- simPower(regSize = regSize[r], startDate = startDate,binMod=binMod,coxMod=coxMod)
    }
    if(coxMod==TRUE){
      bySizeDF$PowerCox[r] <- simPower(regSize = regSize[r], startDate = startDate,binMod=binMod,coxMod=coxMod)
    }
  }
  bySizeDF
}

regSize<-seq(15,50,by=5)

powerPop<-simByPopSize(regSize=regSize, startDate = '2016-01-03',binMod=TRUE,coxMod=TRUE)
plot(powerPop$RegionSize,powerPop$PowerBin,col="black",type='b')
points(powerPop$RegionSize,powerPop$PowerCox,col="red",type='b')


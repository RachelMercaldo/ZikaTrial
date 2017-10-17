library(tidyverse)
library(magrittr)
library(data.table)
library(mefa)
library(gsDesign)
library(survival)
set.seed(628496)

####  Prepping regional weekly Zika incidence data: 

#get the Andersen lab data from https://github.com/andersen-lab/Zika-cases-PAHO
#These files contain #cases by week in PAHO countries from 1/3/2016 thru 5/21/2017 (as of 9/2017)
#The Andersen lab occassionally updates the repository.

fileURLs<-c("https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/caribbean.csv",
            "https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/central_america.csv",
            "https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/south_america.csv")
PAHOdata <- as.data.frame(unlist(lapply(fileURLs,fread), recursive=FALSE)) 

#selecting top-10 places to be to get Zika. I did this in a few minutes in Excel in interest of time, 
#but will be writing some code up to have R do it (yay, reproducibility!)

PAHOdata <- subset(PAHOdata, select = c('Colombia', 'Ecuador', 'Mexico', 'Peru','susp.con.ZIKV.cases','V2')) 

#Andersen data is awesome but is badly formatted (don't tell them I said that...):

PAHOdata <- setnames(PAHOdata[-1,], c('susp.con.ZIKV.cases','V2'), c('year','date'))
PAHOdata[1:52,'year']<-'2016'; PAHOdata[53:nrow(PAHOdata),'year']<-2017  
PAHOdata$date <- as.Date(paste(PAHOdata$date,'-',PAHOdata$year, sep=""),"%d-%b-%Y")

#These are the population sizes for the top-10 countries, taken from PAHO reports, using avg 2016/2017
# http://www.paho.org/hq/index.php?option=com_content&view=article&id=12390&Itemid=42090&lang=en 
popSize<-c(48650000,16506000,128624000,31970000)


#Making it an infection per-capita by region by  multiplying 4.5 and dividing by region pop size
for(reg in 1:(ncol(PAHOdata)-2)){
  PAHOdata[,reg]<-PAHOdata[,reg]*4.5/popSize[reg]
}


PAHOdata<-select(PAHOdata, -year) #getting rid of the extra column for year now that we have date


#8 areas identified by the paper: 'Narino' in Colombia,'Sucumbios' in Ecuador,'Sinaloa' and 'Tamaulipas' in Mexico 
#'Piura','Tumbes','SanManrtin' and 'Ucayali' in Peru


PRF <- c(0.5931, 14.57, 8.73, 10.01, 27.52, 24.46, .3876, .2326) #projectionRateFactor


paho<-data.table(Narino = PAHOdata$Colombia*PRF[1], Sucumbios = PAHOdata$Ecuador*PRF[2], Sinaloa = PAHOdata$Mexico*PRF[3], 
                 Tamaulipas = PAHOdata$Mexico*PRF[4], Piura = PAHOdata$Peru*PRF[5], Tumbes = PAHOdata$Peru*PRF[6], 
                 SanMartin = PAHOdata$Peru*PRF[7], Ucayali = PAHOdata$Peru*PRF[8], date = PAHOdata$date)

####  Parameter list for your favorite trial scenario:
 
makeParms <- function(
  regSize = seq(40,100,by=1),  #regSize = pop size for one region. Can be a vector of possibilities. 8 regions means total of 8*regSize trial participants
  vaccEff = c(0.5,0.8), #Vaccine efficacy. Can be a vector of possibilities.
  symptomRate = 0.225,  
  incubMedian = 5.9,    #Lessler, 2016
  incubDisperse = 1.46, #Lessler, 2016
  startDate = '2016-01-03', #start of data, can change to look at starting trial later, can be a vector of possibilities
  CZS = TRUE, #if you want a trial that looks at CZS (in addition to infection/symptoms), select TRUE!
  CZSTrim1 = 0.13, #Eppes, 2017
  CZSTrim2 = 0.06, #Eppes, 2017
  CZSTrim3 = 0.06, #Eppes, 2017
  femaleOnly = c(TRUE, FALSE), #If CZS, are only females being recruited? Can be a vector of both possibilities.
  TTC = c(TRUE, FALSE), #If CZS, are women recruited Trying To Conceive? Can be a vector of both possibilities.
  pregRate = 0.075, #Taylor, 2003
  contraceptionRate = 0.73,   #UN 2015, Latin America and the Caribbean data for women 15-49
  coxMod = TRUE, 
  binMod = TRUE,
  gs = TRUE,
  iter = 100 
){
  return(as.list(environment()))
}



# Create the trial population with IDs based on region:
makePop <- function(parms) within(parms, {
  regions = colnames(select(paho, -date))
  trial<-data.table(id = paste(substr(rep(regions,each=regSize),1,2),1:regSize,sep=''))
})


#if interested in CZS, the following sets population sex characteristics as desired, 
# sets pregnancy rates for women either all trying to conceive (TTC) or not
CZSandPreg <- function(parms=makePop()) within(parms, {
  if(CZS){
    if(femaleOnly){
      trial$sex <- 'female'
    } else {
      trial$sex <- factor(rbinom(regSize, 1, prob = 0.500), labels = c('male','female'))
    }
    
    if(TTC){
      trial$TTC <- ifelse(trial$sex == 'female','yes','no')
      trial$pregRate<-ifelse(trial$TTC == 'yes',pregRate,0)
    } else {
      trial$TTC <- ifelse(trial$sex == 'female', sample(c('yes','no'), regSize, replace = TRUE, prob = c(1-contraceptionRate,contraceptionRate)),'no')
      trial$pregRate<-ifelse(trial$sex == 'female',ifelse(trial$TTC == 'yes',pregRate,0),0) #assuming contraception is 100% effective
    }
  }
})



#assigning individual RR for everyone.
getIndRR <- function(parms = CZSandPreg()) within(parms, {       
  trial$indRR<-rlnorm(nrow(trial),0,1)
})



#randomize all the people:
randomize <- function(parms = getIndRR()) within(parms, {
  trial$immuneStatus <-factor(rbinom(nrow(trial), 1, prob = 0.500), labels = c('vaccine', 'control'))
})


#### functions to merge with PAHO data and get total rates to simulate infection

#Merge data
mergeData <- function(parms = randomize()) within(parms, {
  trial<-rep(trial,each=nrow(paho))
  paho<-rep(paho,regSize)
  paho<-gather(paho,'region','regRate',c(`Narino`:`Ucayali`))
  trial<-cbind(trial,paho)
  
  immuneDate = as.Date(startDate) + 28  #similar to dengue?
  
  trial$totalRate<-ifelse(trial$date < as.Date(immuneDate),trial$indRR*trial$regRate,
                          trial$indRR*trial$regRate*ifelse(trial$immuneStatus=='vaccine',1-vaccEff,1))
})


#Simulate infection
simInf<-function(parms =mergeData()) within(parms, { 
  trial<-trial[trial$totalRate != 0,] #removing weeks with rates = 0 to avoid NaN with rexp() in next step  
  trial$infectTime<-rexp(nrow(trial),rate=trial$totalRate) 
  
  preImmune<-trial[trial$date < as.Date(immuneDate),]
  preImmune<-preImmune[preImmune$infectTime<=1,]
  
  trial<-trial[!(trial$date < as.Date(immuneDate)) & !(trial$id %in% preImmune$id),] 
  
  infected<-trial[trial$infectTime<=1,]
  infected<-infected[!duplicated(infected$id),] 
  infected$survt<-(infected$date - immuneDate) + 7*infected$infectTime
  
  notInfected<-trial[trial$infectTime>1 & !(trial$id %in% infected$id),]
  notInfected<-notInfected[!duplicated(notInfected$id),] 
  notInfected$survt<-max(trial$date) - immuneDate
  
  trial<-rbind(infected,notInfected)
  trial<-trial[order(trial$date,trial$region),]
  trial$survt<-as.numeric(trial$survt)
  
  trial$status <- ifelse(trial$infectTime <= 1, 1, 0)
})


#simulation for incubation period and symptoms.
symptomatic <- function(parms =simInf()) within(parms, { 
  trial$symptomatic <- ifelse(trial$status == 1, (rbinom(nrow(trial), 1, prob = symptomRate)), 0)
  trial$incubationPeriod <- ifelse(trial$symptomatic == 1, rlnorm(nrow(trial),log(incubMedian),log(incubDisperse)), Inf)
  trial$survtSymptoms <- ifelse(trial$symptomatic == 1, trial$incubationPeriod + trial$survt, Inf)
})



#### functions to simulate pregnancy and get pregnancy time
simPreg <- function(parms = symptomatic()) within(parms, {
  if(CZS){
    trial$pregTime <- ifelse(trial$sex == 'female', 7*(suppressWarnings(rexp(nrow(trial), trial$pregRate))), Inf) #suppress NA warning when rate = 0
    
    birthTime <- trial$pregTime + 280 #assuming pregnancy completed... 
    conceivedPreInf<-trial[(trial$status==1) & (trial$pregTime <= trial$survt),] 
    trial$conceivedPreInf<-ifelse(trial$id %in% conceivedPreInf$id,'yes','no')
  }
})



# Infected in first,second,third trimester?
getTrimesterInf <- function(parms = simPreg()) within(parms,{
  if(CZS){
    trial$trimesterInfected <- 'NotPregInf'
    for(preg in 1:nrow(trial)){
      trial[(trial$pregTime <= trial$survt) & (trial$survt <= (trial$pregTime + 92)), 'trimesterInfected'] <- 'first'
      trial[(trial$pregTime + 93) <= trial$survt & trial$survt <= (trial$pregTime + 184),'trimesterInfected'] <- 'second'
      trial[(trial$pregTime + 185) <= trial$survt & trial$survt <= (trial$pregTime + 300),'trimesterInfected'] <- 'third'
    }
  }
})

# CZS - yes/no, based on trimester infected
CZS <- function(parms = getTrimesterInf()) within(parms,{
  if(CZS){
    trial$czsProb <- NA
    for(preg in 1:nrow(trial)){
      if(trial[preg, 'trimesterInfected'] == 'First'){
        trial[preg, 'czsProb'] <- CZSTrim1
      } else if(trial[preg, 'trimesterInfected'] == 'Second'){
        trial[preg, 'czsProb'] <- CZSTrim2
      } else if(trial[preg, 'trimesterInfected'] == 'Third'){
        trial[preg, 'czsProb'] <- CZSTrim3
      } else {
        trial[preg,'czsProb'] <- 0
      }
    }
    trial$CZS <- suppressWarnings(rbinom(nrow(trial), 1, prob = trial$czsProb))
  }
})




##### Simulate multiple trials with various parameters (edit parm function above for desired parameters. 
##### Any parameter value can become a vector of possibilities you'd like to investigate (of course, that will 
##### make this whole thing take forever))


simTrials <- function(parms = makeParms()) with(parms, { 
  params <- expand.grid(parms)
  params <- params[10:12,]
  data<-data.table(parms = NA, infections = NA, powerCox = NA, powerBin = NA)
    for(p in 1:nrow(params)){
      pvecCox <- rep(NA,params[p,'iter'])
      pvecBin <- rep(NA,params[p,'iter'])
      powerCox <- NA
      powerBin <- NA
      infections <- rep(NA, params[p,'iter'])
      
      for(ii in 1:params[p,'iter']){
        parms = c(params[p,])
        parms = makePop(parms)
        parms = CZSandPreg(parms)
        parms = getIndRR(parms)
        parms = randomize(parms)
        parms = mergeData(parms)
        parms = simInf(parms)
        parms = symptomatic(parms)
        parms = simPreg(parms)
        parms = getTrimesterInf(parms)
        parms = CZS(parms)
        
        infectionMod <- try(coxph(Surv(parms$trial$survt,parms$trial$status==1) ~ parms$trial$immuneStatus, frailty(parms$trial$region, distribution="gamma",sparse=FALSE)),silent=T)
        
        useInfectionMod <-
          !inherits(infectionMod,
                    'try-error')
        
        pvecCox[ii]<-summary(infectionMod)$logtest["pvalue"]
        
        
        powerCox <- mean(pvecCox < 0.05)
        
        
        binMod<-glm(status~immuneStatus,family=binomial(link='logit'),data=parms$trial)
        pvecBin[ii]<-coef(summary(binMod))[2,4]
        
        
        powerBin <- mean(pvecBin < 0.05)
        
        infections[ii] <- sum(parms$trial$status == 1) 
        
      }
      trialOut <- data.table(parms = list(params[p,]), 
                             infections = list(infections), 
                             powerCox = powerCox, 
                             powerBin = powerBin)
      
      data <- rbind(data,trialOut)
    }
    data[-1,]
})

infData<-simTrials(makeParms())

save(infData, file = "simTrialsOutput.Rdata") #save data.table for later plotting.


#GS Design stuff
#try for n.fix:

n.fix<-nBinomial(p1=0.0220623, p2=0.00441246, beta = 0.2)

#n.fix=1316, 658 per arm (for a trial lasting 1 year)

gsArgs <- list(k=4, test.type = 2, n.fix = n.fix, beta=0.2)
gsBounds <- do.call(gsDesign, gsArgs)
plot(gsBounds)

gsBounds

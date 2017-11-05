library(tidyverse)
library(magrittr)
library(data.table)
library(mefa)
library(gsDesign)
library(survival)
library(foreach)
library(iterators)
set.seed(628496)

#Prepping regional weekly Zika incidence data: 

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
popSize<-c(48650000,16506000,130624000,31970000)


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



##### Simulation


makeParms <- function(
  infTrial = TRUE,
  regSize = c(50,100,150,200),  #regSize = pop size for one region. Can be a vector of possibilities. 8 regions means total of 8*regSize trial participants
  vaccEff = c(0.5,0.8), #Vaccine efficacy. Can be a vector of possibilities.
  symptomRate = 0.225,  
  incubMedian = 5.9,    #Lessler, 2016
  incubDisperse = 1.46, #Lessler, 2016
  startDate = '2016-01-03', #start of data, can change to look at starting trial later, can be a vector of possibilities
  CZS = TRUE, #if you want a trial that looks at CZS (in addition to infection/symptoms), select TRUE!
  CZSTrim1 = 0.13, #Eppes, 2017
  CZSTrim2 = 0.06, #Eppes, 2017
  CZSTrim3 = 0.06, #Eppes, 2017
  femaleOnly = TRUE, #If CZS, are only females being recruited? Can be a vector of both possibilities.
  TTC = TRUE, #If CZS, are women recruited Trying To Conceive? Can be a vector of both possibilities.
  pregRate = 0.075, #Taylor, 2003
  contraceptionRate = 0.73,   #UN 2015, Latin America and the Caribbean data for women 15-49
  gs = c(TRUE,FALSE),
  assumedRate = 0.00238,
  iter = 10 
){
  return(as.list(environment())) 
}


#study population fxn:
#make a study population, assign individual RR, and randomize:
makePop <- function(parms = makeParms()) with(parms, {
  regions = colnames(paho[,`Narino`:`Ucayali`])  #get region names from parameters
  trial<-data.table(id = paste(substr(rep(regions,each=regSize),1,2),1:regSize,sep=''))  #assign IDs based on region
  trial$indRR <- rlnorm(nrow(trial),0,1)    #give everyone an individual RR   
  trial$arm <- sample(c('vaccine','control'), nrow(trial), replace = TRUE, prob = c(0.5,0.5)) #randomize to vaccine or control arm. 
  #Should we do this all at once for everyone 
  #in the whole trial or randomize within regions?
  trial
})

#CZS population fxn:
#if interested in CZS, this sets sex characteristics and pregnancy rates for women trying to conceive (TTC) or not:

CZSandPreg <- function(trial,parms) with(parms, {
  if(CZS){
    if(femaleOnly){
      trial$sex <- 'female'
    } else {
      trial$sex <- factor(rbinom(nrow(trial), 1, prob = 0.500), labels = c('male','female'))
    }
    
    if(TTC){
      trial$TTC <- ifelse(trial$sex == 'female','yes','no')
      trial$pregRate<-ifelse(trial$TTC == 'yes',pregRate,0)
    } else {
      trial$TTC <- ifelse(trial$sex == 'female', sample(c('yes','no'), nrow(trial), replace = TRUE, prob = c(1-contraceptionRate,contraceptionRate)),'no')
      trial$pregRate<-ifelse(trial$sex == 'female',ifelse(trial$TTC == 'yes',pregRate,0),0) #assuming contraception is 100% effective
    }
  }
  trial
})

#Merge PAHO data fxn:
mergeData <- function(trial, parms) with(parms, {
  trial<-rep(trial,each=nrow(paho))
  paho<-rep(paho,regSize)
  paho<-gather(paho,'region','regRate',c(`Narino`:`Ucayali`))
  trial<-cbind(trial,paho)
  
  immuneDate = as.Date(startDate) + 30  #similar to dengue?
  
  trial$totalRate<-ifelse(trial$date < as.Date(immuneDate),trial$indRR*trial$regRate,
                          trial$indRR*trial$regRate*ifelse(trial$arm=='vaccine',1-vaccEff,1))
  trial
})

#Simulate infection fxn:
simInf<-function(trial,parms) with(parms, { 
  immuneDate <- as.Date(startDate) + 30
  
  trial<-trial[trial$totalRate != 0,] #removing weeks with rates = 0 to avoid problems with rexp() in next step (does not affect anything. I checked twice) :D  
  trial$infectTime<-rexp(nrow(trial),rate=trial$totalRate) 
  
  preImmune<-trial[trial$date < as.Date(immuneDate),] #if desired, we could look more into this. Not sure if we care. Can always look at total rows of trial at the end and 
  # see how it is different from the total we had at the beginning of this step. That would be the # infected before immunity sets in.
  preImmune<-preImmune[preImmune$infectTime<=1,]
  
  trial<-trial[!(trial$date < as.Date(immuneDate)) & !(trial$id %in% preImmune$id),] #removes all the weeks prior to immunity and all the participants who were infected before immunity
  
  infected<-trial[trial$infectTime<=1,]
  infected<-infected[!duplicated(infected$id),]  #takes first instance of infectTime <= 1 for any id
  infected$survt<-(infected$date - immuneDate) + 7*infected$infectTime
  
  notInfected<-trial[trial$infectTime>1 & !(trial$id %in% infected$id),]
  notInfected<-notInfected[!duplicated(notInfected$id),] 
  notInfected$survt<- Inf
  
  trial<-rbind(infected,notInfected)
  trial<-trial[order(trial$date,trial$region),]
  trial$survt<-as.numeric(trial$survt)
  
  trial$status <- ifelse(trial$infectTime <= 1, 1, 0)
  
  trial
})

#fxn for incubation period and symptoms.
symptomatic <- function(trial,parms) with(parms, { 
  trial$symptomatic <- ifelse(trial$status == 1, (rbinom(nrow(trial), 1, prob = symptomRate)), 0)
  trial$incubationPeriod <- ifelse(trial$symptomatic == 1, rlnorm(nrow(trial),log(incubMedian),log(incubDisperse)), Inf) #I *think* I was supposed to do the log(median), log(dispersion) thing. I *think*
  trial$survtSymptoms <- ifelse(trial$symptomatic == 1, trial$incubationPeriod + trial$survt, Inf)
  trial
})


#fxn to simulate pregnancy and get pregnancy time
simPreg <- function(trial,parms) with(parms, {
  if(CZS){
    trial$pregTime <- ifelse(trial$sex == 'female', 7*(suppressWarnings(rexp(nrow(trial), trial$pregRate))), Inf) #suppress NA warning when rate = 0 (males or not TTC)
    
    trial$birthTime <- trial$pregTime + 280 #assuming pregnancy completed... and that it was completed in a "normal" amount of time.
    
    conceivedPreInf<- trial[(trial$status==1) & (trial$pregTime <= trial$survt),] 
    trial$conceivedPreInf<-ifelse(trial$id %in% conceivedPreInf$id,'yes','no')
    
    trial$trimesterInfected <- 'NotPregInf'
    
    for(preg in 1:nrow(trial)){
      trial[(trial$pregTime <= trial$survt) & (trial$survt <= (trial$pregTime + 92)), 'trimesterInfected'] <- 'First'
      trial[(trial$pregTime + 93) <= trial$survt & trial$survt <= (trial$pregTime + 184),'trimesterInfected'] <- 'Second'
      trial[(trial$pregTime + 185) <= trial$survt & trial$survt <= (trial$pregTime + 280),'trimesterInfected'] <- 'Third'
    }
  }
  trial
})



# CZS - yes/no, based on trimester infected
CZS <- function(trial,parms) with(parms,{
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
  trial
})


analyzeTrial <- function(parms, browse=F) with(parms, {
  if(browse) browser()
  if(infTrial){
    gsDesArgs = list(k=5, 
                     test.type = 2,
                     alpha = 0.025,
                     beta = 0.1,
                     timing = seq(0, 1, l = 6)[-1])
    gsPlan = do.call(gsDesign, gsDesArgs)
    cumRates = assumedRate*c(1,(1-vaccEff))
    fixedSamp = nBinomial(p1 = cumRates[1], p2 = cumRates[2], outtype = 2, beta = 0.1)
    fixedInfs = sum(fixedSamp*cumRates)
    
    contInfs = vaccInfs = date = pvec = parameters = as.numeric(list(NA))
    
    for(i in 1:iter){
      trial<-makePop(parms)
      trial<-CZSandPreg(trial,parms)
      trial<-mergeData(trial,parms)
      trial<-simInf(trial,parms)
      trial<-symptomatic(trial,parms)
      trial<-simPreg(trial,parms)
      trial<-CZS(trial,parms)
      
      infDays = sort(trial$date[trial$status == 1])
      
      if(gs){
        trialIter = data.table(events = round(gsPlan$timing*fixedInfs))
        
        trialIter[,tcal := infDays[events]]
        trialIter <- trialIter[!is.na(tcal)]
        trialIter$trigger = 'events'
        
        if(nrow(trialIter) < gsPlan$k){
          trialIter = rbind(trialIter, data.table(events = NA, tcal = max(trial$date), trigger = 'end time'))
        }
        
        if(nrow(trialIter) == gsPlan$k){
          trialIter = cbind(trialIter, nominalP = pnorm(gsPlan$lower$bound)) 
        }else{
          gsDesArgsAdj = within(gsDesArgs, {
            k = nrow(trialIter)
            timing = seq(0, 1, l = k + 1)[-1]
          })
          if(gsDesArgsAdj$k > 1){
            gsPlanAdj = do.call(gsDesign, gsDesArgsAdj)
            trialIter = cbind(trialIter, nominalP = pnorm(gsPlanAdj$lower$bound))
          }else{
            trialIter = data.table(events = NA, tcal = max(trial$date), trigger = 'end time', nominalP = 0.025)
          }
        }
      }else{
        trialIter = data.table(events = NA, tcal = max(trial$date), trigger = 'end time', nominalP = 0.025)
      }
      
      analysisNum = 0
      trialStopped = FALSE
      
      while(!trialStopped){
        analysisNum = analysisNum + 1
        analysisDate = trialIter[analysisNum, tcal]
        censTrial = as.data.table(trial)
        censTrial[date > analysisDate, date := Inf]
        
        infectionMod <- try(coxph(Surv(rep(0,nrow(censTrial)), censTrial$survt, censTrial$status == 1) ~ censTrial$arm == 'vaccine', 
                                  frailty(censTrial$region, distribution = 'gamma', sparse = FALSE)), silent = TRUE)
        useInfectionMod <- !inherits(infectionMod, 'try-error')
        
        trialIter[analysisNum, rawPinfection := summary(infectionMod)$logtest['pvalue']]
        trialIter[analysisNum, logHRinfection := summary(infectionMod)$coefficients[,1]]
        
        trialIter[analysisNum, vaccGoodInfection := rawPinfection < ifelse(gs, nominalP, .025) & logHRinfection < 0]
        trialIter[analysisNum, vaccBadInfection := rawPinfection < ifelse(gs, nominalP, .025) & logHRinfection > 0]
        
        trialIter$vaccCases[analysisNum] <- nrow(censTrial[arm=='vaccine' & survt != Inf & date < trialIter$tcal[analysisNum],])
        trialIter$contCases[analysisNum] <- nrow(censTrial[arm=='control' & survt != Inf & date < trialIter$tcal[analysisNum],])
        
        earlyStop <- trialIter[analysisNum, vaccGoodInfection | vaccBadInfection]
        
        if(any(earlyStop, na.rm=T)) trialStopped <- T
        if(analysisNum==nrow(trialIter)) trialStopped <- T
      }
      dat<-trialIter[max(nrow(trialIter)), `tcal`]
      
      dat<-as.character(dat[[1]][1])
      
      contInfs[i] <- trialIter[max(nrow(trialIter)),`contCases`]
      vaccInfs[i] <- trialIter[max(nrow(trialIter)), `vaccCases`]
      pvec[i] <- trialIter[max(nrow(trialIter)), `rawPinfection`]
      date[i] <- dat
    }
  }
  out<-data.table(parameters = list(parms),
                   contInfsMean = mean(contInfs),
                   contInfsMedian = median(contInfs),
                   vaccInfsMean = mean(vaccInfs),
                   vaccInfsMedian = median(vaccInfs),
                   power = mean(pvec < 0.05),
                   meanDate = mean(as.Date(date)))
  out
})

params<- makeParms()
params<- expand.grid(params)

results <- foreach(parms = iter(params, by='row')) %do% analyzeTrial(parms) 
trialOut <- do.call('rbind', results)
trialOut <- cbind(params,trialOut)


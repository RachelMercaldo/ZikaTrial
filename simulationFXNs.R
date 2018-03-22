#Simulating trial - functions

makeParms <- function(
  trialType = c('CZStrial','infTrial','symptomTrial'),
  regSize = seq(10,1500, by = 5), #regSize = pop size for one region. 8 regions means total of 8*regSize trial participants
  vaccEff = c(.5,.7,.8,.9), #Vaccine efficacy. 
  startDate = c('2016-01-03','2016-04-03','2016-04-17', '2016-05-15','2016-06-19','2016-07-03', '2016-08-14', '2016-09-18', '2016-10-09', '2016-11-13', '2016-11-27'), #startdate of data
  symptomRate = 0.225,
  incubK = 3.132, #Lessler
  incubC = 6.632, #Lessler
  persistBloodK = 2.007, #Lessler
  persistBloodC = 11.171, #Lessler
  maxEnddate = '2017-05-21', #last week of infection rates
  maxCZSdate = '2018-02-25', #last possible birthdate for women pregnant by maxEndDate.
  testInterval = c(7,14,28), #how many days separate lab tests
  CZSTrim1 = 0.11, #Eppes, 2017
  CZSTrim2 = 0.06, #Eppes, 2017
  CZSTrim3 = 0.06, #Eppes, 2017
  TTC = c(TRUE,FALSE), #If CZS, are women recruited Trying To Conceive (TTC)? If FALSE, assumes 1-contraceptionRate for prevalence of women TTC
  lastPregRate = 0.001, #Taylor, 2003, gave .05 for end of 12 months. Assumption that it would be much lower, but not 0, for 18 months 
  firstPregRate = 0.30, #Taylor, 2003
  contraceptionRate = 0.73, #UN 2015, Latin America and the Caribbean data for women 15-49
  assumedRate = 0.00238, #Assumed rate of infection, for group sequential trial design
  assumedRateSymptoms = 0.00238*0.225, #Assumed rate of symptoms, for group sequential trial design
  assumedRateCZS = 0.00238*0.06, #Assumed rate of CZS, for group sequential trial design
  iter = 500 
){
  parms <- expand.grid(as.list(environment()))
  parms<-parms[!(parms$trialType != 'CZStrial' & (parms$TTC == TRUE)),]
  parms<-parms[!(parms$trialType !='infTrial' & (parms$testInterval %in% c(14,28))),]
  parms<-parms[!(parms$trialType != 'CZStrial' & (parms$startDate %in% c('2016-04-03','2016-07-03','2016-04-17', '2016-05-15','2016-06-19', '2016-08-14', '2016-09-18', '2016-10-09', '2016-11-13', '2016-11-27'))),]
  parms<-parms[!(parms$trialType != 'CZStrial' & (parms$regSize > 1000)),]
}
  
# parms<-makeParms()
# parms<-parms[1,]
# parms$iter<-50
# parms

#Make a study population, assign individual risks, and randomize to vaccine or control:
makePop <- function(paho, parms = makeParms()) with(parms, {
  regions<-colnames(paho)
  regions<-regions[1:8]  #get region names from paho data
  trial<-data.table(id = paste(substr(rep(regions,each=regSize),1,2),1:regSize,sep=''))  #assign unique IDs based on region
  trial$indRR <- rlnorm(nrow(trial),0,1)    #give everyone an individual risk  
  trial$arm <- sample(c('vaccine','control'), nrow(trial), replace = TRUE, prob = c(0.5,0.5)) #randomize to vaccine or control arm. 
  trial
})

# trial<-makePop(paho,parms)
# head(trial)

#Merge PAHO data:
#Following this function, each participant will have n rows, with n determined by the number of weeks in the trial (dependent on startDate)
#For each week, the participant is assigned a risk (totalRate) that is the product of their individual risk and their regional risk that week
mergeData <- function(trial, paho, parms) with(parms, {
  trial<-rep(trial,each=nrow(paho))
  paho<-rep(paho,regSize)
  paho<-gather(paho,'region','regRate',c(1:8))
  trial<-cbind(trial,paho)
  
  immuneDate = as.Date(startDate) + 30  #assuming 1 month until vaccine is protective
  
  trial$totalRate<-ifelse(trial$date < as.Date(immuneDate),trial$indRR*trial$regRate,
                          trial$indRR*trial$regRate*ifelse(trial$arm=='vaccine',1-vaccEff,1))
  trial
})

# trial<-mergeData(trial,paho,parms)
# head(trial)

#If a CZS trial, this gets conception time by assigning each woman monthly conception 'risk' and simulating time-to-conception
simPreg <- function(trial,parms, browse=F) with(parms, {
  if(browse) browser()
  if(trialType == 'CZStrial'){
    immuneDate <- as.Date(startDate) + 30
    mos<-as.numeric(as.Date(maxEnddate)-as.Date(startDate))
    mosLengthOut<-mos/28
    cycleProbs<-rev(seq(lastPregRate,firstPregRate,length.out = 19)) 
    cycleProbs<-cycleProbs[1:mosLengthOut]
    
    pregTrial<-trial[!duplicated(trial$id),]
    pregTrial<-rep(pregTrial, each =mosLengthOut)
    pregTrial$month <- 1:mosLengthOut
    pregTrial$cycleProbs<-cycleProbs
    
    if(TTC){
      pregTrial$TTC<-1
      pregTrial$pregTime<-ifelse(pregTrial$TTC==1, rexp(nrow(pregTrial),rate=pregTrial$cycleProbs),Inf)
    } else {
      pregTrial$TTC <- rbinom(nrow(pregTrial),1,prob=(1-contraceptionRate))  #assuming contraception is 100% effective
      pregTrial$pregTime<-ifelse(pregTrial$TTC==1, rexp(nrow(pregTrial),rate=pregTrial$cycleProbs),Inf)
    }
    
    pregTrial<-pregTrial[pregTrial$pregTime<=1,]
    pregTrial<-pregTrial[!duplicated(pregTrial$id),] #takes first instance of pregTime <= 1 for any id.
    pregTrial$conceptionTime<-pregTrial$pregTime*28 + pregTrial$month*28
    
    pregTrial<-pregTrial[,c('id','conceptionTime')]
    trial<-merge(trial,pregTrial, by = "id", all.x = TRUE)
    trial<-trial[order(trial$id, trial$date),]
  }
  trial
})

# trial<-simPreg(trial,parms)
# head(trial)


#Simulate infection. Using weekly total risk, simulate time-to-infection. 
#Following this, individuals infected before the immune date (30 days post start) are removed,
#  and all remaining individuals have survival times until first week infectTime <= 1, or Inf for those uninfected
simInf<-function(trial,parms) with(parms, { 
  immuneDate <- as.Date(startDate) + 30 #assuming 1 month until vaccine is protective
  
  trial<-trial[trial$totalRate != 0,] #removing weeks with rates = 0 to avoid problems with rexp() in next step 
  trial$infectTime<-rexp(nrow(trial),rate=trial$totalRate) 
  
  #identify infections prior to protective immunity.
  preImmune<-trial[trial$date < as.Date(immuneDate),] 
  preImmune<-preImmune[preImmune$infectTime<=1,]
  
  trial<-trial[!(trial$date < as.Date(immuneDate)) & !(trial$id %in% preImmune$id),] 
  #removes all the weeks prior to immunity and all the participants who were infected before immunity
  
  infected<-trial[trial$infectTime<=1,]
  infected<-infected[!duplicated(infected$id),]  #takes first instance of infectTime <= 1 for any id
  infected$survt<-(infected$date - immuneDate) + 7*infected$infectTime
  
  notInfected<-trial[trial$infectTime>1 & !(trial$id %in% infected$id),]
  notInfected<-notInfected[!duplicated(notInfected$id),] 
  notInfected$survt<- Inf
  
  trial<-rbind(infected,notInfected)
  trial<-trial[order(trial$id,trial$date),]
  trial$survt<-as.numeric(trial$survt)
  
  trial$status <- ifelse(trial$infectTime <= 1, 1, 0)
  
  trial
})

# trial<-simInf(trial,parms)
# head(trial)


#CZS, based on trimester infected. 
getCZSoutcome <- function(trial,parms,browse=F) with(parms, {
  if(browse) browser()
  if(trialType == 'CZStrial'){
    
    #If conception time is NA, replace with INF, then calculate birthTime and transform to date:
    trial$conceptionTime <- replace(trial$conceptionTime, is.na(trial$conceptionTime), Inf)
    trial$birthTime <- trial$conceptionTime + 280
    trial$birthDate <- as.Date(startDate) + trial$birthTime
    
    trial<-as.data.table(trial)
    trial$trimesterInfected <- as.character(Inf)  #keeping a record of trimester infected. Could skip, but useful for debugging
    
    setDT(trial)[conceptionTime <= survt & survt <= (conceptionTime + 93), trimesterInfected := 'First']
    setDT(trial)[(conceptionTime + 94) <= survt & survt <= (conceptionTime + 184), trimesterInfected := 'Second'] 
    setDT(trial)[(conceptionTime + 185) <= survt & survt <= (conceptionTime + 280),trimesterInfected := 'Third']
    suppressWarnings(setDT(trial)[conceptionTime == Inf | survt == Inf, trimesterInfected := Inf])
    
    trial$czsProb <- 0
    
    trial[trimesterInfected == 'First', czsProb := CZSTrim1]
    trial[trimesterInfected == 'Second', czsProb := CZSTrim2]
    trial[trimesterInfected == 'Third', czsProb := CZSTrim3]
    
    trial$CZS <- suppressWarnings(rbinom(nrow(trial), 1, prob = trial$czsProb)) #CZS yes/no
  }
  trial
})

# trial<-getCZSoutcome(trial,parms)
# head(trial)

#get symptomatic cases: generate incubation period and calculate time-to-symptoms
symptomatic <- function(trial,parms) with(parms, { 
  trial$symptomatic <- ifelse(trial$status == 1, (rbinom(nrow(trial), 1, prob = symptomRate)), 0) 
  trial$incubationPeriod <- ifelse(trial$symptomatic == 1, rweibull(nrow(trial),incubK, incubC), Inf) 
  trial$survtSymptoms <- ifelse(trial$symptomatic == 1, trial$incubationPeriod + trial$survt, Inf) 
  trial
})

trial<-symptomatic(trial,parms)
head(trial)

#quick interval function for persistence sim:
in_interval <- function(x, lower, upper){
  lower <= x & x <= upper
}

#test results for zikv, given viral persistence
persistence <- function(trial,parms,browse = F) with(parms, {
  if(browse) browser()
  Dates <- seq.Date(as.Date(startDate),as.Date(maxEnddate), by = testInterval) #which dates will be lab test dates
  testDates <- list()
  for(dat in 1:length(Dates)){   #in terms of number of days since trial start:
    testDates[dat] <- Dates[dat] - as.Date(startDate)
  }
  
  trial$firstDetectableBlood<-ifelse(trial$symptomatic == 1 & trial$survtSymptoms < (trial$survt + 2), trial$survtSymptoms, trial$survt+2) #assume detectable after 2 days
  trial$lastDetectableBlood<-rweibull(nrow(trial), persistBloodK, persistBloodC) + trial$survt
  
  trial$testResult <- 0
  
  for(ind in 1:nrow(trial)){   #evil for-loop. If any of the test dates lies within the RNA detection days, testResult = 1.
    if(trial[ind,'status'] == 1){
    trial[ind,'testResult']<-ifelse(any(in_interval(unlist(testDates),as.numeric(trial[ind,'firstDetectableBlood']),as.numeric(trial[ind,'lastDetectableBlood']))),1,0)
    }else{trial[ind,'testResult'] <- 0}
  }
  #the above uses symptoms only to determine beginning of detection interval. If it's a symptom trial, however, symptoms must be present for testResult=1.
  if(trialType == 'symptomTrial'){
    trial$testResult<-ifelse(trial$testResult == 1 & trial$symptomatic == 1, 1, 0) 
  }
  trial
})

# trial<-persistence(trial,parms)
# trial[trial$status==1,]   #take a look at everyone infected, regardless of testResult

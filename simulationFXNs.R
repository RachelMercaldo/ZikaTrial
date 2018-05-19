#Simulating trial - functions

makeParms <- function(
  trialType = c('CZStrial','infTrial','symptomTrial'),
  regSize = seq(10,1500, by = 5), #regSize = pop size for one region. 8 regions means total of 8*regSize trial participants
  vaccEff = c(.5,.7,.8,.9), #Vaccine efficacy. 
  startDate = c('2016-01-03',
                '2016-04-10',
                '2016-05-15',
                '2016-06-19',
                '2016-07-17',
                '2016-08-14',
                '2016-09-18'), #startdate of data
  symptomRate = 0.225, #Petersen, 2016
  incubK = 3.132, #Lessler
  incubC = 6.632, #Lessler
  persistBloodK = 2.007, #Lessler, 2016
  persistBloodC = 11.171, #Lessler, 2016
  maxEnddate = '2017-05-21', #last week of infection rates
  maxCZSdate = '2018-02-25', #last possible birthdate for women pregnant by maxEndDate.
  testInterval = c(7,14,28), #how many days separate lab tests
  CZSTrim1 = 0.15, #Reynolds, 2017
  CZSTrim2 = 0.0227, #Reynolds, 2017 
  CZSTrim3 = 0.0227, #Reynolds, 2017
  multiCZSrate = 0.21, #Reynolds, 2017
  TTC = c(TRUE,FALSE), #If CZS, are women recruited Trying To Conceive (TTC)? If FALSE, assumes 1-contraceptionRate for prevalence of women TTC
  startPregRate = 0.2126546, #Taylor, 2003, estimated 0.3 for first cycle, 0.05 for 12th. alpha = 0.3, beta = 0.861299. Mean of first 6 months for start.
  contraceptionRate = 0.73, #UN 2015, Latin America and the Caribbean data for women 15-49
  assumedRate = 0.0005, #Assumed rate, for group sequential trial design
  iter = 250 
){
  parms <- expand.grid(as.list(environment()))
  parms<-parms[!(parms$trialType !='infTrial' & (parms$testInterval %in% c(14,28))),]
  parms<-parms[!(parms$trialType != 'CZStrial' & (parms$startDate %in% c('2016-04-10','2016-05-15','2016-06-19', '2016-07-17','2016-08-14', '2016-09-18'))),]
  parms<-parms[!(parms$trialType != 'CZStrial' & (parms$regSize > 1000)),]
}


#create study population, assign individual risks, and randomize to vaccine or control:
makePop <- function(paho, parms = makeParms()) with(parms, {
  regions<-colnames(paho)
  regions<-regions[1:8]  #get region names from paho data
  trial<-data.table(id = paste(substr(rep(regions,each=regSize),1,2),1:regSize,sep=''))  #assign unique IDs based on region
  trial$indRR <- rlnorm(nrow(trial),0,1)    #give everyone an individual risk  
  trial$arm <- sample(c('vaccine','control'), nrow(trial), replace = TRUE, prob = c(0.5,0.5)) #randomize to vaccine or control arm. 
  trial
})


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


#If a CZS trial, this gets conception time by assigning each woman monthly conception 'risk' and simulating time-to-conception
simPreg <- function(trial,parms) with(parms, {
  if(trialType == 'CZStrial'){
    immuneDate <- as.Date(startDate) + 30
    mos<-as.numeric(as.Date(maxEnddate)-as.Date(startDate))
    mosLengthOut<-mos/28
    p<-cyclePs(parms)
    cycleProbs<-p[1:mosLengthOut]
    
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


#Using weekly total risk, simulate time-to-infection. 
#Individuals infected before the immune date (30 days post start) are removed, and all remaining 
#individuals have survival times until first week infectTime <= 1, or Inf for those uninfected
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
  infected$dup<-ifelse(duplicated(infected$id),1,0) #identify multiple infections
  infected$dup<-as.numeric(infected$dup) 
  dups<-infected[infected$dup==1,] #pull out participants with multiple infections
  multiInfect<-infected[infected$id %in% dups$id,]
  
  if(nrow(multiInfect) > 0) multiInfect$dup<-1
  multiInfect$survt<-as.numeric((multiInfect$date-immuneDate) + 7*multiInfect$infectTime)
  
  infected<-infected[!duplicated(infected$id) & !(infected$id %in% multiInfect$id),] #takes first instance of infectTime <= 1 for any id
  infected$survt<-(infected$date - immuneDate) + 7*infected$infectTime
  
  notInfected<-trial[trial$infectTime>1 & !(trial$id %in% infected$id) & !(trial$id %in% multiInfect$id),]
  notInfected<-notInfected[!duplicated(notInfected$id),] 
  notInfected$dup<-0
  notInfected$survt<- Inf
  
  trial<-rbind(infected,notInfected)
  trial<-rbind(trial,multiInfect)
  
  trial<-trial[order(trial$id,trial$date),]
  trial$survt<-as.numeric(trial$survt)
  
  trial$status <- ifelse(trial$infectTime <= 1, 1, 0)
  
  trial
})


#CZS, based on trimester infected. 
getCZSoutcome <- function(trial,parms) with(parms, {
  if(trialType == 'CZStrial'){
    
    #If conception time is NA, replace with INF, then calculate birthTime and transform to date:
    trial$conceptionTime <- replace(trial$conceptionTime, is.na(trial$conceptionTime), Inf)
    trial$birthTime <- trial$conceptionTime + 280
    trial$birthDate <- as.Date(startDate) + trial$birthTime
    
    trial<-as.data.table(trial)
    trial$trimesterInfected <- Inf #keeping a record of trimester infected. Could skip, but useful for debugging
    
    setDT(trial)[conceptionTime <= survt & survt <= (conceptionTime + 93), trimesterInfected := 1]
    setDT(trial)[(conceptionTime + 94) <= survt & survt <= (conceptionTime + 184), trimesterInfected := 2] 
    setDT(trial)[(conceptionTime + 185) <= survt & survt <= (conceptionTime + 280),trimesterInfected := 3]
    suppressWarnings(setDT(trial)[conceptionTime == Inf | survt == Inf, trimesterInfected := Inf])
    
    trial$czsProb <- 0
    
    multiInfect<-trial[trial$dup==1,]
    multiInfect<-multiInfect[,{x=unique(trimesterInfected)
    nnn=x[is.finite(x)]
    nn=length(nnn)
    list(n=nn)},by=id]
    
    multiInfect<-multiInfect[multiInfect$n > 1,]
    
    trial<-trial[!duplicated(trial$id),]
    trial[trimesterInfected == 1, czsProb := CZSTrim1]
    trial[trimesterInfected == 2, czsProb := CZSTrim2]
    trial[trimesterInfected == 3, czsProb := CZSTrim3]
    
    trial$czsProb<-ifelse(trial$trimesterInfected == 1 & trial$id %in% multiInfect$id,multiCZSrate,trial$czsProb)
    
    trial$CZS <- suppressWarnings(rbinom(nrow(trial), 1, prob = trial$czsProb)) #CZS yes/no
    
  }
  trial
})


#get symptomatic cases: generate incubation period and calculate time-to-symptoms
symptomatic <- function(trial,parms) with(parms, { 
  trial$symptomatic <- ifelse(trial$status == 1, (rbinom(nrow(trial), 1, prob = symptomRate)), 0) 
  trial$incubationPeriod <- ifelse(trial$symptomatic == 1, rweibull(nrow(trial),incubK, incubC), Inf) 
  trial$survtSymptoms <- ifelse(trial$symptomatic == 1, trial$incubationPeriod + trial$survt, Inf) 
  trial
})


#test results for zikv, given viral persistence in blood
persistence <- function(trial,parms) with(parms, {
  Dates <- seq.Date(as.Date(startDate),as.Date(maxEnddate), by = testInterval) #which dates will be lab test dates
  testDates <- list()
  
  testDates <- Dates - as.Date(startDate)
  
  trial$firstDetectableBlood<-ifelse(trial$symptomatic == 1 & trial$survtSymptoms < (trial$survt + 2), trial$survtSymptoms, trial$survt+2) #assume detectable after 2 days
  trial$lastDetectableBlood<-rweibull(nrow(trial), persistBloodK, persistBloodC) + trial$survt
  
  trial$testResult <- 0
  trial<-as.data.table(trial)
  temp <- trial[status==1, .(id, status, firstDetectableBlood, lastDetectableBlood)]
  temp2 <- temp[rep(1:nrow(temp), each=length(testDates))]
  temp2$testDate <- rep(testDates, nrow(temp))
  temp2[,testResult := in_interval(testDate, firstDetectableBlood, lastDetectableBlood)] 
  temp3 <- temp2[,.(testResult=any(testResult>0)), id]
  posIDs <- temp3[testResult==T, id]
  trial[id %in% posIDs, testResult:=1]
  ## trial[testResult==1,.(id, firstDetectableBlood, lastDetectableBlood)]
  
  if(trialType == 'symptomTrial'){
    trial$testResult<-ifelse(trial$testResult == 1 & trial$symptomatic == 1, 1, 0) 
  }
  trial
})
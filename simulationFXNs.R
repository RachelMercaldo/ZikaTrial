##Simulation functions
## Code for all functions used to simulate trials, called in analyzefxn()
#########################################################################
## Code accompanies:
## 
## Mercaldo, RA, Bellan, SE. Evaluation of Alternative Endpoints for Zika Virus Vaccine Efficacy Trials. 2019.
##
## Rachel Mercaldo, 2019
## License at bottom.

makeParms <- function(
   trialType = c('infTrial','symptomTrial','CZStrial'),
   regSize = seq(50,28000,by=50), #regSize = pop size for one region. 8 regions means total of 8*regSize trial participants
   vaccEff = c(.5,.7,.9), #Vaccine efficacy. 
   startRule = c(1,2,3), 
   symptomRate = 0.225, #Petersen, 2016
   incubK = 3.132, #Lessler, 2016, converted to Weibull parameters
   incubC = 6.632, #Lessler, 2016, converted to Weibull parameter
   persistBloodK = 2.007, #Lessler, 2016, converted to Weibull parameter
   persistBloodC = 11.171, #Lessler, 2016, converted to Weibull parameter
   maxEnddate = '2017-07-23', #last week of infection rate data
   maxCZSdate = '2018-02-25', #last possible birthdate for women pregnant by maxEndDate.
   testInterval = c(7,14,28), #how many days separate lab tests
   CZSTrim1 = 0.15, #Reynolds, 2016
   CZSTrim2 = 0.0227, #Reynolds, 2016
   CZSTrim3 = 0.0227, #Reynolds, 2016
   multiCZSrate = 0.21, #Reynolds, 2016
   TTC = c(TRUE, FALSE), #If CZS, are women recruited Trying To Conceive (TTC)? If FALSE, assumes 1-contraceptionRate for prevalence of women TTC
   startPregRate = 0.20747, #Taylor, 2003, estimated 0.3 for first cycle, 0.05 for 12th. alpha = 0.3, beta = 0.861299. Mean of first 6 months for start.
   contraceptionRate = 0.73, #UN 2015, Latin America and the Caribbean data for women 15-49
   assumedRate = 0.0005, #Low assumed rate of any outcome, for group sequential trial design
   iter = 250 
){
   parms <- expand.grid(as.list(environment()))
   parms<-parms[!(parms$trialType !='infTrial' & (parms$testInterval %in% c(14,28))),] #unless the scenario is an infection trial, remove options for 14/28 days between tests
   parms<-parms[!(parms$trialType != 'CZStrial' & (parms$regSize > 9000)),] #unless the scenario is a CZS trial, remove options for > 9000 regsize 
}


nestPaho <- function(paho){
   temp<-cbind(paho[,7:8],sample(paho[,1:6],4, replace = FALSE)) #choose 4 of 6 countries for this iteration
   temp<-melt(temp, id.vars = c(1,2),
              measure.vars = c(3:6),
              variable.name = "country",
              value.name = "rate")
   temp<-nest(temp,data = c(date, week, rate))
}


#create study population, assign individual risks, and randomize to vaccine or control:
makePop <- function(dat, parms) with(parms, {
   
   if(startRule == 1){
      start_date<-dat$date[dat$rate >= 0.00025][1] 
   }else if(startRule == 2){
      start_date<-dat$date[dat$rate >= 0.0005][1]
   }else if(startRule == 3){
      dat<-dat[1:40,]
      max_date<-dat$date[dat$rate == max(dat$rate)]
      start_date<-max_date + 7
   }
   
   trial<-data.table(startDate = start_date,
                     id = 1:regSize,
                     indRR = rlnorm(regSize, 0, 1),
                     arm = sample(c('vaccine','control'), regSize, replace = TRUE, prob = c(0.5,0.5)))
   
})



#Merge PAHO data:
#Following this function, each participant will have n rows, with n determined by the number of weeks in the trial (dependent on startDate)
#For each week, the participant is assigned a risk (totalRate) that is the product of their individual risk and their regional risk that week
mergeData <- function(temp, paho, parms) with(parms, {
   trial<-unnest(temp, trial)
   trial$id<-paste(substr(trial$country,1,3),trial$id, sep="")
   
   trial<-trial[rep(seq_len(nrow(trial)), each = nrow(paho)),]
   countries<-as.character(temp$country)
   paho<-subset(paho, select = c("date", "week", countries[1:4]))
   paho<-paho[rep(seq_len(nrow(paho)), regSize),]
   
   paho<-gather(paho, 'country','rate', c(3:6))
   
   trial<-cbind(trial,paho)
   
   trial$immuneDate <- as.Date(trial$startDate) + 30
   trial$totalRate <- ifelse(trial$date < as.Date(trial$immuneDate),trial$indRR*trial$rate,
                           trial$indRR*trial$rate*ifelse(trial$arm=='vaccine',1-vaccEff,1))
   
   trial<-trial[trial$date >= trial$startDate,]
   trial<-trial[,c(3:12)]
   
   trial
})

#Helper fxn:
cyclePs<-function(parms,mos = 21) with(parms, {
   probs<-rep(startPregRate,mos)
   for(i in 2:mos){
      probs[i]=probs[i-1]-(probs[i-1]*.1503)
   }
   probs
})


#If a CZS trial, this gets conception time by assigning each woman monthly conception 'risk' and simulating time-to-conception
simPreg <- function(trial,parms) with(parms, {
   
   if(trialType == 'CZStrial'){
      trial$mos<-floor((as.numeric(as.Date(maxEnddate)-as.Date(trial$startDate)))/28)
      
      #pulling out out each country. co1, co2, co3, co4, as trial length will differ for each:
      cos<-unique(trial$country)
      
      pregTrial <- data.frame()
      
      for(i in 1:4){
         co<-trial[trial$country==cos[i],]
         co<-co[!duplicated(co$id),]
         co<-as.data.frame(co[rep(seq_len(nrow(co)),each=co$mos[1]),])
         co$cycleProbs<-cyclePs(parms, co$mos[1])
         co$month<-1:co$mos[1]
         
         pregTrial <- bind_rows(pregTrial,co)
      }
      
      pregTrial$pregTime<-rexp(nrow(pregTrial),rate = pregTrial$cycleProbs)
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
simInf<-function(trial,parms,browse = F) with(parms, { 
   if(browse) browser()
   trial$infectTime <- NA
   trial$infectTime<-suppressWarnings(ifelse(trial$totalRate == 0, Inf, rexp(nrow(trial),rate=trial$totalRate)))
   
   #identify infections prior to protective immunity.
   preImmune<-trial[trial$date < as.Date(trial$immuneDate),] 
   preImmune<-preImmune[preImmune$infectTime<=1,]
   preImmune<-preImmune[!duplicated(preImmune$id),]
   
   trial<-trial[!(trial$date < as.Date(trial$immuneDate)) & !(trial$id %in% preImmune$id),] 
   
   #removes all the weeks prior to immunity and all the participants who were infected before immunity
   
   infected<-trial[trial$infectTime<=1,]
   infected$dup<-ifelse(duplicated(infected$id),1,0)
   infected$dup<-as.numeric(infected$dup)
   dups<-infected[infected$dup==1,]
   
   multiInfect<-infected[infected$id %in% dups$id,]
   if(nrow(multiInfect) > 0) multiInfect$dup<-1
   multiInfect$survt<-as.numeric((multiInfect$date-multiInfect$immuneDate) + 7*multiInfect$infectTime)
   
   infected<-infected[!duplicated(infected$id) & !(infected$id %in% multiInfect$id),] #takes first instance of infectTime <= 1 for any id
   infected$survt<-as.numeric((infected$date - infected$immuneDate) + 7*infected$infectTime)
   
   notInfected<-trial[trial$infectTime>1 & !(trial$id %in% infected$id) & !(trial$id %in% multiInfect$id),]
   notInfected<-notInfected[!duplicated(notInfected$id),] 
   notInfected$dup<-0
   notInfected$survt<- Inf
   
   trial<-rbind(infected,notInfected)
   trial<-rbind(trial, multiInfect)
   
   trial<-trial[order(trial$id,trial$date),]
   trial$survt<-as.numeric(trial$survt)
   
   trial$status <- ifelse(trial$infectTime <= 1, 1, 0)
   
   if(trialType != 'CZStrial'){
      trial<-trial[!duplicated(trial$id),]
   }
   
   out<-list(trial,preImmune) #will separate trial and preImmune dataframes in AnalyzeFXN before applying getCZSoutcome()
   out
})


#CZS, based on trimester infected. 
getCZSoutcome <- function(trial,parms) with(parms, {
   if(trialType == 'CZStrial'){
      
      #If conception time is NA, replace with INF, then calculate birthTime and transform to date:
      trial$conceptionTime <- replace(trial$conceptionTime, is.na(trial$conceptionTime), Inf)
      trial$birthTime <- trial$conceptionTime + 280
      trial$birthDate <- as.Date(trial$startDate) + trial$birthTime
      
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
      
      trial$multiInfected<-ifelse(trial$trimesterInfected==1 & trial$id %in% multiInfect$id, 1, 0)
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

#in_interval helper function for peristence fxn, below
in_interval <- function(x, lower, upper){
   lower <= x & x <= upper
}


#test results for zikv, given viral persistence in blood
persistence <- function(trial,parms, browse = F) with(parms, {
   if(browse) browser()
   cos<-unique(trial$country)
   
   coList <- list()
   
   for(i in 1:4){
      co<-trial[trial$country==cos[i],]
      dat<-co$startDate[1]
      
      codats<-seq.Date(as.Date(dat), as.Date(maxEnddate), by = testInterval)
      testDates<-list()
      testDates<-codats - as.Date(co$startDate[1])
      
      co$firstDetectableBlood<-ifelse(co$symptomatic == 1 & co$survtSymptoms < (co$survt + 2), co$survtSymptoms, co$survt+2) #assume detectable after 2 days
      co$lastDetectableBlood<-rweibull(nrow(co), persistBloodK, persistBloodC) + co$survt
      
      co$testResult <- 0
      co<-as.data.table(co)
      temp <- co[status==1, .(id, status, firstDetectableBlood, lastDetectableBlood)]
      temp2 <- temp[rep(1:nrow(temp), each=length(testDates))]
      temp2$testDate <- rep(testDates, nrow(temp))
      temp2<-temp2[,testResult := in_interval(testDate, firstDetectableBlood, lastDetectableBlood)] 
      temp3 <- temp2[,.(testResult=any(testResult>0)), id]
      temp3 <- temp3[!is.na(temp3$id),]
      posIDs <- temp3[testResult==T, id]
      co$testResult<-ifelse(co$id %in% posIDs, 1, 0)
      
      coList[i] <- list(co)
   }
   
   trial<-do.call('rbind', coList)
   
   
   if(trialType == 'symptomTrial'){
      trial$testResult<-ifelse(trial$testResult == 1 & trial$symptomatic == 1, 1, 0) 
   }
   trial
})

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
##Analysis FXN
## Overall function to run simulation fxns and analyze trial data
################################################################
## Code accompanies:
## 
## Mercaldo, RA, Bellan, SE. Evaluation of Alternative Endpoints for Zika Virus Vaccine Efficacy Trials. 2019.
##
## Rachel Mercaldo, 2019
## License at bottom.

analyzeTrial <- function(parms, browse = F) with(parms, {
   if(browse) browser()
   set.seed(628496)
   
   gsDesArgs = list(k=5, 
                    test.type = 2,
                    alpha = 0.025,
                    beta = 0.1,
                    timing = seq(0, 1, l = 6)[-1])
   gsPlan = do.call(gsDesign, gsDesArgs)
   
   cumRates = assumedRate*c(1,(1-vaccEff))
   fixedSamp = nBinomial(p1 = cumRates[1], p2 = cumRates[2], outtype = 2, beta = 0.1)
   fixed = sum(fixedSamp*cumRates)
   
   contEvents = vaccEvents = date = pvec = as.numeric(list())
   ColFirst = CosFirst = EcuFirst = MexFirst = PanFirst = PerFirst = as.numeric(list())
   ColPreI = CosPreI = EcuPreI = MexPreI = PanPreI = PerPreI = as.numeric(list())
   
   for(i in 1:iter){
      temp<-nestPaho(paho)
      temp$trial<-map(temp$data, makePop, parms)
      trial<-mergeData(temp, paho, parms)
      trial<-simPreg(trial,parms)
      trial<-simInf(trial,parms)
      preImmune<-trial[2][[1]] #all preImmune infections
      trial<-trial[1][[1]]
      trial<-getCZSoutcome(trial, parms)
      trial<-symptomatic(trial,parms)
      
      
      if(trialType == 'infTrial'){
         trial<-persistence(trial,parms, F)
      }else if(trialType == 'symptomTrial'){
         trial<-persistence(trial,parms)
      }else(trial<-trial)
      
      if(trialType == 'infTrial'){
         trial$outcome <- trial$testResult
         trial$time <- as.Date(trial$date)
         maxDate <- as.Date(max(trial$date))
      }else if(trialType == 'symptomTrial'){
         trial$outcome <- trial$testResult
         trial$time <- as.Date(trial$date)
         maxDate <- as.Date(max(trial$date))
         trial$survt<-trial$survtSymptoms
      }else if(trialType == 'CZStrial'){
         trial$outcome <- trial$CZS
         trial$time <- as.Date(trial$birthDate)
         maxDate<-as.Date(maxCZSdate)
      }
      
      days = sort(trial$time[trial$outcome == 1])
      
      trialIter = data.table(events = round(gsPlan$timing*fixed))
      
      trialIter[,tcal := days[events]]
      trialIter <- trialIter[!is.na(tcal)]
      trialIter$trigger = 'events'
      
      if(nrow(trialIter) < gsPlan$k){
         trialIter = rbind(trialIter, data.table(events = sum(trial$outcome), tcal = as.Date(maxDate), trigger = 'end time'))
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
            trialIter = data.table(events = sum(trial$outcome), tcal = as.Date(maxDate), trigger = 'end time', nominalP = 0.025)
         }
      }
      
      analysisNum = 0
      trialStopped = FALSE
      
      while(!trialStopped){
         analysisNum = analysisNum + 1
         analysisDate = trialIter[analysisNum, tcal]
         censTrial = as.data.table(trial)
         cens = censTrial[time > analysisDate, survt := Inf]
         cens = cens[time > analysisDate, outcome := 0]
         
         if(sum(cens$outcome[cens$arm=='vaccine']) == 0){
            standInV <- data.table(time = trialIter$tcal[analysisNum], arm = 'vaccine', survt = as.numeric(trialIter$tcal[analysisNum]-min(cens$immuneDate)), outcome = 1)
            cens <- rbind(cens, standInV, fill=TRUE)
            standInC <- data.table(time = trialIter$tcal[analysisNum], arm = 'control', survt = as.numeric(trialIter$tcal[analysisNum]-min(cens$immuneDate[!is.na(cens$immuneDate)])), outcome = 1)
            cens <- rbind(cens, standInC, fill=TRUE)
         } 
         if(sum(cens$outcome[cens$arm=='control']) == 0){
            standInV <- data.table(time = trialIter$tcal[analysisNum], arm = 'vaccine', survt = as.numeric(trialIter$tcal[analysisNum]-min(cens$immuneDate)), outcome = 1)
            cens <- rbind(cens, standInV, fill=TRUE)
            standInC <- data.table(time = trialIter$tcal[analysisNum], arm = 'control', survt = as.numeric(trialIter$tcal[analysisNum]-min(cens$immuneDate[!is.na(cens$immuneDate)])), outcome = 1)
            cens <- rbind(cens, standInC, fill=TRUE)
         } 
         
         if(trialType == 'infTrial'){
            mod <- try(coxph(Surv(rep(0,nrow(cens)), cens$survt, cens$outcome == 1) ~ cens$arm == 'vaccine', 
                             frailty(cens$country, distribution = 'gamma', sparse = FALSE)), silent = TRUE)
            useMod <- !inherits(mod, 'try-error')
            
            
            trialIter[analysisNum, rawP := summary(mod)$logtest['pvalue']]
            trialIter[analysisNum, logHR := summary(mod)$coefficients[,1]]
         }else if(trialType == 'symptomTrial'){
            mod <- try(coxph(Surv(rep(0,nrow(cens)), cens$survt, cens$outcome == 1) ~ cens$arm == 'vaccine', 
                             frailty(cens$country, distribution = 'gamma', sparse = FALSE)), silent = TRUE)
            useMod <- !inherits(mod, 'try-error')
            
            
            trialIter[analysisNum, rawP := summary(mod)$logtest['pvalue']]
            trialIter[analysisNum, logHR := summary(mod)$coefficients[,1]]
         }else if(trialType == 'CZStrial'){
            mod <- try(glm(outcome ~ arm, family=binomial(link=logit),data=cens))
            useMod <- !inherits(mod, 'try-error')
            
            
            trialIter[analysisNum, rawP := coef(summary(mod))[2,4]]
            trialIter[analysisNum, logHR := coef(summary(mod))[2,1]]
         }
         
         trialIter[analysisNum, vaccGood := rawP <  nominalP & logHR < 0]
         trialIter[analysisNum, vaccBad := rawP < nominalP & logHR > 0]
         
         trialIter$vaccCases[analysisNum] <- nrow(cens[arm=='vaccine' & outcome == 1 & time <= trialIter$tcal[analysisNum],])
         trialIter$contCases[analysisNum] <- nrow(cens[arm=='control' & outcome == 1 & time <= trialIter$tcal[analysisNum],])
         
         earlyStop <- trialIter[analysisNum, vaccGood | vaccBad]
         
         if(any(earlyStop, na.rm=T)) trialStopped <- T
         if(analysisNum==nrow(trialIter)) trialStopped <- T
      }
      trialIter<-trialIter[complete.cases(trialIter),]
      
      dat<-trialIter[max(nrow(trialIter)), `tcal`]
      
      dat<-as.character(dat[[1]][1])
      
      contEvents[i] <- trialIter[max(nrow(trialIter)),`contCases`]
      vaccEvents[i] <- trialIter[max(nrow(trialIter)), `vaccCases`]
      pvec[i] <- trialIter[max(nrow(trialIter)), `rawP`]
      date[i] <- dat
      
      
      ColFirst[i] = nrow(cens[cens$country == 'Colombia' & cens$trimesterInfected == '1',])
      ColPreI[i] = nrow(preImmune[preImmune$country == 'Colombia',])
      
      CosFirst[i] = nrow(cens[cens$country == 'Costa.Rica' & cens$trimesterInfected == '1',])
      CosPreI[i] = nrow(preImmune[preImmune$country == 'Costa.Rica',])
      
      EcuFirst[i] = nrow(cens[cens$country == 'Ecuador' & cens$trimesterInfected == '1',])
      EcuPreI[i] = nrow(preImmune[preImmune$country == 'Ecuador',])
      
      MexFirst[i] = nrow(cens[cens$country == 'Mexico' & cens$trimesterInfected == '1',])
      MexPreI[i] = nrow(preImmune[preImmune$country == 'Mexico',])
      
      PanFirst[i] = nrow(cens[cens$country == 'Panama' & cens$trimesterInfected == '1',])
      PanPreI[i] = nrow(preImmune[preImmune$country == 'Panama',])
      
      PerFirst[i] = nrow(cens[cens$country == 'Peru' & cens$trimesterInfected == '1',])
      PerPreI[i] = nrow(preImmune[preImmune$country == 'Peru',])
      
   }
   out<-data.table(contMean = mean(contEvents),
                   contMedian = median(contEvents),
                   vaccMean = mean(vaccEvents),
                   vaccMedian = median(vaccEvents),
                   power = mean(pvec < 0.05),
                   meanDate = mean(as.Date(date)),
                   medianDate = median(as.Date(date)),
                   ColFirst = mean(ColFirst),
                   ColPreI = mean(ColPreI, na.rm = T),
                   CosFirst = mean(CosFirst),
                   CosPreI = mean(CosPreI,na.rm = T),
                   EcuFirst = mean(EcuFirst),
                   EcuPreI = mean(EcuPreI, na.rm = T),
                   MexFirst = mean(MexFirst),
                   MexPreI = mean(MexPreI, na.rm = T),
                   PanFirst = mean(PanFirst),
                   PanPreI = mean(PanPreI, na.rm = T),
                   PerFirst = mean(PerFirst),
                   PerPreI = mean(PerPreI, na.rm = T)
   )
   
   out
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
require(gsDesign); require(data.table)
seq(0,1, l = 4)[-1] ## i.e. 3 analyses spaced evenly


gsDesArgs = list(k=3,## total number of analyses
                 test.type=2, ## two-tailed symmetric
                 alpha=0.025, ## type I error (always given as one tailed)
                 beta=0.1, ## 1-power (90% power)
                 timing = seq(0,1, l = 4)[-1]) ## spacing of the interim analyses in terms of amount
## of information as a proportion of full informations the [-1] is to exclude 1, since the input is
## the timing of interim (not final) analyses


gsPlan <- do.call(gsDesign, gsDesArgs)
## runs gsDesign with the arguments specified above. do.call is useful because you can change the
## argument list flexibly elsewhere than the call to the function itself
gsPlan

## Then we need to figure out the sample size for the trial that would be necessary to detect a
## vaccine efficacy of X with Y power. 

## Let say our assumptions are below
assumedVaccEff = .7 ## assumed vaccine efficacy
assumedCumHazard = 10^-4 ## average cumulative hazard of infection (or whatever)
beta = .1               ## desired power

?nBinomial ##
## a function for calculating sample size/power etc for non-sequential trial designs (i.e. no
## interim analyses). We use this function to figure out the amount of information/sample size we'd
## want without interim analyses and then have our final analysis happen at that amount of
## information so we're appropriately powered. By adding interim analyses, we reduce our power ever
## so slightly from the amount specified here, but we gain the potential to end the trial earlier.

## cumulative hazard for those vaccinated & not
cumHazs <- assumedCumHazard*c(1, (1-assumedVaccEff)) #
## sample sizes (in terms of number of people) if this were a fixed analysis trial
fixedSampSize_participants  <- nBinomial(p1 = cumHazs[1], p2 = cumHazs[2], outtype=2, beta = beta) 
## sample sizes (in terms of number of *infections*) if this were a fixed analysis trial
fixedSampSize_infections <- sum(sampSizes*cumHazs) 
fixedSampSize_infections ## about 36 infections would be enough.

## Note that the above calculation is fairly insensitive to assumedCumHazard, unless we assume a
## very large cumHazard. Also note that I've been al ittle sloopy up here with cumulative hazards
## (which are proportions) and multiplying them by the vaccine efficacy. Actually we should be
## multiplying the hazard rate itself times the vaccine efficacy. But the main goal here was to get
## an approximate sample size estimate for this vaccine efficacy and power (there are many other
## power calculators out there that we could have used too).

## Now that we know that 36 is the # of infections for the final analysis, we take the
gsPlan

## Pretend that the following gives the event times (infection or symptoms, or whatever our primary endpoint is)
N <- 1000
infrate <- .001 ## infection rate
maxDurationDay <- 30*18 ## end trial after 18 months no matter what
dat <- data.table(id = 1:N, infectDay = rexp(N, rate = infrate), vacc=sample(c(T,F), N, replace=T))
dat[infectDay > maxDurationDay, infectDay := Inf] ## just to make it look more like your simulation
dat$infectDay

## figur out when to do the analyses in terms of now many events have yet happened
intTab <- data.table(events = round(gsPlan$timing * fixedSampSize_infections))
## Get vector of event (infection) timings
infDays <- dat$infectDay[dat$infectDay!=Inf] 
infDays <- sort(infDays)

infDays[36] ## the 36th infection is when we have the last analysis
ceiling(infDays[36]) ## round up

gs <- T

if(gs) { ## if doing a group-sequential design
    ## Calculate interim analyses times
    intTab[, tcal:= ceiling(infDays[events])] ## calendar time (as opposed to information time); *CEILING MEANS MIGHT HAVE MORE INFO*
    intTab <- intTab[!is.na(tcal)]
    ## what is triggering the analysis? reaching a threshold number of events, or the end of the trial?
    intTab$trigger <- 'events' 
    intTab <- intTab[tcal <= maxDurationDay] ## no analyses after maxDurationDay
    ## If didn't do all analyses because we never got to the fixedSampSize_infections, add one more analysis at maximum trial duration
    if(nrow(intTab) < gsPlan$k) intTab <- rbind(intTab, data.table(events=NA, tcal=maxDurationDay, trigger='end time'))
    intTab

    if(nrow(intTab)==gsPlan$k) { ## if we manage to have enough events such that we could do all interim
        ## & final analyses.
        intTab <- cbind(intTab, nominalP = pnorm(gsPlan$lower$bound))
        ## for instance if we need to cross a Z of 3.01, that means that we need the observed P to be < pnorm(-3.01)
        ## In theory it's the same as nominal P in gsPlan, but I had trouble extracting that from the gsDesign object
    }else{ ## if don't have full # of analyses, must readjust design to spend all
        ## remaining alpha at maximum trial duration. We'll do this by rerunning gsDesign where we
        ## specify that we have different timings as originally which includes the analyses we have
        ## already done, and one final analysis.
        gsDesArgsAdj <- within(gsDesArgs, {
            k <- nrow(intTab)
            timing <- c(timing[k-1],1)
        })
        if(gsDesArgsAdj$k>1) { ## if planning a sequential design (otherwise no need to worry about this)
            gsPlanAdj <- do.call(gsDesign, gsDesArgsAdj) ## rerun gsDesign with different timings
            intTab <- cbind(intTab, nominalP = pnorm(gsPlanAdj$lower$bound))
        }else{ ## gsDesARgsAdj$k==1 means we only do one final analysis due to too few events. It's
            ## as if we're doing non-sequential design and we just use regular z quantiles for p
            ## value calculation
            intTab <- data.table(events = NA, tcal = maxDurationDay, trigger = 'end time', nominalP = .025)
        }
    }
}else{ ## gs==F so doing non-sequential design
    intTab <- data.table(events = NA, tcal = maxDurationDay, trigger = 'end time', nominalP = .025)
}
## Add columns to keep track of how many cases there are in each arm and whether we think the
## vaccine is statistically significantly good or bad at any interim/final analysis
intTab$contCases <- intTab$vaccCases <- intTab$rawP <- as.numeric(NA)
intTab$vaccGood <- intTab$vaccBad <- as.logical(NA)

## Do the analyses
analysisNum <- 0
trialStopped <- F
intStats <- list() ## to store results
while(!trialStopped) { 
    analysisNum <- analysisNum+1 ## iterate
    analysisDay <- intTab[analysisNum, tcal]
    censdat <- dat ## copy dat
    censdat[infectDay > analysisDay, infectDay := Inf] ## censor by analysis time
    ## Note: in case of 0-event arms you may need to add 1 event to each arm if this breaks your survival analysis function
    ## 
    ## Then run statistical analysis on censdat. I know you're using a survival model but I'm not going
    ## to write everything out here since I think you're using package survival and I used package
    ## coxme last. but the key here is to record the p value for the vaccine efficacy estimate of
    ## your model.
    pvalFromModel <- .004 ## pretend you've done the analysis here & this is it
    logHRFromModel <- log(.6) ## pretend this is your estimate of the logHR[vacc vs unvacc]
    intTab[analysisNum, rawP:= pvalFromModel]
    intTab[analysisNum, logHR:= logHRFromModel]
    intTab$vaccCases[analysisNum] <- censdat[vacc==T, sum(infectDay<analysisDay)]
    intTab$contCases[analysisNum] <- censdat[vacc==F, sum(infectDay<analysisDay)]
    intTab[analysisNum, vaccGood :=  rawP < ifelse(gs, nominalP, .025) & logHR >0]
    intTab[analysisNum, vaccBad :=  rawP < ifelse(gs, nominalP, .025) & logHR <0]
    earlyStop <- intTab[analysisNum, vaccGood | vaccBad]
    ## Determine whether trial stopped for early stopping or last analysis
    if(any(earlyStop, na.rm=T))  ## returns T if T, F if NA or F
        trialStopped <- T
    if(analysisNum==nrow(intTab)) trialStopped <- T
}
intTab
## in this example the rawP crosses significance in the second interim analysis and you would record
## the trial as having ended on the 26th day of the trial. Note that the number of events aren't
## really mirrored in vaccCases & contCases since I made up that data set and the p values & logHR
## separately to write this quickly.






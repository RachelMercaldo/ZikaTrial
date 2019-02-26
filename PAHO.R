## Code to download and process PAHO data
################################################################
## Code accompanies:
## 
## Mercaldo, RA, Bellan, SE. Evaluation of Alternative Endpoints for Zika Virus Vaccine Efficacy Trials. 2019.
##
## Rachel Mercaldo, 2019
## License at bottom.

library(tidyverse)
library(data.table)
library(mefa)
library(ggplot2)
library(scales)

# Obtain digitized PAHO Zika case numbers, by epidemiological week, from Andersen Github repo:

fileURLs<-c("https://raw.githubusercontent.com/andersen-lab/zika-epidemiology/master/paho_case_numbers/caribbean.csv",
            "https://raw.githubusercontent.com/andersen-lab/zika-epidemiology/master/paho_case_numbers/central_america.csv",
            "https://raw.githubusercontent.com/andersen-lab/zika-epidemiology/master/paho_case_numbers/south_america.csv")
PAHOdata <- as.data.frame(unlist(lapply(fileURLs,fread), recursive=FALSE)) 

# selecting from PAHOdata the 4 best places to be to get Zika, based on projection. 
# 'susp.con.ZIKV.cases' downloaded as the year column, 'V2' is the date (mo/day)  

PAHOdata <- subset(PAHOdata, select = c('Colombia', 'Costa.Rica',  
                                        'Ecuador', 'Mexico', 
                                        'Panama', 'Peru',
                                        'susp.con.ZIKV.cases','V2')) 

#Reformating data to remove extra row in header and to include a single date column

PAHOdata <- setnames(PAHOdata[-1,], c('susp.con.ZIKV.cases','V2'), c('year','date'))

PAHOdata[1:52,'year']<-'2016'; PAHOdata[53:nrow(PAHOdata),'year']<-2017  #if Andersen lab continues updating repo, 
#may need to update this line
PAHOdata$date <- as.Date(paste(PAHOdata$date,'-',PAHOdata$year, sep=""),"%d-%b-%Y")

PAHOdata<-select(PAHOdata, -year) #get rid of the extra column for year
PAHOdata$Week<-1:nrow(PAHOdata)


#Only 2016 rates (first 52 weeks of data):
PAHOdata16<-PAHOdata[1:52,]
PAHOdata17<-PAHOdata[53:82,]

Rates2016<-c(sum(PAHOdata16$Colombia),
             sum(PAHOdata16$Costa.Rica),
             sum(PAHOdata16$Ecuador),
             sum(PAHOdata16$Mexico),
             sum(PAHOdata16$Panama),
             sum(PAHOdata16$Peru))

Rates2017<-c(sum(PAHOdata17$Colombia),
             sum(PAHOdata17$Costa.Rica),
             sum(PAHOdata17$Ecuador),
             sum(PAHOdata17$Mexico),
             sum(PAHOdata17$Panama),
             sum(PAHOdata17$Peru))

# PRF is the scale to multiply 2016 Andersen data to represent the projected incidence rate while 
# maintaining seasonal variation seen country-wide in the first wave of the epidemic
PRF16 <- 0.01/Rates2016 
PRF17 <- (0.01*(nrow(PAHOdata17)/52))/Rates2017


paho16<-data.table(Colombia = PAHOdata16$Colombia*PRF16[1], 
                   Costa.Rica = PAHOdata16$Costa.Rica*PRF16[2], 
                   Ecuador = PAHOdata16$Ecuador*PRF16[3], Mexico = PAHOdata16$Mexico*PRF16[4], 
                   Panama = PAHOdata16$Panama*PRF16[5], Peru = PAHOdata16$Peru*PRF16[6])

paho17<-data.table(Colombia = PAHOdata17$Colombia*PRF17[1], 
                   Costa.Rica = PAHOdata17$Costa.Rica*PRF17[2], 
                   Ecuador = PAHOdata17$Ecuador*PRF17[3], Mexico = PAHOdata17$Mexico*PRF17[4], 
                   Panama = PAHOdata17$Panama*PRF17[5], Peru = PAHOdata17$Peru*PRF17[6])
# Use PRF to scale PAHO data
# Only 2016 was used to make scale factor, 2016 and 2017 both scaled


paho<-rbind(paho16,paho17)
paho$date <- PAHOdata[1:82,'date']
paho$week <- 1:nrow(paho)

save(paho, file = 'paho.Rdata')

##visualize 
#Scaled:
pahoPlot<-melt.data.table(paho,id.vars = 8, measure.vars = 1:6,variable.name = "Country")
ggplot(pahoPlot, aes(x=week, y=value, group=Country, color=Country))+geom_line(size = 1) + labs(x = 'Week', y = 'Rate') + 
   theme_classic() 

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
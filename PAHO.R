library(tidyverse)
library(data.table)
library(mefa)
library(ggplot2)
library(scales)

# Obtain digitized PAHO Zika case numbers, by epidemiological week, from Andersen Github repo:

fileURLs<-c("https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/caribbean.csv",
            "https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/central_america.csv",
            "https://raw.githubusercontent.com/andersen-lab/Zika-cases-PAHO/master/south_america.csv")
PAHOdata <- as.data.frame(unlist(lapply(fileURLs,fread), recursive=FALSE)) 

# selecting from PAHOdata the 4 best places to be to get Zika, based on projection. 
# 'susp.con.ZIKV.cases' downloaded as the year column, 'V2' is the date (mo/day)  

PAHOdata <- subset(PAHOdata, select = c('Colombia', 'Ecuador', 'Mexico', 'Peru','susp.con.ZIKV.cases','V2')) 

#Reformating data to remove extra row in header and to include a single date column

PAHOdata <- setnames(PAHOdata[-1,], c('susp.con.ZIKV.cases','V2'), c('year','date'))
PAHOdata[1:52,'year']<-'2016'; PAHOdata[53:nrow(PAHOdata),'year']<-2017  #if Andersen lab continues updating repo, 
#may need to update this line
PAHOdata$date <- as.Date(paste(PAHOdata$date,'-',PAHOdata$year, sep=""),"%d-%b-%Y")

#These are the population sizes for the four countries identified in projections:
popSize<-c(48650000,16506000,130624000,31970000) #PAHO, 2016


#Making it an infection per-capita, by region, by multiplying 4.5 (due to asymptomatic rate) and dividing by pop size
for(reg in 1:(ncol(PAHOdata)-2)){
  PAHOdata[,reg]<-PAHOdata[,reg]*4.5/popSize[reg]
}


PAHOdata<-select(PAHOdata, -year) #get rid of the extra column for year


# 8 areas identified by projections within these 4 countries: 'Narino' in Colombia,'Sucumbios' in Ecuador,'Sinaloa' 
# and 'Tamaulipas' in Mexico, 'Piura','Tumbes','SanMartin' and 'Ucayali' in Peru

#Only 2016 rates (first 52 weeks of data):
PAHOdata16<-PAHOdata[1:52,]

ExpectedRates<-c(0.005, 0.0755, 0.0244, 0.0028, 0.071, 0.0631, 0.001, 0.0006) #ZIKAVAT collaboration projected median rate for 2017

Rates2016<-c(sum(PAHOdata16$Colombia), 
             sum(PAHOdata16$Ecuador),
             sum(PAHOdata16$Mexico),
             sum(PAHOdata16$Mexico),
             sum(PAHOdata16$Peru),
             sum(PAHOdata16$Peru),
             sum(PAHOdata16$Peru),
             sum(PAHOdata16$Peru))

# PRF is the scale to multiply 2016 Andersen data to represent the projected incidence rate while 
# maintaining seasonal variation seen country-wide in the first wave of the epidemic
PRF <- ExpectedRates/Rates2016 


# Use PRF to scale PAHO data
# Only 2016 was used to make scale factor, 2016 and 2017 both scaled
paho<-data.table(Narino = PAHOdata$Colombia*PRF[1], Sucumbios = PAHOdata$Ecuador*PRF[2], Sinaloa = PAHOdata$Mexico*PRF[3], 
                 Tamaulipas = PAHOdata$Mexico*PRF[4], Piura = PAHOdata$Peru*PRF[5], Tumbes = PAHOdata$Peru*PRF[6], 
                 SanMartin = PAHOdata$Peru*PRF[7], Ucayali = PAHOdata$Peru*PRF[8], date = PAHOdata$date)
paho$Week <- 1:nrow(paho)

##save
#save(paho, file='paho.Rdata')


##visualize

#2016 original:
pahoPlot<-paho[1:52,]
pahoPlot<-melt.data.table(pahoPlot,id.vars = 10, measure.vars = 1:8,variable.name = "Region")
ggplot(pahoPlot, aes(x=Week, y=value, group=Region, color=Region))+geom_line() + labs(y='Rate') + 
  theme_classic() + scale_y_continuous(limits = c(0,.01))

#2016 scaled:
PAHOdata16$Week<-1:nrow(PAHOdata16)
paho16Plot<-melt.data.table(as.data.table(PAHOdata16), id.vars = 6, measure.vars = 1:4, variable.name = 'Countries')
ggplot(paho16Plot, aes(x=Week,y=value, group=Countries, color=Countries)) + geom_line() + labs(y='Rate') + 
  theme_classic() + scale_y_continuous(labels = comma)

#All data, scaled: 
allpaho<-melt.data.table(paho,id.vars = 10, measure.vars = 1:8,variable.name = "Region")
ggplot(allpaho, aes(x=Week, y=value, group=Region, color=Region))+geom_line() + labs(y='Rate') + 
  theme_classic()

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
PAHOdata$Week<-1:nrow(PAHOdata)

# 8 areas identified by projections within these 4 countries: 'Narino' in Colombia,'Sucumbios' in Ecuador,'Sinaloa' 
# and 'Tamaulipas' in Mexico, 'Piura','Tumbes','SanMartin' and 'Ucayali' in Peru

#Only 2016 rates (first 52 weeks of data):
PAHOdata16<-PAHOdata[1:52,]
PAHOdata17<-PAHOdata[53:82,]

Rates2016<-c(sum(PAHOdata16$Colombia), 
             sum(PAHOdata16$Ecuador),
             sum(PAHOdata16$Mexico),
             sum(PAHOdata16$Peru))

Rates2017<-c(sum(PAHOdata17$Colombia),
             sum(PAHOdata17$Ecuador),
             sum(PAHOdata17$Mexico),
             sum(PAHOdata17$Peru))

# PRF is the scale to multiply 2016 Andersen data to represent the projected incidence rate while 
# maintaining seasonal variation seen country-wide in the first wave of the epidemic
PRF16 <- 0.01/Rates2016 
PRF17 <- (0.01*(nrow(PAHOdata17)/52))/Rates2017


paho16<-data.table(Colombia = PAHOdata16$Colombia*PRF16[1], Ecuador = PAHOdata16$Ecuador*PRF16[2],
                   Mexico = PAHOdata16$Mexico*PRF16[3], Peru = PAHOdata16$Peru*PRF16[4])

paho17<-data.table(Colombia = PAHOdata17$Colombia*PRF17[1], Ecuador = PAHOdata17$Ecuador*PRF17[2],
                   Mexico = PAHOdata17$Mexico*PRF17[3], Peru = PAHOdata17$Peru*PRF17[4])
# Use PRF to scale PAHO data
# Only 2016 was used to make scale factor, 2016 and 2017 both scaled


paho<-rbind(paho16,paho17)
paho$date <- PAHOdata[1:82,'date']
paho$week <- 1:nrow(paho)

save(paho, file = 'paho.Rdata')

##visualize (shown for only paho, not paho2 or paho50)

#Original:
pahoDataPlot<-melt.data.table(as.data.table(PAHOdata), id.vars = 6, measure.vars = 1:4, variable.name = 'Countries')
ggplot(pahoDataPlot, aes(x=Week,y=value, group=Countries, color=Countries)) + geom_line() + labs(y='Rate') + 
  theme_classic() + scale_y_continuous(labels = comma)

#Scaled:
pahoPlot<-melt.data.table(paho,id.vars = 6, measure.vars = 1:4,variable.name = "Country")
ggplot(pahoPlot, aes(x=week, y=value, group=Country, color=Country))+geom_line() + labs(x = 'Week', y = 'Rate') + 
  theme_classic() 


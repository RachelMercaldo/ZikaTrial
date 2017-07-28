library(data.table)
library(ggplot2)

set.seed(628496)

#for simplicity, just putting the shape parameters here. They
#  are based on the mean/variance of incidence rates per 100k peeps per year in 
#  five PAHO regions.


alpha <- 0.4065
beta <- 172.65  #is this huge? This seems huge.


#Creating our population and keeping them in a df for future activities:

createPop <- function(nPerClus, nClus){
  data.frame(cluster = rep(1:nClus, each=nPerClus),
             clusterID = rep(1:nPerClus, nClus),
             studyID = 1:(nClus*nPerClus))
}


#Assigning treatment like last time, calling it immuneStatus so we can add in a 
#  vaccinatedNotProtected group sometime soon as discussed. 
#  For now, splitting them 50/50 within clusters

getImmuneStatus <- function(df){
  immunity <- setkey(setDT(df), cluster)[,list(immuneStatus=rbinom(clusterID, 1, prob=0.500)), by=cluster]
  
  factor(immunity$immuneStatus, labels = c('control','vacProtected'))
}


#The big one! The one that I'm most unsure of, anyway. I used the beta distribution
#  for cluster-level risk as discussed, and then shoved the resulting rate into 
#  log-normal, which I just wrote out since I am not certain if there is some 
#  smart way of doing it - working harder, not smarter.

getProb <- function(df, a = alpha, b = beta, time, eff){
  df$zikaProb <- NA
  df$temp<-NA
  for(clus in unique(df$cluster)){
    df$temp[df$cluster==clus] <- rbeta(1, a, b)
    df$zikaProb <- ifelse(df$immuneStatus=='control', 
                                          1-exp(-(df$temp*time)), 
                                          1-exp(-(((1-eff)*(df$temp))*time)))
  }
  df$zikaProb
}


#Actually calculating power for this pop: pretty much the same as the last time. 
#  I did attempt to add a way to track the number of cases in each arm.

powerFxn <- function(df,iter=100){
  p<-rep(NA,iter)
  nCaseVac <- rep(NA, iter)
  nCaseCon <- rep(NA, iter)
  
  for(i in 1:iter){
    df$zika <- rbinom(df$studyID, 1, prob = df$zikaProb)
    nCaseVac[i] <- nrow(df[df$zika == 1 & df$immuneStatus == 'vacProtected'])
    nCaseCon[i] <- nrow(df[df$zika == 1 & df$immuneStatus == 'control'])
    
    mod <- glm(zika ~ immuneStatus, family = binomial(link = 'logit'), data = df, control=glm.control(maxit=50))   
    p[i] <- coef(summary(mod))[2,4]
  }
  list(mean(p<0.05), p, nCaseVac,nCaseCon)
}


#Example:

nPerClus<-data.frame('NumberInCluster' = seq(100,500,25), 'power' = NA)

for(i in nPerClus$NumberInCluster){
  pop <- createPop(i,nClus=20)
  pop$immuneStatus <- getImmuneStatus(pop)
  pop$zikaProb <- getProb(pop, alpha, beta, 1, eff=0.80)
  power <- powerFxn(pop, iter = 1000)
  nPerClus$power[nPerClus$NumberInCluster==i] <- unlist(power[1])
}

ggplot(data=nPerClus, aes(x=NumberInCluster, y=power)) + geom_line() 

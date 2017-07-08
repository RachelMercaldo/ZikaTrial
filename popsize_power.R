## First two mean/sd  risk values based on location and stolen from PAHO data 
## Others are czs means, from CDC and the articles sitting on my desk in Athens
## (I can fill in the inevitable gaps when I'm back)

## symptoms (leading to "suspected" and "confirmed"):  0.000839, sd 0.000732  (6 mo)
## Infection (ends up basically being 4 to 5* the above): 0.003776, sd 0.0033  (6 mo)
## CZS, any trimester: 0.06 
## CZS, first trimester: 0.11

library(dplyr) 
library(truncnorm) 
library(ggplot2)

pop<-data.frame(npop=seq(1000,15000,by=100)) #numbers that actually ran on stupid laptop

set.seed(628496)

#this fxn requires one to plug in stuff:
fxn <- function(pop, mean = 0.003776, sd = 0.0033, eff = 0.8, i = 2000){
  df <- data.frame(id=1:pop)
  
  df$treatment <- factor(rbinom(df$id, 1, prob = 0.500), labels = c('vaccine','control'))
  
  df$risk <- ifelse(df$treatment=='control',rtruncnorm(pop, a=0, b=1, mean = mean, sd = sd) 
                                           , (1-eff)*rtruncnorm(pop, a=0, b=1, mean = mean, sd = sd))
  
  p<-list(rep(NA,i))
  
  for(j in 1:i){
    df$zika <- rbinom(df$id, 1, prob = df$risk) 
    mod <- glm(zika ~ treatment, family = binomial(link = 'logit'), data = df, control=glm.control(maxit=50))   
    p[j] <- coef(summary(mod))[2,4]
  }
  return(length(p[p < 0.05])/length(p))
}

pop$power<-apply(pop,1,fxn) #probably a million faster ways of doing this?


pop %>% ggplot(aes(x=npop,y=power)) + geom_line() 



#A terribly slow but slightly understandable way to run all possible values of pop#, 
## eff, etc, that you may want

pops <- seq(10,20000,by=1) #Very intense. 
effs <- c(seq(0.5, 0.95, by = .15), .99)
means <- c(0.000839, 0.003776, 0.06, 0.11) 
sds <- c(0.000732, 0.0033, 0, 0) 

dat <- expand.grid(npop = pops, eff = effs, mean = means)

#"why am I up at 3am?" approach to getting the SDs in there:
dat$sd<-ifelse(dat$mean==0.000839,0.000732,NA) 
dat$sd<-ifelse(dat$mean==0.003776,0.0033,dat$sd)
dat$sd<-ifelse(is.na(dat$sd),0,dat$sd)

head(dat)

#Now dat is a data frame of 399,820 obs: all possible combinations of pop #, mean/sd, eff
#The fxn below is the same as above, except it loops through dat instead of through pop # 
## with fixed eff/inc/sd values like the version above. 


set.seed(628496)

bigFxn <- function(dat, i = 99){      #set i to 99 because I have no patience
  df <- data.frame(id=1:dat[1])
  
  df$treatment <- factor(rbinom(df$id, 1, prob = 0.500), labels = c('vaccine','control'))
  
  df$risk <- ifelse(df$treatment=='control',rtruncnorm(dat[1], a=0, b=1, mean = dat[3], sd = dat[4]) 
                    , (1-dat[2])*rtruncnorm(dat[1], a=0, b=1, mean = dat[3], sd = dat[4]))
  
  p<-list(rep(NA,i))
  
  for(j in 1:i){
    df$zika <- rbinom(df$id, 1, prob = df$risk) 
    mod <- glm(zika ~ treatment, family = binomial(link = 'logit'), data = df, control=glm.control(maxit=50))   
    p[j] <- coef(summary(mod))[2,4]
  }
  return(length(p[p < 0.05])/length(p))
}

#dat$power<-apply(dat,1,bigFxn)  ##don't run it unless you want to wait forever


#dat %>% ggplot(aes(x=npop,y=power)) + geom_line() + facet_wrap(mean ~ eff)


##wee example with only 2 effs/means and teeny # of pop values

p<-seq(1000,10000,by=1000)
e<-c(0.5,0.8)
m<-c(0.0008,0.004)

dat2<-expand.grid(npop = p, eff = e, mean = m)

dat2$sd<-ifelse(dat2$mean==0.0008,0.0007,0.003)

dat2 #40 obs

dat2$power<-apply(dat2,1,bigFxn)

dat2 %>% ggplot(aes(x=npop,y=power)) + geom_line() + facet_wrap(mean ~ eff)

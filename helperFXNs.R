in_interval <- function(x, lower, upper){
  lower <= x & x <= upper
}

cycleProbs<-function(start=startPregRate,mos=19){
  probs<-rep(start,mos)
  #probs[1]<-start
  for(i in 2:mos){
    probs[i]=probs[i-1]-(probs[i-1]*.1387)
  }
  probs
}

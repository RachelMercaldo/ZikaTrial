in_interval <- function(x, lower, upper){
  lower <= x & x <= upper
}


cycleProbs<-function(parms,mos=19,browse=F) with(parms,{
  if(browse) browser()
  probs<-rep(startPregRate,mos)
  for(i in 2:mos){
    probs[i]=probs[i-1]-(probs[i-1]*.1387)
  }
  probs
})

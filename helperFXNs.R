in_interval <- function(x, lower, upper){
  lower <= x & x <= upper
}


cyclePs<-function(parms,browse=F) with(parms, {
  if(browse) browser()
  probs<-rep(startPregRate,mosLengthOut)
  for(i in 2:mosLengthOut){
    probs[i]=probs[i-1]-(probs[i-1]*.1387)
  }
  probs
})

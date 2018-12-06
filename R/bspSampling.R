##BSPSampling.R

bspSamples<-function(bsp, reps=10000){
  draws<-matrix(0,ncol = reps, nrow=length(bsp$support))
  alpha1=bsp$precisionAt[1]*bsp$centeringMeasure[1]
  beta1=bsp$precisionAt[1]*(1-bsp$centeringMeasure[1])
  draws[1,]<-rbeta(reps, alpha1, beta1)
  for (i in 2:length(bsp$support)){
    alphak<-bsp$precisionAt[i]*(bsp$centeringMeasure[i]-bsp$centeringMeasure[i-1])
    betak=bsp$precisionAt[i]*(1-(bsp$centeringMeasure[i]))
    epsk<-1-draws[i-1,]
    #if (any(epsk>1 | epsk<0))print(i)
    draws[i,]<-draws[i-1,]+rbeta(reps, alphak, betak)*epsk
  }
  draws
}
a=bspSamples(bsp, reps=15)
matplot(a, type="l")


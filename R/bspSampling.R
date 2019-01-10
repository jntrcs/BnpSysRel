##BSPSampling.R


#'Creates a matrix where each column represents a sampled CDF from the provided BSP
#'
#'@param bsp Any bsp object
#'@param reps Number of samples to draw, default 10000
#'
#'@return A matrix. Each column can be plotted against the BSP's support to view a sampled CDF from the BSP
#'@export
#'
#' @examples
#' bsp=bsp(c(1:3), centeringMeasure = c(.1,.9, .98), precision = 2)
#' samples=bspSamples(bsp)
#' #To visualize
#' matplot(samples, bsp$support, type='l')
#'
bspSampling<-function(bsp, reps=10000){
  draws<-matrix(0,ncol = reps, nrow=length(bsp$support))
  firstTime<-bsp$support[1]
  alpha1=evaluate_precision(bsp, firstTime)*bsp$centeringMeasure[1]
  beta1=evaluate_precision(bsp, firstTime)*(1-bsp$centeringMeasure[1])
  draws[1,]<-rbeta(reps, alpha1, beta1)
  for (i in 2:length(bsp$support)){
    prec<-evaluate_precision(bsp, bsp$support[i])
    alphak<-prec*(bsp$centeringMeasure[i]-bsp$centeringMeasure[i-1])
    betak=prec*(1-(bsp$centeringMeasure[i]))
    epsk<-1-draws[i-1,]
    #if (any(epsk>1 | epsk<0))print(i)
    draws[i,]<-draws[i-1,]+rbeta(reps, alphak, betak)*epsk
  }
  draws
}



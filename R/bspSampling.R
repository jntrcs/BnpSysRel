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
  draws<-matrix(0,ncol = reps, nrow=length(bsp$support)-1)
  ##skip the first spot where it's 0
  firstTime<-bsp$support[2]
  alpha1=evaluate_precision(bsp, firstTime)*bsp$centeringMeasure[2]
  beta1=evaluate_precision(bsp, firstTime)*(1-bsp$centeringMeasure[2])
  draws[1,]<-rbeta(reps, alpha1, beta1)
  for (i in 3:(nrow(draws)+1)){
    prec<-evaluate_precision(bsp, bsp$support[i])
    alphak<-prec*(bsp$centeringMeasure[i]-bsp$centeringMeasure[i-1])
    betak=prec*(1-(bsp$centeringMeasure[i]))
    epsk<-1-draws[i-2,]
    #if (any(epsk>1 | epsk<0))print(i)
    draws[i-1,]<-draws[i-2,]+rbeta(reps, alphak, betak)*epsk
  }
  draws
}


#'Finds approximate (simulated) confidence intervals for the percent failed
#'at a given time
#'
#'@param bsp The bsp object
#'@param times Place where you would like confidence intervals
#'
#'@note If you intend to call this function multiple times on the same object,
#'calling bsp$Samples<-bspSampling(bsp) before bspConfint
#' will increase performance.
#'
#'@return A 2xlength(times) matrix with upper and lower bounds of the confidence interval
#'@export
#'
#' @examples
#' bsp=bsp(c(1:3), centeringMeasure = c(.1,.9, .98), precision = 2)
#' bspConfint(bsp, 1:3)
#'
bspConfint<-function(bsp, times, alpha=.05){
  if(is.null(bsp$Samples))bsp$Samples<-bspSampling(bsp, reps=1000)
  rows<-sapply(times, FUN=function(time)sum(time>=bsp$support)-1)
  intervals= sapply(rows, FUN=
        function(row) quantile(bsp$Samples[row, ], c(alpha/2, 1-alpha/2)))
  intervals[is.na(intervals)]<-0
  intervals
}

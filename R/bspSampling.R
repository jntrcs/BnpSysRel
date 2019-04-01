##BSPSampling.R

#just a test change
#'Creates a matrix where each column represents a sampled CDF from the provided BSP
#'
#'@param bsp Any bsp object
#'@param reps Number of samples to draw, default 10000
#'
#'@return A matrix. Each column can be plotted against
#'the BSP's support to view a sampled CDF from the BSP.
#'Rownames are the support point for that row. Use samplesAt to
#'get extract sampled distribution at single point
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
  rownames(draws)<-bsp$support[-1]
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
#'@param conf.level Alpha level for confidence (default .05)
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
bspConfint<-function(bsp, times, conf.level=.05){
  if(is.null(bsp$Samples))bsp$Samples<-bspSampling(bsp, reps=1000)
  if (conf.level<0 |conf.level>1)stop("conf.level must be between 0 and 1")
  if (conf.level>.6)warning(paste0("Alpha level set to ", conf.level,
                                   ". Are you sure you didn't mean ",
                                   1-conf.level,"?"))
  rows<-sapply(times, FUN=function(time)sum(time>=bsp$support)-1)
  intervals= sapply(rows, FUN=
        function(row) quantile(bsp$Samples[row, ], c(conf.level/2, 1-conf.level/2)))
  intervals[is.na(intervals)]<-0
  colnames(intervals)<-times
  intervals
}

#'Extract distribution at select point
#'
#'@param samples The samples produced by bspSampling()
#'@param time The spot on the support at which you would like to examine the distribution
#'
#'@return A vector of samples drawn from that point on the support
#'
#'@export
#'
#'@examples
#'bsp=bsp(1:3, c(.2,.5,.8), 3)
#'samples=bspSampling(bsp)
#'marginal = sampleAt(samples, 2.5)
samplesAt<-function(samples, time){
  times<-as.numeric(rownames(samples))
  index<-sum(times<=time)
  if(index==0)stop("Time less than all of support")
  return(samples[index,])
}


#' #'Finds approximate confidence interval by fitting a beta distribution to the BSP moemnts
#' #'
#' #'
#' bspConfint2<-function(bsp,times, conf.level=.05){
#'   E1<-evaluate_centering_measure(bsp, times)
#'   E2<-evaluate_second_moment(bsp, times)
#'   v=E2-E1^2
#'   alpha<-E1*(E1*(1-E1)/v-1)
#'   beta<-alpha*(1/E1-1)
#'   uppers<-qbeta(1-conf.level/2, alpha, beta)
#'   lowers<-qbeta(conf.level/2, alpha, beta)
#'   return(rbind(lowers,uppers))
#'
#' }
